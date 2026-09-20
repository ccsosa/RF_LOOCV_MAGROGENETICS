#' Spatial Prediction Workflow for Genetic Heterozygosity (Ho) using Support Vector Regression
#'
#' @description
#' Executes an end-to-end spatial modeling and prediction pipeline for observed genetic 
#' heterozygosity (\emph{Ho}) across geographic landscapes. The workflow integrates spatial 
#' aggregation of genetic point data, environmental raster decorrelation, hyperparameter-tuned 
#' Support Vector Regression (SVR), cross-validation, extrapolation filtering via Multivariate 
#' Environmental Similarity Surface (MESS), and automated cartographic generation.
#'
#' @details
#' The pipeline performs the following sequential steps:
#' \enumerate{
#'   \item \strong{Data Loading & Spatial Alignment:} Reads bioclimatic rasters and genetic sample records, aligning off-coverage coordinates to the nearest valid raster grid cell.
#'   \item \strong{Grid Aggregation:} Aggregates genetic sampling points located within the same raster cell, computing mean observed heterozygosity (\emph{Ho}) and sample weights.
#'   \item \strong{Feature Selection:} Removes collinear predictors using a absolute correlation cutoff (\code{cor_cutoff}) applied to bioclimatic variables.
#'   \item \strong{Model Training & Tuning:} Fits an SVR model (\code{svmRadial} or \code{svmLinear}) via \code{caret::train} using sample counts as fitting weights, optimized over a hyperparameter grid.
#'   \item \strong{Cross-Validation:} Evaluates generalization performance using either Leave-One-Out Cross-Validation (LOOCV) or Repeated $K$-Fold CV ($5 \times 10$).
#'   \item \strong{Raster Prediction & MESS Masking:} Predicts \emph{Ho} across the target species distribution area (SDM) in memory-safe chunks and masks non-analogous environmental regions ($MESS \le 0$).
#'   \item \strong{Export:} Writes clean datasets, model performance metrics, scatterplots, spatial prediction rasters (\code{.tif}), thematic maps (\code{.pdf}), and complete workspace states (\code{.RData}).
#' }
#'
#' @param outdir \code{character}. Output directory path where all metrics, figures, rasters, and workspace files will be saved.
#' @param sp_name \code{character}. Exact target species name used to filter records in the input Excel dataset.
#' @param raster_dir \code{character}. Directory path containing individual predictor raster layers in \code{.tif} format. Must include \code{"lon.tif"} and \code{"lat.tif"}.
#' @param data_path \code{character}. File path to the Excel file (\code{.xlsx}) containing genetic data, species names (\code{sp}), observed heterozygosity (\code{Ho}), and geographic coordinates (\code{lon}, \code{lat}).
#' @param sdm_path \code{character}. File path to the binary Species Distribution Model / Ecological Niche Model raster (\code{.tif}) used to constrain geographic predictions.
#' @param n_cores \code{integer}. Number of processing cores allocated for parallel model training via \code{doParallel}. Default is \code{6L}.
#' @param cor_cutoff \code{numeric}. Absolute pairwise Pearson correlation threshold used for removing collinear predictors (e.g., \code{0.5}). Default is \code{0.5}.
#' @param addLonLat \code{logical}. If \code{TRUE}, explicit geographic coordinates (\code{lon} and \code{lat}) are retained as model predictors alongside bioclimatic variables. Default is \code{FALSE}.
#' @param use_loocv \code{logical}. Cross-validation scheme controller. If \code{TRUE}, executes Leave-One-Out Cross-Validation (LOOCV). If \code{FALSE}, runs a 5-fold CV repeated 10 times. Default is \code{FALSE}.
#' @param kernel_type \code{character}. Support Vector Machine kernel specification passed to \code{caret}. Supported options are \code{"svmRadial"} (Radial Basis Function / RBF) and \code{"svmLinear"}. Default is \code{"svmRadial"}.
#'
#' @return Invisibly returns \code{NULL}. All outputs (metrics CSVs, figures, GeoTIFFs, PDF maps, and \code{.RData} environments) are written to \code{outdir}.
#'
#' @author Jorge (Tesis Macrogenética)
#' @keywords spatial-prediction macrogenetics SVR caret terra tmap
#'
#' @examples
#' \dontrun{
#' # Run pipeline using RBF Kernel with Repeated Cross-Validation
#' svr_ho_pipeline(
#'   outdir      = "C:/Research/Output",
#'   sp_name     = "Crocodylus moreletii",
#'   raster_dir  = "C:/Research/Rasters",
#'   data_path   = "C:/Research/Data/Genetic_Data.xlsx",
#'   sdm_path    = "C:/Research/SDM/C_moreletii_P10.tif",
#'   n_cores     = 8L,
#'   cor_cutoff  = 0.5,
#'   addLonLat   = FALSE,
#'   use_loocv   = FALSE,
#'   kernel_type = "svmLinear"
#' )
#' }
#' @export
# ==============================================================================
# REQUIRED LIBRARIES / SCRIPT DEPENDENCIES
# ==============================================================================
library(dplyr)        # Data manipulation and %>% pipe operators
library(readxl)       # Reading Excel files (.xlsx)
library(sf)           # Handling spatial vector data
library(terra)        # High-performance spatial raster operations and processing
library(caret)        # SVR model training, tuneGrid, and feature selection
library(doParallel)   # Parallel processing execution
library(kernlab)      # Backend engine executing 'svmLinear' in caret
library(yardstick)    # Model performance metrics calculation (RMSE, MAE, R2)
library(ggplot2)      # Data visualization and plotting
library(ggpmisc)      # Annotating equations and R2 in plots (stat_poly_eq)
library(patchwork)    # Combining multiple plots (p_train + p_test)
library(raster)       # Converting 'terra' objects to 'raster' for dismo support
library(dismo)        # MESS (Multivariate Environmental Similarity Surface) masking
library(RColorBrewer) # Color palettes for spatial maps
library(tmap)         # Thematic map generation and export

svr_ho_pipeline <- function(outdir, 
                            sp_name, 
                            raster_dir, 
                            data_path, 
                            sdm_path, 
                            n_cores = 6L,
                            cor_cutoff = 0.5,
                            addLonLat = FALSE,
                            use_loocv = FALSE,
                            kernel_type = "svmLinear") {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  if(isTRUE(addLonLat)){
    message("Using lon and lat as predictors")
  }
  
  # ------------------------------------------------------------------------------
  # 1. RASTER & GENETIC DATA LOADING
  # ------------------------------------------------------------------------------
  bios_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)
  bios_names <- sub("\\.tif$", "", list.files(raster_dir, pattern = "\\.tif$"))
  
  bios <- terra::rast(bios_files)
  names(bios) <- bios_names
  
  if (!all(c("lon", "lat") %in% names(bios))) {
    stop("Error: 'lon' and/or 'lat' layers are missing from the raster directory.")
  }
  
  # Filtrado básico: solo valores validos y Ho > 0
  data <- readxl::read_xlsx(data_path, sheet = "data") %>%
    dplyr::filter(sp == sp_name, !is.na(Ho), Ho > 0, !is.na(lat), !is.na(lon))
  
  my_sf_object <- sf::st_as_sf(data, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  points_vect  <- terra::vect(my_sf_object)
  
  # ------------------------------------------------------------------------------
  # 2. CELL ASSIGNMENT & SPATIAL AGGREGATION
  # ------------------------------------------------------------------------------
  cell_ids <- terra::cells(bios[[1]], points_vect)[, "cell"]
  na_cells <- which(is.na(cell_ids))
  
  if (length(na_cells) > 0) {
    message(sprintf("Assigning %d off-coverage points to nearest valid raster cell...", length(na_cells)))
    extracted_cells <- terra::extract(
      bios[[1]],
      points_vect[na_cells, ],
      buffer = 5000,
      fun = "first",
      cells = TRUE,
      exact = FALSE
    )
    cell_ids[na_cells] <- extracted_cells$cell
  }
  
  my_sf_object$cell_id  <- cell_ids
  coords_cells          <- terra::xyFromCell(bios[[1]], cell_ids)
  my_sf_object$cell_lon <- coords_cells[, 1]
  my_sf_object$cell_lat <- coords_cells[, 2]
  my_sf_object          <- my_sf_object[!is.na(my_sf_object$cell_id), ]
  
  # Agregación por celda (sin filtro por n)
  data_aggregated_sf <- my_sf_object %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      n_samples_in_cell = dplyr::n(),
      Ho  = mean(Ho, na.rm = TRUE),
      lon = mean(cell_lon, na.rm = TRUE),
      lat = mean(cell_lat, na.rm = TRUE),
      .groups = "drop"
    )
  
  message("Modelling with ", nrow(data_aggregated_sf), " aggregated grid cells.")
  
  bios_focal <- terra::focal(bios, w = 15, fun = mean, na.policy = "only", na.rm = TRUE)
  gc()
  
  unique_points_vect <- terra::vect(data_aggregated_sf)
  ext_direct         <- terra::extract(bios_focal, unique_points_vect)
  
  if (any(!complete.cases(ext_direct))) {
    message("Coastal points with NA values detected. Extracting values from nearest valid pixel...")
    na_rows <- which(!complete.cases(ext_direct))
    for (i in na_rows) {
      ext_nearest <- terra::extract(bios_focal, unique_points_vect[i, ], nearest = TRUE)
      ext_direct[i, names(bios)] <- ext_nearest[, names(bios)]
    }
  }
  
  ext_direct_clean <- ext_direct %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE))) %>%
    dplyr::select(-ID)
  
  data_aggregated_sp <- cbind(data_aggregated_sf, ext_direct_clean)
  rm(ext_direct_clean); gc()
  
  # ------------------------------------------------------------------------------
  # 3. FEATURE SELECTION & DECORRELATION
  # ------------------------------------------------------------------------------
  data_sel <- data_aggregated_sp %>%
    dplyr::select(Ho, n_samples_in_cell, dplyr::all_of(names(bios))) %>%
    sf::st_drop_geometry()
  
  valid_rows <- complete.cases(data_sel)
  data_sel   <- data_sel[valid_rows, ]
  n_samples  <- nrow(data_sel)
  
  bio_vars_only <- setdiff(names(bios), c("lon", "lat"))
  
  cor_matrix <- cor(as.data.frame(data_sel[, bio_vars_only]))
  to_remove  <- caret::findCorrelation(cor_matrix, cutoff = cor_cutoff)
  
  selected_bios <- if (length(to_remove) > 0) bio_vars_only[-to_remove] else bio_vars_only
  
  if(isTRUE(addLonLat)){
    predictors_list <- c("lon", "lat", selected_bios)
  } else {
    predictors_list <- selected_bios
  }
  
  data_sel_model <- as.data.frame(data_sel[, c("Ho", "n_samples_in_cell", predictors_list)])
  
  # ------------------------------------------------------------------------------
  # 4. SVR MODEL TRAINING
  # ------------------------------------------------------------------------------
  weights_vec    <- log(data_sel_model$n_samples_in_cell + 1)
  data_sel_model <- data_sel_model %>% dplyr::select(-n_samples_in_cell)
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  if(isTRUE(use_loocv)){
    message("Validation scheme: LOOCV")
    train_ctrl <- caret::trainControl(
      method = "LOOCV",
      savePredictions = "final",
      allowParallel = TRUE
    )  
  } else {
    message("Validation scheme: Repeated 5-Fold CV (10 repeats)")
    train_ctrl <- caret::trainControl(
      method = "repeatedcv",
      number = 5,
      repeats = 10,
      savePredictions = "final",
      allowParallel = TRUE
    )
  }
  
  set.seed(123)
  
  if (kernel_type == "svmRadial") {
    svm_grid <- expand.grid(
      # Valores bajos/medios de sigma aseguran gradientes ambientales suaves
      sigma = c(0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1),
      # Regularización moderada para evitar que C > 50 memorice los datos de entrenamiento
      C     = c(0.05, 0.1, 0.5, 1, 2, 5, 10, 20)
    )
    final_svr <- caret::train(
      Ho ~ .,
      data       = data_sel_model,
      method     = "svmRadial",
      weights    = weights_vec,
      preProcess = c("center", "scale"),
      trControl  = train_ctrl,
      tuneGrid   = svm_grid,
      metric     = "RMSE"
    )
  } else {
    svm_grid <- expand.grid(
      # Búsqueda logarítmica de 30 puntos en el rango [10^-3, 10^2]
      C = 10^seq(-3, 2, length.out = 30)
    )
    final_svr <- caret::train(
      Ho ~ .,
      data       = data_sel_model,
      method     = "svmLinear",
      weights    = weights_vec,
      preProcess = c("center", "scale"),
      trControl  = train_ctrl,
      tuneGrid   = svm_grid,
      metric     = "RMSE"
    )
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  # ------------------------------------------------------------------------------
  # 5. MODEL METRICS & PREDICTION CONSOLIDATION
  # ------------------------------------------------------------------------------
  # Entrenamiento
  train_pred_vals <- predict(final_svr, newdata = data_sel_model)
  train_df <- data.frame(
    Observed  = data_sel_model$Ho,
    Predicted = train_pred_vals
  )
  
  train_rmse <- yardstick::rmse_vec(train_df$Observed, train_df$Predicted)
  train_mae  <- yardstick::mae_vec(train_df$Observed, train_df$Predicted)
  train_rsq  <- yardstick::rsq_vec(train_df$Observed, train_df$Predicted)
  
  # Filtrar por el hiperparámetro óptimo
  if (kernel_type == "svmRadial") {
    cv_preds_raw <- final_svr$pred %>%
      dplyr::filter(sigma == final_svr$bestTune$sigma, C == final_svr$bestTune$C)
  } else {
    cv_preds_raw <- final_svr$pred %>%
      dplyr::filter(C == final_svr$bestTune$C)
  }
  
  # PROMEDIADO OUT-OF-FOLD (Elimina columnas verticales en las 10 repeticiones)
  cv_df <- cv_preds_raw %>%
    dplyr::group_by(rowIndex) %>%
    dplyr::summarise(
      Observed  = mean(obs, na.rm = TRUE),
      Predicted = mean(pred, na.rm = TRUE),
      .groups   = "drop"
    )
  
  # Métricas de validación reales sobre datos promediados
  test_rmse <- yardstick::rmse_vec(cv_df$Observed, cv_df$Predicted)
  test_mae  <- yardstick::mae_vec(cv_df$Observed, cv_df$Predicted)
  test_rsq  <- yardstick::rsq_vec(cv_df$Observed, cv_df$Predicted)
  
  metrics_out <- data.frame(
    .metric = c("rmse", "mae", "rsq", "rmse", "mae", "rsq"),
    mean    = c(test_rmse, test_mae, test_rsq, train_rmse, train_mae, train_rsq),
    dataset = c("validation_test", "validation_test", "validation_test", "train", "train", "train")
  )
  
  suffix <- paste0("_", kernel_type, if(addLonLat) "_LONLAT" else "")
  write.csv(metrics_out, file.path(outdir, paste0(sp_name, "_SVM_Metrics", suffix, ".csv")), row.names = FALSE)
  
  # ------------------------------------------------------------------------------
  # 6. EVALUATIVE PLOTS
  # ------------------------------------------------------------------------------
  p_train <- ggplot(train_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "darkgreen") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "darkgreen") +
    labs(title = "A) Training",
         subtitle = paste("n =", n_samples, "cells"),
         x = "Observed Ho", y = "Predicted Ho") +
    theme_bw(13)
  
  p_test <- ggplot(cv_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "blue") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "blue") +
    labs(title = "B) Validation (Test)",
         subtitle = ifelse(use_loocv, "Leave-One-Out CV", "Repeated 5-Fold CV (Averaged)"),
         x = "Observed Ho", y = "Predicted Ho") +
    theme_bw(13)
  
  ggsave(
    file.path(outdir, paste0(sp_name, "_SVM_Performance", suffix, ".png")),
    p_train + p_test, width = 12, height = 5.5, dpi = 600
  )
  
  # ------------------------------------------------------------------------------
  # 7. SPATIAL PREDICTION (RASTER)
  # ------------------------------------------------------------------------------
  sdm <- terra::rast(sdm_path) * 1
  sdm[sdm == 0] <- NA
  
  bios_selected <- bios[[predictors_list]]
  bios_selected <- terra::resample(bios_selected, sdm, method = "bilinear")
  bios_selected <- bios_selected * sdm
  
  valid_cells <- terra::cells(bios_selected[[1]])
  
  temp_dt <- as.data.frame(terra::extract(bios_selected, valid_cells))
  if ("ID" %in% colnames(temp_dt)) temp_dt$ID <- NULL
  temp_dt$cell_idx <- valid_cells
  
  valid_rows <- complete.cases(temp_dt[, predictors_list, drop = FALSE])
  temp_dt    <- temp_dt[valid_rows, ]
  temp_dt$prediction <- NA_real_
  
  chunk_size <- 50000
  n_rows     <- nrow(temp_dt)
  n_chunks   <- ceiling(n_rows / chunk_size)
  
  for (i in seq_len(n_chunks)) {
    s_idx <- ((i - 1) * chunk_size) + 1
    e_idx <- min(i * chunk_size, n_rows)
    
    newdata_chunk <- temp_dt[s_idx:e_idx, predictors_list, drop = FALSE]
    pred_res <- predict(final_svr, newdata = newdata_chunk)
    temp_dt$prediction[s_idx:e_idx] <- pred_res
  }
  
  pred_raster <- bios_selected[[1]] * NA
  names(pred_raster) <- "Predicted_Ho"
  pred_raster[temp_dt$cell_idx] <- temp_dt$prediction
  
  # ------------------------------------------------------------------------------
  # 8. MESS MASKING
  # ------------------------------------------------------------------------------
  bios_selected_stack <- raster::stack(bios_selected)
  ref_data <- as.data.frame(data_sel_model[, predictors_list, drop = FALSE])
  
  mess_rn <- terra::rast(dismo::mess(x = bios_selected_stack, v = ref_data))
  mess_rn <- terra::resample(mess_rn, sdm, method = "bilinear")
  
  mess_mask <- mess_rn
  mess_mask[mess_mask <= 0 | is.infinite(mess_mask)] <- NA
  mess_mask[!is.na(mess_mask)] <- 1
  
  pred_raster_no_interpolated <- pred_raster * mess_mask
  
  terra::writeRaster(
    pred_raster_no_interpolated, 
    file.path(outdir, paste0(sp_name, "_Ho_Map_SVM", suffix, ".tif")), 
    overwrite = TRUE
  )
  
  # ------------------------------------------------------------------------------
  # 9. MAP EXPORT & WORKSPACE SAVING
  # ------------------------------------------------------------------------------
  data("World", package = "tmap")
  
  map_out <- tm_shape(pred_raster_no_interpolated) + 
    tm_raster(
      col.legend = tm_legend(title = "Ho"),
      palette = RColorBrewer::brewer.pal(7, "YlGnBu")
    ) +
    tm_shape(World) + 
    tm_borders(col = "grey40") +
    tm_shape(my_sf_object) +
    tm_symbols(col = "grey20", size = 0.15, alpha = 0.6) +
    tm_shape(data_aggregated_sf) +
    tm_symbols(col = "red", size = 0.35, shape = 21, border.col = "black", border.lwd = 0.5) +
    tm_graticules(labels.size = 0.7) +
    tm_layout(inner.margins = 0, legend.outside = TRUE, legend.outside.position = "right")
  
  tmap::tmap_save(
    tm       = map_out, 
    filename = file.path(outdir, paste0(sp_name, "_Ho_Map_SVM", suffix, ".pdf")), 
    width    = 25, 
    height   = 15, 
    units    = "cm", 
    dpi      = 600
  )
  
  save.image(file.path(outdir, paste0(sp_name, "_Ho_Workspace", suffix, ".RData")))
  
  message("Pipeline completed successfully!")
  invisible(NULL)
}



# ------------------------------------------------------------------------------
# Crocodylus acutus
# ------------------------------------------------------------------------------
# With geographic coordinates (LON/LAT) + LOOCV
CA <- svr_ho_pipeline(
  outdir      = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics",
  sp_name     = "Crocodylus acutus",
  raster_dir  = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
  data_path   = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
  sdm_path    = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif",
  n_cores     = 12,
  cor_cutoff  = 0.5,
  addLonLat   = TRUE,
  use_loocv   = TRUE,
  kernel_type = "svmLinear"
)

# Bioclimatic variables only + LOOCV
CA2 <- svr_ho_pipeline(
  outdir      = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics",
  sp_name     = "Crocodylus acutus",
  raster_dir  = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
  data_path   = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
  sdm_path    = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif",
  n_cores     = 8,
  cor_cutoff  = 0.5,
  addLonLat   = FALSE,
  use_loocv   = TRUE,
  kernel_type = "svmLinear"
)


# ------------------------------------------------------------------------------
# Crocodylus intermedius (Commented)
# ------------------------------------------------------------------------------
# CI <- svr_ho_pipeline(
#   outdir      = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics",
#   sp_name     = "Crocodylus intermedius",
#   raster_dir  = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
#   data_path   = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
#   sdm_path    = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_intermedius/Crocodylus_intermedius_Binario_P10.tif",
#   n_cores     = 8,
#   cor_cutoff  = 0.6,
#   addLonLat   = FALSE,
#   use_loocv   = TRUE,
#   kernel_type = "svmRadial"
# )


# ------------------------------------------------------------------------------
# Crocodylus moreletii
# ------------------------------------------------------------------------------
# With geographic coordinates (LON/LAT) + Repeated 5-Fold CV
CM <- svr_ho_pipeline(
  outdir      = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics",
  sp_name     = "Crocodylus moreletii",
  raster_dir  = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
  data_path   = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
  sdm_path    = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif",
  n_cores     = 8,
  cor_cutoff  = 0.5,
  addLonLat   = TRUE,
  use_loocv   = FALSE,
  kernel_type = "svmLinear" # Recomended over svmLinear for non-linear response
)

# Bioclimatic variables only + Repeated 5-Fold CV
CM2 <- svr_ho_pipeline(
  outdir      = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics",
  sp_name     = "Crocodylus moreletii",
  raster_dir  = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
  data_path   = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
  sdm_path    = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif",
  n_cores     = 8,
  cor_cutoff  = 0.5,
  addLonLat   = FALSE,
  use_loocv   = FALSE,
  kernel_type = "svmLinear"
)
