#' Spatial Prediction Workflow for Genetic Heterozygosity (Ho) using Support Vector Regression
#'
#' @param outdir \code{character}. Output directory path.
#' @param sp_name \code{character}. Exact target species name.
#' @param raster_dir \code{character}. Directory path containing individual predictor raster layers.
#' @param data_path \code{character}. File path to the Excel file containing genetic data.
#' @param sdm_path \code{character}. File path to the binary SDM raster.
#' @param n_cores \code{integer}. Number of processing cores. Default is \code{6L}.
#' @param cor_cutoff \code{numeric}. Absolute pairwise Pearson correlation threshold. Default is \code{0.5}.
#' @param addLonLat \code{logical}. Retain explicit lon/lat coordinates. Default is \code{FALSE}.
#' @param use_loocv \code{logical}. If \code{TRUE}, executes LOOCV (ignored if \code{use_blockcv = TRUE}). Default is \code{FALSE}.
#' @param use_blockcv \code{logical}. If \code{TRUE} and points >= 50, executes Spatial Block CV via \code{blockCV}. Default is \code{FALSE}.
#' @param kernel_type \code{character}. SVM kernel (\code{"svmRadial"} or \code{"svmLinear"}). Default is \code{"svmRadial"}.
#' @param aggregate_occs_cells \code{logical}. Aggregate occurrences by raster cell. Default is \code{TRUE}.
#'
#' @return Invisibly returns \code{NULL}.
#' @export

library(dplyr)
library(readxl)
library(sf)
library(terra)
library(caret)
library(doParallel)
library(kernlab)
library(yardstick)
library(ggplot2)
library(ggpmisc)
library(patchwork)
library(raster)
library(dismo)
library(spThin)
library(RColorBrewer)
library(tmap)
library(blockCV)

svr_ho_pipeline <- function(outdir, 
                            sp_name, 
                            raster_dir, 
                            data_path, 
                            sdm_path, 
                            n_cores = 6L,
                            cor_cutoff = 0.5,
                            addLonLat = FALSE,
                            use_loocv = FALSE,
                            use_blockcv = FALSE,
                            kernel_type = "svmRadial",
                            aggregate_occs_cells = TRUE) {
  
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
  
  data <- readxl::read_xlsx(data_path, sheet = "data") %>%
    dplyr::filter(sp == sp_name, !is.na(Ho), Ho > 0, !is.na(lat), !is.na(lon))
  
  my_sf_object <- sf::st_as_sf(data, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  points_vect  <- terra::vect(my_sf_object)
  
  if(isTRUE(aggregate_occs_cells)){
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
    
    data_aggregated_sf <- my_sf_object %>%
      dplyr::group_by(cell_id) %>%
      dplyr::summarise(
        n_samples_in_cell = dplyr::n(),
        Ho  = median(Ho, na.rm = TRUE),
        Ho_mad  = mad(Ho, na.rm = TRUE),
        lon = unique(cell_lon, na.rm = TRUE),
        lat = unique(cell_lat, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data_aggregated_sf <- my_sf_object
  }
  
  message("Modelling with ", nrow(data_aggregated_sf), " aggregated grid cells/points.")
  
  bios_focal <- terra::focal(bios, w = 15, fun = mean, na.policy = "only", na.rm = TRUE)
  gc()
  
  unique_points_vect <- terra::vect(data_aggregated_sf)
  ext_direct         <- terra::extract(bios, unique_points_vect,"mean")
  
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
  
  if(nrow(data_aggregated_sp)>50){
    message("more than 50 occurrences are availables, performing spatial thin at 5 km")
    data_aggregated_sp <- terra::thin(vect(data_aggregated_sp),5000)
    data_aggregated_sp<- st_as_sf(data_aggregated_sp)
  }  # ------------------------------------------------------------------------------
  # 3. FEATURE SELECTION & DECORRELATION
  # ------------------------------------------------------------------------------
  if(isTRUE(aggregate_occs_cells)){
    data_sel <- data_aggregated_sp %>%
      dplyr::select(Ho, n_samples_in_cell, dplyr::all_of(names(bios))) %>%
      sf::st_drop_geometry()
  } else {
    data_sel <- data_aggregated_sp %>%
      dplyr::select(Ho, dplyr::all_of(names(bios))) %>%
      sf::st_drop_geometry()
  }
  
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
  
  if(isTRUE(aggregate_occs_cells)){
    data_sel_model <- as.data.frame(data_sel[, c("Ho", "n_samples_in_cell", predictors_list)])
  } else {
    data_sel_model <- as.data.frame(data_sel[, c("Ho", predictors_list)])
  }
  
  
  #data used for model in CSV
  if(isTRUE(aggregate_occs_cells)){
    file_to_save <- file.path(outdir, paste0(sp_name, "_data_used_for_model_CELLS", ".csv"))
  } else {
    file_to_save <- file.path(outdir, paste0(sp_name, "_data_used_for_model_POINTS",  ".csv"))
  }
  write.csv(data_sel_model, file_to_save, row.names = FALSE)
  
  #predictor list
  if(isTRUE(aggregate_occs_cells)){
    file_to_save <- file.path(outdir, paste0(sp_name, "predictor_CELLS", ".csv"))
  } else {
    file_to_save <- file.path(outdir, paste0(sp_name, "predictor_points",  ".csv"))
  }
  write.csv(predictors_list, file_to_save, row.names = FALSE)
  
  # ------------------------------------------------------------------------------
  # 4. SVR MODEL TRAINING & VALIDATION SCHEME EVALUATION
  # ------------------------------------------------------------------------------
  if(isTRUE(aggregate_occs_cells)){
    weights_vec    <- sqrt(data_sel_model$n_samples_in_cell)
    data_sel_model <- data_sel_model %>% dplyr::select(-n_samples_in_cell)
  } else {
    weights_vec    <- NULL
  }
  
  applied_blockcv <- FALSE
  
  if (isTRUE(use_blockcv)) {
    if (n_samples >= 50) {
      message(sprintf("Validation scheme: Spatial Block Cross-Validation (blockCV) [n = %d >= 50]", n_samples))
      
      data_sf_model <- data_aggregated_sp[valid_rows, ]
      
      set.seed(123)
      sb <- blockCV::cv_spatial(
        x = data_sf_model,
        k = 4,
        selection = "random",
        progress = FALSE
      )
      
      train_index <- lapply(sb$folds_list, function(fold) fold[[1]])
      test_index <- lapply(sb$folds_list, function(fold) fold[[2]])
      
      names(train_index) <- paste0("Fold", seq_along(train_index))
      names(test_index) <- names(train_index)
      
      
      train_ctrl <- caret::trainControl(
        method = "cv",
        number = length(train_index),
        index = train_index,
        search = "random",
        indexOut = test_index,
        savePredictions = "final",
        allowParallel = TRUE
      )
      applied_blockcv <- TRUE
      
    } else {
      message(sprintf("[ESCENARIO DESCARTADO] 'use_blockcv = TRUE' solicitado, pero n = %d (< 7).", n_samples))
      message("--> Finalizando ejecucion de este escenario sin generar resultados.\n")
      return(NULL)
    }
  }
  
  if (!applied_blockcv) {
    if (isTRUE(use_loocv)) {
      message(sprintf("Validation scheme: Leave-One-Out CV (LOOCV) [n = %d]", n_samples))
      train_ctrl <- caret::trainControl(
        method = "LOOCV",
        savePredictions = "final",
        allowParallel = TRUE
      )  
    } else {
      message(sprintf("Validation scheme: Repeated 5-Fold CV (10 repeats) [n = %d]", n_samples))
      train_ctrl <- caret::trainControl(
        method = "repeatedcv",
        number = 5,
        repeats = 10,
        savePredictions = "final",
        allowParallel = TRUE
      )
    }
  }
  
  # Etiqueta estandarizada para identificar la validación en los archivos
  val_scheme_str <- if (applied_blockcv) {
    "BLOCKCV"
  } else if (isTRUE(use_loocv)) {
    "LOOCV"
  } else {
    "REPEATEDCV"
  }
  
  # ------------------------------------------------------------------------------
  # ENTRENAMIENTO SVM
  # ------------------------------------------------------------------------------
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  set.seed(123)
  
  if (kernel_type == "svmRadial") {
    svm_grid <- expand.grid(
      sigma = c(0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1,0.5,1),
      C     = c(0.05, 0.1, 0.5, 1, 2, 5, 10, 20,50,100)
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
  # A) PREDICCIONES DE ENTRENAMIENTO
  train_pred_vals <- predict(final_svr, newdata = data_sel_model)
  train_df <- data.frame(
    rowIndex  = 1:nrow(data_sel_model),
    Observed  = data_sel_model$Ho,
    Predicted = train_pred_vals
  )
  
  train_rmse <- yardstick::rmse_vec(train_df$Observed, train_df$Predicted)
  train_mae  <- yardstick::mae_vec(train_df$Observed, train_df$Predicted)
  train_rsq  <- yardstick::rsq_vec(train_df$Observed, train_df$Predicted)
  
  # B) PREDICCIONES OUT-OF-FOLD (VALIDACIÓN CV)
  if (kernel_type == "svmRadial") {
    best_sig <- final_svr$bestTune$sigma
    best_c   <- final_svr$bestTune$C
    cv_preds_raw <- dplyr::inner_join(
      final_svr$pred, 
      final_svr$bestTune, 
      by = names(final_svr$bestTune)
    )  } else {
    best_c   <- final_svr$bestTune$C
    cv_preds_raw <- dplyr::inner_join(
      final_svr$pred, 
      final_svr$bestTune, 
      by = names(final_svr$bestTune)
    )
  }
  
  # Si el esquema incluye repeticiones (ej. repeatedcv), promediamos por observacion (rowIndex)
  cv_df <- cv_preds_raw %>%
    dplyr::group_by(rowIndex) %>%
    dplyr::summarise(
      Observed  = mean(obs, na.rm = TRUE),             # Es idéntico en todas las repeticiones
      Predicted = mean(pred, na.rm = TRUE), # Promedio de predicción fuera de pliegue
      .groups   = "drop"
    ) %>%
    dplyr::arrange(rowIndex)
  
  test_rmse <- yardstick::rmse_vec(cv_df$Observed, cv_df$Predicted)
  test_mae  <- yardstick::mae_vec(cv_df$Observed, cv_df$Predicted)
  test_rsq  <- yardstick::rsq_vec(cv_df$Observed, cv_df$Predicted)
  
  # Inclusión explícita del esquema de validación y parámetros en las métricas
  metrics_out <- data.frame(
    .metric           = c("rmse", "mae", "rsq", "rmse", "mae", "rsq"),
    mean              = c(test_rmse, test_mae, test_rsq, train_rmse, train_mae, train_rsq),
    dataset           = c("validation_test", "validation_test", "validation_test", "train", "train", "train"),
    validation_scheme = val_scheme_str,
    kernel_type       = kernel_type,
    add_lon_lat       = addLonLat
  )
  
  # Construcción del sufijo único para nombres de archivo
  suffix <- paste0("_", kernel_type, "_", val_scheme_str, if(addLonLat) "_LONLAT" else "")
  
  if(isTRUE(aggregate_occs_cells)){
    file_to_save <- file.path(outdir, paste0(sp_name, "_SVM_Metrics_CELLS", suffix, ".csv"))
  } else {
    file_to_save <- file.path(outdir, paste0(sp_name, "_SVM_Metrics", suffix, ".csv"))
  }
  write.csv(metrics_out, file_to_save, row.names = FALSE)
  
  # ------------------------------------------------------------------------------
  # 6. EVALUATIVE PLOTS
  # ------------------------------------------------------------------------------
  if(isTRUE(aggregate_occs_cells)){
    file_to_save <- file.path(outdir, paste0(sp_name, "_SVM_Performance_CELLS", suffix, ".png"))
  } else {
    file_to_save <- file.path(outdir, paste0(sp_name, "_SVM_Performance", suffix, ".png"))  
  }
  
  p_train <- ggplot(train_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "darkgreen") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "darkgreen") +
    labs(title = "A) Training Fit",
         subtitle = paste("n =", n_samples, "cells (Resubstitution)"),
         x = "Observed Ho", y = "Predicted Ho") +
    theme_bw(13)
  
  test_subtitle <- if (applied_blockcv) {
    "Spatial Block Cross-Validation (blockCV)"
  } else if (use_loocv) {
    "Leave-One-Out CV (LOOCV)"
  } else {
    "Repeated 5-Fold CV (Averaged)"
  }
  
  p_test <- ggplot(cv_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "blue") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "blue") +
    labs(title = "B) Cross-Validation (Test)",
         subtitle = test_subtitle,
         x = "Observed Ho", y = "Predicted Ho") +
    theme_bw(13)
  
  ggsave(
    file_to_save,
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
  
  valid_rows_pred <- complete.cases(temp_dt[, predictors_list, drop = FALSE])
  temp_dt         <- temp_dt[valid_rows_pred, ]
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
  
  if(isTRUE(aggregate_occs_cells)){
    file_to_save <- file.path(outdir, paste0(sp_name, "_Ho_Map_SVM_CELLS", suffix, ".tif"))
  } else {
    file_to_save <- file.path(outdir, paste0(sp_name, "_Ho_Map_SVM", suffix, ".tif"))
  }
  
  terra::writeRaster(
    pred_raster_no_interpolated, 
    file_to_save, 
    overwrite = TRUE
  )
  
  # ------------------------------------------------------------------------------
  # 9. MAP EXPORT & WORKSPACE SAVING
  # ------------------------------------------------------------------------------
  data("World", package = "tmap")
  
  if(isTRUE(aggregate_occs_cells)){
    file_to_save <- file.path(outdir, paste0(sp_name, "_Ho_Map_SVM_CELLS", suffix, ".pdf"))
  } else {
    file_to_save <- file.path(outdir, paste0(sp_name, "_Ho_Map_SVM", suffix, ".pdf"))
  }
  
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
    filename = file_to_save, 
    width    = 25, 
    height   = 15, 
    units    = "cm", 
    dpi      = 600
  )
  
  if(isTRUE(aggregate_occs_cells)){
    file_to_save <- file.path(outdir, paste0(sp_name, "_Ho_Workspace_CELLS", suffix, ".RData"))
  } else {
    file_to_save <- file.path(outdir, paste0(sp_name, "_Ho_Workspace", suffix, ".RData"))
  }
  
  save.image(file_to_save)
  
  message(sprintf("Pipeline completed successfully for '%s' using [%s]!", sp_name, val_scheme_str))
  invisible(NULL)
}

# ==============================================================================
# CONFIGURACIÓN GENERAL DE RUTAS Y PARÁMETROS
# ==============================================================================
out_dir    <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics" #donde se guardan
r_dir      <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s" #ruta de las capas
d_path     <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx" #tabla de entrada (Ho)
sdm_base   <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval" #irectorio donde estan los modelos binarizados
# =============================================================================
# 1. DEFINICIÓN DE RUTAS BASE
# =============================================================================
sdm_ca <- file.path(sdm_base, "Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif") #C. acutus (SDM 1)
sdm_cm <- file.path(sdm_base, "Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif") #C. moreletti (SDM 1)

# Mapas de información por especie y tipo de validación cruzada
species_info <- list(
  CA = list(name = "Crocodylus acutus", sdm = sdm_ca),
  CM = list(name = "Crocodylus moreletii", sdm = sdm_cm)
)

cv_configs <- list(
  LOOCV   = list(use_loocv = TRUE,  use_blockcv = FALSE),
  BLOCKCV = list(use_loocv = FALSE, use_blockcv = TRUE),
  KFOLD   = list(use_loocv = FALSE, use_blockcv = FALSE) # Ajusta banderas si KFOLD requiere otra op.
)

# =============================================================================
# 2. GENERACIÓN DE LA GRILLA DE ESCENARIOS
# =============================================================================
scenarios_grid <- expand.grid(
  sp_code              = c("CA", 
                           "CM"),
  aggregate_occs_cells = c(TRUE, FALSE),
  addLonLat            = c(TRUE, FALSE),
  cv_type              = c(
                          "LOOCV",
                           "BLOCKCV",
                           "KFOLD"
                           ),
  stringsAsFactors     = FALSE
)

cat("Total de escenarios a ejecutar:", nrow(scenarios_grid), "\n")

# =============================================================================
# 3. EJECUCIÓN AUTOMATIZADA
# =============================================================================
results <- list()

for (i in 1:nrow(scenarios_grid)) {
  
  row <- scenarios_grid[i, ]
  
  # Extraer configuración de especie y validación cruzada
  sp_data <- species_info[[row$sp_code]]
  cv_data <- cv_configs[[row$cv_type]]
  
  # Construir sufijos dinámicos para el directorio
  cell_tag <- if (isTRUE(row$aggregate_occs_cells)) "CELLS" else "POINTS"
  ll_tag   <- if (isTRUE(row$addLonLat)) "WITH_LONLAT" else "NO_LONLAT"
  cv_tag   <- row$cv_type
  
  
  # Crear nombre único de carpeta de salida
  folder_name <- paste(row$sp_code, cell_tag, ll_tag, cv_tag, sep = "_")
  outdir_path <- file.path(out_dir, folder_name)
  
  cat(sprintf("\n------------------------------------------------------------\n"))
  cat(sprintf("[%d/%d] Ejecutando: %s\n", i, nrow(scenarios_grid), folder_name))
  cat(sprintf("------------------------------------------------------------\n"))
  
  # Llamada a la función del pipeline
  results[[folder_name]] <- svr_ho_pipeline(
    outdir               = outdir_path,
    sp_name              = sp_data$name,
    raster_dir           = r_dir,
    data_path            = d_path,
    sdm_path             = sp_data$sdm,
    n_cores              = 8,
    cor_cutoff           = 0.5,
    addLonLat            = row$addLonLat,
    use_blockcv          = cv_data$use_blockcv,
    use_loocv            = cv_data$use_loocv,
    kernel_type          = "svmRadial",
    aggregate_occs_cells = row$aggregate_occs_cells
  )
  tmpFiles(current=TRUE, orphan=TRUE, old=TRUE, remove=TRUE)
  
}

cat("\n¡Procesamiento finalizado exitosamente para todos los escenarios!\n")