################################################################################
# MACROGENETICS PREDICTION WORKFLOW (RANDOM FOREST VIA RANGER)
# Continuous geographic prediction of heterozygosity (Ho) over a raster.
################################################################################

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ranger)
  library(patchwork)
  library(dplyr)
  library(ggplot2)
  library(ggpmisc)
  library(tmap)
  library(Boruta)
  library(ape)
  library(RColorBrewer)
  library(readxl)
  library(caret)
  library(dismo)
  library(raster)
})

#' Nested LOOCV Random Forest Workflow for Spatial Prediction of Heterozygosity (Ho)
#'
#' Fits a Random Forest (\code{ranger}) model to predict continuous genetic 
#' observed heterozygosity (\eqn{H_o}) from environmental predictors (e.g., bioclimatic rasters) 
#' and optional geographic coordinates. The pipeline implements fully nested Leave-One-Out 
#' Cross-Validation (LOOCV), spatial aggregation/thinning, raster focal smoothing for coastal 
#' point extraction, a final production model, spatial prediction masked by a binary Species 
#' Distribution Model (SDM) footprint, Multivariate Environmental Similarity Surface (MESS) 
#' masking to prevent extrapolation, and spatial autocorrelation assessment of LOOCV residuals 
#' via Moran's I.
#'
#' @param outdir Character. Output directory path where all resulting figures (PNG), 
#'   tables (CSV), GeoTIFF rasters, and the saved \code{.RData} workspace image will be stored.
#' @param sp_name Character. Target species name matching the \code{"sp"} column in the input Excel 
#'   file. Used as a prefix for saved output files.
#' @param raster_dir Character. Directory path containing environmental/bioclimatic raster layers 
#'   in GeoTIFF (\code{.tif}) format used as model predictors.
#' @param data_path Character. File path to the Excel workbook (\code{.xlsx}) containing genetic 
#'   and spatial point data in a sheet named \code{"data"}. Must include \code{"sp"}, \code{"Ho"}, 
#'   \code{"lon"}, and \code{"lat"} columns.
#' @param sdm_path Character. File path to a binary SDM raster (\code{.tif}) used to restrict the 
#'   spatial prediction footprint (1 = suitable/predicted region, 0/NA = masked out).
#' @param N_CORES Integer. Number of CPU threads passed to \code{\link[ranger]{ranger}} for model 
#'   fitting and parallel spatial raster prediction. Default is \code{6L}.
#' @param boruta_threshold Integer. Consensus threshold defining the minimum number of independent 
#'   Boruta runs (out of \code{n_seeds}) in which a variable must be classified as "Confirmed" 
#'   to be retained for modeling. Default is \code{30}.
#' @param n_seeds Integer. Total number of random seed iterations used to evaluate variable 
#'   selection stability via \code{\link[Boruta]{Boruta}}. Default is \code{30}.
#' @param cor_cutoff Numeric. Absolute Pearson correlation coefficient threshold (between 0 and 1) 
#'   passed to \code{\link[caret]{findCorrelation}} to remove redundant bioclimatic predictors 
#'   prior to feature selection. Default is \code{0.5}.
#' @param addLonLat Logical. If \code{TRUE}, explicit longitude (\code{"lon"}) and latitude 
#'   (\code{"lat"}) variables are appended to the predictor pool alongside environmental layers. 
#'   Default is \code{FALSE}.
#' @param aggregate_occs_cells Logical. If \code{TRUE}, multiple occurrence samples falling within 
#'   the same raster grid cell are spatially aggregated into a single observation (using median 
#'   \eqn{H_o} and cell centroid coordinates) weighted by sample size. If \code{FALSE}, raw individual 
#'   points are modeled directly. Default is \code{TRUE}.
#'
#' @details
#' The function executes a full macrogenetic spatial predictive pipeline comprising the following steps:
#' \enumerate{
#'   \item \strong{Data Ingestion & Filtering}: Bioclimatic rasters and genetic Excel data are loaded. 
#'     Records are filtered for the target species (\code{sp_name}) and non-missing positive \eqn{H_o} values.
#'   \item \strong{Grid Cell Aggregation & Coastal Extraction}: When \code{aggregate_occs_cells = TRUE}, 
#'     points inside the same raster cell are collapsed to the cell centroid, calculating median \eqn{H_o} 
#'     and sample counts (\code{n_samples_in_cell}). Points offset from raster coverage are snapped to 
#'     the nearest valid cell. Coastal points returning \code{NA} environmental values are extracted 
#'     using a 15-cell focal-averaged raster stack (\code{terra::focal}).
#'   \item \strong{Spatial Thinning}: If the dataset exceeds 99 spatial observations, spatial thinning 
#'     at a 5 km minimum distance is applied via \code{terra::thin}.
#'   \item \strong{Decorrelation & Feature Selection}: Highly collinear predictors above \code{cor_cutoff} 
#'     are removed using pairwise correlation matrix analysis (\code{caret::findCorrelation}).
#'   \item \strong{Nested LOOCV Validation}: Evaluates model generalization capacity without spatial leakage. 
#'     For each fold \eqn{i \in \{1, \dots, N\}}:
#'     \itemize{
#'       \item Boruta feature selection is executed across \code{n_seeds} on the \eqn{N-1} training set.
#'       \item Hyperparameter tuning (\code{num.trees}, \code{mtry}, \code{min.node.size}) is performed via 
#'         a grid search optimizing Out-Of-Bag (OOB) \eqn{R^2}.
#'       \item A Random Forest model is fitted on training fold \eqn{i} and evaluates the single held-out fold \eqn{i}.
#'     }
#'   \item \strong{Production Model & Spatial Prediction}: A final production Random Forest model is trained 
#'     on the full dataset. Spatial prediction is generated chunk-wise across the raster stack within the 
#'     binary SDM footprint.
#'   \item \strong{MESS Extrapolation Masking & Autocorrelation Check}: Multivariate Environmental Similarity 
#'     Surface (\code{dismo::mess}) masking removes spatial regions requiring environmental extrapolation 
#'     (\eqn{MESS \le 0}). Finally, Moran's I (\code{ape::Moran.I}) measures spatial autocorrelation in the 
#'     LOOCV residuals using inverse-distance spatial weights with units dropped to avoid S3 class conflicts.
#' }
#'
#' @return Invisibly returns \code{NULL}. The primary function outputs are written directly to \code{outdir}:
#' \itemize{
#'   \item \code{*_LOOCV_Performance_RF_VALUES_nested.csv}: Observed vs. predicted \eqn{H_o} values under nested LOOCV.
#'   \item \code{*_LOOCV_nested_vars_per_fold.csv}: Predictors selected by Boruta per cross-validation fold.
#'   \item \code{*_LOOCV_nested_hyperparams_per_fold.csv}: Optimal hyperparameters selected per fold.
#'   \item \code{*_Variable_Importance_RF.png}: Permutation variable importance plot for the final model.
#'   \item \code{*_Model_Performance_Train_vs_LOOCV.png}: Comparative regression scatter plots (Train fit vs. LOOCV test).
#'   \item \code{*_Ho_Macrogenetics_Map_RF_final.tif}: Final continuous predicted \eqn{H_o} surface masked by SDM and MESS.
#'   \item \code{*_Moran_I_LOOCV_residuals.csv}: Spatial autocorrelation test statistics for LOOCV residuals.
#'   \item \code{*_Ho_Macrogenetics_Map_RF_final.RData}: Complete saved R environment workspace image.
#' }
#'
#' @author Chrystian Camilo Sosa Arango
#'
#' @export
RF_LOOCV <- function(outdir, 
                     sp_name, 
                     raster_dir,
                     data_path, 
                     sdm_path,
                     N_CORES = 6L,
                     boruta_threshold = 30, 
                     n_seeds = 30,
                     cor_cutoff = 0.5,
                     addLonLat = FALSE,
                     aggregate_occs_cells = TRUE
) {
  
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
  ext_direct         <- terra::extract(bios, unique_points_vect, "mean")
  
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
  
  if(nrow(data_aggregated_sp) > 99){
    message("more than 100 occurrences are availables, performing spatial thin at 5 km")
    data_aggregated_sp <- terra::thin(vect(data_aggregated_sp), 5000)
    data_aggregated_sp <- st_as_sf(data_aggregated_sp)
  }
  
  # ------------------------------------------------------------------------------
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
  
  # ------------------------------------------------------------------------------
  # HELPER FUNCTIONS FOR NESTED LOOCV
  # ------------------------------------------------------------------------------
  run_boruta_selection <- function(train_data, all_vars, threshold, n_seeds) {
    formula_full <- as.formula(paste("Ho ~", paste(all_vars, collapse = " + ")))
    resultados <- vector("list", n_seeds)
    for (s in seq_len(n_seeds)) {
      set.seed(s)
      boruta_run <- Boruta::Boruta(formula_full, data = train_data, doTrace = 0)
      boruta_run <- Boruta::TentativeRoughFix(boruta_run)
      resultados[[s]] <- data.frame(
        seed = s,
        variable = names(boruta_run$finalDecision),
        decision = as.character(boruta_run$finalDecision)
      )
    }
    tabla <- do.call(rbind, resultados)
    resumen <- tabla %>%
      dplyr::filter(decision == "Confirmed") %>%
      dplyr::count(variable, name = "n_confirmed") %>%
      dplyr::arrange(desc(n_confirmed))
    
    vars_sel <- resumen$variable[resumen$n_confirmed >= threshold]
    if (length(vars_sel) == 0) vars_sel <- all_vars
    vars_sel
  }
  
  run_grid_search <- function(train_data, vars_sel) {
    formula_opt <- as.formula(paste("Ho ~", paste(vars_sel, collapse = " + ")))
    n_tr <- nrow(train_data)
    
    # Dynamic adjustment of min.node.size and mtry
    min_nodes <- unique(pmax(1, floor(c(n_tr * 0.10, n_tr * 0.20, n_tr * 0.30))))
    mtry_max  <- length(vars_sel)
    mtry_vals <- unique(pmax(1, pmin(mtry_max, c(2, floor(sqrt(mtry_max)), floor(mtry_max / 2)))))
    
    tuning_grid <- expand.grid(
      num.trees = c(500, 1000, 1500),
      mtry = mtry_vals,
      min.node.size = min_nodes
    )
    
    res <- vector("list", nrow(tuning_grid))
    for (g in seq_len(nrow(tuning_grid))) {
      fit <- ranger::ranger(
        formula = formula_opt,
        data = train_data,
        case.weights = if("n_samples_in_cell" %in% names(train_data)) sqrt(train_data$n_samples_in_cell) else NULL,
        num.trees = tuning_grid$num.trees[g],
        mtry = tuning_grid$mtry[g],
        min.node.size = tuning_grid$min.node.size[g],
        importance = "permutation",
        replace = FALSE,
        seed = 100,
        num.threads = N_CORES
      )
      res[[g]] <- data.frame(
        num.trees = tuning_grid$num.trees[g],
        mtry = tuning_grid$mtry[g],
        min.node.size = tuning_grid$min.node.size[g],
        OOB_R2 = fit$r.squared,
        MSE = fit$prediction.error
      )
    }
    tuning_df <- do.call(rbind, res)
    tuning_df[order(tuning_df$OOB_R2, decreasing = TRUE), ][1, ]
  }
  
  # ------------------------------------------------------------------------------
  # 4. NESTED LOOCV
  # ------------------------------------------------------------------------------
  loocv_preds  <- numeric(n_samples)
  loocv_params <- vector("list", n_samples)
  loocv_vars   <- vector("list", n_samples)
  
  message("Running Nested LOOCV...")
  for (i in seq_len(n_samples)) {
    train_i <- data_sel_model[-i, ]
    test_i  <- data_sel_model[i, , drop = FALSE]
    
    vars_i <- run_boruta_selection(train_i, predictors_list, boruta_threshold, n_seeds)
    loocv_vars[[i]] <- data.frame(fold = i, variable = vars_i)
    
    best_params_i <- run_grid_search(train_i, vars_i)
    loocv_params[[i]] <- best_params_i
    
    formula_i <- as.formula(paste("Ho ~", paste(vars_i, collapse = " + ")))
    model_i <- ranger::ranger(
      formula = formula_i,
      data = train_i,
      case.weights = if("n_samples_in_cell" %in% names(train_i)) sqrt(train_i$n_samples_in_cell) else NULL,
      num.trees = best_params_i$num.trees,
      mtry = best_params_i$mtry,
      min.node.size = best_params_i$min.node.size,
      replace = FALSE,
      seed = 100,
      importance = "permutation",
      num.threads = N_CORES
    )
    loocv_preds[i] <- predict(model_i, data = test_i)$predictions
    if (i %% 5 == 0) message(sprintf("  LOOCV Iteration %d/%d completed", i, n_samples))
  }
  
  loocv_df <- data.frame(Observed = data_sel_model$Ho, Predicted = loocv_preds)
  write.csv(loocv_df, file.path(outdir, paste0(sp_name, "_LOOCV_Performance_RF_VALUES_nested.csv")), row.names = FALSE)
  
  vars_df <- do.call(rbind, loocv_vars)
  write.csv(vars_df, file.path(outdir, paste0(sp_name, "_LOOCV_nested_vars_per_fold.csv")), row.names = FALSE)
  
  params_df <- do.call(rbind, loocv_params)
  params_df$fold <- seq_len(n_samples)
  write.csv(params_df, file.path(outdir, paste0(sp_name, "_LOOCV_nested_hyperparams_per_fold.csv")), row.names = FALSE)
  
  res_loocv   <- caret::postResample(loocv_df$Predicted, loocv_df$Observed)
  sse_test    <- sum((loocv_df$Observed - loocv_df$Predicted)^2)
  sst_test    <- sum((loocv_df$Observed - mean(loocv_df$Observed))^2)
  r2_nse_test <- 1 - (sse_test / sst_test)
  
  # ------------------------------------------------------------------------------
  # 5. FINAL PRODUCTION MODEL
  # ------------------------------------------------------------------------------
  vars_seleccionadas <- run_boruta_selection(data_sel, predictors_list, boruta_threshold, n_seeds)
  formula_rf_optimized <- as.formula(paste("Ho ~", paste(vars_seleccionadas, collapse = " + ")))
  best_params <- run_grid_search(data_sel_model, vars_seleccionadas)
  
  best_model <- ranger::ranger(
    formula = formula_rf_optimized,
    data = data_sel_model,
    case.weights = if("n_samples_in_cell" %in% names(data_sel_model)) sqrt(data_sel_model$n_samples_in_cell) else NULL,
    num.trees = best_params$num.trees,
    mtry = best_params$mtry,
    min.node.size = best_params$min.node.size,
    importance = "permutation",
    replace = FALSE,
    seed = 100,
    num.threads = N_CORES
  )
  
  # Variable Importance
  imps <- data.frame(
    var = names(best_model$variable.importance),
    imps = best_model$variable.importance / max(best_model$variable.importance)
  )
  
  p_imp <- ggplot(imps, aes(x = imps, y = reorder(var, imps))) +
    geom_point(size = 5, colour = "#ff6767") +
    labs(x = "Relative Importance", y = "Ecological Predictors") +
    theme_bw(14)
  
  ggsave(file.path(outdir, paste0(sp_name, "_Variable_Importance_RF.png")), p_imp, width = 7, height = 5, dpi = 300)
  
  # ------------------------------------------------------------------------------
  # 6. SPATIAL PREDICTION AND MESS MASK
  # ------------------------------------------------------------------------------
  sdm <- terra::rast(sdm_path) * 1
  sdm[sdm == 0] <- NA
  
  bios_selected <- bios[[vars_seleccionadas]]
  bios_selected <- terra::resample(bios_selected, sdm, method = "bilinear")
  bios_selected <- bios_selected * sdm
  valid_cells   <- terra::cells(bios_selected[[1]])
  
  temp.dt <- as.data.frame(terra::extract(bios_selected, valid_cells))
  temp.dt$cell_idx <- valid_cells
  
  valid_rows <- complete.cases(temp.dt[, vars_seleccionadas, drop = FALSE])
  temp.dt <- temp.dt[valid_rows, ]
  temp.dt$prediction <- NA_real_
  
  chunk_size <- 50000
  n_rows     <- nrow(temp.dt)
  n_chunks   <- ceiling(n_rows / chunk_size)
  
  for (i in seq_len(n_chunks)) {
    s_idx <- ((i - 1) * chunk_size) + 1
    e_idx <- min(i * chunk_size, n_rows)
    
    newdata_chunk <- temp.dt[s_idx:e_idx, vars_seleccionadas, drop = FALSE]
    pred_res <- predict(best_model, data = newdata_chunk, num.threads = N_CORES)
    temp.dt$prediction[s_idx:e_idx] <- pred_res$predictions
  }
  
  pred_raster <- bios_selected[[1]] * NA
  names(pred_raster) <- "Predicted_Ho"
  pred_raster[temp.dt$cell_idx] <- temp.dt$prediction
  
  # MESS calculation
  bios_selected_stack <- raster::stack(bios_selected)
  mess_RN <- terra::rast(dismo::mess(x = bios_selected_stack, v = data_sel[, vars_seleccionadas]))
  mess_RN <- terra::resample(mess_RN, sdm, method = "bilinear")
  
  mess_RN2 <- mess_RN
  mess_RN2[mess_RN2 <= 0 | is.infinite(mess_RN2)] <- NA
  mess_RN2[!is.na(mess_RN2)] <- 1
  
  pred_raster_no_interpolated <- pred_raster * mess_RN2
  terra::writeRaster(pred_raster_no_interpolated, file.path(outdir, paste0(sp_name, "_Ho_Macrogenetics_Map_RF_final.tif")), overwrite = TRUE)
  
  # ------------------------------------------------------------------------------
  # 7. COMPARISON AND SPATIAL AUTOCORRELATION (MORAN'S I)
  # ------------------------------------------------------------------------------
  train_preds <- predict(best_model, data = data_sel)$predictions
  train_df    <- data.frame(Observed = data_sel$Ho, Predicted = train_preds)
  res_train   <- caret::postResample(train_df$Predicted, train_df$Observed)
  
  p_train <- ggplot(train_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "darkgreen") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "darkgreen") +
    labs(title = "A) Final Model Fit (Train)", subtitle = paste("n =", n_samples, "cells"), x = "Observed Ho", y = "Predicted Ho") +
    theme_bw(13)
  
  p_test <- ggplot(loocv_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "blue") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "blue") +
    labs(title = "B) Nested LOOCV Validation (Test)", subtitle = "Real predictive capacity assessment", x = "Observed Ho", y = "Predicted Ho") +
    theme_bw(13)
  
  ggsave(file.path(outdir, paste0(sp_name, "_Model_Performance_Train_vs_LOOCV.png")), p_train + p_test, width = 12, height = 5.5, dpi = 300)
  
  # # Spatial autocorrelation with ape::Moran.I
  # residuos_loocv <- data_sel_model$Ho - loocv_preds
  # 
  # # Ensure single POINT geometries and strip spatial units
  # pts_sf <- sf::st_centroid(data_aggregated_sp)
  # coords_dist_mat <- units::drop_units(sf::st_distance(pts_sf)) / 1000
  # 
  # # Build spatial weights matrix
  # diag(coords_dist_mat) <- NA
  # pesos_espaciales <- 1 / coords_dist_mat
  # diag(pesos_espaciales) <- 0
  # pesos_espaciales[is.na(pesos_espaciales) | is.infinite(pesos_espaciales)] <- 0
  # 
  # # Row-standardize weights
  # row_sums <- rowSums(pesos_espaciales)
  # pesos_espaciales <- pesos_espaciales / ifelse(row_sums == 0, 1, row_sums)
  # 
  # # Calculate Moran's I
  # moran_result <- ape::Moran.I(residuos_loocv, pesos_espaciales)
  # 
  # write.csv(
  #   data.frame(observed = moran_result$observed, expected = moran_result$expected, sd = moran_result$sd, p_value = moran_result$p.value),
  #   file.path(outdir, paste0(sp_name, "_Moran_I_LOOCV_residuals.csv")),
  #   row.names = FALSE
  # )
  
  message("Process successfully finished for: ", sp_name)
  save.image(file.path(outdir, paste0(sp_name, "_Ho_Macrogenetics_Map_RF_final.RData")))
}

CA <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/CA_RF_LONLAT_COORDS",
               sp_name = "Crocodylus acutus",
               raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
               data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
               sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif",
               N_CORES = 8,
               boruta_threshold = 30,
               n_seeds = 30,
               cor_cutoff = 0.5,
               addLonLat = T,
               aggregate_occs_cells = TRUE)


CA2 <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/CA_RF_LONLAT_POINTS",
               sp_name = "Crocodylus acutus",
               raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
               data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
               sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif",
               N_CORES = 8,
               boruta_threshold = 30,
               n_seeds = 30,
               cor_cutoff = 0.5,
               addLonLat = T,
               aggregate_occs_cells = F)
CA3 <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/CA_RF_NO_LONLAT_COORDS",
               sp_name = "Crocodylus acutus",
               raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
               data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
               sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif",
               N_CORES = 8,
               boruta_threshold = 30,
               n_seeds = 30,
               cor_cutoff = 0.5,
               addLonLat = F,
               aggregate_occs_cells = TRUE)


CA4 <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/CA_RF_NO_LONLAT_POINTS",
                sp_name = "Crocodylus acutus",
                raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
                data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
                sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif",
                N_CORES = 8,
                boruta_threshold = 30,
                n_seeds = 30,
                cor_cutoff = 0.5,
                addLonLat = F,
                aggregate_occs_cells = F)

################################################################################
CM <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/CM_RF_LONLAT_COORDS",
               sp_name = "Crocodylus moreletii",
               raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
               data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
               sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif",
               boruta_threshold = 30,
               n_seeds = 30,
               cor_cutoff = 0.5,
               addLonLat = T,
               aggregate_occs_cells = TRUE)


CM2 <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/CM_RF_LONLAT_POINTS",
                sp_name = "Crocodylus moreletii",
                raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
                data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
                sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif",
                N_CORES = 8,
                boruta_threshold = 30,
                n_seeds = 30,
                cor_cutoff = 0.5,
                addLonLat = T,
                aggregate_occs_cells = F)
CM3 <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/CM_RF_NO_LONLAT_COORDS",
                sp_name = "Crocodylus moreletii",
                raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
                data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
                sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif",
                N_CORES = 8,
                boruta_threshold = 30,
                n_seeds = 30,
                cor_cutoff = 0.5,
                addLonLat = F,
                aggregate_occs_cells = TRUE)


CM4 <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/CM_RF_NO_LONLAT_POINTS",
                sp_name = "Crocodylus moreletii",
                raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
                data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
                sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif",
                N_CORES = 8,
                boruta_threshold = 30,
                n_seeds = 30,
                cor_cutoff = 0.5,
                addLonLat = F,
                aggregate_occs_cells = F)