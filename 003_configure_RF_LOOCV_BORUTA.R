################################################################################
# MACROGENETICS PREDICTION WORKFLOW (RANDOM FOREST VIA RANGER)
# Continuous geographic prediction of heterozygosity (Ho) over a raster.
################################################################################

suppressPackageStartupMessages({
  library(terra)
  library(fields)
  library(adespatial)
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
})

outdir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# species name, used for output file naming
sp_name <- "Crocodylus acutus"
# directory containing the environmental raster layers
raster_dir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s" 
# path to the Ho (heterozygosity) data table
data_path <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx"
# path to the binary SDM (species distribution model) raster
sdm_path <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif"
# number of cores used by ranger
N_CORES <- 6L  # modifiable parameter, never detectCores()-1
# Boruta consensus threshold (out of n_seeds runs)
boruta_threshold <- 30   # change to 8L if you'd rather require 8/10 instead of full consensus
# number of random seeds used to assess Boruta stability
n_seeds <- 30

#' Nested LOOCV Random Forest workflow for spatial prediction of Ho
#'
#' Fits a Random Forest (ranger) model to predict a continuous genetic
#' heterozygosity value (Ho) from environmental (bioclimatic) covariates,
#' validates it with a fully nested Leave-One-Out Cross-Validation (both
#' Boruta variable selection AND hyperparameter grid search are redone
#' independently within every fold, so there is no information leakage
#' from the held-out cell into either step), fits a final production model
#' on the complete dataset, predicts across the species' raster extent
#' (restricted to the SDM footprint), masks the prediction to the region of
#' environmental interpolation using a standard MESS mask, and checks the
#' LOOCV residuals for spatial autocorrelation (Moran's I).
#'
#' @param outdir Output directory where all results (CSVs, plots, rasters,
#'   the saved .RData workspace) are written.
#' @param sp_name Species name, used as a prefix for every output file.
#' @param raster_dir Directory containing the bioclimatic/environmental
#'   raster layers (.tif files) used as predictors.
#' @param data_path Path to the Excel file with the genetic data table
#'   (must contain a "data" sheet with columns sp, Ho, lon, lat).
#' @param sdm_path A SpatRaster filename: The filename must be a binary
#'  species distribution model used to  restrict the prediction/extrapolation area
#' @param N_CORES Integer. Number of CPU threads used by ranger::ranger()
#'   and ranger's predict(). Kept as an explicit, user-modifiable parameter
#'   instead of parallel::detectCores() - 1.
#' @param boruta_threshold Integer. Minimum number of Boruta runs (out of
#'   n_seeds) in which a variable must be "Confirmed" for it to be included
#'   in the model. 10 = full consensus across all seeds.
#' @param n_seeds Integer. Number of random seeds used to repeat the Boruta
#'   variable-selection procedure, so its stability can be assessed.
#'
#' @details
#' Workflow steps performed inside the function:
#' \enumerate{
#'   \item Load bioclimatic rasters and the genetic (Ho) point data, filter
#'     to the target species, and build an sf/SpatVector point object.
#'   \item Snap each point to its nearest valid raster cell (buffering
#'     out-of-coverage points by 5 km when needed) and aggregate multiple
#'     points falling in the same cell (mean Ho, sample count per cell).
#'   \item Extract environmental values using a focal-filled raster (to
#'     avoid NAs at coastal/marine edges), with a nearest-pixel fallback
#'     for any remaining NA (coastal) cells.
#'   \item Run a fully nested Leave-One-Out Cross-Validation: for every
#'     held-out cell, Boruta variable selection and a ranger hyperparameter
#'     grid search (num.trees x mtry x min.node.size, chosen by OOB R2) are
#'     both re-run using only the n-1 training cells, then a model is fit
#'     and used to predict the held-out cell. This avoids leaking any
#'     information about a test cell into variable selection or tuning.
#'   \item Report LOOCV performance (RMSE, Pearson R2, Nash-Sutcliffe R2),
#'     and the stability of selected variables and hyperparameters across
#'     folds.
#'   \item Fit a final production Random Forest on the complete dataset
#'     (Boruta selection + grid search re-run on all n_samples cells — this
#'     step is not validated against those same cells, so it does not leak
#'     into the LOOCV estimate above).
#'   \item Predict the final model across the raster extent in chunks (for
#'     memory efficiency), restricted to the SDM (species distribution
#'     model) footprint.
#'   \item Compute a standard MESS (Multivariate Environmental Similarity
#'     Surface) mask (Elith et al. 2010): cells with MESS <= 0 indicate
#'     environmental extrapolation and are masked out (NA); MESS > 0 cells
#'     (interpolation) are kept.
#'   \item Export prediction maps (raster + tmap PDF), a variable
#'     importance plot (permutation importance from the final model), and
#'     a two-panel comparison of in-sample (train) vs. LOOCV (test)
#'     observed-vs-predicted performance.
#'   \item Check the LOOCV residuals for spatial autocorrelation using
#'     Moran's I (inverse-distance spatial weights) and save the result.
#' }
#'
#' All outputs are written to \code{outdir}, prefixed with \code{sp_name}.
#'
#' @return Invisibly nothing; called for its side effects (files written to
#'   \code{outdir} and an .RData workspace image saved at the end).
#'   
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

RF_LOOCV <- function(outdir, sp_name, raster_dir, data_path, sdm_path, N_CORES = 6L, boruta_threshold = 30, n_seeds = 30) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # ------------------------------------------------------------------------------
  # 1. RASTER LAYERS AND DATA LOADING & PREPROCESSING
  # ------------------------------------------------------------------------------
  bios_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)
  bios_names <- sub("\\.tif$", "", list.files(raster_dir, pattern = "\\.tif$"))
  
  bios <- terra::rast(bios_files)
  names(bios) <- bios_names
  
  data <- readxl::read_xlsx(data_path, sheet = "data") %>% 
    dplyr::filter(sp == sp_name, !is.na(Ho), !is.na(lat), !is.na(lon))
  
  my_sf_object <- sf::st_as_sf(data, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  puntos_vect  <- terra::vect(my_sf_object)
  
  # ------------------------------------------------------------------------------
  # 2. NEAREST CELL ASSIGNMENT AND SPATIAL AGGREGATION PER CELL
  # ------------------------------------------------------------------------------
  cell_ids <- terra::cells(bios[[1]], puntos_vect)[, "cell"]
  na_cells <- which(is.na(cell_ids))
  
  if (length(na_cells) > 0) {
    message(sprintf("Assigning %d out-of-coverage points to the nearest valid raster cell...", length(na_cells)))
    extracted_cells <- terra::extract(
      bios[[1]],
      puntos_vect[na_cells, ],
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
      Ho  = mean(Ho, na.rm = TRUE),
      lon = mean(cell_lon, na.rm = TRUE),
      lat = mean(cell_lat, na.rm = TRUE),
      .groups = "drop"
    )
  
  bios_focal <- terra::focal(bios, w = 15, fun = mean, na.policy = "only", na.rm = TRUE)
  gc()
  
  puntos_unicos_vect <- terra::vect(data_aggregated_sf)
  ext_direct <- terra::extract(bios_focal, puntos_unicos_vect)
  
  if (any(!complete.cases(ext_direct))) {
    message("Coastal points with NAs detected. Extracting from the nearest coastal pixel...")
    na_rows <- which(!complete.cases(ext_direct))
    for (i in na_rows) {
      ext_nearest <- terra::extract(bios_focal, puntos_unicos_vect[i, ], nearest = TRUE)
      ext_direct[i, names(bios)] <- ext_nearest[, names(bios)]
    }
  }
  
  ext_direct_clean <- ext_direct %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE))) %>%
    dplyr::select(-ID)
  
  data_aggregated_sf_sp <- cbind(data_aggregated_sf, ext_direct_clean)
  rm(ext_direct_clean); gc()
  
  # ------------------------------------------------------------------------------
  # 3. VARIABLE SELECTION AND DECOLLINEARIZATION
  # ------------------------------------------------------------------------------
  data_sel <- data_aggregated_sf_sp %>%
    dplyr::select(Ho, n_samples_in_cell, lon, lat, dplyr::all_of(names(bios))) %>%
    sf::st_drop_geometry()
  
  valid_rows <- complete.cases(data_sel)
  data_sel   <- data_sel[valid_rows, ]
  data_aggregated_sf_sp_sel <- data_aggregated_sf_sp[valid_rows, ]
  
  n_samples  <- nrow(data_sel)
  cor_matrix <- cor(data_sel[, names(bios)])
  to_remove  <- caret::findCorrelation(cor_matrix, cutoff = 0.6)
  
  vars <- if (length(to_remove) > 0) colnames(cor_matrix)[-to_remove] else colnames(cor_matrix)
  
  data_sel <- data_sel[c("Ho", "n_samples_in_cell", vars)]
  data_aggregated_sf_sp_sel <- data_aggregated_sf_sp_sel[, c("Ho", "n_samples_in_cell", vars)]
  
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
        case.weights = sqrt(train_data$n_samples_in_cell),
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
    train_i <- data_sel[-i, ]
    test_i  <- data_sel[i, , drop = FALSE]
    
    vars_i <- run_boruta_selection(train_i, vars, boruta_threshold, n_seeds)
    loocv_vars[[i]] <- data.frame(fold = i, variable = vars_i)
    
    best_params_i <- run_grid_search(train_i, vars_i)
    loocv_params[[i]] <- best_params_i
    
    formula_i <- as.formula(paste("Ho ~", paste(vars_i, collapse = " + ")))
    model_i <- ranger::ranger(
      formula = formula_i,
      data = train_i,
      case.weights = sqrt(train_i$n_samples_in_cell),
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
  
  loocv_df <- data.frame(Observed = data_sel$Ho, Predicted = loocv_preds)
  write.csv(loocv_df, file.path(outdir, paste0(sp_name, "_LOOCV_Performance_RF_VALUES_nested.csv")), row.names = FALSE)
  
  vars_df <- do.call(rbind, loocv_vars)
  write.csv(vars_df, file.path(outdir, paste0(sp_name, "_LOOCV_nested_vars_per_fold.csv")), row.names = FALSE)
  
  params_df <- do.call(rbind, loocv_params)
  params_df$fold <- seq_len(n_samples)
  write.csv(params_df, file.path(outdir, paste0(sp_name, "_LOOCV_nested_hyperparams_per_fold.csv")), row.names = FALSE)
  
  res_loocv <- caret::postResample(loocv_df$Predicted, loocv_df$Observed)
  sse_test  <- sum((loocv_df$Observed - loocv_df$Predicted)^2)
  sst_test  <- sum((loocv_df$Observed - mean(loocv_df$Observed))^2)
  r2_nse_test <- 1 - (sse_test / sst_test)
  
  # ------------------------------------------------------------------------------
  # 5. FINAL PRODUCTION MODEL
  # ------------------------------------------------------------------------------
  vars_seleccionadas <- run_boruta_selection(data_sel, vars, boruta_threshold, n_seeds)
  formula_rf_optimized <- as.formula(paste("Ho ~", paste(vars_seleccionadas, collapse = " + ")))
  best_params <- run_grid_search(data_sel, vars_seleccionadas)
  
  best_model <- ranger::ranger(
    formula = formula_rf_optimized,
    data = data_sel,
    case.weights = sqrt(data_sel$n_samples_in_cell),
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
  
  # Spatial autocorrelation with ape::Moran.I
  residuos_loocv <- data_sel$Ho - loocv_preds
  pts_sf <- data_aggregated_sf_sp_sel %>% sf::st_cast("POINT")
  
  coords_dist_mat <- as.matrix(sf::st_distance(pts_sf)) / 1000
  diag(coords_dist_mat) <- NA
  pesos_espaciales <- 1 / coords_dist_mat
  diag(pesos_espaciales) <- 0
  pesos_espaciales[is.na(pesos_espaciales) | is.infinite(pesos_espaciales)] <- 0
  
  row_sums <- rowSums(pesos_espaciales)
  pesos_espaciales <- pesos_espaciales / ifelse(row_sums == 0, 1, row_sums)
  
  moran_result <- ape::Moran.I(residuos_loocv, pesos_espaciales)
  
  write.csv(
    data.frame(observed = moran_result$observed, expected = moran_result$expected, sd = moran_result$sd, p_value = moran_result$p.value),
    file.path(outdir, paste0(sp_name, "_Moran_I_LOOCV_residuals.csv")),
    row.names = FALSE
  )
  
  message("Process successfully finished for: ", sp_name)
  save.image(file.path(outdir, paste0(sp_name, "_Ho_Macrogenetics_Map_RF_final.RData")))
}

CA <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics",
               sp_name = "Crocodylus acutus",
               raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
               data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
               sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif",
               N_CORES = 8,
               boruta_threshold = 30,
               n_seeds = 30)

CI <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics",
               sp_name = "Crocodylus intermedius",
               raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
               data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
               sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_intermedius/Crocodylus_intermedius_Binario_P10.tif",
               N_CORES = 8,
               boruta_threshold = 30,
               n_seeds = 30)

CM <- RF_LOOCV(outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics",
               sp_name = "Crocodylus moreletii",
               raster_dir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s",
               data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
               sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif",
               N_CORES = 8,
               boruta_threshold = 30,
               n_seeds = 30)
