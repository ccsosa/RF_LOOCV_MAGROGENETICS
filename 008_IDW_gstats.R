#' @title Spatial Interpolation, Cross-Validation, and Mapping of Genetic Diversity (Ho)
#' @description Encapsulated function to process genetic data, aggregate by SDM cells,
#' optimize IDW power (idp) via LOOCV using gstat, and generate high-resolution maps and GeoTiffs.
#' @author Chrystian Sosa (Modified via AI collaborator) -- gstat version
#' @date 2026-09-22

# 0. Global Parameters -----------------------------------------------------
N_CORES <- 6  # Fixed number of cores for parallel idp search. Modify here if needed.
#https://rpubs.com/maquiroga/631984
#https://www.paulamoraga.com/book-spatial/spatial-interpolation-methods.html
# 1. Load Required Libraries ----------------------------------------------
library(dplyr)
library(terra)
library(gstat)
library(sp)
library(sf)
library(readxl)
library(Metrics)
library(ggplot2)
library(ggpmisc)
library(patchwork)
library(tmap)
library(parallel)

# 2. Define Function ------------------------------------------------------
run_idw_analysis <- function(
    sp_name,
    data_path,
    sdm_path,
    outdir,
    suffix = "_IDW_gstat",
    idp_seq = seq(0.001, 10, 0.1),  # coarser than spatstat version (gstat.cv is heavier); tune as needed
    n_cores = N_CORES
) {
  message(paste(">>> Starting gstat IDW analysis for:", sp_name))
  
  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # Load SDM Raster
  sdm <- terra::rast(sdm_path) * 1
  sdm[sdm == 0] <- NA
  
  # Load, Filter, and Aggregate Data by SDM Raster Cell
  raw_data <- readxl::read_xlsx(data_path, sheet = "data") %>%
    dplyr::filter(sp == sp_name, !is.na(Ho), Ho > 0, !is.na(lat), !is.na(lon))
  
  coords_matrix <- as.matrix(raw_data[, c("lon", "lat")])
  raw_data$cell_id <- terra::cellFromXY(sdm, coords_matrix)
  raw_data <- raw_data %>% dplyr::filter(!is.na(cell_id))
  
  my_sf_object <- sf::st_as_sf(raw_data, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  data_aggregated <- raw_data %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      Ho = mean(Ho, na.rm = TRUE),
      Ho_sd = sd(Ho, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  data_aggregated$Ho_sd[is.na(data_aggregated$Ho_sd)] <- 0
  
  cell_coords <- terra::xyFromCell(sdm, data_aggregated$cell_id)
  data_aggregated$lon <- cell_coords[, 1]
  data_aggregated$lat <- cell_coords[, 2]
  
  data_aggregated_sf <- sf::st_as_sf(data_aggregated, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  # gstat works on sp (Spatial*) objects
  data_sp <- sf::as_Spatial(data_aggregated_sf)
  
  # Cross-Validation to Optimize IDW Power (idp) -- parallelized over n_cores
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, { library(gstat); library(sp); library(Metrics) })
  parallel::clusterExport(cl, varlist = c("data_sp"), envir = environment())
  
  mse_result <- parallel::parSapply(cl, idp_seq, function(p) {
    idw_model <- gstat::gstat(formula = Ho ~ 1, locations = data_sp, set = list(idp = p))
    cv_i <- gstat::gstat.cv(idw_model, nfold = nrow(data_sp), verbose = FALSE)
    Metrics::mse(cv_i$observed, cv_i$var1.pred)
  })
  
  parallel::stopCluster(cl)
  on.exit(NULL)  # cluster already stopped, avoid double stop
  
  optimal_power <- idp_seq[which.min(mse_result)]
  
  # Final LOOCV at optimal power
  idw_model_opt <- gstat::gstat(formula = Ho ~ 1, locations = data_sp, set = list(idp = optimal_power))
  cv_opt <- gstat::gstat.cv(idw_model_opt, nfold = nrow(data_sp), verbose = FALSE)
  rmse_val <- Metrics::rmse(cv_opt$observed, cv_opt$var1.pred)
  
  # In-sample (training) predictions: predict at the same locations used to fit the model.
  # Note: IDW is an exact interpolator, so this fit will be near-perfect (RMSE ~ 0) by
  # construction -- it's shown only as a contrast against the honest LOOCV/test performance.
  train_pred <- predict(idw_model_opt, newdata = data_sp, debug.level = 0)
  rmse_train <- Metrics::rmse(data_sp$Ho, train_pred$var1.pred)
  
  # Print key metrics directly to console
  message(sprintf("--- Results for %s ---", sp_name))
  message(sprintf("Optimal Power (idp): %.4f", optimal_power))
  message(sprintf("Training RMSE (in-sample, exact-interp): %.4f", rmse_train))
  message(sprintf("LOOCV RMSE (test):                       %.4f", rmse_val))
  
  # Performance Plots: A) Training (in-sample) vs B) LOOCV (test)
  train_df <- data.frame(Observed = data_sp$Ho, Predicted = train_pred$var1.pred)
  test_df  <- data.frame(Observed = cv_opt$observed, Predicted = cv_opt$var1.pred)
  n_samples <- nrow(data_aggregated)
  
  file_to_save_plot <- file.path(outdir, paste0(sp_name, "_IDW_Performance_CELLS", suffix, ".png"))
  
  p_train <- ggplot(train_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "darkgreen") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "darkgreen") +
    labs(
      title = paste("A) Training (in-sample) -", sp_name),
      subtitle = sprintf("n = %d cells | Opt. idp: %.3f | RMSE: %.4f", n_samples, optimal_power, rmse_train),
      x = "Observed Ho", y = "Predicted Ho"
    ) +
    theme_bw(13)
  
  p_test <- ggplot(test_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "steelblue4") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "steelblue4") +
    labs(
      title = paste("B) LOOCV (test) -", sp_name),
      subtitle = sprintf("n = %d cells | Opt. idp: %.3f | RMSE: %.4f", n_samples, optimal_power, rmse_val),
      x = "Observed Ho", y = "Predicted Ho"
    ) +
    theme_bw(13)
  
  p_combined <- p_train + p_test
  
  ggsave(file_to_save_plot, p_combined, width = 14, height = 5.5, dpi = 600)
  
  # Generate Final Raster Output -------------------------------------------
  # Build prediction grid from all non-NA SDM cell centers
  grid_df <- as.data.frame(sdm, xy = TRUE, na.rm = TRUE)[, c("x", "y")]
  sp::coordinates(grid_df) <- ~x + y
  sp::proj4string(grid_df) <- sp::CRS("EPSG:4326")
  
  idw_pred <- predict(idw_model_opt, newdata = grid_df)
  
  idw_pred_df <- as.data.frame(idw_pred)[, c("x", "y", "var1.pred")]
  idw_raster <- terra::rast(idw_pred_df, type = "xyz", crs = "EPSG:4326")
  
  idw_raster <- terra::resample(idw_raster, sdm)
  idw_raster <- terra::mask(idw_raster, sdm)
  names(idw_raster) <- "Predicted_Ho"
  
  # Save Raster as GeoTiff
  raster_output_path <- file.path(outdir, paste0(sp_name, "_Ho_raster_CELLS", suffix, ".tif"))
  writeRaster(idw_raster, raster_output_path, overwrite = TRUE)
  
  # Map Generation with tmap
  data("World", package = "tmap")
  
  file_to_save_pdf <- file.path(outdir, paste0(sp_name, "_Ho_Map_IDW_CELLS", suffix, ".pdf"))
  file_to_save_rdata <- file.path(outdir, paste0(sp_name, "_Ho_Workspace_CELLS", suffix, ".RData"))
  
  map_out <- tm_shape(idw_raster) +
    tm_raster(
      col.legend = tm_legend(title = "Predicted Ho"),
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
  
  tmap::tmap_save(tm = map_out, filename = file_to_save_pdf, width = 25, height = 15, units = "cm", dpi = 600)
  
  # Save Workspace Environment
  save.image(file_to_save_rdata)
  message(paste(">>> Successfully completed analysis for:", sp_name))
}

# 3. Example Execution ----------------------------------------------------
run_idw_analysis(
  sp_name = "Crocodylus acutus",
  data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
  sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif",
  outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/IDW_CA_gstat",
  suffix = "_IDW_gstat"
)

run_idw_analysis(
  sp_name = "Crocodylus moreletii",
  data_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx",
  sdm_path = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif",
  outdir = "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/IDW_CM_gstat",
  suffix = "_IDW_gstat"
)