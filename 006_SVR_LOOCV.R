#' @title Spatial Interpolation, Cross-Validation, and Mapping of Genetic Diversity (Ho)
#' @description Processes genetic data, aggregates points by SDM raster cells, 
#' performs IDW cross-validation (LOOCV at points), builds diagnostic plots 
#' using ggpmisc, and exports high-resolution maps via tmap and workspace images.
#' @author Jorge (Modified via AI collaborator)
#' @date 2026-09-22

# 1. Load Required Libraries ----------------------------------------------
library(dplyr)
library(terra)
library(spatstat)
library(sf)
library(readxl)
library(Metrics)
library(ggplot2)
library(ggpmisc)
library(patchwork)
library(tmap)

# 2. Define Parameters and Paths ------------------------------------------
outdir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/CA_RF_NO_LONLAT_COORDS"
sp_name <- "Crocodylus acutus"
raster_dir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s"
data_path <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx"
sdm_path <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/Crocodylus_acutus/Crocodylus_acutus_Binario_P10.tif"

aggregate_occs_cells <- TRUE
suffix <- "_IDW"

# Create output directory if it doesn't exist
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# 3. Load SDM Raster First (needed for cell-level aggregation) -----------
sdm <- terra::rast(sdm_path) * 1
sdm[sdm == 0] <- NA

# 4. Load, Filter, and Aggregate Data by SDM Raster Cell ------------------
raw_data <- readxl::read_xlsx(data_path, sheet = "data") %>%
  dplyr::filter(sp == sp_name, !is.na(Ho), Ho > 0, !is.na(lat), !is.na(lon))

# Assign each point to its corresponding SDM raster cell ID
coords_matrix <- as.matrix(raw_data[, c("lon", "lat")])
raw_data$cell_id <- terra::cellFromXY(sdm, coords_matrix)

# Filter out points that fall outside the SDM valid area
raw_data <- raw_data %>% dplyr::filter(!is.na(cell_id))

# Keep raw points as sf for mapping purposes later
my_sf_object <- sf::st_as_sf(raw_data, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Aggregate Ho by raster cell (mean, standard deviation, and sample count n)
data_aggregated <- raw_data %>%
  dplyr::group_by(cell_id) %>%
  dplyr::summarise(
    Ho = mean(Ho, na.rm = TRUE),
    Ho_sd = sd(Ho, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

# Replace NA standard deviations (occurring when n = 1) with 0
data_aggregated$Ho_sd[is.na(data_aggregated$Ho_sd)] <- 0

# Get the exact center coordinates of those SDM raster cells
cell_coords <- terra::xyFromCell(sdm, data_aggregated$cell_id)
data_aggregated$lon <- cell_coords[, 1]
data_aggregated$lat <- cell_coords[, 2]

# Create sf object for aggregated cells
data_aggregated_sf <- sf::st_as_sf(data_aggregated, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# 5. Prepare Spatial Objects and Window -----------------------------------
points_vect <- terra::vect(data_aggregated_sf)
coords <- terra::crds(points_vect)

sdm_ext <- terra::ext(sdm)

win <- spatstat.geom::owin(
  xrange = c(sdm_ext[1], sdm_ext[2]), 
  yrange = c(sdm_ext[3], sdm_ext[4])
)

# Create Marked Point Pattern (ppp) object
ppp_obj <- spatstat.geom::ppp(
  x = coords[, 1], 
  y = coords[, 2], 
  window = win
)
marks(ppp_obj) <- data_aggregated$Ho

# 6. Cross-Validation to Optimize IDW Power -------------------------------
powers <- seq(0.001, 10, 0.01)
mse_result <- NULL

for(power in powers){
  print(paste("Testing power:", power))
  CV_idw <- spatstat.explore::idw(ppp_obj, power = power, at = "points")
  mse_result <- c(mse_result, Metrics::mse(ppp_obj$marks, CV_idw))
}

optimal_power <- powers[which.min(mse_result)]
plot(powers, mse_result, type = "l", main = "MSE vs. IDW Power (Cell-aggregated)", xlab = "Power", ylab = "MSE")

# 7. Predicted vs. Observed Evaluation Plot (Train vs CV Test Style) --------
# Resubstitution (Training Fit: predicting the points back with IDW including themselves)
train_preds <- as.numeric(spatstat.explore::idw(ppp_obj, power = optimal_power, at = "points", leaveoneout = FALSE))
train_df <- data.frame(Observed = ppp_obj$marks, Predicted = train_preds)

# Leave-One-Out Cross-Validation (Test Fit)
cv_preds <- spatstat.explore::idw(ppp_obj, power = optimal_power, at = "points", leaveoneout = TRUE)
cv_df <- data.frame(Observed = ppp_obj$marks, Predicted = as.numeric(cv_preds))

n_samples <- nrow(data_aggregated)

if(isTRUE(aggregate_occs_cells)){
  file_to_save_plot <- file.path(outdir, paste0(sp_name, "_IDW_Performance_CELLS", suffix, ".png"))
} else {
  file_to_save_plot <- file.path(outdir, paste0(sp_name, "_IDW_Performance", suffix, ".png"))  
}

p_train <- ggplot(train_df, aes(x = Observed, y = Predicted)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  stat_poly_line(color = "darkgreen") +
  stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
  geom_point(size = 3, alpha = 0.8, color = "darkgreen") +
  labs(title = "A) IDW",
       subtitle = paste("n =", n_samples, "cells (Resubstitution)"),
       x = "Observed Ho", y = "Predicted Ho") +
  theme_bw(13)

ggsave(
  file_to_save_plot,
  p_train, width = 12, height = 5.5, dpi = 600
)

# 8. Generate Final Raster Output -----------------------------------------
idw_raster <- terra::rast(
  spatstat.explore::idw(ppp_obj, power = optimal_power, at = "pixels"),
  crs = "EPSG:4326"
)

idw_raster <- terra::resample(idw_raster, sdm)
idw_raster <- terra::mask(idw_raster, sdm)

# Assign a clean name to avoid "lyr.1" in tmap
names(idw_raster) <- "Predicted_Ho"
pred_raster_no_interpolated <- idw_raster 

# 9. Map Generation with tmap ---------------------------------------------
data("World", package = "tmap")

if(isTRUE(aggregate_occs_cells)){
  file_to_save_pdf <- file.path(outdir, paste0(sp_name, "_Ho_Map_IDW_CELLS", suffix, ".pdf"))
  file_to_save_rdata <- file.path(outdir, paste0(sp_name, "_Ho_Workspace_CELLS", suffix, ".RData"))
} else {
  file_to_save_pdf <- file.path(outdir, paste0(sp_name, "_Ho_Map_IDW", suffix, ".pdf"))
  file_to_save_rdata <- file.path(outdir, paste0(sp_name, "_Ho_Workspace", suffix, ".RData"))
}

map_out <- tm_shape(pred_raster_no_interpolated) + 
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

tmap::tmap_save(
  tm       = map_out, 
  filename = file_to_save_pdf, 
  width    = 25, 
  height   = 15, 
  units    = "cm", 
  dpi      = 600
)

# 10. Save Workspace Environment ------------------------------------------
save.image(file_to_save_rdata)