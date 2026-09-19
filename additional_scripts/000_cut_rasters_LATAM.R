library(terra)

# 1. Define paths
dir_in  <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval/climate/wc2.1_2.5m"
dir_out <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER"
tmpl_path <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/Latin America/template.tif"

# 2. Get files and load as a multi-layer SpatRaster directly
files <- list.files(dir_in, pattern = "\\.tif$", full.names = TRUE)
data  <- terra::rast(files)

# 3. Load template and crop (mask optional if you need exact shape match)
template <- terra::rast(tmpl_path)
data2    <- terra::crop(data, template)

# 4. Write all layers in a single, optimized command
terra::writeRaster(
  data2, 
  filename  = file.path(dir_out, paste0(names(data2), ".tif")),
  overwrite = TRUE
)