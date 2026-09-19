################################################################################
# MACROGENETICS PREDICTION WORKFLOW (RANDOM FOREST VIA RANGER)
# obtain lon lat layers
################################################################################

suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

# directory containing the environmental raster layers
raster_dir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s" 
x <- rast(list.files(raster_dir[[1]],full.names = T))
# Crear capas de coordenadas XY con la misma extensión y resolución que bios
lon_layer <- terra::init(x[[1]], "x")
lat_layer <- terra::init(x[[1]], "y")

names(lon_layer) <- "lon"
names(lat_layer) <- "lat"

writeRaster(lon_layer,"D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s/lon.tif",overwrite=T)
writeRaster(lat_layer,"D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s/lat.tif",overwrite=T)
