library(geodata)

outdir <- "D:/LAND_USE_GEODATA"
# Descarga capas de fracción de cobertura Copernicus (100m nativo)
# var: "trees", "grassland", "shrubs", "cropland", "built up", 
#      "bare", "snow", "water permanent", "water seasonal", "moss lichen"
path_lc <- file.path(outdir, "landcover")
dir.create(path_lc, showWarnings = FALSE)

# Ejemplo: descargar fracción de agua y bosque (ajusta 'var' a lo que necesites)
lc_water <- geodata::landcover(var = "water", path = path_lc)
lc_mangroves <- geodata::landcover(var = "mangroves", path = path_lc)
lc_wetland <- geodata::landcover(var = "wetland", path = path_lc)
lc_trees <- geodata::landcover(var = "trees", path = path_lc)


template <- rast("D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/Latin America/Latin America 30S/template_30s.tif")

lc_water <- terra::resample(lc_water,template)
lc_mangroves <- terra::resample(lc_mangroves,template)
lc_wetland <- terra::resample(lc_wetland,template)
lc_trees <- terra::resample(lc_trees,template)

dir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s"
terra::writeRaster(lc_water,paste0(dir,"/","lc_water.tif"))
terra::writeRaster(lc_mangroves,paste0(dir,"/","lc_mangroves.tif"))
terra::writeRaster(lc_wetland,paste0(dir,"/","lc_wetland.tif"))
terra::writeRaster(lc_trees,paste0(dir,"/","lc_trees.tif"))
