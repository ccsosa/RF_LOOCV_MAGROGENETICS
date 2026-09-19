library(terra)
library(sf)

dir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/NaturalEarthData"
outdir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/Latin America"
template <- terra::rast("D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/Latin America/Latin America 30S/template_30s.tif")
vars <- c("ne_10m_coastline", "ne_10m_rivers_lake_centerlines", "ne_10m_urban_areas","coast_rivers")

crs_planar <- "+proj=aea +lat_1=19 +lat_2=-41 +lat_0=-11 +lon_0=-76 +datum=WGS84 +units=m"

for (i in 1:length(vars)) {
  template2 <- terra::init(template, fun = NA)
  feat <- vect(paste0(dir, "/", vars[[i]], ".shp"))
  feat_crop <- terra::crop(feat, template)
  feat_lines <- terra::rasterize(feat_crop, template, touches = TRUE)

  feat_lines_planar <- terra::project(feat_lines, crs_planar, method = "near")

  x_dist_planar <- terra::distance(feat_lines_planar, unit = "km")
  # plot(x_dist_planar)

  x_dist <- terra::project(x_dist_planar, template2)
  x_dist <- terra::mask(x_dist, template)
  crs(x_dist) <- crs(template)
  par(mfrow=c(1,2))
  plot(feat_crop)
  plot(x_dist)
  terra::writeRaster(x_dist, paste0(outdir, "/", vars[[i]], "_30s_dist.tif"), overwrite = TRUE)
}

# 
# library(terra)
# library(sf)
# 
# dir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/NaturalEarthData"
# outdir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/Latin America"
# template <- terra::rast("D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/climate/wc2.1_2.5m_bio_1.tif")
# vars <- c("ne_10m_coastline", "ne_10m_rivers_lake_centerlines", "ne_10m_urban_areas")
# 
# for (i in seq_along(vars)) {
#   print(i)
#   template2 <- terra::init(template, fun = NA)
#   feat <- vect(paste0(dir, "/", vars[[i]], ".shp"))
#   feat_crop <- terra::crop(feat, template)
#   feat_lines <- terra::rasterize(feat_crop, template, touches = TRUE)
#   
#   x_dist <- terra::distance(feat_lines, unit = "km", method = "haversine")
#   x_dist <- terra::mask(x_dist, template)
#   
#   terra::writeRaster(x_dist, paste0(outdir, "/", vars[[i]], "_dist.tif"), overwrite = TRUE)
# }
