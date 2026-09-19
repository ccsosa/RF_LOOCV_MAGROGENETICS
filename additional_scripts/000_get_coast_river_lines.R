library(sf)
library(terra)
CL <- sf::st_read(
"D:/PROGRAMAS/Dropbox/TESIS_JORGE/NaturalEarthData/ne_10m_rivers_lake_centerlines.shp")
COAST <- sf::st_read("D:/PROGRAMAS/Dropbox/TESIS_JORGE/NaturalEarthData/ne_10m_coastline.shp")

# Load raster to define the extent (adjust path)

# Cargar capas bioclimáticas
raster_dir <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/RASTER/test_layers_30s"
bios_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)
bios_names <- sub("\\.tif$", "", list.files(raster_dir, pattern = "\\.tif$"))
bios <- terra::rast(bios_files)
names(bios) <- bios_names
r <- bios[[1]]
# Build extent polygon from the raster, matching CRS of each vector layer
ext_poly_CL <- st_as_sfc(st_bbox(ext(r), crs = crs(r))) |> st_transform(st_crs(CL))
ext_poly_COAST <- st_as_sfc(st_bbox(ext(r), crs = crs(r))) |> st_transform(st_crs(COAST))

# Crop (st_intersection keeps only geometry inside the box, clipping cut lines)
CL_crop    <- st_intersection(st_make_valid(CL), ext_poly_CL)
COAST_crop <- st_intersection(st_make_valid(COAST), ext_poly_COAST)

# Union each cropped layer into a single geometry, then combine both
CL_union    <- st_union(CL_crop)
COAST_union <- st_union(COAST_crop)

combined <- st_union(CL_union, COAST_union)
plot(st_geometry(combined))

sf::write_sf(combined,"D:/PROGRAMAS/Dropbox/TESIS_JORGE/NaturalEarthData/coast_rivers.shp")
