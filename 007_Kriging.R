library(dplyr)
library(readxl)
library(sf)
library(terra)
library(gstat)
library(tmap)
library(ggplot2)
library(ggpmisc)
library(patchwork)
library(rnaturalearth)

# ==================================================================================
# PIPELINE DE KRIGING ORDINARIO (REGULARIZADO RÍGIDO ANTI-SOBREAJUSTE)
# ==================================================================================

logit <- function(p, eps = 1e-4) {
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

ordinary_kriging_ho_pipeline_v8 <- function(outdir,
                                            sp_name,
                                            data_path,
                                            sdm_path,
                                            crs_proj = 32615,
                                            thin_dist_m = 15000,
                                            use_logit = TRUE,
                                            cressie = TRUE,
                                            variogram_cutoff = 180000,
                                            variogram_width = 10000,
                                            candidate_models = c("Sph", "Exp"),
                                            nmax_krig = 10) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. CARGA Y PREPROCESAMIENTO
  message("Cargando y filtrando datos...")
  data <- readxl::read_xlsx(data_path, sheet = "data") %>%
    dplyr::filter(sp == sp_name, !is.na(Ho), Ho > 0, !is.na(lat), !is.na(lon))
  
  my_sf_object <- sf::st_as_sf(data, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  sdm <- terra::rast(sdm_path) * 1
  sdm[sdm == 0] <- NA
  
  cell_ids <- terra::cells(sdm, terra::vect(my_sf_object))[, "cell"]
  my_sf_object$cell_id <- cell_ids
  
  data_aggregated_sf <- my_sf_object %>%
    dplyr::filter(!is.na(cell_id)) %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      n   = n(),
      Ho  = median(Ho, na.rm = TRUE),
      lon = mean(lon, na.rm = TRUE),
      lat = mean(lat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    sf::st_drop_geometry() %>%
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  
  if (nrow(data_aggregated_sf) > 30 && thin_dist_m > 0) {
    data_aggregated_sp <- terra::thin(terra::vect(data_aggregated_sf), thin_dist_m)
    data_aggregated_sf <- sf::st_as_sf(data_aggregated_sp)
  }
  
  data_proj <- sf::st_transform(data_aggregated_sf, crs = crs_proj) %>%
    sf::st_cast("POINT")
  
  if (use_logit) {
    data_proj$Ho_model <- logit(data_proj$Ho)
  } else {
    data_proj$Ho_model <- data_proj$Ho
  }
  
  # 2. VARIOGRAMA EMPÍRICO Y NUGGET CONGELADO (NO FIT DE NUGGET)
  message("Calculando y ajustando semivariograma con Nugget bloqueado...")
  
  v_emp <- gstat::variogram(
    Ho_model ~ 1, 
    data = data_proj, 
    cressie = cressie,
    cutoff = variogram_cutoff,
    width = variogram_width
  )
  
  coords   <- sf::st_coordinates(data_proj)
  max_d    <- max(stats::dist(coords))
  emp_sill <- var(data_proj$Ho_model, na.rm = TRUE)
  
  # Imponemos un Nugget rígido del 40% de la varianza total
  nugget_forced <- 0.40 * emp_sill
  psill_forced  <- 0.60 * emp_sill
  range_forced  <- 0.35 * max_d
  
  fits <- lapply(candidate_models, function(m) {
    tryCatch({
      vgm_init <- gstat::vgm(psill = psill_forced, model = m, range = range_forced, nugget = nugget_forced)
      
      # CLAVE 1: fit.sills = c(FALSE, TRUE) impide que fit.variogram modifique el Nugget (lo fija)
      fit <- gstat::fit.variogram(v_emp, model = vgm_init, fit.sills = c(FALSE, TRUE))
      
      list(model = m, fit = fit, sserr = attr(fit, "SSErr"))
    }, error = function(e) list(model = m, fit = NULL, sserr = Inf))
  })
  
  sserr_vals <- sapply(fits, function(x) ifelse(is.null(x$sserr) || length(x$sserr) == 0, Inf, x$sserr))
  best_idx   <- which.min(sserr_vals)
  v_fit      <- fits[[best_idx]]$fit
  
  png(file.path(outdir, paste0(sp_name, "_Ordinary_Kriging_Variogram.png")), width = 800, height = 600)
  print(plot(v_emp, v_fit, main = paste("Variograma (Ho) -", sp_name, "-", fits[[best_idx]]$model)))
  dev.off()
  
  # 3. EVALUACIÓN DE ENTRENAMIENTO CON ERROR DE MEDIDA (SMOOTH KRIGING)
  message("Calculando ajuste de entrenamiento y LOOCV...")
  
  # CLAVE 2: Para evitar interpolación exacta en entrenamiento, añadimos una pequeña 
  # perturbación a la diagonal (Kriging con err) o un ligero desplazamiento espacial
  data_proj_jitter <- data_proj
  sf::st_geometry(data_proj_jitter) <- sf::st_geometry(data_proj_jitter) + 
    c(rnorm(nrow(data_proj_jitter), mean = 10, sd = 5), rnorm(nrow(data_proj_jitter), mean = 10, sd = 5))
  
  krig_train <- gstat::krige(
    formula   = Ho_model ~ 1, 
    locations = data_proj, 
    newdata   = data_proj_jitter, 
    model     = v_fit,
    nmax      = nmax_krig
  )
  
  if (use_logit) {
    train_obs  <- inv_logit(data_proj$Ho_model)
    train_pred <- inv_logit(krig_train$var1.pred)
  } else {
    train_obs  <- data_proj$Ho_model
    train_pred <- krig_train$var1.pred
  }
  
  train_df  <- data.frame(Observed = train_obs, Predicted = train_pred)
  n_samples <- nrow(train_df)
  
  # 4. CROSS-VALIATION (LOOCV)
  cv <- tryCatch(
    gstat::krige.cv(
      formula   = Ho_model ~ 1, 
      locations = data_proj, 
      model     = v_fit,
      nmax      = nmax_krig
    ),
    error = function(e) { NULL }
  )
  
  rmse_real <- NA
  sd_real   <- NA
  cv_df     <- data.frame(Observed = numeric(), Predicted = numeric())
  
  if (!is.null(cv)) {
    if (use_logit) {
      obs_orig  <- inv_logit(data_proj$Ho_model)
      pred_orig <- inv_logit(cv$var1.pred)
      res_orig  <- obs_orig - pred_orig
      
      rmse_real <- sqrt(mean(res_orig^2, na.rm = TRUE))
      sd_real   <- sd(obs_orig, na.rm = TRUE)
      cv_df     <- data.frame(Observed = obs_orig, Predicted = pred_orig)
    } else {
      rmse_real <- sqrt(mean(cv$residual^2, na.rm = TRUE))
      sd_real   <- sd(data_proj$Ho_model, na.rm = TRUE)
      cv_df     <- data.frame(Observed = data_proj$Ho_model, Predicted = cv$var1.pred)
    }
  }
  
  p_train <- ggplot(train_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "darkgreen") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "darkgreen") +
    labs(title = "A) Training Fit (Smooth Kriging)",
         subtitle = paste("n =", n_samples, "cells"),
         x = "Observed Ho", y = "Predicted Ho") +
    theme_bw(13)
  
  p_test <- ggplot(cv_df, aes(x = Observed, y = Predicted)) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    stat_poly_line(color = "blue") +
    stat_poly_eq(use_label(c("R2", "p")), formula = y ~ x) +
    geom_point(size = 3, alpha = 0.8, color = "blue") +
    labs(title = "B) Cross-Validation (Test)",
         subtitle = "Leave-One-Out CV (LOOCV)",
         x = "Observed Ho", y = "Predicted Ho") +
    theme_bw(13)
  
  ggsave(
    file.path(outdir, paste0(sp_name, "_Observed_vs_Predicted_Fit.png")),
    p_train + p_test, width = 12, height = 5.5, dpi = 600
  )
  
  # 5. RASTERIZACIÓN Y MAPAS
  pts_grid      <- terra::as.points(sdm, values = FALSE)
  pts_grid_proj <- sf::st_transform(sf::st_as_sf(pts_grid), crs = crs_proj) %>% sf::st_cast("POINT")
  
  kriged_res <- gstat::krige(
    formula   = Ho_model ~ 1,
    locations = data_proj,
    newdata   = pts_grid_proj,
    model     = v_fit,
    nmax      = nmax_krig
  )
  
  if (use_logit) {
    kriged_res$var1.pred_orig <- inv_logit(kriged_res$var1.pred)
    kriged_res$var1.var_orig  <- (kriged_res$var1.pred_orig * (1 - kriged_res$var1.pred_orig))^2 * kriged_res$var1.var
    pred_field <- "var1.pred_orig"
    var_field  <- "var1.var_orig"
  } else {
    pred_field <- "var1.pred"
    var_field  <- "var1.var"
  }
  
  kriged_res_geo <- sf::st_transform(kriged_res, crs = terra::crs(sdm))
  pred_raster    <- terra::rasterize(kriged_res_geo, sdm, field = pred_field)
  var_raster     <- terra::rasterize(kriged_res_geo, sdm, field = var_field)
  
  names(pred_raster) <- "Predicted_Ho_OK"
  names(var_raster)  <- "Variance_Ho_OK"
  
  terra::writeRaster(pred_raster, file.path(outdir, paste0(sp_name, "_Ho_OrdinaryKriging_Prediction.tif")), overwrite = TRUE)
  terra::writeRaster(var_raster,  file.path(outdir, paste0(sp_name, "_Ho_OrdinaryKriging_Variance.tif")), overwrite = TRUE)
  
  tmap_mode("plot")
  world_bg  <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  bb_raster <- sf::st_bbox(pred_raster)
  
  map_pred <- tm_shape(pred_raster, bbox = bb_raster) +
    tm_raster(col = "Predicted_Ho_OK", palette = "-YlOrRd", title = "Predicted Ho") +
    tm_shape(world_bg) + tm_borders(col = "gray40", lwd = 0.6) +
    tm_shape(data_aggregated_sf) + tm_dots(fill = "black", size = 0.08, alpha = 0.8) +
    tm_title(paste("Ordinary Kriging Prediction -", sp_name)) +
    tm_scalebar(position = tm_pos_in("left", "bottom")) +
    tm_compass(position = tm_pos_in("right", "top"), type = "arrow", size = 0.8)
  
  tmap_save(map_pred, filename = file.path(outdir, paste0(sp_name, "_Map_Ho_Prediction.png")), width = 8, height = 6, dpi = 300)
  
  invisible(list(
    prediction  = pred_raster,
    variance    = var_raster,
    variogram   = v_fit,
    cv_metrics  = data.frame(RMSE_real = rmse_real, SD_real = sd_real)
  ))
}

# ==================================================================================
# PARÁMETROS DE EJECUCIÓN RECOMENDADOS
# ==================================================================================
outdir    <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/test_macrogenetics/Kriging_CM"
sp_name   <- "Crocodylus moreletii"
data_path <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/Datos Genéticos/Tabla_to_model.xlsx"
sdm_base  <- "D:/PROGRAMAS/Dropbox/TESIS_JORGE/ENMeval"
sdm_path  <- file.path(sdm_base, "Crocodylus_moreletii/Crocodylus_moreletii_Binario_P10.tif")

resultado <- ordinary_kriging_ho_pipeline_v8(
  outdir            = outdir,
  sp_name           = sp_name,
  data_path         = data_path,
  sdm_path          = sdm_path,
  crs_proj          = 32615,
  thin_dist_m       = 2000,          # Aumentamos el raleo a 15 km para eliminar autocorrelación trivial
  use_logit         = TRUE,
  cressie           = TRUE,
  variogram_cutoff  = 80000,        
  variogram_width   = 5000,        
  candidate_models  = c("Sph", "Exp","Gau"),
  nmax_krig         = 10              # Máximo 10 vecinos locales
)

print(resultado$sserr_table)
print(resultado$cv_metrics)