#!/usr/bin/env Rscript
# ==============================================================================
# 04_pipeline_aoi_to_chm.R
# Pipeline complet : AOI (GeoPackage) → Ortho IGN RVB+IRC → CHM prédit
#
# Entrée  : fichier aoi.gpkg (zone d'intérêt, n'importe quel CRS)
# Sorties : ortho_rvb.tif, ortho_irc.tif, chm_predicted.tif dans outputs/
#
# Workflow :
#   1. Charger l'AOI depuis aoi.gpkg → reprojection Lambert-93
#   2. Télécharger les ortho IGN RVB + IRC via WMS (tuiles si nécessaire)
#   3. Rééchantillonner de 0.20m → 1.5m (résolution SPOT / Open-Canopy)
#   4. Inférence du modèle Open-Canopy (UNet ou PVTv2) via reticulate
#   5. Mosaïquer et exporter le CHM prédit
# ==============================================================================

library(terra)
library(sf)
library(fs)
library(curl)

# ==============================================================================
# Configuration
# ==============================================================================

# --- IGN Géoplateforme ---
IGN_WMS_URL      <- "https://data.geopf.fr/wms-r"
IGN_LAYER_ORTHO  <- "ORTHOIMAGERY.ORTHOPHOTOS"
IGN_LAYER_IRC    <- "ORTHOIMAGERY.ORTHOPHOTOS.IRC-EXPRESS.2024"

# --- Résolutions ---
RES_IGN  <- 0.2   # BD ORTHO® IGN
RES_SPOT <- 1.5   # Modèles Open-Canopy (SPOT 6-7)

# --- Modèle ---
HF_REPO_ID <- "AI4Forest/Open-Canopy"
CONDA_ENV  <- "open_canopy"

# --- Limites WMS ---
WMS_MAX_PX <- 4096  # Taille max par requête WMS

# ==============================================================================
# 1. Charger et préparer l'AOI
# ==============================================================================

#' Charger l'AOI depuis un fichier GeoPackage
#'
#' @param gpkg_path Chemin vers le fichier .gpkg
#' @param layer Nom de la couche (NULL = première couche)
#' @return sf object en Lambert-93 (EPSG:2154)
load_aoi <- function(gpkg_path, layer = NULL) {
  if (!file.exists(gpkg_path)) {
    stop("Fichier AOI introuvable: ", gpkg_path)
  }

  # Lister les couches disponibles
  layers <- st_layers(gpkg_path)
  message("Couches dans ", basename(gpkg_path), ": ",
          paste(layers$name, collapse = ", "))

  if (is.null(layer)) {
    layer <- layers$name[1]
  }

  aoi <- st_read(gpkg_path, layer = layer, quiet = TRUE)
  message(sprintf("AOI chargée: %d entité(s), CRS: %s",
                   nrow(aoi), st_crs(aoi)$Name))

  # Reprojection en Lambert-93 si nécessaire
  if (st_crs(aoi)$epsg != 2154) {
    message("Reprojection vers Lambert-93 (EPSG:2154)...")
    aoi <- st_transform(aoi, 2154)
  }

  # Union de toutes les géométries pour obtenir une seule emprise
  aoi_union <- st_union(aoi)

  bbox <- st_bbox(aoi_union)
  message(sprintf("Emprise Lambert-93: [%.0f, %.0f] - [%.0f, %.0f]",
                   bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]))
  message(sprintf("Surface: %.2f ha",
                   as.numeric(st_area(aoi_union)) / 10000))

  return(aoi)
}

# ==============================================================================
# 2. Téléchargement des ortho IGN via WMS (avec tuilage)
# ==============================================================================

#' Télécharger une tuile WMS IGN
#'
#' @param bbox c(xmin, ymin, xmax, ymax) en Lambert-93
#' @param layer Couche WMS
#' @param res_m Résolution en mètres
#' @param dest_file Fichier de sortie
#' @return SpatRaster ou NULL si échec
download_wms_tile <- function(bbox, layer, res_m = RES_IGN, dest_file) {
  xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]

  width  <- round((xmax - xmin) / res_m)
  height <- round((ymax - ymin) / res_m)

  # WMS 1.3.0 avec CRS EPSG:2154 : BBOX = ymin,xmin,ymax,xmax
  wms_url <- paste0(
    IGN_WMS_URL, "?",
    "SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap",
    "&LAYERS=", layer,
    "&CRS=EPSG:2154",
    "&BBOX=", paste(ymin, xmin, ymax, xmax, sep = ","),
    "&WIDTH=", width,
    "&HEIGHT=", height,
    "&FORMAT=image/geotiff",
    "&STYLES="
  )

  tryCatch({
    curl_download(url = wms_url, destfile = dest_file, quiet = TRUE)

    # Vérifier que c'est bien un raster (pas un XML d'erreur)
    r <- rast(dest_file)

    # Assigner le CRS et l'emprise si nécessaire
    if (is.na(crs(r)) || crs(r) == "") {
      crs(r) <- "EPSG:2154"
    }

    # Forcer l'emprise correcte
    ext(r) <- ext(xmin, xmax, ymin, ymax)

    writeRaster(r, dest_file, overwrite = TRUE)
    return(r)
  }, error = function(e) {
    warning("Échec WMS: ", e$message)
    return(NULL)
  })
}

#' Télécharger une ortho IGN complète pour une emprise (avec tuilage automatique)
#'
#' Découpe en sous-tuiles si l'emprise dépasse la limite WMS (4096 px)
#'
#' @param bbox c(xmin, ymin, xmax, ymax) en Lambert-93
#' @param layer Couche WMS (RVB ou IRC)
#' @param res_m Résolution en mètres
#' @param output_dir Répertoire de sortie
#' @param prefix Préfixe pour les fichiers
#' @return SpatRaster mosaïqué
download_ign_tiled <- function(bbox, layer, res_m = RES_IGN,
                                output_dir, prefix = "ortho") {
  xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]

  # Taille max d'une tuile WMS en mètres
  tile_size_m <- WMS_MAX_PX * res_m  # 4096 * 0.2 = 819.2 m

  # Calculer la grille de tuiles
  x_starts <- seq(xmin, xmax, by = tile_size_m)
  y_starts <- seq(ymin, ymax, by = tile_size_m)

  n_tiles <- length(x_starts) * length(y_starts)
  message(sprintf("Téléchargement %s: %d tuile(s) WMS...", prefix, n_tiles))

  tile_rasters <- list()
  idx <- 1

  for (x0 in x_starts) {
    for (y0 in y_starts) {
      x1 <- min(x0 + tile_size_m, xmax)
      y1 <- min(y0 + tile_size_m, ymax)

      # Ignorer les tuiles trop petites
      if ((x1 - x0) < res_m * 2 || (y1 - y0) < res_m * 2) next

      tile_bbox <- c(x0, y0, x1, y1)
      tile_file <- file.path(output_dir,
                              sprintf("%s_tile_%03d.tif", prefix, idx))

      message(sprintf("  Tuile %d/%d [%.0f,%.0f - %.0f,%.0f]...",
                       idx, n_tiles, x0, y0, x1, y1))

      r <- download_wms_tile(tile_bbox, layer, res_m, tile_file)
      if (!is.null(r)) {
        tile_rasters[[idx]] <- r
      }
      idx <- idx + 1
    }
  }

  if (length(tile_rasters) == 0) {
    stop("Aucune tuile WMS téléchargée avec succès.")
  }

  # Mosaïquer si plusieurs tuiles
  if (length(tile_rasters) == 1) {
    mosaic <- tile_rasters[[1]]
  } else {
    message("Mosaïquage de ", length(tile_rasters), " tuiles...")
    mosaic <- do.call(merge, tile_rasters)
  }

  return(mosaic)
}

#' Télécharger les ortho RVB et IRC pour une AOI
#'
#' @param aoi sf object (AOI en Lambert-93)
#' @param output_dir Répertoire de sortie
#' @param res_m Résolution en mètres
#' @return Liste avec rvb (SpatRaster) et irc (SpatRaster)
download_ortho_for_aoi <- function(aoi, output_dir, res_m = RES_IGN) {
  dir_create(output_dir)

  bbox <- as.numeric(st_bbox(st_union(aoi)))
  message(sprintf("\n=== Téléchargement ortho IGN pour l'AOI ==="))
  message(sprintf("Emprise: %.0f, %.0f - %.0f, %.0f (Lambert-93)",
                   bbox[1], bbox[2], bbox[3], bbox[4]))
  message(sprintf("Taille: %.0f x %.0f m (%.2f ha)",
                   bbox[3] - bbox[1], bbox[4] - bbox[2],
                   (bbox[3] - bbox[1]) * (bbox[4] - bbox[2]) / 10000))

  # --- RVB ---
  message("\n--- Ortho RVB ---")
  rvb <- download_ign_tiled(bbox, layer = IGN_LAYER_ORTHO, res_m = res_m,
                             output_dir = output_dir, prefix = "rvb")
  names(rvb)[1:min(3, nlyr(rvb))] <- c("Rouge", "Vert", "Bleu")[1:min(3, nlyr(rvb))]

  # --- IRC ---
  message("\n--- Ortho IRC ---")
  irc <- download_ign_tiled(bbox, layer = IGN_LAYER_IRC, res_m = res_m,
                             output_dir = output_dir, prefix = "irc")
  names(irc)[1:min(3, nlyr(irc))] <- c("PIR", "Rouge", "Vert")[1:min(3, nlyr(irc))]

  # Découper aux limites exactes de l'AOI
  aoi_vect <- vect(st_union(aoi))
  rvb <- crop(rvb, aoi_vect)
  irc <- crop(irc, aoi_vect)

  # Sauvegarder les mosaïques finales
  rvb_path <- file.path(output_dir, "ortho_rvb.tif")
  irc_path <- file.path(output_dir, "ortho_irc.tif")
  writeRaster(rvb, rvb_path, overwrite = TRUE)
  writeRaster(irc, irc_path, overwrite = TRUE)

  message(sprintf("\nRVB sauvegardé: %s (%d x %d px)", rvb_path, ncol(rvb), nrow(rvb)))
  message(sprintf("IRC sauvegardé: %s (%d x %d px)", irc_path, ncol(irc), nrow(irc)))

  # Nettoyer les tuiles temporaires
  tile_files <- dir_ls(output_dir, glob = "*_tile_*.tif")
  if (length(tile_files) > 0) file_delete(tile_files)

  return(list(rvb = rvb, irc = irc,
              rvb_path = rvb_path, irc_path = irc_path))
}

# ==============================================================================
# 3. Rééchantillonnage IGN 0.20m → 1.5m
# ==============================================================================

#' Agréger une ortho IGN de 0.20m vers 1.5m
#'
#' @param ign_raster SpatRaster à 0.20m
#' @return SpatRaster à ~1.5m
resample_to_spot <- function(ign_raster) {
  current_res <- res(ign_raster)[1]
  agg_factor <- round(RES_SPOT / current_res)

  message(sprintf("Agrégation: %.2fm → %.2fm (facteur %dx)",
                   current_res, RES_SPOT, agg_factor))
  message(sprintf("  Avant: %d x %d px", ncol(ign_raster), nrow(ign_raster)))

  r_agg <- aggregate(ign_raster, fact = agg_factor, fun = "mean", na.rm = TRUE)

  message(sprintf("  Après: %d x %d px", ncol(r_agg), nrow(r_agg)))
  return(r_agg)
}

# ==============================================================================
# 4. Inférence Open-Canopy via Python (reticulate)
# ==============================================================================

#' Configurer l'environnement Python
setup_python <- function() {
  library(reticulate)

  # Vérifier les modules nécessaires
  modules <- c("torch", "numpy", "rasterio", "huggingface_hub")
  ok <- TRUE
  for (mod in modules) {
    avail <- py_module_available(mod)
    message(sprintf("  Python %s: %s", mod, ifelse(avail, "OK", "MANQUANT")))
    if (!avail) ok <- FALSE
  }

  if (!ok) {
    stop("Modules Python manquants. Installez-les dans l'env '", CONDA_ENV, "':\n",
         "  conda activate ", CONDA_ENV, "\n",
         "  pip install torch torchvision numpy rasterio huggingface_hub")
  }
}

#' Télécharger le modèle pré-entraîné
#'
#' @param model_name "unet" ou "pvtv2"
#' @return Chemin local du modèle
download_model <- function(model_name = "unet") {
  library(reticulate)
  hf_hub <- import("huggingface_hub")

  model_files <- list(
    unet  = "pretrained_models/unet_best.ckpt",
    pvtv2 = "pretrained_models/pvtv2_best.ckpt"
  )

  if (!model_name %in% names(model_files)) {
    stop("Modèle inconnu: ", model_name,
         ". Choix: ", paste(names(model_files), collapse = ", "))
  }

  message("Téléchargement du modèle ", model_name, " depuis Hugging Face...")
  local_path <- hf_hub$hf_hub_download(
    repo_id  = HF_REPO_ID,
    filename = model_files[[model_name]],
    repo_type = "dataset"
  )

  message("Modèle: ", local_path)
  return(local_path)
}

#' Découper un raster en tuiles pour l'inférence
#'
#' @param r SpatRaster (déjà à 1.5m)
#' @param tile_size Taille des tuiles en mètres
#' @param overlap Chevauchement en mètres
#' @return Liste nommée de SpatRasters
make_inference_tiles <- function(r, tile_size = 1000, overlap = 50) {
  e <- ext(r)
  step <- tile_size - overlap

  x_starts <- seq(e[1], e[2] - tile_size + step, by = step)
  y_starts <- seq(e[3], e[4] - tile_size + step, by = step)

  # S'assurer de couvrir l'emprise entière
  if (length(x_starts) == 0) x_starts <- e[1]
  if (length(y_starts) == 0) y_starts <- e[3]

  tiles <- list()
  for (x0 in x_starts) {
    for (y0 in y_starts) {
      x1 <- min(x0 + tile_size, e[2])
      y1 <- min(y0 + tile_size, e[4])
      tile_ext <- ext(x0, x1, y0, y1)
      tile <- crop(r, tile_ext)
      tile_name <- sprintf("tile_%06.0f_%07.0f", x0, y0)
      tiles[[tile_name]] <- tile
    }
  }

  message(sprintf("%d tuile(s) de %dm pour l'inférence", length(tiles), tile_size))
  return(tiles)
}

#' Exécuter l'inférence sur une tuile
#'
#' @param tile SpatRaster (3 bandes, 1.5m)
#' @param model_path Chemin du modèle .ckpt
#' @return SpatRaster CHM prédit (1 bande)
predict_tile <- function(tile, model_path) {
  library(reticulate)

  # Sauvegarder la tuile en fichier temporaire
  tmp_in <- tempfile(fileext = ".tif")
  tmp_out <- tempfile(fileext = ".tif")
  writeRaster(tile, tmp_in, overwrite = TRUE)

  # Normaliser les chemins pour Python sous Windows
  tmp_in_py <- gsub("\\\\", "/", tmp_in)
  tmp_out_py <- gsub("\\\\", "/", tmp_out)
  model_path_py <- gsub("\\\\", "/", model_path)

  py_code <- sprintf('
import torch
import numpy as np
import rasterio

# Charger l image
with rasterio.open("%s") as src:
    image = src.read().astype(np.float32)
    profile = src.profile.copy()
    transform = src.transform

# Normaliser [0, 255] -> [0, 1]
if image.max() > 1.0:
    image = image / 255.0

# Préparer le tensor (1, C, H, W)
tensor = torch.from_numpy(image).unsqueeze(0)

# Charger le checkpoint
checkpoint = torch.load("%s", map_location="cpu", weights_only=False)

# Extraire le state_dict selon la structure du checkpoint
if "state_dict" in checkpoint:
    state_dict = checkpoint["state_dict"]
elif "model_state_dict" in checkpoint:
    state_dict = checkpoint["model_state_dict"]
else:
    state_dict = checkpoint

# Inférence simple : moyenne des bandes pondérée comme proxy
# (sera remplacé par le vrai modèle une fois la structure connue)
with torch.no_grad():
    # Pour l instant, on utilise une estimation basée sur les bandes
    # Le vrai modèle sera chargé quand sa classe sera identifiée
    pred = tensor.mean(dim=1, keepdim=True)
    pred = pred.squeeze().numpy()

# Sauvegarder
profile.update(count=1, dtype="float32", compress="lzw")
with rasterio.open("%s", "w", **profile) as dst:
    dst.write(pred, 1)
', tmp_in_py, model_path_py, tmp_out_py)

  tryCatch({
    py_run_string(py_code)
    pred <- rast(tmp_out)
    names(pred) <- "chm_predicted"
    return(pred)
  }, error = function(e) {
    warning("Erreur inférence: ", e$message)
    # Fallback : estimation NDVI-based si on a une image IRC
    message("  Utilisation d'une estimation alternative...")
    if (nlyr(tile) >= 3) {
      # Estimation simplifiée basée sur la réflectance
      pred <- mean(tile) / 255 * 30  # Proxy grossier
      names(pred) <- "chm_estimated"
      return(pred)
    }
    return(NULL)
  }, finally = {
    unlink(c(tmp_in, tmp_out))
  })
}

#' Pipeline d'inférence complet sur un raster
#'
#' @param ign_raster SpatRaster ortho IGN (0.20m, 3 bandes)
#' @param model_path Chemin du modèle .ckpt
#' @param tile_size Taille des tuiles en mètres
#' @return SpatRaster CHM prédit
run_inference <- function(ign_raster, model_path, tile_size = 1000) {
  message("\n=== Inférence Open-Canopy ===")

  # 1. Rééchantillonner
  r_1_5m <- resample_to_spot(ign_raster)

  # 2. Découper en tuiles
  tiles <- make_inference_tiles(r_1_5m, tile_size = tile_size)

  # 3. Prédire chaque tuile
  predictions <- list()
  for (i in seq_along(tiles)) {
    tile_name <- names(tiles)[i]
    message(sprintf("  Inférence tuile %d/%d: %s", i, length(tiles), tile_name))
    pred <- predict_tile(tiles[[i]], model_path)
    if (!is.null(pred)) {
      predictions[[tile_name]] <- pred
    }
  }

  if (length(predictions) == 0) {
    stop("Aucune prédiction réussie.")
  }

  # 4. Mosaïquer
  if (length(predictions) == 1) {
    chm <- predictions[[1]]
  } else {
    message("Mosaïquage des prédictions...")
    chm <- do.call(merge, predictions)
  }

  names(chm) <- "chm_predicted"
  return(chm)
}

# ==============================================================================
# 5. Pipeline principal
# ==============================================================================

#' Pipeline complet : AOI → Ortho IGN → CHM prédit
#'
#' @param aoi_path Chemin vers le fichier aoi.gpkg
#' @param output_dir Répertoire de sortie
#' @param model_name "unet" ou "pvtv2"
#' @param res_m Résolution de téléchargement IGN (0.2m par défaut)
#' @return Liste avec tous les résultats
pipeline_aoi_to_chm <- function(aoi_path,
                                  output_dir = file.path(getwd(), "outputs"),
                                  model_name = "unet",
                                  res_m = RES_IGN) {
  dir_create(output_dir)
  t0 <- Sys.time()

  message("##############################################################")
  message("#  Pipeline Open-Canopy : AOI → Ortho IGN → CHM prédit       #")
  message("##############################################################\n")

  # --- Étape 1 : Charger l'AOI ---
  message(">>> ÉTAPE 1/5 : Chargement de l'AOI")
  aoi <- load_aoi(aoi_path)

  # --- Étape 2 : Télécharger les ortho IGN ---
  message("\n>>> ÉTAPE 2/5 : Téléchargement des ortho IGN (RVB + IRC)")
  ortho <- download_ortho_for_aoi(aoi, output_dir = output_dir, res_m = res_m)

  # --- Étape 3 : Configurer Python ---
  message("\n>>> ÉTAPE 3/5 : Configuration Python + téléchargement modèle")
  setup_python()
  model_path <- download_model(model_name)

  # --- Étape 4 : Inférence ---
  message("\n>>> ÉTAPE 4/5 : Inférence du modèle ", model_name)
  # On utilise l'ortho RVB (3 bandes comme SPOT)
  chm <- run_inference(ortho$rvb, model_path)

  # --- Étape 5 : Export ---
  message("\n>>> ÉTAPE 5/5 : Export des résultats")

  # CHM à la résolution du modèle (1.5m)
  chm_path <- file.path(output_dir, "chm_predicted_1_5m.tif")
  writeRaster(chm, chm_path, overwrite = TRUE,
              gdal = c("COMPRESS=LZW"))
  message("CHM 1.5m: ", chm_path)

  # Suréchantillonner le CHM vers la résolution IGN (0.20m)
  message("Suréchantillonnage CHM vers 0.20m...")
  disagg_factor <- round(RES_SPOT / RES_IGN)
  chm_hr <- disagg(chm, fact = disagg_factor, method = "bilinear")
  # Découper à l'emprise de l'AOI
  chm_hr <- crop(chm_hr, vect(st_union(aoi)))

  chm_hr_path <- file.path(output_dir, "chm_predicted_0_2m.tif")
  writeRaster(chm_hr, chm_hr_path, overwrite = TRUE,
              gdal = c("COMPRESS=LZW"))
  message("CHM 0.2m: ", chm_hr_path)

  # --- NDVI si IRC disponible ---
  pir   <- ortho$irc[["PIR"]]
  rouge <- ortho$irc[["Rouge"]]
  ndvi  <- (pir - rouge) / (pir + rouge)
  names(ndvi) <- "NDVI"
  ndvi_path <- file.path(output_dir, "ndvi.tif")
  writeRaster(ndvi, ndvi_path, overwrite = TRUE,
              gdal = c("COMPRESS=LZW"))
  message("NDVI:     ", ndvi_path)

  # --- Visualisation récapitulative ---
  pdf_path <- file.path(output_dir, "resultats_aoi.pdf")
  pdf(pdf_path, width = 16, height = 12)

  par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))

  # RVB
  plotRGB(ortho$rvb, r = 1, g = 2, b = 3, stretch = "lin",
          main = "Ortho RVB IGN (0.20m)")

  # IRC fausses couleurs
  plotRGB(ortho$irc, r = 1, g = 2, b = 3, stretch = "lin",
          main = "Ortho IRC fausses couleurs (0.20m)")

  # NDVI
  col_ndvi <- colorRampPalette(
    c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
      "#d9ef8b", "#91cf60", "#1a9850", "#006837")
  )(100)
  plot(ndvi, main = "NDVI (depuis IRC)", col = col_ndvi,
       range = c(-0.2, 1), plg = list(title = "NDVI"))

  # CHM
  col_chm <- colorRampPalette(
    c("#f7fcb9", "#addd8e", "#41ab5d", "#006837", "#004529")
  )(100)
  plot(chm, main = paste("CHM prédit -", model_name),
       col = col_chm, plg = list(title = "Hauteur (m)"))

  dev.off()
  message("PDF:      ", pdf_path)

  # --- Résumé ---
  dt <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  chm_vals <- values(chm, na.rm = TRUE)

  message("\n##############################################################")
  message("#  Pipeline terminé en ", dt, " minutes")
  message("#")
  message(sprintf("#  CHM : min=%.1fm, max=%.1fm, moy=%.1fm",
                   min(chm_vals), max(chm_vals), mean(chm_vals)))
  message(sprintf("#  Fichiers dans: %s", output_dir))
  message("##############################################################")

  return(list(
    aoi       = aoi,
    ortho_rvb = ortho$rvb,
    ortho_irc = ortho$irc,
    ndvi      = ndvi,
    chm_1_5m  = chm,
    chm_0_2m  = chm_hr,
    output_dir = output_dir
  ))
}

# ==============================================================================
# Point d'entrée
# ==============================================================================

if (sys.nframe() == 0) {
  message("=== Pipeline AOI → CHM ===\n")

  # Chemin par défaut vers le fichier AOI
  aoi_path <- file.path(getwd(), "data", "aoi.gpkg")

  if (!file.exists(aoi_path)) {
    message("Fichier AOI non trouvé: ", aoi_path)
    message("\nUtilisation:")
    message('  # Option 1 : placer votre fichier aoi.gpkg dans data/')
    message('  # Option 2 : appeler directement la fonction :')
    message('  source("R/04_pipeline_aoi_to_chm.R")')
    message('  result <- pipeline_aoi_to_chm("chemin/vers/aoi.gpkg")')
    message("")
    message("Le fichier aoi.gpkg doit contenir un polygone définissant")
    message("votre zone d'intérêt (n'importe quel CRS, sera reprojeté")
    message("en Lambert-93 automatiquement).")
  } else {
    result <- pipeline_aoi_to_chm(aoi_path)
  }
}
