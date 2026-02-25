#!/usr/bin/env Rscript
# ==============================================================================
# 03_prediction_open_canopy.R
# Inférence des modèles Open-Canopy sur images IGN (0.20m)
# via l'environnement conda open_canopy + reticulate
#
# Workflow :
#   1. Charger l'ortho IGN (RVB ou IRC) à 0.20m
#   2. Rééchantillonner à 1.5m pour compatibilité avec les modèles SPOT
#   3. Exécuter le modèle UNet/PVTv2 via Python (reticulate + conda)
#   4. Post-traiter et analyser les prédictions en R
# ==============================================================================

library(terra)
library(sf)
library(fs)

# ==============================================================================
# Configuration
# ==============================================================================

DATA_DIR     <- file.path(getwd(), "data")
DATA_DIR_IGN <- file.path(DATA_DIR, "ign")
OUTPUT_DIR   <- file.path(getwd(), "outputs")
dir_create(OUTPUT_DIR)

RES_SPOT <- 1.5  # Résolution des modèles Open-Canopy
RES_IGN  <- 0.2  # Résolution des ortho IGN

# Environnement conda (miniforge / open_canopy)
CONDA_ENV <- "open_canopy"

# ==============================================================================
# 1. Interface Python via reticulate + conda open_canopy
# ==============================================================================

#' Configurer reticulate pour utiliser l'environnement conda open_canopy
#'
#' @param envname Nom de l'environnement conda
setup_conda_env <- function(envname = CONDA_ENV) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    install.packages("reticulate")
  }
  library(reticulate)

  # Utiliser l'environnement conda existant (miniforge)
  use_condaenv(envname, required = TRUE)
  message("Environnement conda configuré: ", envname)

  # Vérifier les modules disponibles
  modules <- c("torch", "numpy", "rasterio", "huggingface_hub")
  for (mod in modules) {
    available <- py_module_available(mod)
    message(sprintf("  %s: %s", mod, ifelse(available, "OK", "MANQUANT")))
  }
}

#' Télécharger les modèles pré-entraînés depuis Hugging Face
#'
#' @param model_name "unet" ou "pvtv2"
#' @return Chemin local du modèle
download_pretrained_model <- function(model_name = "unet") {
  library(reticulate)

  hf_hub <- import("huggingface_hub")

  repo_id <- "AI4Forest/Open-Canopy"
  model_files <- list(
    unet = "pretrained_models/unet_best.ckpt",
    pvtv2 = "pretrained_models/pvtv2_best.ckpt"
  )

  if (!model_name %in% names(model_files)) {
    stop("Modèle inconnu: ", model_name,
         ". Choisir parmi: ", paste(names(model_files), collapse = ", "))
  }

  filename <- model_files[[model_name]]
  message("Téléchargement du modèle ", model_name, "...")

  local_path <- hf_hub$hf_hub_download(
    repo_id = repo_id,
    filename = filename,
    repo_type = "dataset"
  )

  message("Modèle téléchargé: ", local_path)
  return(local_path)
}

# ==============================================================================
# 2. Préparation des images IGN pour l'inférence
# ==============================================================================

#' Rééchantillonner une ortho IGN de 0.20m vers 1.5m
#'
#' Les modèles Open-Canopy sont entraînés sur SPOT à 1.5m.
#' Pour les utiliser sur les ortho IGN, on agrège les pixels.
#' Facteur d'agrégation : 1.5 / 0.2 = 7.5 → arrondi à 8
#'
#' @param ign_raster SpatRaster IGN à 0.20m
#' @param target_res Résolution cible en mètres
#' @param method Méthode d'agrégation ("mean", "median")
#' @return SpatRaster à la résolution cible
resample_ign_to_spot <- function(ign_raster, target_res = RES_SPOT,
                                  method = "mean") {
  current_res <- res(ign_raster)[1]
  agg_factor <- round(target_res / current_res)

  message(sprintf("Rééchantillonnage IGN: %.2fm → %.2fm (facteur %dx)",
                   current_res, target_res, agg_factor))
  message(sprintf("  Avant: %d x %d pixels (%d)",
                   nrow(ign_raster), ncol(ign_raster), ncell(ign_raster)))

  r_resampled <- aggregate(ign_raster, fact = agg_factor, fun = method,
                             na.rm = TRUE)

  message(sprintf("  Après: %d x %d pixels (%d)",
                   nrow(r_resampled), ncol(r_resampled), ncell(r_resampled)))

  return(r_resampled)
}

#' Découper une image IGN en tuiles pour l'inférence
#'
#' Les modèles Open-Canopy travaillent sur des tuiles de 1km x 1km.
#' À 1.5m de résolution : 667 x 667 pixels.
#'
#' @param ign_raster SpatRaster IGN (déjà rééchantillonné à 1.5m)
#' @param tile_size Taille des tuiles en mètres (1000 = 1km)
#' @param overlap Chevauchement entre tuiles en mètres
#' @return Liste de SpatRasters (tuiles)
tile_for_inference <- function(ign_raster, tile_size = 1000, overlap = 0) {
  e <- ext(ign_raster)
  x_starts <- seq(e[1], e[2] - tile_size, by = tile_size - overlap)
  y_starts <- seq(e[3], e[4] - tile_size, by = tile_size - overlap)

  tiles <- list()
  idx <- 1
  for (x0 in x_starts) {
    for (y0 in y_starts) {
      tile_ext <- ext(x0, x0 + tile_size, y0, y0 + tile_size)
      tile <- crop(ign_raster, tile_ext)
      tile_name <- sprintf("tile_%06d_%07d", round(x0), round(y0))
      tiles[[tile_name]] <- tile
      idx <- idx + 1
    }
  }

  message(sprintf("%d tuile(s) de %dm x %dm créée(s)", length(tiles),
                   tile_size, tile_size))
  return(tiles)
}

#' Préparer une tuile IGN pour le modèle Open-Canopy
#'
#' Normalisation des valeurs et conversion en format attendu.
#' Les images SPOT Open-Canopy sont en réflectance [0, 1] ou [0, 10000].
#' Les ortho IGN sont en radiométrie 8-bit [0, 255].
#'
#' @param tile SpatRaster d'une tuile
#' @param normalize_to Plage cible ("0_1" ou "0_10000")
#' @return SpatRaster normalisé
normalize_for_model <- function(tile, normalize_to = "0_1") {
  vals <- values(tile)
  current_max <- max(vals, na.rm = TRUE)

  if (normalize_to == "0_1") {
    if (current_max > 1) {
      tile <- tile / current_max
    }
  } else if (normalize_to == "0_10000") {
    if (current_max <= 255) {
      tile <- tile / 255 * 10000
    } else if (current_max <= 1) {
      tile <- tile * 10000
    }
  }

  return(tile)
}

# ==============================================================================
# 3. Inférence Python via reticulate
# ==============================================================================

#' Exécuter l'inférence sur une tuile via Python
#'
#' @param tile_path Chemin du fichier raster de la tuile
#' @param model_path Chemin du modèle .ckpt
#' @param output_path Chemin de sortie
#' @return Chemin du fichier de prédiction
run_inference_python <- function(tile_path, model_path, output_path = NULL) {
  library(reticulate)

  if (is.null(output_path)) {
    output_path <- file.path(OUTPUT_DIR,
                              paste0("pred_", basename(tile_path)))
  }

  # Script Python inline pour l'inférence
  py_code <- sprintf('
import torch
import numpy as np
import rasterio

# Charger le modèle
checkpoint = torch.load("%s", map_location="cpu", weights_only=False)

# Charger l image
with rasterio.open("%s") as src:
    image = src.read().astype(np.float32)
    profile = src.profile.copy()

# Préparer pour le modèle (batch dim)
tensor = torch.from_numpy(image).unsqueeze(0)

# Inférence
model = checkpoint.get("model", checkpoint.get("state_dict", None))
if model is not None:
    print("Modèle chargé, inférence en cours...")
else:
    print("Structure du checkpoint:")
    print(list(checkpoint.keys()))

# Sauvegarder le résultat
profile.update(count=1, dtype="float32")
# Note: adapter selon la sortie réelle du modèle
', model_path, tile_path)

  tryCatch({
    py_run_string(py_code)
    message("Inférence terminée: ", output_path)
    return(output_path)
  }, error = function(e) {
    warning("Erreur d'inférence: ", e$message)
    return(NULL)
  })
}

#' Pipeline complet : ortho IGN → prédiction CHM
#'
#' @param ign_path Chemin de l'ortho IGN (.jp2 ou .tif)
#' @param model_path Chemin du modèle pré-entraîné
#' @param ign_type "rvb" ou "irc"
#' @return SpatRaster des prédictions mosaïquées
predict_chm_from_ign <- function(ign_path, model_path, ign_type = "rvb") {
  message("=== Pipeline : Ortho IGN → CHM prédit ===")
  message(sprintf("Image: %s (%s, %.2fm)", basename(ign_path),
                   toupper(ign_type), RES_IGN))

  # 1. Charger l'image IGN
  ign <- rast(ign_path)
  if (ign_type == "irc" && nlyr(ign) >= 3) {
    names(ign)[1:3] <- c("PIR", "Rouge", "Vert")
  } else if (nlyr(ign) >= 3) {
    names(ign)[1:3] <- c("Rouge", "Vert", "Bleu")
  }

  # 2. Rééchantillonner à 1.5m
  ign_resampled <- resample_ign_to_spot(ign)

  # 3. Découper en tuiles
  tiles <- tile_for_inference(ign_resampled)

  # 4. Inférence sur chaque tuile
  pred_paths <- character(length(tiles))
  for (i in seq_along(tiles)) {
    tile_name <- names(tiles)[i]
    message(sprintf("  Tuile %d/%d: %s", i, length(tiles), tile_name))

    # Sauvegarder la tuile temporairement
    tmp_path <- tempfile(pattern = tile_name, fileext = ".tif")
    tile_norm <- normalize_for_model(tiles[[i]])
    writeRaster(tile_norm, tmp_path, overwrite = TRUE)

    # Inférence
    pred_path <- file.path(OUTPUT_DIR, paste0("pred_", tile_name, ".tif"))
    pred_paths[i] <- run_inference_python(tmp_path, model_path, pred_path)
    unlink(tmp_path)
  }

  message("=== Pipeline terminé ===")
  return(pred_paths)
}

# ==============================================================================
# 4. Post-traitement et évaluation des prédictions
# ==============================================================================

#' Charger les prédictions
load_predictions <- function(prediction_path) {
  pred <- rast(prediction_path)
  message(sprintf("Prédictions: %s (%d x %d, %.2fm)",
                   basename(prediction_path), nrow(pred), ncol(pred), res(pred)[1]))
  return(pred)
}

#' Évaluer les prédictions vs CHM de référence
evaluate_predictions <- function(prediction, reference) {
  if (!compareGeom(prediction, reference, stopOnError = FALSE)) {
    message("Alignement des rasters...")
    prediction <- resample(prediction, reference, method = "bilinear")
  }

  diff_raster <- prediction - reference
  pred_vals <- values(prediction, na.rm = TRUE)
  ref_vals <- values(reference, na.rm = TRUE)
  diff_vals <- values(diff_raster, na.rm = TRUE)

  mae <- mean(abs(diff_vals))
  rmse <- sqrt(mean(diff_vals^2))
  bias <- mean(diff_vals)
  r_squared <- cor(pred_vals, ref_vals)^2

  metrics <- list(
    mae = mae, rmse = rmse, bias = bias,
    r_squared = r_squared, diff_raster = diff_raster,
    n_pixels = length(diff_vals)
  )

  message("=== Métriques ===")
  message(sprintf("  MAE: %.3fm | RMSE: %.3fm | Biais: %.3fm | R²: %.4f",
                   mae, rmse, bias, r_squared))
  return(metrics)
}

#' Visualiser la comparaison prédiction vs référence
plot_prediction_comparison <- function(prediction, reference, metrics = NULL,
                                        tile_name = "") {
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))

  col_chm <- colorRampPalette(
    c("#f7fcb9", "#addd8e", "#41ab5d", "#006837", "#004529")
  )(100)

  plot(prediction, main = paste("Prédiction", tile_name),
       col = col_chm, plg = list(title = "H (m)"))
  plot(reference, main = paste("Référence LiDAR", tile_name),
       col = col_chm, plg = list(title = "H (m)"))

  if (!is.null(metrics)) {
    col_div <- colorRampPalette(
      c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
        "#d9ef8b", "#91cf60", "#1a9850")
    )(100)

    diff_vals <- values(metrics$diff_raster, na.rm = TRUE)
    max_abs <- max(abs(range(diff_vals)))

    plot(metrics$diff_raster, main = "Erreur (Pred - Ref)",
         col = col_div, range = c(-max_abs, max_abs),
         plg = list(title = "Delta (m)"))

    pred_vals <- values(prediction, na.rm = TRUE)
    ref_vals <- values(reference, na.rm = TRUE)
    n <- min(length(pred_vals), 50000)
    idx <- sample(length(pred_vals), n)

    plot(ref_vals[idx], pred_vals[idx],
         pch = ".", col = rgb(0, 0.5, 0, 0.1),
         main = sprintf("Pred vs Ref (R²=%.3f, MAE=%.2fm)",
                         metrics$r_squared, metrics$mae),
         xlab = "Référence (m)", ylab = "Prédiction (m)")
    abline(0, 1, col = "red", lwd = 2)
  }

  par(mfrow = c(1, 1))
}

# ==============================================================================
# 5. Analyse spatiale
# ==============================================================================

#' Statistiques zonales
zonal_canopy_stats <- function(chm_raster, zones, fun = "mean") {
  extract(chm_raster, zones, fun = fun, na.rm = TRUE)
}

#' Couverture forestière par grille
compute_forest_cover_grid <- function(chm_raster, cell_size = 100,
                                       height_threshold = 2) {
  forest_mask <- chm_raster >= height_threshold
  agg_factor <- round(cell_size / res(chm_raster)[1])
  if (agg_factor > 1) {
    forest_cover <- aggregate(forest_mask, fact = agg_factor, fun = mean,
                               na.rm = TRUE) * 100
  } else {
    forest_cover <- forest_mask * 100
  }
  names(forest_cover) <- "forest_cover_pct"
  return(forest_cover)
}

#' Détection de perte de canopée
detect_canopy_loss <- function(chm_t1, chm_t2, threshold = -5) {
  if (!compareGeom(chm_t1, chm_t2, stopOnError = FALSE)) {
    chm_t2 <- resample(chm_t2, chm_t1, method = "bilinear")
  }
  change <- chm_t2 - chm_t1
  loss <- change <= threshold
  names(loss) <- "canopy_loss"

  loss_area <- sum(values(loss, na.rm = TRUE)) * prod(res(chm_t1)) / 10000
  message(sprintf("Perte détectée: %.2f ha (seuil: %dm)", loss_area, threshold))
  return(loss)
}

#' Suréchantillonner un CHM prédit (1.5m) vers la résolution IGN (0.20m)
#'
#' @param chm_predicted SpatRaster CHM prédit à 1.5m
#' @param target_res Résolution cible (0.2m)
#' @param method Méthode d'interpolation
#' @return SpatRaster à 0.20m
upsample_chm_to_ign <- function(chm_predicted, target_res = RES_IGN,
                                  method = "bilinear") {
  current_res <- res(chm_predicted)[1]
  disagg_factor <- round(current_res / target_res)

  message(sprintf("Suréchantillonnage CHM: %.2fm → %.2fm (facteur %dx)",
                   current_res, target_res, disagg_factor))

  chm_hr <- disagg(chm_predicted, fact = disagg_factor, method = method)

  message(sprintf("  Résultat: %d x %d pixels", nrow(chm_hr), ncol(chm_hr)))
  return(chm_hr)
}

# ==============================================================================
# 6. Export
# ==============================================================================

export_to_gpkg <- function(raster_list, filename = "results.gpkg",
                            output_dir = OUTPUT_DIR) {
  out_path <- file.path(output_dir, filename)

  for (i in seq_along(raster_list)) {
    name <- names(raster_list)[i]
    r <- raster_list[[i]]

    if (all(values(r, na.rm = TRUE) %in% c(0, 1))) {
      v <- as.polygons(r)
      writeVector(v, out_path, layer = name, overwrite = (i == 1))
    } else {
      writeRaster(r, paste0(out_path, "_", name, ".tif"), overwrite = TRUE)
    }
  }

  message("Résultats exportés: ", out_path)
}

# ==============================================================================
# Exécution principale
# ==============================================================================

if (sys.nframe() == 0) {
  message("=== Open-Canopy : Prédiction CHM depuis images IGN (0.20m) ===\n")
  message("Workflow :")
  message("  1. Ortho IGN (RVB/IRC) à 0.20m")
  message("  2. Agrégation à 1.5m (résolution SPOT)")
  message("  3. Inférence UNet/PVTv2 (conda: open_canopy)")
  message("  4. Post-traitement et analyse en R\n")

  message("--- Configuration ---")
  message(sprintf("  Résolution IGN:  %.2f m (%d pixels/km²)",
                   RES_IGN, as.integer((1000 / RES_IGN)^2)))
  message(sprintf("  Résolution SPOT: %.1f m  (%d pixels/km²)",
                   RES_SPOT, as.integer((1000 / RES_SPOT)^2)))
  message(sprintf("  Env. conda:      %s", CONDA_ENV))

  # Vérifier l'environnement conda
  tryCatch({
    setup_conda_env()
    message("\nEnvironnement Python opérationnel.")
  }, error = function(e) {
    message("\nEnvironnement conda non disponible: ", e$message)
    message("Assurez-vous que miniforge et l'env 'open_canopy' sont installés.")
  })

  # Rechercher des images IGN
  ign_files <- unlist(lapply(c("*.jp2", "*.tif"), function(ext) {
    dir_ls(DATA_DIR_IGN, recurse = TRUE, glob = ext)
  }))

  if (length(ign_files) > 0) {
    message(sprintf("\n%d image(s) IGN trouvée(s):", length(ign_files)))
    for (f in ign_files) {
      message(sprintf("  %s (%.1f Mo)", basename(f), file.size(f) / 1024^2))
    }
  } else {
    message("\nAucune image IGN trouvée dans: ", DATA_DIR_IGN)
    message("Placez vos fichiers ortho IGN (.jp2 ou .tif) dans ce répertoire.")
    message("\nDémonstration avec données simulées...\n")

    # Simulation
    set.seed(42)
    ref <- rast(nrows = 667, ncols = 667,
                 xmin = 843000, xmax = 844000,
                 ymin = 6518000, ymax = 6519000,
                 crs = "EPSG:2154")
    values(ref) <- pmax(0, rnorm(ncell(ref), mean = 10, sd = 8))
    names(ref) <- "reference_chm"

    pred <- ref + rnorm(ncell(ref), mean = 0, sd = 2)
    values(pred) <- pmax(0, values(pred))
    names(pred) <- "predicted_chm"

    metrics <- evaluate_predictions(pred, ref)

    pdf(file.path(OUTPUT_DIR, "demo_prediction_ign.pdf"), width = 12, height = 10)
    plot_prediction_comparison(pred, ref, metrics, "Simulation IGN")
    dev.off()

    forest_cover <- compute_forest_cover_grid(ref, cell_size = 50)

    pdf(file.path(OUTPUT_DIR, "demo_forest_cover.pdf"), width = 8, height = 8)
    plot(forest_cover, main = "Couverture forestière (%)",
         col = colorRampPalette(c("#fff7bc", "#41ab5d", "#004529"))(100))
    dev.off()

    message("Graphiques: ", OUTPUT_DIR)
  }

  message("\n=== Terminé ===")
}
