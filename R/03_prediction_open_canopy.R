#!/usr/bin/env Rscript
# ==============================================================================
# 03_prediction_open_canopy.R
# Utilisation des modèles pré-entraînés Open-Canopy via Python (reticulate)
# et post-traitement des prédictions en R
# ==============================================================================

library(terra)
library(sf)
library(fs)

# ==============================================================================
# Configuration
# ==============================================================================

DATA_DIR <- file.path(getwd(), "data")
OUTPUT_DIR <- file.path(getwd(), "outputs")
dir_create(OUTPUT_DIR)

# ==============================================================================
# 1. Interface avec les modèles Python via reticulate
# ==============================================================================

#' Configurer l'environnement Python pour Open-Canopy
#'
#' Installe les dépendances Python nécessaires via reticulate
#' @param envname Nom de l'environnement conda/virtualenv
setup_python_env <- function(envname = "open-canopy") {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    install.packages("reticulate")
  }
  library(reticulate)

  # Créer un environnement virtuel si nécessaire
  if (!virtualenv_exists(envname)) {
    message("Création de l'environnement Python...")
    virtualenv_create(envname)
    virtualenv_install(envname, packages = c(
      "torch", "torchvision", "numpy", "rasterio",
      "huggingface_hub", "safetensors"
    ))
  }

  use_virtualenv(envname, required = TRUE)
  message("Environnement Python configuré: ", envname)
}

#' Télécharger les modèles pré-entraînés depuis Hugging Face
#'
#' @param model_name Nom du modèle: "unet" ou "pvtv2"
#' @param dest_dir Répertoire de destination
#' @return Chemin local du modèle
download_pretrained_model <- function(model_name = "unet",
                                       dest_dir = DATA_DIR) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Le package reticulate est nécessaire. ",
         "Installez-le avec: install.packages('reticulate')")
  }
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
# 2. Post-traitement des prédictions en R
# ==============================================================================

#' Charger les prédictions de hauteur de canopée
#'
#' @param prediction_path Chemin vers le raster de prédiction (.tif)
#' @return SpatRaster
load_predictions <- function(prediction_path) {
  pred <- rast(prediction_path)
  message(sprintf("Prédictions chargées: %s", basename(prediction_path)))
  message(sprintf("  Dimensions: %d x %d", nrow(pred), ncol(pred)))
  message(sprintf("  Résolution: %.2f m", res(pred)[1]))
  return(pred)
}

#' Comparer les prédictions avec les CHM de référence (ground truth)
#'
#' @param prediction SpatRaster des prédictions
#' @param reference SpatRaster du CHM de référence (LiDAR)
#' @return Liste avec les métriques d'évaluation
evaluate_predictions <- function(prediction, reference) {
  # Aligner les rasters si nécessaire
  if (!compareGeom(prediction, reference, stopOnError = FALSE)) {
    message("Alignement des rasters...")
    prediction <- resample(prediction, reference, method = "bilinear")
  }

  # Calculer les différences
  diff_raster <- prediction - reference

  pred_vals <- values(prediction, na.rm = TRUE)
  ref_vals <- values(reference, na.rm = TRUE)
  diff_vals <- values(diff_raster, na.rm = TRUE)

  # Métriques
  mae <- mean(abs(diff_vals))
  rmse <- sqrt(mean(diff_vals^2))
  bias <- mean(diff_vals)
  r_squared <- cor(pred_vals, ref_vals)^2

  metrics <- list(
    mae = mae,
    rmse = rmse,
    bias = bias,
    r_squared = r_squared,
    diff_raster = diff_raster,
    n_pixels = length(diff_vals)
  )

  message(sprintf("=== Métriques d'évaluation ==="))
  message(sprintf("  MAE:  %.3f m", mae))
  message(sprintf("  RMSE: %.3f m", rmse))
  message(sprintf("  Biais: %.3f m", bias))
  message(sprintf("  R²:   %.4f", r_squared))
  message(sprintf("  Pixels: %d", length(diff_vals)))

  return(metrics)
}

#' Visualiser la comparaison prédiction vs référence
#'
#' @param prediction SpatRaster des prédictions
#' @param reference SpatRaster de référence
#' @param metrics Liste de métriques (sortie de evaluate_predictions)
#' @param tile_name Nom de la tuile
plot_prediction_comparison <- function(prediction, reference, metrics = NULL,
                                        tile_name = "") {
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))

  # Palette commune
  col_chm <- colorRampPalette(
    c("#f7fcb9", "#addd8e", "#41ab5d", "#006837", "#004529")
  )(100)

  # Prédiction
  plot(prediction, main = paste("Prédiction", tile_name),
       col = col_chm, plg = list(title = "H (m)"))

  # Référence
  plot(reference, main = paste("Référence LiDAR", tile_name),
       col = col_chm, plg = list(title = "H (m)"))

  # Erreur
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

    # Scatterplot
    pred_vals <- values(prediction, na.rm = TRUE)
    ref_vals <- values(reference, na.rm = TRUE)

    # Sous-échantillonner pour le graphique
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
# 3. Agrégation et analyse spatiale
# ==============================================================================

#' Calculer des statistiques zonales de hauteur de canopée
#'
#' @param chm_raster SpatRaster du CHM
#' @param zones SpatVector des zones (polygones)
#' @param fun Fonction d'agrégation
#' @return data.frame avec les statistiques par zone
zonal_canopy_stats <- function(chm_raster, zones, fun = "mean") {
  stats <- extract(chm_raster, zones, fun = fun, na.rm = TRUE)
  return(stats)
}

#' Calculer la couverture forestière par cellule de grille
#'
#' @param chm_raster SpatRaster du CHM
#' @param cell_size Taille des cellules de grille (en mètres)
#' @param height_threshold Seuil de hauteur pour considérer un arbre (m)
#' @return SpatRaster de couverture forestière (%)
compute_forest_cover_grid <- function(chm_raster, cell_size = 100,
                                       height_threshold = 2) {
  # Créer un masque binaire forêt/non-forêt
  forest_mask <- chm_raster >= height_threshold

  # Agréger à la résolution de la grille
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

#' Détecter les zones de perte de canopée
#'
#' @param chm_t1 CHM à la date 1
#' @param chm_t2 CHM à la date 2
#' @param threshold Seuil de perte en mètres
#' @return SpatRaster binaire (1 = perte détectée)
detect_canopy_loss <- function(chm_t1, chm_t2, threshold = -5) {
  if (!compareGeom(chm_t1, chm_t2, stopOnError = FALSE)) {
    chm_t2 <- resample(chm_t2, chm_t1, method = "bilinear")
  }

  change <- chm_t2 - chm_t1
  loss <- change <= threshold
  names(loss) <- "canopy_loss"

  loss_area <- sum(values(loss, na.rm = TRUE)) * prod(res(chm_t1)) / 10000
  message(sprintf("Surface de perte détectée: %.2f ha (seuil: %d m)",
                   loss_area, threshold))

  return(loss)
}

# ==============================================================================
# 4. Export pour SIG
# ==============================================================================

#' Exporter les résultats en GeoPackage
#'
#' @param raster_list Liste nommée de SpatRasters
#' @param filename Nom du fichier de sortie
#' @param output_dir Répertoire de sortie
export_to_gpkg <- function(raster_list, filename = "open_canopy_results.gpkg",
                            output_dir = OUTPUT_DIR) {
  out_path <- file.path(output_dir, filename)

  for (i in seq_along(raster_list)) {
    name <- names(raster_list)[i]
    r <- raster_list[[i]]

    # Convertir en polygones pour les données catégorielles
    if (is.logical(values(r)[1]) || all(values(r, na.rm = TRUE) %in% c(0, 1))) {
      v <- as.polygons(r)
      writeVector(v, out_path, layer = name,
                   overwrite = (i == 1))
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
  message("=== Open-Canopy: Prédiction et Post-traitement ===\n")
  message("Ce script nécessite soit:")
  message("  1. Des prédictions pré-calculées (.tif)")
  message("  2. Python + PyTorch pour exécuter les modèles\n")

  # Rechercher les fichiers de prédiction
  pred_files <- dir_ls(DATA_DIR, recurse = TRUE,
                        glob = "*predict*|*pred*|*output*")

  if (length(pred_files) > 0) {
    message("Fichiers de prédiction trouvés:")
    for (f in pred_files) message("  ", f)
  } else {
    message("Aucun fichier de prédiction trouvé.")
    message("Démonstration avec données simulées...\n")

    # Simulation
    set.seed(42)
    ref <- rast(nrows = 667, ncols = 667, xmin = 0, xmax = 1000,
                 ymin = 0, ymax = 1000)
    values(ref) <- pmax(0, rnorm(ncell(ref), mean = 10, sd = 8))
    names(ref) <- "reference_chm"

    pred <- ref + rnorm(ncell(ref), mean = 0, sd = 2)
    values(pred) <- pmax(0, values(pred))
    names(pred) <- "predicted_chm"

    # Évaluation
    metrics <- evaluate_predictions(pred, ref)

    # Visualisation
    pdf(file.path(OUTPUT_DIR, "demo_prediction.pdf"), width = 12, height = 10)
    plot_prediction_comparison(pred, ref, metrics, "Simulation")
    dev.off()

    # Couverture forestière
    forest_cover <- compute_forest_cover_grid(ref, cell_size = 50)

    pdf(file.path(OUTPUT_DIR, "demo_forest_cover.pdf"), width = 8, height = 8)
    plot(forest_cover, main = "Couverture forestière (%)",
         col = colorRampPalette(c("#fff7bc", "#41ab5d", "#004529"))(100))
    dev.off()

    message("\nGraphiques sauvegardés dans: ", OUTPUT_DIR)
  }

  message("\n=== Terminé ===")
}
