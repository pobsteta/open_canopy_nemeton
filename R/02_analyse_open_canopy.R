#!/usr/bin/env Rscript
# ==============================================================================
# 02_analyse_open_canopy.R
# Analyse et visualisation du dataset Open-Canopy
# Chargement des images SPOT 6-7 et des CHM dérivés du LiDAR
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
# 1. Chargement des données
# ==============================================================================

#' Charger une image SPOT 6-7
#'
#' Les images SPOT sont des rasters multispectraux à 1.5 m de résolution
#'
#' @param tile_path Chemin vers le fichier .tif
#' @return SpatRaster
load_spot_image <- function(tile_path) {
  if (!file.exists(tile_path)) {
    stop("Fichier introuvable: ", tile_path)
  }

  r <- rast(tile_path)
  message(sprintf("Image SPOT chargée: %s", basename(tile_path)))
  message(sprintf("  Dimensions: %d x %d pixels", nrow(r), ncol(r)))
  message(sprintf("  Bandes: %d", nlyr(r)))
  message(sprintf("  Résolution: %.2f x %.2f m", res(r)[1], res(r)[2]))
  message(sprintf("  CRS: %s", crs(r, describe = TRUE)$name))
  return(r)
}

#' Charger un Canopy Height Model (CHM) dérivé du LiDAR
#'
#' @param tile_path Chemin vers le fichier .tif
#' @return SpatRaster
load_chm <- function(tile_path) {
  if (!file.exists(tile_path)) {
    stop("Fichier introuvable: ", tile_path)
  }

  r <- rast(tile_path)
  message(sprintf("CHM chargé: %s", basename(tile_path)))
  message(sprintf("  Dimensions: %d x %d pixels", nrow(r), ncol(r)))
  message(sprintf("  Résolution: %.2f x %.2f m", res(r)[1], res(r)[2]))

  vals <- values(r, na.rm = TRUE)
  if (length(vals) > 0) {
    message(sprintf("  Hauteur min: %.2f m", min(vals)))
    message(sprintf("  Hauteur max: %.2f m", max(vals)))
    message(sprintf("  Hauteur moyenne: %.2f m", mean(vals)))
  }

  return(r)
}

#' Lister et charger toutes les tuiles d'un split
#'
#' @param split Split du dataset ("train", "val", "test")
#' @param data_type "images" ou "lidar"
#' @param data_dir Répertoire racine des données
#' @return Liste de SpatRasters
load_tiles <- function(split = "test", data_type = "images",
                        data_dir = DATA_DIR) {
  tile_dir <- file.path(data_dir, split, data_type)

  if (!dir.exists(tile_dir)) {
    message("Répertoire introuvable: ", tile_dir)
    message("Recherche de fichiers .tif dans: ", data_dir)
    tif_files <- dir_ls(data_dir, recurse = TRUE, glob = "*.tif")
    if (length(tif_files) > 0) {
      message("Fichiers .tif trouvés:")
      for (f in tif_files) message("  ", f)
    }
    return(list())
  }

  tif_files <- dir_ls(tile_dir, glob = "*.tif")
  message(sprintf("%d tuile(s) trouvée(s) dans %s/%s",
                   length(tif_files), split, data_type))

  tiles <- lapply(tif_files, function(f) {
    tryCatch(rast(f), error = function(e) {
      warning("Impossible de charger: ", f, " - ", e$message)
      NULL
    })
  })

  names(tiles) <- tools::file_path_sans_ext(basename(tif_files))
  tiles <- Filter(Negate(is.null), tiles)
  return(tiles)
}

# ==============================================================================
# 2. Visualisation
# ==============================================================================

#' Visualiser une image SPOT en couleurs naturelles (RGB)
#'
#' @param spot_raster SpatRaster de l'image SPOT
#' @param title Titre du graphique
#' @param bands Indices des bandes RGB (défaut: 1, 2, 3)
plot_spot_rgb <- function(spot_raster, title = "Image SPOT 6-7",
                           bands = c(1, 2, 3)) {
  if (nlyr(spot_raster) >= 3) {
    plotRGB(spot_raster, r = bands[1], g = bands[2], b = bands[3],
            stretch = "lin", main = title)
  } else {
    plot(spot_raster, main = title)
  }
}

#' Visualiser le Canopy Height Model
#'
#' @param chm_raster SpatRaster du CHM
#' @param title Titre du graphique
#' @param col_palette Palette de couleurs
plot_chm <- function(chm_raster, title = "Canopy Height Model (m)",
                      col_palette = NULL) {
  if (is.null(col_palette)) {
    col_palette <- colorRampPalette(
      c("#f7fcb9", "#addd8e", "#41ab5d", "#006837", "#004529")
    )(100)
  }

  plot(chm_raster, main = title, col = col_palette,
       plg = list(title = "Hauteur (m)"))
}

#' Visualiser côte à côte l'image SPOT et le CHM
#'
#' @param spot_raster SpatRaster de l'image SPOT
#' @param chm_raster SpatRaster du CHM
#' @param tile_name Nom de la tuile
plot_spot_chm_comparison <- function(spot_raster, chm_raster,
                                       tile_name = "") {
  par(mfrow = c(1, 2), mar = c(2, 2, 3, 4))

  plot_spot_rgb(spot_raster,
                title = paste("SPOT 6-7", tile_name))
  plot_chm(chm_raster,
            title = paste("CHM LiDAR", tile_name))

  par(mfrow = c(1, 1))
}

#' Créer une carte de classification de la canopée
#'
#' @param chm_raster SpatRaster du CHM
#' @param breaks Seuils de hauteur pour la classification
#' @param labels Étiquettes des classes
#' @param title Titre
plot_canopy_classes <- function(chm_raster,
                                 breaks = c(0, 2, 5, 10, 20, Inf),
                                 labels = c("Sol/herbe (<2m)",
                                            "Arbustes (2-5m)",
                                            "Petits arbres (5-10m)",
                                            "Arbres moyens (10-20m)",
                                            "Grands arbres (>20m)"),
                                 title = "Classes de hauteur de canopée") {
  classes <- classify(chm_raster, rcl = breaks, include.lowest = TRUE)

  colors <- c("#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#253494")
  plot(classes, main = title, col = colors,
       levels = labels, type = "classes",
       plg = list(legend = labels))
}

# ==============================================================================
# 3. Statistiques et analyse
# ==============================================================================

#' Calculer les statistiques de hauteur de canopée
#'
#' @param chm_raster SpatRaster du CHM
#' @return data.frame avec les statistiques
compute_chm_stats <- function(chm_raster) {
  vals <- values(chm_raster, na.rm = TRUE)

  stats <- data.frame(
    n_pixels = length(vals),
    n_na = sum(is.na(values(chm_raster))),
    min_height = min(vals),
    max_height = max(vals),
    mean_height = mean(vals),
    median_height = median(vals),
    sd_height = sd(vals),
    q25 = quantile(vals, 0.25),
    q75 = quantile(vals, 0.75),
    pct_forest = sum(vals >= 2) / length(vals) * 100,
    pct_tall_trees = sum(vals >= 20) / length(vals) * 100,
    stringsAsFactors = FALSE
  )

  return(stats)
}

#' Histogramme des hauteurs de canopée
#'
#' @param chm_raster SpatRaster du CHM
#' @param title Titre du graphique
#' @param n_breaks Nombre de classes
plot_chm_histogram <- function(chm_raster,
                                title = "Distribution des hauteurs de canopée",
                                n_breaks = 50) {
  vals <- values(chm_raster, na.rm = TRUE)

  hist(vals, breaks = n_breaks, main = title,
       xlab = "Hauteur (m)", ylab = "Fréquence",
       col = "#41ab5d", border = "white")

  abline(v = mean(vals), col = "red", lwd = 2, lty = 2)
  legend("topright",
         legend = sprintf("Moyenne: %.1f m", mean(vals)),
         col = "red", lty = 2, lwd = 2, bty = "n")
}

#' Calculer la différence de canopée entre deux dates
#'
#' @param chm_t1 CHM à la date 1
#' @param chm_t2 CHM à la date 2
#' @return SpatRaster de la différence (t2 - t1)
compute_canopy_change <- function(chm_t1, chm_t2) {
  if (!compareGeom(chm_t1, chm_t2, stopOnError = FALSE)) {
    message("Alignement des rasters...")
    chm_t2 <- resample(chm_t2, chm_t1, method = "bilinear")
  }

  change <- chm_t2 - chm_t1
  names(change) <- "height_change"
  return(change)
}

#' Visualiser les changements de canopée
#'
#' @param change_raster SpatRaster de la différence de hauteur
#' @param title Titre
plot_canopy_change <- function(change_raster,
                                title = "Changement de hauteur de canopée (m)") {
  # Palette divergente: rouge (perte) -> blanc (stable) -> bleu (gain)
  col_palette <- colorRampPalette(
    c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
      "#d9ef8b", "#91cf60", "#1a9850")
  )(100)

  vals <- values(change_raster, na.rm = TRUE)
  max_abs <- max(abs(range(vals)))

  plot(change_raster, main = title, col = col_palette,
       range = c(-max_abs, max_abs),
       plg = list(title = "Delta H (m)"))
}

# ==============================================================================
# 4. Export des résultats
# ==============================================================================

#' Exporter un raster en GeoTIFF
#'
#' @param raster_obj SpatRaster à exporter
#' @param filename Nom du fichier de sortie
#' @param output_dir Répertoire de sortie
export_raster <- function(raster_obj, filename, output_dir = OUTPUT_DIR) {
  out_path <- file.path(output_dir, filename)
  writeRaster(raster_obj, out_path, overwrite = TRUE)
  message("Raster exporté: ", out_path)
  return(out_path)
}

#' Exporter les statistiques en CSV
#'
#' @param stats_df data.frame de statistiques
#' @param filename Nom du fichier de sortie
#' @param output_dir Répertoire de sortie
export_stats <- function(stats_df, filename = "chm_statistics.csv",
                          output_dir = OUTPUT_DIR) {
  out_path <- file.path(output_dir, filename)
  write.csv(stats_df, out_path, row.names = FALSE)
  message("Statistiques exportées: ", out_path)
  return(out_path)
}

# ==============================================================================
# Exécution principale
# ==============================================================================

if (sys.nframe() == 0) {
  message("=== Open-Canopy: Analyse et Visualisation ===\n")

  # Rechercher les fichiers .tif téléchargés
  tif_files <- dir_ls(DATA_DIR, recurse = TRUE, glob = "*.tif")

  if (length(tif_files) == 0) {
    message("Aucun fichier .tif trouvé dans ", DATA_DIR)
    message("Exécutez d'abord: Rscript R/01_download_open_canopy.R")
    message("\nDémonstration avec des données simulées...")

    # Données simulées pour démonstration
    set.seed(42)
    demo_chm <- rast(nrows = 667, ncols = 667, xmin = 0, xmax = 1000,
                      ymin = 0, ymax = 1000)
    values(demo_chm) <- pmax(0, rnorm(ncell(demo_chm), mean = 8, sd = 6))
    names(demo_chm) <- "canopy_height"

    message("\n--- Statistiques du CHM simulé ---")
    stats <- compute_chm_stats(demo_chm)
    print(stats)

    # Visualisation
    pdf(file.path(OUTPUT_DIR, "demo_analysis.pdf"), width = 12, height = 8)

    par(mfrow = c(2, 2))
    plot_chm(demo_chm, title = "CHM simulé")
    plot_chm_histogram(demo_chm, title = "Distribution des hauteurs")
    plot_canopy_classes(demo_chm, title = "Classes de végétation")

    # Simulation de changement
    demo_chm_t2 <- demo_chm - abs(rnorm(ncell(demo_chm), mean = 1, sd = 2))
    change <- compute_canopy_change(demo_chm, demo_chm_t2)
    plot_canopy_change(change, title = "Changement simulé")

    dev.off()
    message("\nGraphiques sauvegardés: ", file.path(OUTPUT_DIR, "demo_analysis.pdf"))
  } else {
    message(sprintf("%d fichier(s) .tif trouvé(s)", length(tif_files)))

    # Séparer images et CHM
    spot_files <- tif_files[grep("images|spot", tif_files, ignore.case = TRUE)]
    chm_files <- tif_files[grep("lidar|chm", tif_files, ignore.case = TRUE)]

    message(sprintf("  Images SPOT: %d", length(spot_files)))
    message(sprintf("  CHM LiDAR: %d", length(chm_files)))

    # Analyser chaque tuile CHM
    all_stats <- data.frame()
    for (chm_path in chm_files) {
      tile_name <- tools::file_path_sans_ext(basename(chm_path))
      message(sprintf("\n--- Analyse de la tuile: %s ---", tile_name))

      chm <- load_chm(chm_path)
      stats <- compute_chm_stats(chm)
      stats$tile <- tile_name
      all_stats <- rbind(all_stats, stats)

      # Visualisation
      pdf_path <- file.path(OUTPUT_DIR, paste0(tile_name, "_analysis.pdf"))
      pdf(pdf_path, width = 10, height = 8)

      par(mfrow = c(2, 2))
      plot_chm(chm, title = paste("CHM -", tile_name))
      plot_chm_histogram(chm)
      plot_canopy_classes(chm)

      # Si image SPOT correspondante disponible
      spot_match <- spot_files[grep(tile_name, spot_files)]
      if (length(spot_match) > 0) {
        spot <- load_spot_image(spot_match[1])
        plot_spot_rgb(spot, title = paste("SPOT -", tile_name))
      }

      dev.off()
      message("  Graphiques: ", pdf_path)
    }

    # Exporter les statistiques
    if (nrow(all_stats) > 0) {
      export_stats(all_stats)
      message("\n=== Résumé global ===")
      message(sprintf("Hauteur moyenne: %.2f m", mean(all_stats$mean_height)))
      message(sprintf("Couvert forestier moyen: %.1f%%",
                       mean(all_stats$pct_forest)))
    }
  }

  message("\n=== Analyse terminée ===")
}
