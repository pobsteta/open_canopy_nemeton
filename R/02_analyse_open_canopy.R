#!/usr/bin/env Rscript
# ==============================================================================
# 02_analyse_open_canopy.R
# Analyse et visualisation :
#   - Images SPOT 6-7 (Open-Canopy, 1.5m)
#   - Ortho IGN RVB + IRC (BD ORTHO®, 0.20m)
#   - CHM dérivés du LiDAR
# ==============================================================================

library(terra)
library(sf)
library(fs)

# ==============================================================================
# Configuration
# ==============================================================================

DATA_DIR     <- file.path(getwd(), "data")
DATA_DIR_HF  <- file.path(DATA_DIR, "open_canopy")
DATA_DIR_IGN <- file.path(DATA_DIR, "ign")
OUTPUT_DIR   <- file.path(getwd(), "outputs")
dir_create(OUTPUT_DIR)

# Résolutions
RES_SPOT <- 1.5
RES_IGN  <- 0.2

# ==============================================================================
# 1. Chargement des données
# ==============================================================================

#' Charger une image SPOT 6-7 (Open-Canopy)
#'
#' @param tile_path Chemin vers le fichier .tif
#' @return SpatRaster
load_spot_image <- function(tile_path) {
  if (!file.exists(tile_path)) stop("Fichier introuvable: ", tile_path)

  r <- rast(tile_path)
  message(sprintf("Image SPOT chargée: %s", basename(tile_path)))
  message(sprintf("  Dimensions: %d x %d | Bandes: %d | Résolution: %.2f m",
                   nrow(r), ncol(r), nlyr(r), res(r)[1]))
  return(r)
}

#' Charger une ortho IGN (RVB ou IRC, 0.20m)
#'
#' Supporte JPEG2000 (.jp2) et GeoTIFF (.tif)
#' RVB : bandes Rouge, Vert, Bleu
#' IRC : bandes Proche Infrarouge (PIR), Rouge, Vert
#'
#' @param file_path Chemin vers le fichier .jp2 ou .tif
#' @param type "rvb" ou "irc"
#' @return SpatRaster avec bandes nommées
load_ign_ortho <- function(file_path, type = "rvb") {
  if (!file.exists(file_path)) stop("Fichier introuvable: ", file_path)

  r <- rast(file_path)

  if (type == "rvb" && nlyr(r) >= 3) {
    names(r)[1:3] <- c("Rouge", "Vert", "Bleu")
  } else if (type == "irc" && nlyr(r) >= 3) {
    names(r)[1:3] <- c("PIR", "Rouge", "Vert")
  }

  message(sprintf("Ortho IGN %s chargée: %s", toupper(type), basename(file_path)))
  message(sprintf("  Dimensions: %d x %d | Bandes: %d (%s)",
                   nrow(r), ncol(r), nlyr(r), paste(names(r), collapse = ", ")))
  message(sprintf("  Résolution: %.2f m | CRS: %s",
                   res(r)[1], crs(r, describe = TRUE)$name))
  return(r)
}

#' Charger un CHM dérivé du LiDAR
#'
#' @param tile_path Chemin vers le fichier .tif
#' @return SpatRaster
load_chm <- function(tile_path) {
  if (!file.exists(tile_path)) stop("Fichier introuvable: ", tile_path)

  r <- rast(tile_path)
  message(sprintf("CHM chargé: %s", basename(tile_path)))
  message(sprintf("  Dimensions: %d x %d | Résolution: %.2f m",
                   nrow(r), ncol(r), res(r)[1]))

  vals <- values(r, na.rm = TRUE)
  if (length(vals) > 0) {
    message(sprintf("  Hauteur: min=%.2f, max=%.2f, moy=%.2f m",
                     min(vals), max(vals), mean(vals)))
  }
  return(r)
}

#' Lister et charger toutes les tuiles d'un split
load_tiles <- function(split = "test", data_type = "images",
                        data_dir = DATA_DIR_HF) {
  tile_dir <- file.path(data_dir, split, data_type)

  if (!dir.exists(tile_dir)) {
    message("Répertoire introuvable: ", tile_dir)
    return(list())
  }

  tif_files <- dir_ls(tile_dir, glob = "*.tif")
  message(sprintf("%d tuile(s) dans %s/%s", length(tif_files), split, data_type))

  tiles <- lapply(tif_files, function(f) {
    tryCatch(rast(f), error = function(e) {
      warning("Impossible de charger: ", f); NULL
    })
  })
  names(tiles) <- tools::file_path_sans_ext(basename(tif_files))
  Filter(Negate(is.null), tiles)
}

# ==============================================================================
# 2. Indices spectraux depuis les ortho IRC IGN
# ==============================================================================

#' Calculer le NDVI depuis une ortho IRC IGN
#'
#' NDVI = (PIR - Rouge) / (PIR + Rouge)
#' Nécessite une image IRC avec bandes PIR et Rouge
#'
#' @param irc_raster SpatRaster IRC IGN (bandes: PIR, Rouge, Vert)
#' @return SpatRaster du NDVI (valeurs entre -1 et 1)
compute_ndvi <- function(irc_raster) {
  pir <- irc_raster[["PIR"]]
  rouge <- irc_raster[["Rouge"]]

  ndvi <- (pir - rouge) / (pir + rouge)
  names(ndvi) <- "NDVI"

  vals <- values(ndvi, na.rm = TRUE)
  message(sprintf("NDVI calculé: min=%.3f, max=%.3f, moy=%.3f",
                   min(vals), max(vals), mean(vals)))
  return(ndvi)
}

#' Calculer le GNDVI (Green NDVI) depuis une ortho IRC IGN
#'
#' GNDVI = (PIR - Vert) / (PIR + Vert)
#'
#' @param irc_raster SpatRaster IRC IGN
#' @return SpatRaster du GNDVI
compute_gndvi <- function(irc_raster) {
  pir <- irc_raster[["PIR"]]
  vert <- irc_raster[["Vert"]]

  gndvi <- (pir - vert) / (pir + vert)
  names(gndvi) <- "GNDVI"
  return(gndvi)
}

#' Calculer le NDRE (Normalized Difference Red Edge) approximé
#'
#' Approximation basée sur la différence PIR - Vert normalisée.
#' Utile pour différencier les espèces et l'état de santé.
#'
#' @param irc_raster SpatRaster IRC IGN
#' @return SpatRaster
compute_savi <- function(irc_raster, L = 0.5) {
  pir <- irc_raster[["PIR"]]
  rouge <- irc_raster[["Rouge"]]

  savi <- ((pir - rouge) / (pir + rouge + L)) * (1 + L)
  names(savi) <- "SAVI"
  return(savi)
}

#' Créer un masque de végétation à partir du NDVI
#'
#' @param ndvi_raster SpatRaster du NDVI
#' @param threshold Seuil NDVI pour considérer de la végétation
#' @return SpatRaster binaire (1 = végétation)
mask_vegetation <- function(ndvi_raster, threshold = 0.3) {
  veg_mask <- ndvi_raster >= threshold
  names(veg_mask) <- "vegetation"

  pct <- sum(values(veg_mask, na.rm = TRUE)) / sum(!is.na(values(veg_mask))) * 100
  message(sprintf("Végétation détectée (NDVI >= %.2f): %.1f%%", threshold, pct))
  return(veg_mask)
}

# ==============================================================================
# 3. Visualisation
# ==============================================================================

#' Visualiser une image en couleurs naturelles (RGB)
#'
#' Fonctionne avec SPOT (Open-Canopy) et ortho IGN RVB
#'
#' @param raster_rgb SpatRaster avec au moins 3 bandes
#' @param title Titre du graphique
#' @param bands Indices des bandes RGB
plot_rgb <- function(raster_rgb, title = "Image RGB", bands = c(1, 2, 3)) {
  if (nlyr(raster_rgb) >= 3) {
    plotRGB(raster_rgb, r = bands[1], g = bands[2], b = bands[3],
            stretch = "lin", main = title)
  } else {
    plot(raster_rgb, main = title)
  }
}

#' Visualiser une ortho IRC en fausses couleurs
#'
#' Affichage IRC classique : PIR en rouge, Rouge en vert, Vert en bleu
#' La végétation active apparaît en rouge vif.
#'
#' @param irc_raster SpatRaster IRC IGN (PIR, Rouge, Vert)
#' @param title Titre
plot_irc <- function(irc_raster, title = "Ortho IRC IGN (0.20m)") {
  # Affichage fausses couleurs : PIR=R, R=G, V=B
  plotRGB(irc_raster, r = 1, g = 2, b = 3,
          stretch = "lin", main = title)
}

#' Visualiser le NDVI
#'
#' @param ndvi_raster SpatRaster du NDVI
#' @param title Titre
plot_ndvi <- function(ndvi_raster, title = "NDVI") {
  col_ndvi <- colorRampPalette(
    c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
      "#d9ef8b", "#91cf60", "#1a9850", "#006837")
  )(100)

  plot(ndvi_raster, main = title, col = col_ndvi, range = c(-0.2, 1),
       plg = list(title = "NDVI"))
}

#' Visualiser le CHM
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

#' Comparaison multi-panneaux : IGN RVB + IRC + NDVI + CHM
#'
#' @param rvb_raster SpatRaster ortho RVB IGN
#' @param irc_raster SpatRaster ortho IRC IGN
#' @param chm_raster SpatRaster CHM (optionnel)
#' @param tile_name Nom de la tuile
plot_ign_full_comparison <- function(rvb_raster, irc_raster,
                                      chm_raster = NULL,
                                      tile_name = "") {
  n_panels <- ifelse(is.null(chm_raster), 3, 4)
  ncols <- ifelse(n_panels == 4, 2, 3)
  nrows <- ifelse(n_panels == 4, 2, 1)

  par(mfrow = c(nrows, ncols), mar = c(2, 2, 3, 4))

  # RVB
  plot_rgb(rvb_raster, title = paste("Ortho RVB 0.20m", tile_name))

  # IRC fausses couleurs
  plot_irc(irc_raster, title = paste("IRC fausses couleurs", tile_name))

  # NDVI
  ndvi <- compute_ndvi(irc_raster)
  plot_ndvi(ndvi, title = paste("NDVI", tile_name))

  # CHM si disponible
  if (!is.null(chm_raster)) {
    plot_chm(chm_raster, title = paste("CHM LiDAR", tile_name))
  }

  par(mfrow = c(1, 1))
}

#' Comparaison SPOT (1.5m) vs IGN (0.20m) sur la même zone
#'
#' @param spot_raster SpatRaster SPOT
#' @param ign_raster SpatRaster ortho IGN RVB
#' @param tile_name Nom de la tuile
plot_resolution_comparison <- function(spot_raster, ign_raster,
                                        tile_name = "") {
  par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))

  plot_rgb(spot_raster,
           title = sprintf("SPOT 6-7 (%.1fm) %s", RES_SPOT, tile_name))
  plot_rgb(ign_raster,
           title = sprintf("Ortho IGN (%.1fm) %s", RES_IGN, tile_name))

  par(mfrow = c(1, 1))
}

#' Classification de la canopée
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
# 4. Statistiques et analyse
# ==============================================================================

#' Statistiques de hauteur de canopée
compute_chm_stats <- function(chm_raster) {
  vals <- values(chm_raster, na.rm = TRUE)

  data.frame(
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
}

#' Statistiques spectrales d'une ortho IRC IGN
#'
#' @param irc_raster SpatRaster IRC IGN
#' @return data.frame avec statistiques par bande + NDVI
compute_irc_stats <- function(irc_raster) {
  ndvi <- compute_ndvi(irc_raster)
  ndvi_vals <- values(ndvi, na.rm = TRUE)

  stats <- data.frame(
    pir_mean = mean(values(irc_raster[["PIR"]], na.rm = TRUE)),
    rouge_mean = mean(values(irc_raster[["Rouge"]], na.rm = TRUE)),
    vert_mean = mean(values(irc_raster[["Vert"]], na.rm = TRUE)),
    ndvi_mean = mean(ndvi_vals),
    ndvi_median = median(ndvi_vals),
    ndvi_sd = sd(ndvi_vals),
    pct_vegetation = sum(ndvi_vals >= 0.3) / length(ndvi_vals) * 100,
    pct_dense_veg = sum(ndvi_vals >= 0.6) / length(ndvi_vals) * 100,
    pct_bare_soil = sum(ndvi_vals < 0.1) / length(ndvi_vals) * 100,
    stringsAsFactors = FALSE
  )

  return(stats)
}

#' Croiser NDVI (IGN IRC) et CHM pour caractériser la végétation
#'
#' @param ndvi_raster SpatRaster NDVI (depuis ortho IRC IGN)
#' @param chm_raster SpatRaster CHM LiDAR
#' @return SpatRaster avec classes croisées
cross_ndvi_chm <- function(ndvi_raster, chm_raster) {
  # Aligner les résolutions si nécessaire
  if (!compareGeom(ndvi_raster, chm_raster, stopOnError = FALSE)) {
    message("Rééchantillonnage du NDVI vers la résolution du CHM...")
    ndvi_raster <- resample(ndvi_raster, chm_raster, method = "bilinear")
  }

  # Classes croisées :
  # 1 = Sol nu (NDVI < 0.2, CHM < 2m)
  # 2 = Herbe/culture (NDVI >= 0.2, CHM < 2m)
  # 3 = Arbuste sans feuilles (NDVI < 0.3, CHM 2-5m)
  # 4 = Arbuste feuillu (NDVI >= 0.3, CHM 2-5m)
  # 5 = Arbre peu vigoureux (NDVI < 0.4, CHM >= 5m)
  # 6 = Arbre vigoureux (NDVI >= 0.4, CHM >= 5m)

  cross <- chm_raster * 0  # template

  cross[ndvi_raster < 0.2  & chm_raster < 2]   <- 1
  cross[ndvi_raster >= 0.2 & chm_raster < 2]   <- 2
  cross[ndvi_raster < 0.3  & chm_raster >= 2 & chm_raster < 5] <- 3
  cross[ndvi_raster >= 0.3 & chm_raster >= 2 & chm_raster < 5] <- 4
  cross[ndvi_raster < 0.4  & chm_raster >= 5]  <- 5
  cross[ndvi_raster >= 0.4 & chm_raster >= 5]  <- 6

  names(cross) <- "ndvi_chm_class"
  return(cross)
}

#' Histogramme des hauteurs de canopée
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

#' Différence de canopée entre deux dates
compute_canopy_change <- function(chm_t1, chm_t2) {
  if (!compareGeom(chm_t1, chm_t2, stopOnError = FALSE)) {
    chm_t2 <- resample(chm_t2, chm_t1, method = "bilinear")
  }
  change <- chm_t2 - chm_t1
  names(change) <- "height_change"
  return(change)
}

#' Visualiser les changements de canopée
plot_canopy_change <- function(change_raster,
                                title = "Changement de hauteur de canopée (m)") {
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
# 5. Export
# ==============================================================================

export_raster <- function(raster_obj, filename, output_dir = OUTPUT_DIR) {
  out_path <- file.path(output_dir, filename)
  writeRaster(raster_obj, out_path, overwrite = TRUE)
  message("Raster exporté: ", out_path)
  return(out_path)
}

export_stats <- function(stats_df, filename = "statistics.csv",
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
  message("=== Analyse : Open-Canopy (SPOT 1.5m) + IGN Ortho (0.20m) ===\n")

  # --- Rechercher les fichiers disponibles ---
  all_images <- unlist(lapply(c("*.tif", "*.tiff", "*.jp2"), function(ext) {
    dir_ls(DATA_DIR, recurse = TRUE, glob = ext)
  }))

  if (length(all_images) == 0) {
    message("Aucun fichier image trouvé dans ", DATA_DIR)
    message("Exécutez d'abord: Rscript R/01_download_open_canopy.R")
    message("\nDémonstration avec des données simulées...\n")

    # --- Simulation IGN IRC + CHM ---
    set.seed(42)

    # Simuler une ortho IRC IGN à 0.20m (1km² = 5000x5000 px)
    # On réduit à 500x500 pour la démo
    demo_irc <- rast(nrows = 500, ncols = 500, nlyrs = 3,
                      xmin = 843000, xmax = 843100,
                      ymin = 6518000, ymax = 6518100,
                      crs = "EPSG:2154")
    names(demo_irc) <- c("PIR", "Rouge", "Vert")
    values(demo_irc) <- cbind(
      PIR = pmax(0, pmin(255, rnorm(ncell(demo_irc), 180, 40))),
      Rouge = pmax(0, pmin(255, rnorm(ncell(demo_irc), 90, 30))),
      Vert = pmax(0, pmin(255, rnorm(ncell(demo_irc), 100, 30)))
    )

    # Simuler un CHM
    demo_chm <- rast(nrows = 500, ncols = 500,
                      xmin = 843000, xmax = 843100,
                      ymin = 6518000, ymax = 6518100,
                      crs = "EPSG:2154")
    values(demo_chm) <- pmax(0, rnorm(ncell(demo_chm), mean = 8, sd = 6))
    names(demo_chm) <- "canopy_height"

    # NDVI
    ndvi <- compute_ndvi(demo_irc)

    # Stats
    message("\n--- Statistiques IRC ---")
    irc_stats <- compute_irc_stats(demo_irc)
    print(irc_stats)

    message("\n--- Statistiques CHM ---")
    chm_stats <- compute_chm_stats(demo_chm)
    print(chm_stats)

    # Visualisation
    pdf(file.path(OUTPUT_DIR, "demo_ign_analysis.pdf"), width = 14, height = 10)

    par(mfrow = c(2, 3), mar = c(2, 2, 3, 4))
    plot_irc(demo_irc, title = "IRC fausses couleurs (0.20m)")
    plot_ndvi(ndvi, title = "NDVI")
    plot_chm(demo_chm, title = "CHM LiDAR")
    plot_chm_histogram(demo_chm)
    plot_canopy_classes(demo_chm)

    # Croisement NDVI x CHM
    cross <- cross_ndvi_chm(ndvi, demo_chm)
    cross_labels <- c("Sol nu", "Herbe/culture", "Arbuste nu",
                       "Arbuste feuillu", "Arbre stress", "Arbre vigoureux")
    cross_colors <- c("#d7191c", "#fdae61", "#ffffbf",
                       "#abdda4", "#a6611a", "#018571")
    plot(cross, main = "Croisement NDVI x CHM", col = cross_colors,
         levels = cross_labels, type = "classes")

    dev.off()
    message("\nGraphiques: ", file.path(OUTPUT_DIR, "demo_ign_analysis.pdf"))

  } else {
    message(sprintf("%d fichier(s) image trouvé(s)\n", length(all_images)))

    # Séparer par type
    ign_irc <- all_images[grep("IRC|irc|infrarouge", all_images, ignore.case = TRUE)]
    ign_rvb <- all_images[grep("ign", all_images, ignore.case = TRUE)]
    ign_rvb <- setdiff(ign_rvb, ign_irc)
    spot_files <- all_images[grep("spot|open_canopy", all_images, ignore.case = TRUE)]
    chm_files <- all_images[grep("lidar|chm", all_images, ignore.case = TRUE)]

    message(sprintf("  Ortho IGN IRC: %d", length(ign_irc)))
    message(sprintf("  Ortho IGN RVB: %d", length(ign_rvb)))
    message(sprintf("  Images SPOT:   %d", length(spot_files)))
    message(sprintf("  CHM LiDAR:     %d", length(chm_files)))

    # Analyser les fichiers IRC IGN
    for (irc_path in ign_irc) {
      tile_name <- tools::file_path_sans_ext(basename(irc_path))
      message(sprintf("\n--- Analyse IRC IGN: %s ---", tile_name))

      irc <- load_ign_ortho(irc_path, type = "irc")
      ndvi <- compute_ndvi(irc)
      stats <- compute_irc_stats(irc)
      print(stats)

      pdf_path <- file.path(OUTPUT_DIR, paste0(tile_name, "_irc_analysis.pdf"))
      pdf(pdf_path, width = 12, height = 6)
      par(mfrow = c(1, 3), mar = c(2, 2, 3, 4))
      plot_irc(irc, title = paste("IRC", tile_name))
      plot_ndvi(ndvi, title = paste("NDVI", tile_name))

      veg <- mask_vegetation(ndvi)
      plot(veg, main = paste("Végétation", tile_name),
           col = c("white", "#1a9850"))
      dev.off()
      message("  Graphiques: ", pdf_path)
    }
  }

  message("\n=== Analyse terminée ===")
}
