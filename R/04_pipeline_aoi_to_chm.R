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
#   3. Combiner RVB + PIR en image 4 bandes (R, V, B, PIR)
#   4. Rééchantillonner de 0.20m → 1.5m (résolution SPOT / Open-Canopy)
#   5. Inférence du modèle Open-Canopy (UNet/SMP, 4 canaux) via reticulate
#   6. Mosaïquer et exporter le CHM prédit
#
# Architecture Open-Canopy :
#   - Modèle UNet (smp) avec encodeur ResNet34, 4 canaux → 1 sortie
#   - Modèle PVTv2 (timm) avec pvt_v2_b3, 4 canaux → 1 sortie
#   - Entrée : 4 bandes (R, G, B, NIR) avec normalisation extraite du checkpoint
#   - Sortie : hauteur de canopée en mètres (targets stockés en dm / 10)
#   - Checkpoint : PyTorch Lightning, clés préfixées "net.seg_model."
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
IGN_LAYER_IRC    <- "ORTHOIMAGERY.ORTHOPHOTOS.IRC"

# --- Millésime (NULL = couche la plus récente disponible) ---
# Exemples : "2024", "2023", "2021"...
# Ortho RVB → couche ORTHOIMAGERY.ORTHOPHOTOS{année}
# IRC       → couche ORTHOIMAGERY.ORTHOPHOTOS.IRC.{année}
MILLESIME_ORTHO <- NULL
MILLESIME_IRC   <- NULL

#' Construire le nom de couche WMS en fonction du millésime
#'
#' @param type "ortho" ou "irc"
#' @param millesime Année (chaîne ou numérique), NULL pour la couche courante
#' @return Nom de couche WMS IGN
ign_layer_name <- function(type = c("ortho", "irc"), millesime = NULL) {
  type <- match.arg(type)
  if (is.null(millesime)) {
    if (type == "ortho") return(IGN_LAYER_ORTHO)
    else                 return(IGN_LAYER_IRC)
  }
  millesime <- as.character(millesime)
  if (type == "ortho") {
    return(paste0("ORTHOIMAGERY.ORTHOPHOTOS", millesime))
  } else {
    return(paste0("ORTHOIMAGERY.ORTHOPHOTOS.IRC.", millesime))
  }
}

# --- Résolutions ---
RES_IGN  <- 0.2   # BD ORTHO® IGN
RES_SPOT <- 1.5   # Modèles Open-Canopy (SPOT 6-7)

# --- Modèle Open-Canopy ---
HF_REPO_ID <- "AI4Forest/Open-Canopy"
CONDA_ENV  <- "open_canopy"
N_INPUT_CHANNELS <- 4  # R, G, B, NIR (ordre SPOT 6-7)
# Chemin vers le code source Open-Canopy (optionnel pour PVTv2)
OPEN_CANOPY_SRC <- NULL  # Auto-détecté, téléchargé, ou module embarqué

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

  # WMS 1.3.0 avec CRS EPSG:2154 : BBOX = xmin,ymin,xmax,ymax (easting, northing)
  wms_url <- paste0(
    IGN_WMS_URL, "?",
    "SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap",
    "&LAYERS=", layer,
    "&CRS=EPSG:2154",
    "&BBOX=", paste(xmin, ymin, xmax, ymax, sep = ","),
    "&WIDTH=", width,
    "&HEIGHT=", height,
    "&FORMAT=image/geotiff",
    "&STYLES="
  )

  tryCatch({
    # Télécharger dans un fichier temporaire
    tmp_file <- tempfile(fileext = ".tif")
    curl_download(url = wms_url, destfile = tmp_file, quiet = TRUE)

    # Vérifier que c'est bien un raster (pas un XML d'erreur)
    r <- rast(tmp_file)

    # Assigner le CRS et l'emprise si nécessaire
    needs_fix <- FALSE
    if (is.na(crs(r)) || crs(r) == "") {
      crs(r) <- "EPSG:2154"
      needs_fix <- TRUE
    }

    # Forcer l'emprise correcte
    ext(r) <- ext(xmin, xmax, ymin, ymax)

    # Écrire le fichier final (depuis le temp, pas de conflit)
    writeRaster(r, dest_file, overwrite = TRUE)

    # Libérer la source temp et nettoyer
    r <- rast(dest_file)
    unlink(tmp_file)

    return(r)
  }, error = function(e) {
    unlink(tmp_file)
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
#' @param millesime_ortho Millésime ortho RVB (NULL = plus récent)
#' @param millesime_irc Millésime IRC (NULL = plus récent)
#' @return Liste avec rvb, irc (SpatRaster) et millésimes utilisés
download_ortho_for_aoi <- function(aoi, output_dir, res_m = RES_IGN,
                                    millesime_ortho = MILLESIME_ORTHO,
                                    millesime_irc = MILLESIME_IRC) {
  dir_create(output_dir)

  # Résoudre les couches WMS selon le millésime
  layer_ortho <- ign_layer_name("ortho", millesime_ortho)
  layer_irc   <- ign_layer_name("irc",   millesime_irc)
  label_ortho <- if (is.null(millesime_ortho)) "plus récent" else millesime_ortho
  label_irc   <- if (is.null(millesime_irc))   "plus récent" else millesime_irc

  # Vérifier si les fichiers existent déjà (cache)
  rvb_path <- file.path(output_dir, "ortho_rvb.tif")
  irc_path <- file.path(output_dir, "ortho_irc.tif")

  if (file.exists(rvb_path) && file.exists(irc_path)) {
    message("\n=== Ortho IGN déjà présentes (cache) ===")
    message(sprintf("  RVB: %s", rvb_path))
    message(sprintf("  IRC: %s", irc_path))
    message("Réutilisation des fichiers existants (supprimez-les pour forcer le re-téléchargement).")

    rvb <- rast(rvb_path)
    irc <- rast(irc_path)
    # Restaurer les noms de bandes (perdus lors de l'écriture GeoTIFF)
    names(rvb)[1:min(3, nlyr(rvb))] <- c("Rouge", "Vert", "Bleu")[1:min(3, nlyr(rvb))]
    names(irc)[1:min(3, nlyr(irc))] <- c("PIR", "Rouge", "Vert")[1:min(3, nlyr(irc))]

    return(list(rvb = rvb, irc = irc,
                rvb_path = rvb_path, irc_path = irc_path,
                millesime_ortho = millesime_ortho,
                millesime_irc = millesime_irc,
                layer_ortho = layer_ortho,
                layer_irc = layer_irc))
  }

  bbox <- as.numeric(st_bbox(st_union(aoi)))

  message(sprintf("\n=== Téléchargement ortho IGN pour l'AOI ==="))
  message(sprintf("Emprise: %.0f, %.0f - %.0f, %.0f (Lambert-93)",
                   bbox[1], bbox[2], bbox[3], bbox[4]))
  message(sprintf("Taille: %.0f x %.0f m (%.2f ha)",
                   bbox[3] - bbox[1], bbox[4] - bbox[2],
                   (bbox[3] - bbox[1]) * (bbox[4] - bbox[2]) / 10000))
  message(sprintf("Millésime RVB: %s (couche: %s)", label_ortho, layer_ortho))
  message(sprintf("Millésime IRC: %s (couche: %s)", label_irc, layer_irc))

  # --- RVB ---
  message("\n--- Ortho RVB ---")
  rvb <- tryCatch(
    download_ign_tiled(bbox, layer = layer_ortho, res_m = res_m,
                       output_dir = output_dir, prefix = "rvb"),
    error = function(e) {
      if (!is.null(millesime_ortho)) {
        message(sprintf("  Couche %s indisponible, fallback sur %s",
                        layer_ortho, IGN_LAYER_ORTHO))
        layer_ortho <<- IGN_LAYER_ORTHO
        label_ortho <<- "plus récent (fallback)"
        download_ign_tiled(bbox, layer = IGN_LAYER_ORTHO, res_m = res_m,
                           output_dir = output_dir, prefix = "rvb")
      } else stop(e)
    }
  )
  names(rvb)[1:min(3, nlyr(rvb))] <- c("Rouge", "Vert", "Bleu")[1:min(3, nlyr(rvb))]

  # --- IRC ---
  message("\n--- Ortho IRC ---")
  irc <- tryCatch(
    download_ign_tiled(bbox, layer = layer_irc, res_m = res_m,
                       output_dir = output_dir, prefix = "irc"),
    error = function(e) {
      if (!is.null(millesime_irc)) {
        message(sprintf("  Couche %s indisponible, fallback sur %s",
                        layer_irc, IGN_LAYER_IRC))
        layer_irc <<- IGN_LAYER_IRC
        label_irc <<- "plus récent (fallback)"
        download_ign_tiled(bbox, layer = IGN_LAYER_IRC, res_m = res_m,
                           output_dir = output_dir, prefix = "irc")
      } else stop(e)
    }
  )
  names(irc)[1:min(3, nlyr(irc))] <- c("PIR", "Rouge", "Vert")[1:min(3, nlyr(irc))]

  # Découper aux limites exactes de l'AOI
  aoi_vect <- vect(st_union(aoi))
  rvb <- crop(rvb, aoi_vect)
  irc <- crop(irc, aoi_vect)

  # Sauvegarder les mosaïques finales
  writeRaster(rvb, rvb_path, overwrite = TRUE)
  writeRaster(irc, irc_path, overwrite = TRUE)

  # Recharger depuis les fichiers consolidés pour que les SpatRaster
  # pointent vers ortho_rvb.tif / ortho_irc.tif (et non les tuiles)
  rvb <- rast(rvb_path)
  irc <- rast(irc_path)
  # Restaurer les noms de bandes (perdus lors de l'écriture GeoTIFF)
  names(rvb)[1:min(3, nlyr(rvb))] <- c("Rouge", "Vert", "Bleu")[1:min(3, nlyr(rvb))]
  names(irc)[1:min(3, nlyr(irc))] <- c("PIR", "Rouge", "Vert")[1:min(3, nlyr(irc))]

  message(sprintf("\nRVB sauvegardé: %s (%d x %d px)", rvb_path, ncol(rvb), nrow(rvb)))
  message(sprintf("IRC sauvegardé: %s (%d x %d px)", irc_path, ncol(irc), nrow(irc)))

  # Nettoyer les tuiles temporaires
  tile_files <- dir_ls(output_dir, glob = "*_tile_*.tif")
  if (length(tile_files) > 0) file_delete(tile_files)

  return(list(rvb = rvb, irc = irc,
              rvb_path = rvb_path, irc_path = irc_path,
              millesime_ortho = millesime_ortho,
              millesime_irc = millesime_irc,
              layer_ortho = layer_ortho,
              layer_irc = layer_irc))
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
#'
#' Vérifie que les modules Python nécessaires sont disponibles.
#' Note : huggingface_hub Python n'est plus requis si le package R hfhub
#' est installé (méthode préférée pour le téléchargement des modèles).
setup_python <- function() {
  library(reticulate)

  # Vérifier si Python est déjà configuré (ex: lancé depuis conda activate)
  py_ok <- tryCatch(nzchar(py_config()$python), error = function(e) FALSE)

  if (!py_ok) {
    # Chercher l'env conda par nom (Miniforge, Miniconda, Anaconda)
    conda_envs <- tryCatch(conda_list(), error = function(e) data.frame())
    if (nrow(conda_envs) > 0) {
      match_idx <- grep(CONDA_ENV, conda_envs$name, ignore.case = TRUE)
      if (length(match_idx) > 0) {
        use_python(conda_envs$python[match_idx[1]], required = TRUE)
        message("  Env conda détecté: ", conda_envs$name[match_idx[1]])
      }
    }
    # Fallback : CONDA_PREFIX
    if (!tryCatch(nzchar(py_config()$python), error = function(e) FALSE)) {
      conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
      if (nzchar(conda_prefix)) {
        py_path <- if (.Platform$OS.type == "windows") {
          file.path(conda_prefix, "python.exe")
        } else {
          file.path(conda_prefix, "bin", "python")
        }
        if (file.exists(py_path)) {
          use_python(py_path, required = TRUE)
          message("  Python (CONDA_PREFIX): ", conda_prefix)
        }
      }
    }
  }

  # Modules Python essentiels (inférence)
  modules_core <- c("torch", "numpy", "rasterio",
                     "segmentation_models_pytorch", "timm")
  # huggingface_hub Python est optionnel si hfhub R est disponible
  has_hfhub_r <- requireNamespace("hfhub", quietly = TRUE)

  ok <- TRUE
  for (mod in modules_core) {
    avail <- py_module_available(mod)
    message(sprintf("  Python %s: %s", mod, ifelse(avail, "OK", "MANQUANT")))
    if (!avail) ok <- FALSE
  }

  # Vérifier huggingface_hub Python seulement si hfhub R n'est pas dispo
  if (has_hfhub_r) {
    message("  hfhub (R natif): OK (téléchargement HF sans Python)")
  } else {
    hf_avail <- py_module_available("huggingface_hub")
    message(sprintf("  Python huggingface_hub: %s",
                     ifelse(hf_avail, "OK", "MANQUANT")))
    if (!hf_avail) ok <- FALSE
  }

  if (!ok) {
    pip_cmd <- "  pip install torch torchvision numpy rasterio segmentation-models-pytorch timm"
    if (!has_hfhub_r) {
      pip_cmd <- paste0(pip_cmd, " huggingface_hub")
    }
    stop("Modules Python manquants. Installez-les dans l'env '", CONDA_ENV, "':\n",
         "  conda activate ", CONDA_ENV, "\n", pip_cmd,
         "\n\nOu installez le package R hfhub pour éviter la dépendance Python huggingface_hub :\n",
         "  install.packages('hfhub')")
  }
}

#' Trouver le nom du fichier checkpoint dans le dataset HF via l'API
#'
#' Interroge l'API Hugging Face pour découvrir les fichiers checkpoint
#' disponibles dans le répertoire pretrained_models/ du dataset.
#'
#' @param repo_id Identifiant du dépôt HF (ex: "AI4Forest/Open-Canopy")
#' @param model_name "unet" ou "pvtv2" pour filtrer
#' @return Chemin du fichier checkpoint dans le dépôt, ou NULL
find_checkpoint_name <- function(repo_id = HF_REPO_ID,
                                  model_name = "pvtv2") {
  url <- paste0("https://huggingface.co/api/datasets/", repo_id)
  resp <- tryCatch(
    httr2::request(url) |> httr2::req_perform(),
    error = function(e) NULL
  )
  if (is.null(resp)) return(NULL)

  info <- jsonlite::fromJSON(httr2::resp_body_string(resp))
  files <- info$siblings$rfilename
  ckpt_files <- files[grepl("^pretrained_models/.*\\.ckpt$", files)]
  if (length(ckpt_files) == 0) return(NULL)

  # Filtrer par nom de modèle
  pattern <- switch(model_name,
    unet  = "unet|smp",
    pvtv2 = "pvt|pvtv2",
    model_name
  )
  matches <- ckpt_files[grep(pattern, ckpt_files, ignore.case = TRUE)]
  if (length(matches) > 0) return(matches[1])
  return(ckpt_files[1])
}

#' Résoudre le chemin réel d'un fichier HuggingFace cache
#'
#' Sur Windows, le cache HF utilise des symlinks (snapshots/ → blobs/)
#' qui ne fonctionnent pas sans le mode développeur. Cette fonction
#' résout le chemin réel ou copie le fichier blob si nécessaire.
#'
#' @param path Chemin retourné par hub_download
#' @return Chemin réel accessible
#' @keywords internal
.resolve_hf_path <- function(path) {
  path <- as.character(path)

  # Si le fichier est directement lisible, tout va bien
  if (file.exists(path) && file.size(path) > 0) {
    return(normalizePath(path, winslash = "/"))
  }

  # Sur Windows : le symlink snapshots/xxx/file → blobs/yyy ne marche pas
  # Chercher le fichier dans le dossier blobs/ du même repo
  if (.Platform$OS.type == "windows" && grepl("snapshots", path)) {
    repo_dir <- sub("/snapshots/.*$", "", path)
    blobs_dir <- file.path(repo_dir, "blobs")

    if (dir.exists(blobs_dir)) {
      # Les blobs sont nommés par leur hash SHA256
      # Chercher le plus gros fichier .ckpt-compatible (le checkpoint)
      blobs <- list.files(blobs_dir, full.names = TRUE)
      if (length(blobs) > 0) {
        sizes <- file.size(blobs)
        # Prendre le plus gros blob (le checkpoint est le plus volumineux)
        best <- blobs[which.max(sizes)]
        if (!is.na(file.size(best)) && file.size(best) > 1e6) {
          message("  Windows: résolution symlink HF → ", best)
          return(normalizePath(best, winslash = "/"))
        }
      }
    }

    # Alternative : lire le fichier symlink comme texte
    # (HF crée parfois un fichier texte contenant le hash du blob)
    if (file.exists(path) && file.size(path) < 200) {
      blob_hash <- trimws(readLines(path, n = 1, warn = FALSE))
      blob_path <- file.path(repo_dir, "blobs", blob_hash)
      if (file.exists(blob_path) && file.size(blob_path) > 1e6) {
        message("  Windows: résolution hash HF → ", blob_path)
        return(normalizePath(blob_path, winslash = "/"))
      }
    }

    warning("Chemin HuggingFace inaccessible: ", path, "\n",
            "  Activez le mode développeur Windows ou téléchargez le modèle manuellement.")
  }

  path
}

#' Télécharger le modèle pré-entraîné depuis Hugging Face
#'
#' Utilise le package R hfhub (natif, sans Python) pour télécharger
#' le fichier checkpoint. Fallback sur Python huggingface_hub si hfhub
#' n'est pas installé.
#'
#' @param model_name "unet" ou "pvtv2"
#' @return Chemin local du modèle
download_model <- function(model_name = "pvtv2") {
  message("Téléchargement du modèle: ", model_name)
  message("Depuis: ", HF_REPO_ID)

  # --- Méthode 1 : hfhub R natif (préféré) ---
  if (requireNamespace("hfhub", quietly = TRUE)) {
    ckpt_name <- find_checkpoint_name(HF_REPO_ID, model_name)
    if (!is.null(ckpt_name)) {
      message("  Téléchargement via hfhub (R natif): ", ckpt_name)
      tryCatch({
        local_path <- hfhub::hub_download(HF_REPO_ID, ckpt_name,
                                            repo_type = "dataset")
        local_path <- .resolve_hf_path(local_path)
        message("  Modèle téléchargé: ", local_path)
        return(local_path)
      }, error = function(e) {
        message("  hfhub échoué: ", e$message, " \u2192 fallback Python")
      })
    }
  }

  # --- Méthode 2 : Python huggingface_hub (fallback) ---
  library(reticulate)
  hf_hub <- import("huggingface_hub")

  # Vérifier que le token HuggingFace est configuré (dataset gated)
  token <- Sys.getenv("HF_TOKEN", unset = "")
  if (token == "") {
    tryCatch({
      stored <- hf_hub$HfFolder$get_token()
      if (is.null(stored) || stored == "") stop("no token")
    }, error = function(e) {
      message("ATTENTION: Aucun token HuggingFace détecté.")
      message("Le dataset AI4Forest/Open-Canopy est en accès restreint.")
      message("Connectez-vous d'abord avec l'une de ces méthodes :")
      message("  1. Dans R:     Sys.setenv(HF_TOKEN = 'hf_votre_token')")
      message("  2. En terminal: huggingface-cli login")
      message("  3. Créez un token sur: https://huggingface.co/settings/tokens")
    })
  }

  message("Recherche des checkpoints via Python huggingface_hub...")

  tryCatch({
    api <- hf_hub$HfApi()

    all_files <- tryCatch({
      files <- api$list_repo_files(HF_REPO_ID, repo_type = "dataset")
      reticulate::iterate(files)
    }, error = function(e) {
      tree <- api$list_repo_tree(
        HF_REPO_ID,
        path_in_repo = "pretrained_models",
        repo_type = "dataset"
      )
      items <- tryCatch(reticulate::iterate(tree), error = function(e2) as.list(tree))
      paths <- character(0)
      for (item in items) {
        p <- tryCatch(item$path, error = function(e3)
          tryCatch(item$rfilename, error = function(e4)
            tryCatch(as.character(item), error = function(e5) "")))
        if (nzchar(p)) paths <- c(paths, p)
      }
      paths
    })

    ckpt_files <- character(0)
    for (f in all_files) {
      fname <- as.character(f)
      if (grepl("^pretrained_models/", fname) && grepl("\\.ckpt$", fname)) {
        ckpt_files <- c(ckpt_files, basename(fname))
      }
    }

    if (length(ckpt_files) == 0) {
      stop("Aucun fichier .ckpt trouvé dans pretrained_models/")
    }

    message("Checkpoints disponibles:")
    for (f in ckpt_files) message("  - ", f)

    pattern <- switch(model_name,
      unet  = "unet|smp",
      pvtv2 = "pvt|pvtv2",
      model_name
    )
    match_idx <- grep(pattern, ckpt_files, ignore.case = TRUE)

    if (length(match_idx) == 0) {
      message("Aucun match pour '", model_name, "', utilisation du premier checkpoint.")
      target_file <- ckpt_files[1]
    } else {
      target_file <- ckpt_files[match_idx[1]]
    }

    message("Téléchargement: ", target_file)
    local_path <- hf_hub$hf_hub_download(
      repo_id  = HF_REPO_ID,
      filename = paste0("pretrained_models/", target_file),
      repo_type = "dataset"
    )
    local_path <- .resolve_hf_path(local_path)
    message("Modèle: ", local_path)
    return(local_path)

  }, error = function(e) {
    message("Impossible de lister le dataset HuggingFace: ", e$message)
    message("Le dataset est peut-être privé (gated). Vérifiez votre HF_TOKEN.")
    message("")
    message("Alternatives :")
    message("  1. Définir votre token : Sys.setenv(HF_TOKEN = 'hf_...')")
    message("  2. Se connecter via CLI : huggingface-cli login")
    message("  3. Installer hfhub R : install.packages('hfhub')")
    message("  4. Fournir le chemin manuellement :")
    message('     result <- pipeline_aoi_to_chm("data/aoi.gpkg",')
    message('       model_path = "chemin/vers/checkpoint.ckpt")')
    stop("Échec du téléchargement du modèle.", call. = FALSE)
  })
}

#' Télécharger le code source Open-Canopy
#'
#' Clone le dépôt GitHub Open-Canopy dans un dossier cache local.
#' Réutilise le clone existant si déjà présent.
#'
#' @param dest Dossier de destination (NULL = cache utilisateur)
#' @param force Forcer le re-téléchargement même si déjà présent
#' @return Chemin vers le dossier Open-Canopy cloné, ou NULL en cas d'échec
#' @export
download_open_canopy_src <- function(dest = NULL, force = FALSE) {
  repo_url <- "https://github.com/fajwel/Open-Canopy.git"

  # Dossier cache par défaut
  if (is.null(dest)) {
    cache_dir <- tools::R_user_dir("opencanopy", "cache")
    dest <- file.path(cache_dir, "Open-Canopy")
  }

  # Vérifier si déjà présent
  if (dir.exists(file.path(dest, "src", "models")) && !force) {
    message("Open-Canopy source déjà en cache: ", dest)
    return(dest)
  }

  # Vérifier que git est disponible
  git_ok <- tryCatch({
    res <- system2("git", "--version", stdout = TRUE, stderr = TRUE)
    length(res) > 0
  }, error = function(e) FALSE)

  if (!git_ok) {
    message("git n'est pas disponible, téléchargement par archive ZIP")
    return(.download_open_canopy_zip(dest))
  }

  # Clone superficiel (seulement le dernier commit)
  message("Clonage d'Open-Canopy (shallow clone)...")
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  # Supprimer le dossier incomplet si force
  if (force && dir.exists(dest)) {
    unlink(dest, recursive = TRUE)
  }

  res <- tryCatch({
    system2("git", c("clone", "--depth", "1", repo_url, dest),
            stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    message("Erreur git clone: ", e$message)
    return(NULL)
  })

  if (dir.exists(file.path(dest, "src", "models"))) {
    message("Open-Canopy source téléchargé: ", dest)
    return(dest)
  }

  message("Le clone git a échoué, tentative par archive ZIP...")
  .download_open_canopy_zip(dest)
}

#' Télécharger Open-Canopy par archive ZIP (fallback sans git)
#' @param dest Dossier de destination
#' @return Chemin vers le dossier Open-Canopy, ou NULL en cas d'échec
#' @keywords internal
.download_open_canopy_zip <- function(dest) {
  zip_url <- "https://github.com/fajwel/Open-Canopy/archive/refs/heads/main.zip"
  tmp_zip <- tempfile(fileext = ".zip")

  tryCatch({
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    message("Téléchargement de l'archive Open-Canopy...")
    download.file(zip_url, tmp_zip, mode = "wb", quiet = TRUE)

    # Extraire dans un dossier temporaire
    tmp_dir <- tempfile("oc_extract_")
    utils::unzip(tmp_zip, exdir = tmp_dir)
    unlink(tmp_zip)

    # GitHub nomme le dossier "Open-Canopy-main"
    extracted <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
    if (length(extracted) == 1) {
      if (dir.exists(dest)) unlink(dest, recursive = TRUE)
      file.rename(extracted[1], dest)
      unlink(tmp_dir, recursive = TRUE)
    }

    if (dir.exists(file.path(dest, "src", "models"))) {
      message("Open-Canopy source téléchargé (ZIP): ", dest)
      return(dest)
    }

    message("Extraction réussie mais structure inattendue dans: ", dest)
    return(NULL)
  }, error = function(e) {
    message("Échec du téléchargement ZIP: ", e$message)
    return(NULL)
  })
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
#' Supporte les deux architectures Open-Canopy :
#' - UNet (SMP, ResNet34) : reconstruction directe via segmentation_models_pytorch
#' - PVTv2 (timm, pvt_v2_b3) : module embarqué ou code source Open-Canopy
#'
#' @param tile SpatRaster (4 bandes : R, G, B, PIR à 1.5m)
#' @param model_path Chemin du modèle .ckpt (PyTorch Lightning)
#' @param model_name "unet" ou "pvtv2" pour la reconstruction
#' @param open_canopy_src Chemin vers le code source Open-Canopy (pour PVTv2)
#' @return SpatRaster CHM prédit (1 bande, en mètres)
predict_tile <- function(tile, model_path, model_name = "pvtv2",
                          open_canopy_src = NULL) {
  library(reticulate)

  # Sauvegarder la tuile en fichier temporaire (chemin long pour Windows)
  tmp_dir <- normalizePath(tempdir(), winslash = "/")
  tmp_in <- tempfile(tmpdir = tmp_dir, fileext = ".tif")
  tmp_out <- tempfile(tmpdir = tmp_dir, fileext = ".tif")
  writeRaster(tile, tmp_in, overwrite = TRUE)

  # Résoudre le chemin modèle (symlinks HF sur Windows)
  model_path <- .resolve_hf_path(model_path)

  # Normaliser les chemins pour Python sous Windows
  tmp_in_py <- gsub("\\\\", "/", tmp_in)
  tmp_out_py <- gsub("\\\\", "/", tmp_out)
  model_path_py <- gsub("\\\\", "/", model_path)

  # Chemin Open-Canopy source (pour PVTv2)
  oc_src_py <- ""
  if (!is.null(open_canopy_src)) {
    oc_src_py <- gsub("\\\\", "/", open_canopy_src)
  }

  py_code <- '
import sys
import os
import torch
import torch.nn as nn
import numpy as np
import rasterio

# ======================================================================
# Charger l image 4 bandes (R, G, B, PIR)
# ======================================================================
with rasterio.open("__INPUT_PATH__") as src:
    image = src.read().astype(np.float32)  # (C, H, W)
    profile = src.profile.copy()

num_bands, H, W = image.shape
print(f"Image chargee: {num_bands} bandes, {H}x{W} px")

# Tensor brut - normalisation appliquee apres lecture du checkpoint
tensor = torch.from_numpy(image).unsqueeze(0)  # (1, C, H, W)
print(f"Tensor: shape={tuple(tensor.shape)}, "
      f"min={tensor.min():.1f}, max={tensor.max():.1f}")

# ======================================================================
# Charger le checkpoint PyTorch Lightning
# ======================================================================
ckpt_path = "__MODEL_PATH__"
oc_src = "__OC_SRC__"
embedded_py = "__EMBEDDED_PY__"

# Ajouter Open-Canopy au path avant chargement (le checkpoint peut
# referencer des classes du module src)
if oc_src and os.path.isdir(oc_src):
    if oc_src not in sys.path:
        sys.path.insert(0, oc_src)
        print(f"Open-Canopy source ajoute au path: {oc_src}")

# Ajouter le module embarque au path (pour PVTv2 sans repo externe)
if embedded_py and os.path.isdir(embedded_py):
    if embedded_py not in sys.path:
        sys.path.insert(0, embedded_py)
        print(f"Module timmNet embarque disponible: {embedded_py}")

import types as _types

class _Dummy:
    """Classe factice pour depickling des references manquantes."""
    def __init__(self, *a, **kw): pass
    def __setstate__(self, state):
        if isinstance(state, dict):
            self.__dict__.update(state)

class _MockMod(_types.ModuleType):
    """Module factice dont les attributs retournent _Dummy."""
    def __getattr__(self, name):
        return _Dummy

class _SrcFinder:
    """Import hook pour simuler le package src (Open-Canopy)."""
    def find_module(self, fullname, path=None):
        if fullname == "src" or fullname.startswith("src."):
            return self
        return None
    def load_module(self, fullname):
        if fullname in sys.modules:
            return sys.modules[fullname]
        mod = _MockMod(fullname)
        mod.__loader__ = self
        mod.__path__ = []
        mod.__package__ = fullname
        sys.modules[fullname] = mod
        return mod

try:
    checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
except (ModuleNotFoundError, ImportError) as e:
    print(f"  Chargement direct echoue: {e}")
    print("  Simulation des modules manquants pour le depickling...")
    _finder = _SrcFinder()
    sys.meta_path.insert(0, _finder)
    try:
        checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    finally:
        sys.meta_path.remove(_finder)
        for _k in list(sys.modules):
            if _k == "src" or _k.startswith("src."):
                del sys.modules[_k]
    print("  Checkpoint charge avec modules simules")

print(f"Checkpoint cles: {list(checkpoint.keys())}")

# ======================================================================
# Extraire mean/std du checkpoint pour normaliser l entree
# ======================================================================
dm_hparams = checkpoint.get("datamodule_hyper_parameters", {})
model_hparams = checkpoint.get("hyper_parameters", {})

print(f"Datamodule hparams: {dm_hparams}")
_nc = model_hparams.get("num_channels", "?")
_mn = dm_hparams.get("mean", "?")
_sd = dm_hparams.get("std", "?")
print(f"Model hparams: num_channels={_nc}, mean={_mn}, std={_sd}")

# Recuperer mean/std pour la normalisation
_ckpt_mean = dm_hparams.get("mean", None)
_ckpt_std = dm_hparams.get("std", None)

# Appliquer la normalisation sur le tensor d entree
if _ckpt_mean is not None and _ckpt_std is not None:
    import ast
    # mean/std peuvent etre des listes (par bande) ou des scalaires
    if isinstance(_ckpt_mean, str):
        try:
            _ckpt_mean = ast.literal_eval(_ckpt_mean)
        except Exception:
            pass
    if isinstance(_ckpt_std, str):
        try:
            _ckpt_std = ast.literal_eval(_ckpt_std)
        except Exception:
            pass

    _mean_t = torch.tensor(_ckpt_mean, dtype=torch.float32)
    _std_t = torch.tensor(_ckpt_std, dtype=torch.float32)

    # Reshape pour broadcasting (1, C, 1, 1)
    if _mean_t.dim() == 1:
        _mean_t = _mean_t.reshape(1, -1, 1, 1)
        _std_t = _std_t.reshape(1, -1, 1, 1)

    # Si les valeurs sont non-triviales (pas 0/1), appliquer la normalisation
    if not (torch.all(_mean_t == 0) and torch.all(_std_t == 1)):
        tensor = (tensor - _mean_t) / _std_t
        print(f"  Normalisation appliquee: mean={_ckpt_mean}, std={_ckpt_std}")
        print(f"  Tensor normalise: min={tensor.min():.3f}, max={tensor.max():.3f}")
    else:
        print("  mean=0, std=1 : pas de normalisation specifique")
        # Verifier si les valeurs sont en 8-bit -> convertir en reflectance
        if tensor.max() > 255:
            tensor = tensor / 10000.0
            print(f"  Conversion reflectance /10000: min={tensor.min():.4f}, max={tensor.max():.4f}")
        elif tensor.max() > 1:
            tensor = tensor / 255.0
            print(f"  Conversion 8-bit /255: min={tensor.min():.4f}, max={tensor.max():.4f}")
else:
    # Pas de mean/std dans le checkpoint : normaliser heuristiquement
    print("  Pas de mean/std dans le checkpoint")
    if tensor.max() > 255:
        tensor = tensor / 10000.0
        print(f"  Conversion reflectance /10000: min={tensor.min():.4f}, max={tensor.max():.4f}")
    elif tensor.max() > 1:
        tensor = tensor / 255.0
        print(f"  Conversion 8-bit /255: min={tensor.min():.4f}, max={tensor.max():.4f}")

state_dict = checkpoint.get("state_dict", checkpoint)
keys = list(state_dict.keys())
print(f"State dict: {len(keys)} parametres")
print(f"  Premieres cles: {keys[:5]}")

# Debug : chercher les cles seg_head dans le checkpoint
_seg_keys = [k for k in keys if "seg_head" in k]
print(f"  Cles seg_head dans checkpoint: {len(_seg_keys)}")
if _seg_keys:
    print(f"    Exemples: {_seg_keys[:5]}")
else:
    print("    AUCUNE cle seg_head trouvee dans le checkpoint!")
    print(f"  Dernieres cles: {keys[-5:]}")

# ======================================================================
# Fonction set_first_layer (reproduction du code Open-Canopy)
# Adapte la premiere couche conv de 3 a N canaux
# ======================================================================
def set_first_layer(model, n_channels):
    if n_channels == 3:
        return
    for module in model.modules():
        if isinstance(module, nn.Conv2d) and module.in_channels == 3:
            break
        if isinstance(module, nn.Linear):
            break
    previous_weight = module.weight.detach()
    if previous_weight.dim() == 4:
        # Conv2d
        n_out = previous_weight.shape[0]
        new_weight = torch.randn(
            n_out, n_channels,
            previous_weight.shape[2], previous_weight.shape[3]
        )
        new_weight[:, :3] = previous_weight
    elif previous_weight.dim() == 2:
        # Linear (ViT-style patch embedding)
        n_out = previous_weight.shape[0]
        n_elem = previous_weight.shape[1] // 3
        new_weight = torch.randn((n_out, n_channels * n_elem))
        new_weight[:, :3 * n_elem] = previous_weight
    module.weight = nn.parameter.Parameter(new_weight)
    if hasattr(module, "in_channels"):
        module.in_channels = n_channels

# ======================================================================
# Nettoyage des cles du state_dict
# ======================================================================
def clean_state_dict(state_dict, prefix_to_strip):
    """Enleve un prefixe des cles du state_dict."""
    clean = {}
    for k, v in state_dict.items():
        new_k = k
        for prefix in prefix_to_strip:
            if new_k.startswith(prefix):
                new_k = new_k[len(prefix):]
                break
        clean[new_k] = v
    return clean

def smart_load_state_dict(model, ckpt_state_dict, prefixes_to_strip):
    """Charge le state_dict en adaptant automatiquement les cles.

    Gere les differences de nommage entre versions de timm :
    - features_only=True peut renommer stages.0 en stages_0
    - Prefixes variables (net., model., net.model.)
    """
    import re

    # Etape 1 : essayer les prefixes dans l ordre
    best_missing = None
    best_dict = None
    best_prefix = None

    # Combinaisons de prefixes a essayer
    prefix_combos = [
        prefixes_to_strip,          # ["net.", "model."] -> strip "net." only
        ["net.model."],             # strip "net.model." together
        ["net."],                   # strip "net." only
        [],                         # no stripping
    ]

    for combo in prefix_combos:
        d = clean_state_dict(ckpt_state_dict, combo) if combo else dict(ckpt_state_dict)
        missing, unexpected = model.load_state_dict(d, strict=False)
        if best_missing is None or len(missing) < len(best_missing):
            best_missing = missing
            best_dict = d
            best_prefix = combo
        if len(missing) == 0:
            print(f"  Prefixe parfait: {combo}")
            return missing, unexpected

    # Etape 2 : construire mapping si beaucoup de cles manquent
    model_keys = set(model.state_dict().keys())
    ckpt_keys = set(best_dict.keys())

    if len(best_missing) > len(model_keys) * 0.1:
        print(f"  Mapping automatique ({len(best_missing)} cles manquantes sur {len(model_keys)})")

        def normalize_key(k):
            """stages_0.blocks_1 -> stages.0.blocks.1"""
            return re.sub(r"_(\\d+)", r".\\1", k)

        # Construire un index inverse : cle normalisee -> cle checkpoint
        ckpt_norm = {normalize_key(k): k for k in ckpt_keys}

        # Combiner matches directs ET normalises
        combined = {}
        matched_direct = 0
        matched_norm = 0

        for model_k in model_keys:
            # 1. Match direct : la cle existe telle quelle dans best_dict
            if model_k in best_dict:
                combined[model_k] = best_dict[model_k]
                matched_direct += 1
            else:
                # 2. Match normalise : stages_0 -> stages.0
                norm_k = normalize_key(model_k)
                if norm_k in ckpt_norm:
                    ckpt_k = ckpt_norm[norm_k]
                    combined[model_k] = best_dict[ckpt_k]
                    matched_norm += 1

        total = matched_direct + matched_norm
        print(f"  {total}/{len(model_keys)} cles trouvees ({matched_direct} directes, {matched_norm} normalisees)")

        if total > len(best_missing) * 0.5:
            missing, unexpected = model.load_state_dict(combined, strict=False)
            print(f"  Apres chargement: {len(missing)} manquantes, {len(unexpected)} inattendues")
            if 0 < len(missing) <= 10:
                for m in missing:
                    print(f"    - {m}")
            return missing, unexpected

    # Etape 3 : utiliser le meilleur resultat
    print(f"  Meilleur prefixe: {best_prefix}, {len(best_missing)} manquantes")
    model.load_state_dict(best_dict, strict=False)
    return best_missing, []

# ======================================================================
# Reconstruire le modele selon l architecture
# ======================================================================
model_name = "__MODEL_NAME__"
model = None

# Detecter le type de modele depuis les cles du checkpoint
has_seg_model = any("seg_model" in k for k in keys)
has_timm_model = any("net.model." in k for k in keys)
has_seg_head = any("net.seg_head." in k for k in keys)

print(f"Detection: seg_model={has_seg_model}, timm_model={has_timm_model}, "
      f"seg_head={has_seg_head}")

# ======================================================================
# UNet (SMP)
# ======================================================================
if has_seg_model or (model_name == "unet" and not has_timm_model):
    import segmentation_models_pytorch as smp

    print("=== Reconstruction: SMP UNet (ResNet34, 4ch, 1 classe) ===")

    seg_model = smp.create_model(
        arch="unet",
        encoder_name="resnet34",
        classes=1,
        in_channels=3,
        encoder_weights=None,
    )
    set_first_layer(seg_model.encoder, num_bands)

    # Cles Lightning : "net.seg_model.encoder.conv1.weight" etc.
    clean_dict = clean_state_dict(state_dict,
        ["net.seg_model.", "model.seg_model.", "net.", "model."])

    missing, unexpected = seg_model.load_state_dict(clean_dict, strict=False)
    if missing:
        print(f"  Cles manquantes: {len(missing)}")
        for m in missing[:3]:
            print(f"    - {m}")
    if unexpected:
        print(f"  Cles inattendues: {len(unexpected)}")

    seg_model.eval()
    # Wrapper pour retourner un tensor (pas un dict)
    model = seg_model
    print("UNet SMP charge avec succes")

# ======================================================================
# PVTv2 (timm + SimpleSegmentationHead)
# ======================================================================
elif has_timm_model or has_seg_head or model_name == "pvtv2":
    print("=== Reconstruction: PVTv2 (timm pvt_v2_b3, 4ch, 1 classe) ===")

    _timmNet = None
    _import_errors = []

    # Strategie 1 : module embarque (aucune dependance externe)
    try:
        from timmnet_standalone import timmNet as _timmNet
        print("Import timmNet depuis module embarque (standalone)")
    except Exception as e:
        _import_errors.append(f"Module embarque: {e}")
        print(f"Import module embarque echoue: {e}")

    # Strategie 2 : code source Open-Canopy (si fourni)
    if _timmNet is None and oc_src and os.path.isdir(oc_src):
        try:
            from src.models.components.timmNet import timmNet as _timmNet
            print("Import timmNet depuis Open-Canopy source")
        except Exception as e:
            _import_errors.append(f"Open-Canopy source: {e}")
            print(f"Import Open-Canopy echoue: {e}")

    if _timmNet is None:
        msg = "Impossible de charger timmNet pour PVTv2.\\n"
        msg += f"  Module embarque: {embedded_py!r}\\n"
        msg += f"  Open-Canopy src: {oc_src!r}\\n"
        for err in _import_errors:
            msg += f"  - {err}\\n"
        msg += "  Verifiez que timm est installe: pip install timm"
        raise RuntimeError(msg)

    # PVTv2 reduit par facteur 32 : img_size doit etre multiple de 32
    _pad_multiple = 32
    _img_size = ((max(H, W) + _pad_multiple - 1) // _pad_multiple) * _pad_multiple
    print(f"  img_size ajuste: {max(H, W)} -> {_img_size} (multiple de {_pad_multiple})")

    pvt_model = _timmNet(
        backbone="pvt_v2_b3.in1k",
        num_classes=1,
        num_channels=num_bands,
        pretrained=False,
        pretrained_path=None,
        img_size=_img_size,
        use_FPN=False,
    )

    # Debug : comparer les cles attendues vs fournies
    model_sample = list(pvt_model.state_dict().keys())[:5]
    ckpt_sample = list(state_dict.keys())[:5]
    print(f"  Cles modele (exemples): {model_sample}")
    print(f"  Cles checkpoint (exemples): {ckpt_sample}")

    missing, unexpected = smart_load_state_dict(
        pvt_model, state_dict, ["net.", "model."])
    if missing:
        print(f"  Cles manquantes finales: {len(missing)}")
        for m in missing[:5]:
            print(f"    - {m}")
    if unexpected:
        print(f"  Cles inattendues finales: {len(unexpected)}")

    pvt_model.eval()
    model = pvt_model
    print("PVTv2 charge avec succes")

# ======================================================================
# Inference
# ======================================================================
if model is None:
    raise RuntimeError("Le modele na pas pu etre charge. Verifiez les logs ci-dessus.")

# Padding a un multiple de 32 pour PVTv2 (ou autre modele par reduction)
_pad = 32
_orig_H, _orig_W = H, W
_pad_H = ((_orig_H + _pad - 1) // _pad) * _pad
_pad_W = ((_orig_W + _pad - 1) // _pad) * _pad

if _pad_H != _orig_H or _pad_W != _orig_W:
    print(f"Padding: {_orig_H}x{_orig_W} -> {_pad_H}x{_pad_W}")
    padded = torch.zeros(1, num_bands, _pad_H, _pad_W, dtype=tensor.dtype)
    padded[:, :, :_orig_H, :_orig_W] = tensor
    tensor = padded

with torch.no_grad():
    output = model(tensor)

    # Le modele retourne {"out": tensor} (Open-Canopy) ou tensor (SMP)
    if isinstance(output, dict):
        pred = output["out"]
    else:
        pred = output

    pred = pred.squeeze().cpu().numpy()

    # Recadrer au dimensions originales (enlever le padding)
    if pred.ndim == 2:
        pred = pred[:_orig_H, :_orig_W]
    elif pred.ndim == 3:
        pred = pred[:, :_orig_H, :_orig_W]

    # Le modele predit en metres (targets = dm / 10)
    pred = np.clip(pred, 0, 50)
    pred = np.round(pred, 1)

    print(f"CHM predit: min={pred.min():.1f}m, max={pred.max():.1f}m, "
          f"mean={pred.mean():.1f}m")

# ======================================================================
# Sauvegarder le resultat
# ======================================================================
profile.update(count=1, dtype="float32", compress="lzw")
with rasterio.open("__OUTPUT_PATH__", "w", **profile) as dst:
    dst.write(pred.astype(np.float32), 1)

print("Prediction sauvegardee")
'
  # Chemin vers le module Python embarqué (timmnet_standalone)
  embedded_py_dir <- system.file("python", package = "opencanopy")
  if (!nzchar(embedded_py_dir)) {
    # Fallback: package non installé, utiliser le chemin source
    embedded_py_dir <- file.path(
      system.file(package = "opencanopy"),
      "python"
    )
    if (!dir.exists(embedded_py_dir)) {
      # Dernier recours : chemin relatif depuis le code source
      pkg_root <- tryCatch(
        rprojroot::find_package_root_file(),
        error = function(e) getwd()
      )
      embedded_py_dir <- file.path(pkg_root, "inst", "python")
    }
  }
  embedded_py_py <- gsub("\\\\", "/", embedded_py_dir)

  # Substitution des placeholders (evite sprintf et ses limites)
  py_code <- gsub("__INPUT_PATH__", tmp_in_py, py_code, fixed = TRUE)
  py_code <- gsub("__MODEL_PATH__", model_path_py, py_code, fixed = TRUE)
  py_code <- gsub("__MODEL_NAME__", model_name, py_code, fixed = TRUE)
  py_code <- gsub("__OC_SRC__", oc_src_py, py_code, fixed = TRUE)
  py_code <- gsub("__EMBEDDED_PY__", embedded_py_py, py_code, fixed = TRUE)
  py_code <- gsub("__OUTPUT_PATH__", tmp_out_py, py_code, fixed = TRUE)

  tryCatch({
    py_run_string(py_code)
    pred_disk <- rast(tmp_out)
    # Forcer la lecture en mémoire AVANT de supprimer le fichier temp
    pred <- setValues(rast(pred_disk), values(pred_disk))
    names(pred) <- "chm_predicted"
    rm(pred_disk)
    return(pred)
  }, error = function(e) {
    stop("Erreur inférence du modèle: ", e$message, call. = FALSE)
  }, finally = {
    unlink(c(tmp_in, tmp_out))
  })
}

#' Combiner les ortho RVB et IRC en image 4 bandes (R, G, B, PIR)
#'
#' Le modèle Open-Canopy attend 4 canaux : Rouge, Vert, Bleu, PIR
#' - RVB fournit les 3 premières bandes
#' - IRC fournit le PIR (bande 1 de l'IRC = Proche Infrarouge)
#'
#' @param rvb SpatRaster ortho RVB (3 bandes : Rouge, Vert, Bleu)
#' @param irc SpatRaster ortho IRC (3 bandes : PIR, Rouge, Vert)
#' @return SpatRaster 4 bandes (Rouge, Vert, Bleu, PIR)
combine_rvb_irc <- function(rvb, irc) {
  message("Combinaison RVB + PIR en image 4 bandes...")

  # Aligner les emprises et résolutions
  if (!compareGeom(rvb, irc, stopOnError = FALSE)) {
    message("  Rééchantillonnage IRC sur la grille RVB...")
    irc <- resample(irc, rvb, method = "bilinear")
  }

  # Extraire PIR (bande 1 de l'IRC)
  pir <- irc[[1]]
  names(pir) <- "PIR"

  # Combiner : R, G, B, PIR
  rgbn <- c(rvb[[1]], rvb[[2]], rvb[[3]], pir)
  names(rgbn) <- c("Rouge", "Vert", "Bleu", "PIR")

  message(sprintf("  Image 4 bandes: %d x %d px, %d bandes",
                   ncol(rgbn), nrow(rgbn), nlyr(rgbn)))
  return(rgbn)
}

#' Pipeline d'inférence complet sur un raster
#'
#' @param rvb SpatRaster ortho RVB (0.20m, 3 bandes)
#' @param irc SpatRaster ortho IRC (0.20m, 3 bandes)
#' @param model_path Chemin du modèle .ckpt
#' @param model_name "unet" ou "pvtv2"
#' @param tile_size Taille des tuiles en mètres
#' @return SpatRaster CHM prédit
run_inference <- function(rvb, irc, model_path, model_name = "pvtv2",
                           tile_size = 1000, open_canopy_src = NULL) {
  message("\n=== Inférence Open-Canopy ===")

  # 1. Combiner RVB + PIR en 4 bandes
  rgbn <- combine_rvb_irc(rvb, irc)

  # 2. Rééchantillonner de 0.20m → 1.5m
  r_1_5m <- resample_to_spot(rgbn)

  # 3. Découper en tuiles
  tiles <- make_inference_tiles(r_1_5m, tile_size = tile_size)

  # 4. Prédire chaque tuile
  predictions <- list()
  for (i in seq_along(tiles)) {
    tile_name <- names(tiles)[i]
    message(sprintf("  Inférence tuile %d/%d: %s", i, length(tiles), tile_name))
    pred <- predict_tile(tiles[[i]], model_path, model_name, open_canopy_src)
    if (!is.null(pred)) {
      predictions[[tile_name]] <- pred
    }
  }

  if (length(predictions) == 0) {
    stop("Aucune prédiction réussie.")
  }

  # 5. Mosaïquer
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
#' @param model_path Chemin local vers un checkpoint .ckpt (optionnel,
#'   sinon téléchargé depuis HuggingFace)
#' @param open_canopy_src Chemin vers le code source Open-Canopy (optionnel,
#'   le module embarqué timmnet_standalone est utilisé par défaut)
#' @param res_m Résolution de téléchargement IGN (0.2m par défaut)
#' @param millesime_ortho Millésime ortho RVB (NULL = plus récent)
#' @param millesime_irc Millésime IRC (NULL = plus récent)
#' @return Liste avec tous les résultats
pipeline_aoi_to_chm <- function(aoi_path,
                                  output_dir = file.path(getwd(), "outputs"),
                                  model_name = "pvtv2",
                                  model_path = NULL,
                                  open_canopy_src = NULL,
                                  res_m = RES_IGN,
                                  millesime_ortho = MILLESIME_ORTHO,
                                  millesime_irc = MILLESIME_IRC) {
  dir_create(output_dir)
  t0 <- Sys.time()

  # Contournement Windows : les chemins courts 8.3 (ex. PASCAL~1.OBS)
  # provoquent des erreurs terra/GDAL. On normalise le tempdir.
  if (.Platform$OS.type == "windows") {
    long_tmpdir <- normalizePath(tempdir(), winslash = "/")
    terraOptions(tempdir = long_tmpdir)
  }

  message("##############################################################")
  message("#  Pipeline Open-Canopy : AOI → Ortho IGN → CHM prédit       #")
  message("##############################################################\n")

  # --- Étape 1 : Charger l'AOI ---
  message(">>> ÉTAPE 1/5 : Chargement de l'AOI")
  aoi <- load_aoi(aoi_path)

  # --- Étape 2 : Télécharger les ortho IGN ---
  message("\n>>> ÉTAPE 2/5 : Téléchargement des ortho IGN (RVB + IRC)")
  ortho <- download_ortho_for_aoi(aoi, output_dir = output_dir, res_m = res_m,
                                   millesime_ortho = millesime_ortho,
                                   millesime_irc = millesime_irc)

  # --- Étape 3 : Configurer Python ---
  message("\n>>> ÉTAPE 3/5 : Configuration Python + téléchargement modèle")
  setup_python()
  if (is.null(model_path)) {
    model_path <- download_model(model_name)
  } else {
    message("Utilisation du modèle local: ", model_path)
    if (!file.exists(model_path)) {
      stop("Fichier modèle introuvable: ", model_path)
    }
  }

  # --- Étape 4 : Inférence ---
  message("\n>>> ÉTAPE 4/5 : Inférence du modèle ", model_name)

  # Auto-détecter ou télécharger le code source Open-Canopy (pour PVTv2)
  if (is.null(open_canopy_src) && model_name == "pvtv2") {
    # 1. Chercher dans les emplacements courants
    candidates <- c(
      file.path(getwd(), "Open-Canopy"),
      file.path(dirname(getwd()), "Open-Canopy"),
      file.path(Sys.getenv("USERPROFILE"), "dev", "Open-Canopy"),
      file.path(Sys.getenv("HOME"), "dev", "Open-Canopy"),
      file.path(tools::R_user_dir("opencanopy", "cache"), "Open-Canopy")
    )
    for (cand in candidates) {
      if (dir.exists(file.path(cand, "src", "models"))) {
        open_canopy_src <- cand
        message("Open-Canopy source détecté: ", open_canopy_src)
        break
      }
    }
    # 2. Sinon, télécharger automatiquement
    if (is.null(open_canopy_src)) {
      message("Open-Canopy source non trouvé localement, téléchargement...")
      open_canopy_src <- download_open_canopy_src()
    }
  }

  # Combiner RVB + IRC (PIR) en 4 bandes, puis inférence
  chm <- run_inference(ortho$rvb, ortho$irc, model_path, model_name,
                        open_canopy_src = open_canopy_src)

  # --- Étape 5 : Export ---
  message("\n>>> ÉTAPE 5/5 : Export des résultats")

  # CHM à la résolution du modèle (1.5m)
  chm_path <- file.path(output_dir, "chm_predicted_1_5m.tif")
  if (file.exists(chm_path)) file.remove(chm_path)
  writeRaster(chm, chm_path, gdal = c("COMPRESS=LZW"))
  message("CHM 1.5m: ", chm_path)

  # Suréchantillonner le CHM vers la résolution IGN (0.20m)
  message("Suréchantillonnage CHM vers 0.20m...")
  disagg_factor <- round(RES_SPOT / RES_IGN)
  chm_hr_path <- file.path(output_dir, "chm_predicted_0_2m.tif")
  if (file.exists(chm_hr_path)) file.remove(chm_hr_path)
  chm_hr <- disagg(chm, fact = disagg_factor, method = "bilinear")
  chm_hr <- crop(chm_hr, vect(st_union(aoi)),
                 filename = chm_hr_path, overwrite = TRUE,
                 gdal = c("COMPRESS=LZW"))
  message("CHM 0.2m: ", chm_hr_path)

  # --- NDVI si IRC disponible ---
  pir   <- ortho$irc[["PIR"]]
  rouge <- ortho$irc[["Rouge"]]
  ndvi  <- (pir - rouge) / (pir + rouge)
  names(ndvi) <- "NDVI"
  ndvi_path <- file.path(output_dir, "ndvi.tif")
  if (file.exists(ndvi_path)) file.remove(ndvi_path)
  writeRaster(ndvi, ndvi_path, overwrite = TRUE,
              gdal = c("COMPRESS=LZW"))
  message("NDVI:     ", ndvi_path)

  # --- Visualisation récapitulative ---
  pdf_path <- file.path(output_dir, "resultats_aoi.pdf")
  pdf(pdf_path, width = 16, height = 12)

  par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))

  # RVB
  label_ortho <- if (is.null(ortho$millesime_ortho)) "plus récent" else ortho$millesime_ortho
  label_irc   <- if (is.null(ortho$millesime_irc))   "plus récent" else ortho$millesime_irc
  plotRGB(ortho$rvb, r = 1, g = 2, b = 3, stretch = "lin",
          main = sprintf("Ortho RVB IGN 0.20m (millésime: %s)", label_ortho))

  # IRC fausses couleurs
  plotRGB(ortho$irc, r = 1, g = 2, b = 3, stretch = "lin",
          main = sprintf("Ortho IRC fausses couleurs 0.20m (millésime: %s)", label_irc))

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

  # --- Visualisation interactive (RStudio) avec ggplot2 + patchwork ---
  p_combined <- NULL
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("patchwork", quietly = TRUE) &&
      requireNamespace("tidyterra", quietly = TRUE)) {

    library(ggplot2)
    library(patchwork)
    library(tidyterra)

    # Panel 1 : Ortho RVB
    p_rvb <- ggplot() +
      geom_spatraster_rgb(data = ortho$rvb, r = 1, g = 2, b = 3,
                          max_col_value = 255) +
      labs(title = sprintf("Ortho RVB 0.20m (%s)", label_ortho)) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))

    # Panel 2 : Ortho IRC fausses couleurs
    p_irc <- ggplot() +
      geom_spatraster_rgb(data = ortho$irc, r = 1, g = 2, b = 3,
                          max_col_value = 255) +
      labs(title = sprintf("Ortho IRC 0.20m (%s)", label_irc)) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))

    # Panel 3 : NDVI
    p_ndvi <- ggplot() +
      geom_spatraster(data = ndvi) +
      scale_fill_gradientn(
        colours = c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
                    "#d9ef8b", "#91cf60", "#1a9850", "#006837"),
        limits = c(-0.2, 1), na.value = "transparent",
        name = "NDVI"
      ) +
      labs(title = "NDVI (depuis IRC)") +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))

    # Panel 4 : CHM
    p_chm <- ggplot() +
      geom_spatraster(data = chm) +
      scale_fill_gradientn(
        colours = c("#f7fcb9", "#addd8e", "#41ab5d", "#006837", "#004529"),
        na.value = "transparent",
        name = "Hauteur (m)"
      ) +
      labs(title = paste("CHM", model_name)) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 10, face = "bold"))

    p_combined <- (p_rvb | p_irc) / (p_ndvi | p_chm) +
      plot_annotation(
        title = "Pipeline Open-Canopy : ortho IGN \u2192 CHM",
        subtitle = sprintf("AOI: %s | CHM: min=%.1fm, max=%.1fm, moy=%.1fm",
                           basename(aoi_path),
                           min(values(chm, na.rm = TRUE)),
                           max(values(chm, na.rm = TRUE)),
                           mean(values(chm, na.rm = TRUE))),
        theme = theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey40")
        )
      )

    print(p_combined)
    message("Plot patchwork affich\u00e9 dans RStudio")
  }

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
    millesime_ortho = ortho$millesime_ortho,
    millesime_irc   = ortho$millesime_irc,
    layer_ortho     = ortho$layer_ortho,
    layer_irc       = ortho$layer_irc,
    ndvi      = ndvi,
    chm_1_5m  = chm,
    chm_0_2m  = chm_hr,
    plot      = p_combined,
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
    message("  # Avec un checkpoint local :")
    message('  result <- pipeline_aoi_to_chm("data/aoi.gpkg",')
    message('    model_path = "chemin/vers/checkpoint.ckpt")')
    message("")
    message("Le fichier aoi.gpkg doit contenir un polygone définissant")
    message("votre zone d'intérêt (n'importe quel CRS, sera reprojeté")
    message("en Lambert-93 automatiquement).")
  } else {
    result <- pipeline_aoi_to_chm(aoi_path)
  }
}
