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
#   - Entrée : 4 bandes (R, G, B, NIR) sans normalisation (mean=0, std=1)
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
# Chemin vers le code source Open-Canopy (nécessaire pour PVTv2)
OPEN_CANOPY_SRC <- NULL  # Sera détecté automatiquement ou passé en paramètre

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
setup_python <- function() {
  library(reticulate)

  # Vérifier les modules nécessaires
  modules <- c("torch", "numpy", "rasterio", "huggingface_hub",
               "segmentation_models_pytorch", "timm")
  ok <- TRUE
  for (mod in modules) {
    avail <- py_module_available(mod)
    message(sprintf("  Python %s: %s", mod, ifelse(avail, "OK", "MANQUANT")))
    if (!avail) ok <- FALSE
  }

  if (!ok) {
    stop("Modules Python manquants. Installez-les dans l'env '", CONDA_ENV, "':\n",
         "  conda activate ", CONDA_ENV, "\n",
         "  pip install torch torchvision numpy rasterio huggingface_hub ",
         "segmentation-models-pytorch timm")
  }
}

#' Télécharger le modèle pré-entraîné depuis Hugging Face
#'
#' Découvre dynamiquement les fichiers .ckpt dans pretrained_models/
#' du dataset AI4Forest/Open-Canopy et télécharge celui correspondant
#' au modèle demandé.
#'
#' @param model_name "unet" ou "pvtv2"
#' @return Chemin local du modèle
download_model <- function(model_name = "unet") {
  library(reticulate)
  hf_hub <- import("huggingface_hub")

  # Vérifier que le token HuggingFace est configuré (dataset gated)
  token <- Sys.getenv("HF_TOKEN", unset = "")
  if (token == "") {
    # Vérifier si huggingface-cli login a été fait
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

  message("Recherche des checkpoints dans le dataset HuggingFace...")

  # Lister dynamiquement les fichiers .ckpt dans pretrained_models/
  tryCatch({
    api <- hf_hub$HfApi()

    # Utiliser list_repo_files (retourne des chaînes) — plus robuste que list_repo_tree
    all_files <- tryCatch({
      files <- api$list_repo_files(
        HF_REPO_ID,
        repo_type = "dataset"
      )
      reticulate::iterate(files)
    }, error = function(e) {
      # Fallback : list_repo_tree avec gestion des attributs
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

    # Filtrer les .ckpt dans pretrained_models/
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

    # Trouver le checkpoint correspondant au modèle demandé
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
    message("Modèle: ", local_path)
    return(local_path)

  }, error = function(e) {
    message("Impossible de lister le dataset HuggingFace: ", e$message)
    message("Le dataset est peut-être privé (gated). Vérifiez votre HF_TOKEN.")
    message("")
    message("Alternatives :")
    message("  1. Définir votre token : Sys.setenv(HF_TOKEN = 'hf_...')")
    message("  2. Se connecter via CLI : huggingface-cli login")
    message("  3. Fournir le chemin manuellement :")
    message('     result <- pipeline_aoi_to_chm("data/aoi.gpkg",')
    message('       model_path = "chemin/vers/checkpoint.ckpt")')
    stop("Échec du téléchargement du modèle.", call. = FALSE)
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
#' - PVTv2 (timm, pvt_v2_b3) : nécessite le code source Open-Canopy
#'
#' @param tile SpatRaster (4 bandes : R, G, B, PIR à 1.5m)
#' @param model_path Chemin du modèle .ckpt (PyTorch Lightning)
#' @param model_name "unet" ou "pvtv2" pour la reconstruction
#' @param open_canopy_src Chemin vers le code source Open-Canopy (pour PVTv2)
#' @return SpatRaster CHM prédit (1 bande, en mètres)
predict_tile <- function(tile, model_path, model_name = "unet",
                          open_canopy_src = NULL) {
  library(reticulate)

  # Sauvegarder la tuile en fichier temporaire (chemin long pour Windows)
  tmp_dir <- normalizePath(tempdir(), winslash = "/")
  tmp_in <- tempfile(tmpdir = tmp_dir, fileext = ".tif")
  tmp_out <- tempfile(tmpdir = tmp_dir, fileext = ".tif")
  writeRaster(tile, tmp_in, overwrite = TRUE)

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

# Open-Canopy utilise mean=0, std=1 : pas de normalisation
tensor = torch.from_numpy(image).unsqueeze(0)  # (1, C, H, W)
print(f"Tensor: shape={tuple(tensor.shape)}, "
      f"min={tensor.min():.1f}, max={tensor.max():.1f}")

# ======================================================================
# Charger le checkpoint PyTorch Lightning
# ======================================================================
ckpt_path = "__MODEL_PATH__"
oc_src = "__OC_SRC__"

# Ajouter Open-Canopy au path avant chargement (le checkpoint peut
# referencer des classes du module src)
if oc_src and os.path.isdir(oc_src):
    if oc_src not in sys.path:
        sys.path.insert(0, oc_src)
        print(f"Open-Canopy source ajoute au path: {oc_src}")

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

state_dict = checkpoint.get("state_dict", checkpoint)
keys = list(state_dict.keys())
print(f"State dict: {len(keys)} parametres")
print(f"  Premieres cles: {keys[:5]}")

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

    if oc_src and os.path.isdir(oc_src):
        try:
            from src.models.components.timmNet import timmNet
            print("Import timmNet depuis Open-Canopy reussi")

            # Reconstruire le modele PVTv2 avec les memes parametres
            # que configs/model/PVTv2_B.yaml + _seg_default.yaml
            pvt_model = timmNet(
                backbone="pvt_v2_b3.in1k",
                num_classes=1,
                num_channels=num_bands,
                pretrained=False,       # pas besoin, on charge le checkpoint
                pretrained_path=None,
                img_size=max(H, W),     # taille de l image d entree
                use_FPN=False,
            )

            # Cles Lightning : "net.model.xxx" et "net.seg_head.xxx"
            clean_dict = clean_state_dict(state_dict, ["net.", "model."])

            missing, unexpected = pvt_model.load_state_dict(
                clean_dict, strict=False)
            if missing:
                print(f"  Cles manquantes: {len(missing)}")
                for m in missing[:3]:
                    print(f"    - {m}")
            if unexpected:
                print(f"  Cles inattendues: {len(unexpected)}")

            pvt_model.eval()
            model = pvt_model
            print("PVTv2 charge avec succes")

        except ImportError as e:
            print(f"Erreur import Open-Canopy: {e}")
            print("Verifiez que le code source est complet.")
        except Exception as e:
            print(f"Erreur reconstruction PVTv2: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("ERREUR: Le modele PVTv2 necessite le code source Open-Canopy.")
        print(f"  Chemin fourni: {oc_src!r}")
        print("  Utilisez le parametre open_canopy_src dans pipeline_aoi_to_chm()")
        print("  Exemple:")
        print(\'    pipeline_aoi_to_chm("aoi.gpkg", model_name="pvtv2",\')
        print(\'      open_canopy_src="C:/Users/.../Open-Canopy")\')

# ======================================================================
# Inference
# ======================================================================
with torch.no_grad():
    if model is not None:
        output = model(tensor)

        # Le modele retourne {"out": tensor} (Open-Canopy) ou tensor (SMP)
        if isinstance(output, dict):
            pred = output["out"]
        else:
            pred = output

        pred = pred.squeeze().cpu().numpy()

        # Le modele predit en metres (targets = dm / 10)
        pred = np.clip(pred, 0, 50)
        pred = np.round(pred, 1)

        print(f"CHM predit: min={pred.min():.1f}m, max={pred.max():.1f}m, "
              f"mean={pred.mean():.1f}m")
    else:
        print("ATTENTION: Modele non charge, fallback estimation NDVI-based")
        img = tensor.squeeze().numpy()
        if num_bands >= 4:
            pir = img[3]
            rouge = img[0]
            ndvi = (pir - rouge) / (pir + rouge + 1e-6)
            pred = np.clip(ndvi * 25, 0, 40)
        else:
            greenness = img[1] / (img.mean(axis=0) + 1e-6)
            pred = np.clip(greenness * 20, 0, 40)

# ======================================================================
# Sauvegarder le resultat
# ======================================================================
profile.update(count=1, dtype="float32", compress="lzw")
with rasterio.open("__OUTPUT_PATH__", "w", **profile) as dst:
    dst.write(pred.astype(np.float32), 1)

print("Prediction sauvegardee")
'
  # Substitution des placeholders (evite sprintf et ses limites)
  py_code <- gsub("__INPUT_PATH__", tmp_in_py, py_code, fixed = TRUE)
  py_code <- gsub("__MODEL_PATH__", model_path_py, py_code, fixed = TRUE)
  py_code <- gsub("__MODEL_NAME__", model_name, py_code, fixed = TRUE)
  py_code <- gsub("__OC_SRC__", oc_src_py, py_code, fixed = TRUE)
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
    warning("Erreur inférence: ", e$message)
    message("  Utilisation d'une estimation NDVI alternative...")
    if (nlyr(tile) >= 4) {
      pir <- tile[[4]]
      rouge <- tile[[1]]
      ndvi <- (pir - rouge) / (pir + rouge + 0.001)
      pred <- clamp(ndvi * 25, lower = 0, upper = 40)
      names(pred) <- "chm_estimated"
      return(pred)
    }
    return(NULL)
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
run_inference <- function(rvb, irc, model_path, model_name = "unet",
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
#' @param open_canopy_src Chemin vers le code source Open-Canopy (nécessaire
#'   pour PVTv2, auto-détecté sinon)
#' @param res_m Résolution de téléchargement IGN (0.2m par défaut)
#' @param millesime_ortho Millésime ortho RVB (NULL = plus récent)
#' @param millesime_irc Millésime IRC (NULL = plus récent)
#' @return Liste avec tous les résultats
pipeline_aoi_to_chm <- function(aoi_path,
                                  output_dir = file.path(getwd(), "outputs"),
                                  model_name = "unet",
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

  # Auto-détecter le code source Open-Canopy si nécessaire (pour PVTv2)
  if (is.null(open_canopy_src) && model_name == "pvtv2") {
    # Chercher dans les emplacements courants
    candidates <- c(
      file.path(getwd(), "Open-Canopy"),
      file.path(dirname(getwd()), "Open-Canopy"),
      file.path(Sys.getenv("USERPROFILE"), "dev", "Open-Canopy"),
      file.path(Sys.getenv("HOME"), "dev", "Open-Canopy")
    )
    for (cand in candidates) {
      if (dir.exists(file.path(cand, "src", "models"))) {
        open_canopy_src <- cand
        message("Open-Canopy source détecté: ", open_canopy_src)
        break
      }
    }
    if (is.null(open_canopy_src)) {
      stop("Le modèle PVTv2 nécessite le code source Open-Canopy.\n",
           "Clonez le dépôt: git clone https://github.com/fajwel/Open-Canopy\n",
           "Puis passez le chemin: open_canopy_src = 'chemin/vers/Open-Canopy'")
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
