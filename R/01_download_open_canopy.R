#!/usr/bin/env Rscript
# ==============================================================================
# 01_download_open_canopy.R
# Téléchargement du dataset Open-Canopy depuis Hugging Face
# + Chargement des données IGN (BD ORTHO® RVB et IRC à 0.20 m)
#
# Dataset HF : AI4Forest/Open-Canopy
# https://huggingface.co/datasets/AI4Forest/Open-Canopy
#
# Données IGN : BD ORTHO® 0.20 m (RVB + IRC)
# https://geoservices.ign.fr/bdortho
# ==============================================================================

# --- Installation des packages nécessaires ---
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

install_if_missing("httr2")
install_if_missing("jsonlite")
install_if_missing("terra")
install_if_missing("sf")
install_if_missing("curl")
install_if_missing("fs")

library(httr2)
library(jsonlite)
library(curl)
library(fs)

# ==============================================================================
# Configuration
# ==============================================================================

# --- Hugging Face ---
HF_REPO_ID <- "AI4Forest/Open-Canopy"
HF_API_URL <- "https://huggingface.co/api/datasets"
HF_TOKEN <- Sys.getenv("HF_TOKEN", unset = "")

# --- IGN Géoplateforme ---
# WMS/WMTS pour les ortho IGN à 0.20 m
IGN_WMTS_URL <- "https://data.geopf.fr/wmts"
IGN_WMS_URL  <- "https://data.geopf.fr/wms-r"
# Couches disponibles
IGN_LAYER_ORTHO <- "ORTHOIMAGERY.ORTHOPHOTOS"
IGN_LAYER_IRC   <- "ORTHOIMAGERY.ORTHOPHOTOS.IRC-EXPRESS.2024"

# --- Répertoires ---
DATA_DIR     <- file.path(getwd(), "data")
DATA_DIR_HF  <- file.path(DATA_DIR, "open_canopy")
DATA_DIR_IGN <- file.path(DATA_DIR, "ign")
dir_create(c(DATA_DIR, DATA_DIR_HF, DATA_DIR_IGN))

# --- Résolutions ---
RES_SPOT <- 1.5  # Résolution SPOT 6-7 (Open-Canopy)
RES_IGN  <- 0.2  # Résolution BD ORTHO® IGN

# ==============================================================================
# A. Fonctions Hugging Face (Open-Canopy)
# ==============================================================================

#' Lister les fichiers du dataset sur Hugging Face
#'
#' @param repo_id Identifiant du dépôt (ex: "AI4Forest/Open-Canopy")
#' @param path Chemin dans le dépôt (ex: "" pour la racine)
#' @param token Token Hugging Face (optionnel)
#' @return data.frame avec les informations des fichiers
hf_list_files <- function(repo_id, path = "", token = "") {
  url <- paste0(HF_API_URL, "/", repo_id, "/tree/main")
  if (nchar(path) > 0) {
    url <- paste0(url, "/", path)
  }

  req <- request(url)
  if (nchar(token) > 0) {
    req <- req |> req_headers(Authorization = paste("Bearer", token))
  }

  resp <- req |>
    req_error(is_error = function(resp) FALSE) |>
    req_perform()

  if (resp_status(resp) != 200) {
    warning("Erreur API HF: ", resp_status(resp))
    return(data.frame())
  }

  content <- resp |> resp_body_json()

  files_df <- do.call(rbind, lapply(content, function(item) {
    data.frame(
      type = item$type %||% NA,
      path = item$path %||% NA,
      size = item$size %||% NA,
      oid = item$oid %||% NA,
      stringsAsFactors = FALSE
    )
  }))

  return(files_df)
}

#' Télécharger un fichier depuis Hugging Face
#'
#' @param repo_id Identifiant du dépôt
#' @param filename Chemin du fichier dans le dépôt
#' @param dest_dir Répertoire de destination local
#' @param token Token Hugging Face (optionnel)
#' @param overwrite Écraser si le fichier existe déjà
#' @return Chemin local du fichier téléchargé
hf_download_file <- function(repo_id, filename, dest_dir, token = "",
                              overwrite = FALSE) {
  local_path <- file.path(dest_dir, filename)

  if (file.exists(local_path) && !overwrite) {
    message("Fichier déjà présent: ", local_path)
    return(local_path)
  }

  dir_create(dirname(local_path))

  url <- paste0(
    "https://huggingface.co/datasets/", repo_id,
    "/resolve/main/", filename
  )

  message("Téléchargement HF: ", filename)

  headers <- list()
  if (nchar(token) > 0) {
    headers[["Authorization"]] <- paste("Bearer", token)
  }

  tryCatch({
    curl_download(
      url = url,
      destfile = local_path,
      handle = new_handle(.list = headers),
      quiet = FALSE
    )
    message("OK: ", local_path)
    return(local_path)
  }, error = function(e) {
    warning("Échec du téléchargement: ", filename, " - ", e$message)
    return(NULL)
  })
}

#' Télécharger un ensemble de fichiers depuis Hugging Face
hf_download_files <- function(repo_id, file_list, dest_dir, token = "",
                               overwrite = FALSE) {
  paths <- character(length(file_list))
  for (i in seq_along(file_list)) {
    paths[i] <- hf_download_file(
      repo_id = repo_id,
      filename = file_list[i],
      dest_dir = dest_dir,
      token = token,
      overwrite = overwrite
    )
  }
  return(paths)
}

#' Télécharger un sous-ensemble du dataset Open-Canopy
#'
#' @param split "train", "val", ou "test"
#' @param n_tiles Nombre de tuiles
#' @param data_type "images" (SPOT), "lidar" (CHM), "lidar_v2", ou "all"
#' @param dest_dir Répertoire de destination
#' @param token Token Hugging Face
download_open_canopy_subset <- function(split = "test",
                                         n_tiles = 5,
                                         data_type = "all",
                                         dest_dir = DATA_DIR_HF,
                                         token = HF_TOKEN) {
  message("=== Téléchargement Open-Canopy (SPOT 1.5m) ===")
  message("Split: ", split, " | Tuiles: ", n_tiles, " | Type: ", data_type)

  files <- hf_list_files(HF_REPO_ID, path = split, token = token)

  if (nrow(files) == 0) {
    for (subdir in c("images", "lidar")) {
      path <- paste0(split, "/", subdir)
      sub_files <- hf_list_files(HF_REPO_ID, path = path, token = token)
      if (nrow(sub_files) > 0) {
        files <- rbind(files, sub_files)
      }
    }
  }

  if (nrow(files) == 0) {
    message("Aucun fichier trouvé pour le split '", split, "'.")
    return(invisible(NULL))
  }

  if (data_type == "images") {
    files <- files[grep("images|spot", files$path, ignore.case = TRUE), ]
  } else if (data_type == "lidar") {
    files <- files[grep("lidar(?!_v2)", files$path, perl = TRUE), ]
  } else if (data_type == "lidar_v2") {
    files <- files[grep("lidar_v2", files$path), ]
  }

  tif_files <- files[grep("\\.tif$", files$path), ]
  if (nrow(tif_files) > n_tiles) {
    tif_files <- tif_files[seq_len(n_tiles), ]
  }

  message(sprintf("\n%d fichier(s) à télécharger.", nrow(tif_files)))

  downloaded <- hf_download_files(
    repo_id = HF_REPO_ID,
    file_list = tif_files$path,
    dest_dir = dest_dir,
    token = token
  )

  message("\n=== Téléchargement terminé ===")
  message(sum(!is.na(downloaded)), " fichier(s) téléchargé(s).")
  return(downloaded)
}

#' Cloner le dataset complet via git (nécessite git-lfs)
hf_git_clone <- function(dest_dir = file.path(DATA_DIR_HF, "Open-Canopy"),
                          token = HF_TOKEN) {
  if (dir.exists(dest_dir)) {
    message("Le répertoire existe déjà: ", dest_dir)
    return(invisible(dest_dir))
  }

  lfs_check <- system("git lfs version", intern = TRUE, ignore.stderr = TRUE)
  if (length(lfs_check) == 0) {
    stop("git-lfs n'est pas installé. ",
         "Installez-le avec: sudo apt install git-lfs")
  }

  message("Clonage du dataset Open-Canopy (~360 Go)...")

  clone_url <- if (nchar(token) > 0) {
    paste0("https://", token, "@huggingface.co/datasets/", HF_REPO_ID)
  } else {
    paste0("https://huggingface.co/datasets/", HF_REPO_ID)
  }

  system2("git", args = c("clone", clone_url, dest_dir))
  message("Clonage terminé: ", dest_dir)
  return(invisible(dest_dir))
}

# ==============================================================================
# B. Fonctions IGN — BD ORTHO® RVB et IRC à 0.20 m
# ==============================================================================

#' Charger une image ortho IGN locale (JPEG2000 ou GeoTIFF)
#'
#' La BD ORTHO® IGN est distribuée en JPEG2000 (.jp2) depuis 2016.
#' Le package terra supporte nativement ce format via GDAL.
#'
#' @param file_path Chemin vers le fichier .jp2 ou .tif
#' @param type "rvb" ou "irc"
#' @return SpatRaster
load_ign_ortho <- function(file_path, type = "rvb") {
  if (!file.exists(file_path)) {
    stop("Fichier introuvable: ", file_path)
  }

  r <- terra::rast(file_path)

  # Nommer les bandes selon le type
  if (type == "rvb" && nlyr(r) >= 3) {
    names(r)[1:3] <- c("Rouge", "Vert", "Bleu")
  } else if (type == "irc" && nlyr(r) >= 3) {
    names(r)[1:3] <- c("PIR", "Rouge", "Vert")
  }

  message(sprintf("Ortho IGN %s chargée: %s", toupper(type), basename(file_path)))
  message(sprintf("  Dimensions: %d x %d pixels", nrow(r), ncol(r)))
  message(sprintf("  Bandes: %d (%s)", nlyr(r), paste(names(r), collapse = ", ")))
  message(sprintf("  Résolution: %.2f x %.2f m", res(r)[1], res(r)[2]))
  message(sprintf("  CRS: %s", crs(r, describe = TRUE)$name))
  message(sprintf("  Emprise: %.0f - %.0f E, %.0f - %.0f N",
                   ext(r)[1], ext(r)[2], ext(r)[3], ext(r)[4]))

  return(r)
}

#' Scanner un répertoire pour trouver les ortho IGN
#'
#' @param dir_path Répertoire à scanner
#' @param recursive Recherche récursive
#' @return data.frame avec les fichiers trouvés et leur type
scan_ign_files <- function(dir_path = DATA_DIR_IGN, recursive = TRUE) {
  extensions <- c("*.jp2", "*.tif", "*.tiff")
  files <- unlist(lapply(extensions, function(ext) {
    dir_ls(dir_path, recurse = recursive, glob = ext)
  }))

  if (length(files) == 0) {
    message("Aucun fichier image trouvé dans: ", dir_path)
    return(data.frame())
  }

  # Détecter le type (RVB ou IRC) à partir du nom de fichier
  df <- data.frame(
    path = files,
    filename = basename(files),
    type = ifelse(grepl("IRC|irc|infrarouge|nir", basename(files),
                         ignore.case = TRUE), "irc", "rvb"),
    size_mb = file.size(files) / 1024^2,
    stringsAsFactors = FALSE
  )

  message(sprintf("Trouvé %d fichier(s) ortho IGN:", nrow(df)))
  message(sprintf("  RVB: %d", sum(df$type == "rvb")))
  message(sprintf("  IRC: %d", sum(df$type == "irc")))

  return(df)
}

#' Télécharger une ortho IGN via WMS (pour une emprise donnée)
#'
#' @param bbox Emprise c(xmin, ymin, xmax, ymax) en Lambert-93 (EPSG:2154)
#' @param layer Couche WMS (RVB ou IRC)
#' @param res_m Résolution en mètres (0.2 par défaut)
#' @param dest_file Chemin du fichier de sortie
#' @return Chemin local du fichier téléchargé
download_ign_wms <- function(bbox, layer = IGN_LAYER_ORTHO, res_m = RES_IGN,
                              dest_file = NULL) {
  xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]

  width  <- round((xmax - xmin) / res_m)
  height <- round((ymax - ymin) / res_m)

  # Limiter à 4096 px (limite WMS typique)
  max_px <- 4096
  if (width > max_px || height > max_px) {
    message("Emprise trop grande pour un seul appel WMS. ",
            "Découpage en tuiles recommandé.")
    scale_factor <- max_px / max(width, height)
    width  <- round(width * scale_factor)
    height <- round(height * scale_factor)
    message(sprintf("  Résolution dégradée: %.2f m",
                     (xmax - xmin) / width))
  }

  if (is.null(dest_file)) {
    layer_short <- ifelse(grepl("IRC", layer), "irc", "rvb")
    dest_file <- file.path(DATA_DIR_IGN,
                            sprintf("ign_%s_%.0f_%.0f.tif",
                                    layer_short, xmin, ymin))
  }

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

  message(sprintf("Téléchargement IGN WMS (%s): %dx%d px...", layer, width, height))

  tryCatch({
    curl_download(url = wms_url, destfile = dest_file, quiet = FALSE)
    message("OK: ", dest_file)

    # Assigner le CRS si nécessaire
    r <- terra::rast(dest_file)
    if (is.na(crs(r)) || crs(r) == "") {
      crs(r) <- "EPSG:2154"
      terra::writeRaster(r, dest_file, overwrite = TRUE)
    }

    return(dest_file)
  }, error = function(e) {
    warning("Échec du téléchargement WMS: ", e$message)
    return(NULL)
  })
}

#' Télécharger les ortho RVB et IRC pour une même emprise
#'
#' @param bbox Emprise c(xmin, ymin, xmax, ymax) en Lambert-93 (EPSG:2154)
#' @param res_m Résolution en mètres
#' @return Liste avec les chemins des fichiers RVB et IRC
download_ign_ortho_pair <- function(bbox, res_m = RES_IGN) {
  message("=== Téléchargement ortho IGN RVB + IRC ===")
  message(sprintf("Emprise: %.0f, %.0f, %.0f, %.0f (Lambert-93)",
                   bbox[1], bbox[2], bbox[3], bbox[4]))
  message(sprintf("Résolution: %.2f m", res_m))

  rvb_path <- download_ign_wms(bbox, layer = IGN_LAYER_ORTHO, res_m = res_m)
  irc_path <- download_ign_wms(bbox, layer = IGN_LAYER_IRC, res_m = res_m)

  return(list(rvb = rvb_path, irc = irc_path))
}

# ==============================================================================
# C. Rééchantillonnage IGN 0.2m → 1.5m (pour compatibilité Open-Canopy)
# ==============================================================================

#' Rééchantillonner une ortho IGN de 0.2m vers 1.5m
#'
#' Pour utiliser les modèles Open-Canopy (entraînés sur SPOT 1.5m),
#' il faut rééchantillonner les images IGN à la même résolution.
#'
#' @param ign_raster SpatRaster de l'ortho IGN (0.2m)
#' @param target_res Résolution cible en mètres (1.5 par défaut)
#' @param method Méthode de rééchantillonnage ("bilinear", "average", "cubic")
#' @return SpatRaster rééchantillonné
resample_ign_to_spot <- function(ign_raster, target_res = RES_SPOT,
                                  method = "average") {
  current_res <- res(ign_raster)[1]
  if (abs(current_res - target_res) < 0.01) {
    message("Résolution déjà à ", target_res, " m")
    return(ign_raster)
  }

  agg_factor <- round(target_res / current_res)
  message(sprintf("Rééchantillonnage: %.2f m → %.2f m (facteur %dx)",
                   current_res, target_res, agg_factor))

  # Agrégation (moyenne des pixels) plutôt que simple rééchantillonnage
  # pour mieux simuler la résolution SPOT
  r_resampled <- terra::aggregate(ign_raster, fact = agg_factor, fun = "mean",
                                   na.rm = TRUE)

  message(sprintf("  Avant: %d x %d pixels", nrow(ign_raster), ncol(ign_raster)))
  message(sprintf("  Après: %d x %d pixels", nrow(r_resampled), ncol(r_resampled)))

  return(r_resampled)
}

# ==============================================================================
# Point d'entrée principal
# ==============================================================================

if (sys.nframe() == 0) {
  message("=== Open-Canopy + IGN Ortho : Chargement des données ===\n")

  # --- A. Hugging Face (Open-Canopy, SPOT 1.5m) ---
  message("--- A. Dataset Hugging Face (Open-Canopy, SPOT 1.5m) ---")
  message("Dataset: ", HF_REPO_ID)
  message("Token HF: ", ifelse(nchar(HF_TOKEN) > 0, "configuré", "non défini"))

  root_files <- hf_list_files(HF_REPO_ID, token = HF_TOKEN)
  if (nrow(root_files) > 0) {
    message("\nStructure du dataset:")
    print(root_files[, c("type", "path")])
  }

  # --- B. Données IGN locales (ortho RVB + IRC, 0.20m) ---
  message("\n--- B. Données IGN locales (BD ORTHO® 0.20m) ---")
  message("Répertoire: ", DATA_DIR_IGN)

  ign_files <- scan_ign_files()
  if (nrow(ign_files) > 0) {
    print(ign_files[, c("filename", "type", "size_mb")])
  } else {
    message("Placez vos fichiers IGN (.jp2 ou .tif) dans: ", DATA_DIR_IGN)
    message("\nPour télécharger via WMS, utilisez:")
    message('  bbox <- c(843000, 6518000, 844000, 6519000)  # Lambert-93')
    message('  download_ign_ortho_pair(bbox)')
  }

  # --- C. Comparaison des résolutions ---
  message("\n--- C. Comparaison SPOT vs IGN ---")
  message(sprintf("SPOT 6-7 : %.1f m/pixel → %d pixels/km²",
                   RES_SPOT, as.integer((1000 / RES_SPOT)^2)))
  message(sprintf("IGN Ortho: %.1f m/pixel → %d pixels/km²",
                   RES_IGN, as.integer((1000 / RES_IGN)^2)))
  message(sprintf("Facteur de résolution: %.1fx", RES_SPOT / RES_IGN))
}
