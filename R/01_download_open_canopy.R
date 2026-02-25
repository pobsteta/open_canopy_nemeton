#!/usr/bin/env Rscript
# ==============================================================================
# 01_download_open_canopy.R
# Téléchargement du dataset Open-Canopy depuis Hugging Face
# Dataset: AI4Forest/Open-Canopy
# https://huggingface.co/datasets/AI4Forest/Open-Canopy
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

# Identifiant du dataset Hugging Face
HF_REPO_ID <- "AI4Forest/Open-Canopy"
HF_API_URL <- "https://huggingface.co/api/datasets"

# Répertoire local de téléchargement
DATA_DIR <- file.path(getwd(), "data")
dir_create(DATA_DIR)

# Token Hugging Face (optionnel pour les datasets publics)
# Définir la variable d'environnement HF_TOKEN ou la renseigner ici
HF_TOKEN <- Sys.getenv("HF_TOKEN", unset = "")

# ==============================================================================
# Fonctions de téléchargement
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
    warning("Erreur lors de la requête API: ", resp_status(resp))
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

  message("Téléchargement: ", filename)

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
    warning("Échec du téléchargement de ", filename, ": ", e$message)
    return(NULL)
  })
}

#' Télécharger un ensemble de fichiers depuis Hugging Face
#'
#' @param repo_id Identifiant du dépôt
#' @param file_list Vecteur de chemins de fichiers
#' @param dest_dir Répertoire de destination
#' @param token Token Hugging Face
#' @param overwrite Écraser les fichiers existants
#' @return Vecteur des chemins locaux
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

# ==============================================================================
# Téléchargement du dataset Open-Canopy (sous-ensemble)
# ==============================================================================

#' Télécharger un sous-ensemble du dataset Open-Canopy
#'
#' Le dataset complet fait ~360 Go. Cette fonction permet de télécharger
#' un nombre limité de tuiles pour exploration et prototypage.
#'
#' @param split Split à télécharger: "train", "val", ou "test"
#' @param n_tiles Nombre de tuiles à télécharger
#' @param data_type Type de données: "images" (SPOT), "lidar" (CHM),
#'   "lidar_v2" (CHM v2), ou "all"
#' @param dest_dir Répertoire de destination
#' @param token Token Hugging Face
download_open_canopy_subset <- function(split = "test",
                                         n_tiles = 5,
                                         data_type = "all",
                                         dest_dir = DATA_DIR,
                                         token = HF_TOKEN) {
  message("=== Téléchargement Open-Canopy ===")
  message("Split: ", split, " | Tuiles: ", n_tiles, " | Type: ", data_type)

  # Lister les fichiers disponibles dans le split
  message("\nListe des fichiers disponibles...")
  files <- hf_list_files(HF_REPO_ID, path = split, token = token)

  if (nrow(files) == 0) {
    # Essayer de lister les sous-dossiers images/lidar
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
    message("Essai de téléchargement direct de tuiles exemples...")
    download_sample_tiles(split, n_tiles, data_type, dest_dir, token)
    return(invisible(NULL))
  }

  # Filtrer par type de données
  if (data_type == "images") {
    files <- files[grep("images|spot", files$path, ignore.case = TRUE), ]
  } else if (data_type == "lidar") {
    files <- files[grep("lidar(?!_v2)", files$path, perl = TRUE), ]
  } else if (data_type == "lidar_v2") {
    files <- files[grep("lidar_v2", files$path), ]
  }

  # Limiter au nombre de tuiles demandé
  tif_files <- files[grep("\\.tif$", files$path), ]
  if (nrow(tif_files) > n_tiles) {
    tif_files <- tif_files[seq_len(n_tiles), ]
  }

  message(sprintf("\n%d fichier(s) à télécharger.", nrow(tif_files)))

  # Télécharger
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

#' Télécharger des tuiles exemples avec des noms connus
#'
#' @param split Split du dataset
#' @param n_tiles Nombre de tuiles
#' @param data_type Type de données
#' @param dest_dir Répertoire de destination
#' @param token Token HF
download_sample_tiles <- function(split = "test", n_tiles = 5,
                                   data_type = "all", dest_dir = DATA_DIR,
                                   token = HF_TOKEN) {
  message("Tentative de téléchargement de tuiles exemples...")

  # Structure typique du dataset:
  # {split}/images/{tile_id}.tif  - Image SPOT 6-7
  # {split}/lidar/{tile_id}.tif   - CHM dérivé du LiDAR
  # {split}/lidar_v2/{tile_id}.tif - CHM v2

  # Essayer de lister les répertoires pour trouver les tuiles
  subdirs <- c("images", "lidar", "lidar_v2")
  if (data_type != "all") {
    subdirs <- data_type
  }

  for (subdir in subdirs) {
    path <- paste0(split, "/", subdir)
    files <- hf_list_files(HF_REPO_ID, path = path, token = token)

    if (nrow(files) > 0) {
      tif_files <- files[grep("\\.tif$", files$path), ]
      if (nrow(tif_files) > n_tiles) {
        tif_files <- tif_files[seq_len(n_tiles), ]
      }

      message(sprintf("Téléchargement de %d fichier(s) depuis %s...",
                       nrow(tif_files), path))

      hf_download_files(
        repo_id = HF_REPO_ID,
        file_list = tif_files$path,
        dest_dir = dest_dir,
        token = token
      )
    }
  }
}

# ==============================================================================
# Téléchargement via git clone (alternative pour le dataset complet)
# ==============================================================================

#' Cloner le dataset complet via git (nécessite git-lfs)
#'
#' @param dest_dir Répertoire de destination
#' @param token Token Hugging Face (optionnel)
hf_git_clone <- function(dest_dir = file.path(DATA_DIR, "Open-Canopy"),
                          token = HF_TOKEN) {
  if (dir.exists(dest_dir)) {
    message("Le répertoire existe déjà: ", dest_dir)
    return(invisible(dest_dir))
  }

  # Vérifier git-lfs
  lfs_check <- system("git lfs version", intern = TRUE, ignore.stderr = TRUE)
  if (length(lfs_check) == 0) {
    stop("git-lfs n'est pas installé. ",
         "Installez-le avec: sudo apt install git-lfs (Linux) ",
         "ou brew install git-lfs (macOS)")
  }

  message("Clonage du dataset Open-Canopy (~360 Go)...")
  message("Ceci peut prendre plusieurs heures selon votre connexion.")

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
# Point d'entrée principal
# ==============================================================================

if (sys.nframe() == 0) {
  message("=== Open-Canopy: Téléchargement depuis Hugging Face ===\n")
  message("Dataset: ", HF_REPO_ID)
  message("Répertoire: ", DATA_DIR)
  message("Token HF: ", ifelse(nchar(HF_TOKEN) > 0, "configuré", "non défini"))
  message("")

  # Lister les fichiers à la racine du dataset
  message("Exploration de la structure du dataset...")
  root_files <- hf_list_files(HF_REPO_ID, token = HF_TOKEN)
  if (nrow(root_files) > 0) {
    message("\nStructure à la racine:")
    print(root_files[, c("type", "path")])
  }

  # Télécharger un sous-ensemble de tuiles de test
  message("\n--- Téléchargement d'un sous-ensemble de test ---")
  downloaded <- download_open_canopy_subset(
    split = "test",
    n_tiles = 3,
    data_type = "all"
  )

  message("\nPour télécharger le dataset complet, utilisez:")
  message('  hf_git_clone()')
  message("\nPour télécharger plus de tuiles:")
  message('  download_open_canopy_subset(split = "train", n_tiles = 100)')
}
