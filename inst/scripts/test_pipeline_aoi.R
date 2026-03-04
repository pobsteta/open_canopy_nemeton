#!/usr/bin/env Rscript
# ==============================================================================
# test_pipeline_aoi.R
# Script de test du package open_canopy_nemeton sur un fichier aoi.gpkg
#
# Usage :
#   Rscript inst/scripts/test_pipeline_aoi.R
#   Rscript inst/scripts/test_pipeline_aoi.R chemin/vers/mon_aoi.gpkg
#   Rscript inst/scripts/test_pipeline_aoi.R data/aoi.gpkg pvtv2
#
# Arguments (optionnels) :
#   1. Chemin vers le fichier .gpkg (défaut: data/aoi.gpkg)
#   2. Nom du modèle : "unet" ou "pvtv2" (défaut: pvtv2)
# ==============================================================================

# --- Paramètres CLI ---
args <- commandArgs(trailingOnly = TRUE)
aoi_path   <- if (length(args) >= 1) args[1] else "data/aoi.gpkg"
model_name <- if (length(args) >= 2) args[2] else "pvtv2"

cat("==============================================================\n")
cat("  Test du package open_canopy_nemeton\n")
cat("==============================================================\n\n")

# --- 1. Vérifier les dépendances R ---
cat(">>> Étape 1 : Vérification des dépendances R\n")

pkgs_required <- c("terra", "sf", "httr2", "jsonlite", "curl", "fs", "reticulate")
pkgs_optional <- c("hfhub", "ggplot2", "patchwork", "tidyterra")

ok <- TRUE
for (pkg in pkgs_required) {
  avail <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  [%s] %s (requis)\n", ifelse(avail, "OK", "!!"), pkg))
  if (!avail) ok <- FALSE
}
for (pkg in pkgs_optional) {
  avail <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  [%s] %s (optionnel)\n", ifelse(avail, "OK", "--"), pkg))
}

if (!ok) {
  stop("Packages R manquants. Installez-les avec :\n",
       "  install.packages(c('",
       paste(pkgs_required, collapse = "', '"), "'))",
       call. = FALSE)
}

has_hfhub <- requireNamespace("hfhub", quietly = TRUE)
cat(sprintf("\n  Méthode de téléchargement HF : %s\n",
            ifelse(has_hfhub,
                   "hfhub (R natif, recommandé)",
                   "Python huggingface_hub (fallback)")))

# --- 2. Vérifier le fichier AOI ---
cat("\n>>> Étape 2 : Vérification du fichier AOI\n")

if (!file.exists(aoi_path)) {
  stop("Fichier AOI introuvable : ", aoi_path, "\n",
       "  Placez un fichier .gpkg dans data/ ou passez le chemin en argument.\n",
       "  Usage : Rscript inst/scripts/test_pipeline_aoi.R chemin/vers/aoi.gpkg",
       call. = FALSE)
}

library(sf)
aoi <- st_read(aoi_path, quiet = TRUE)
cat(sprintf("  Fichier  : %s\n", aoi_path))
cat(sprintf("  Couches  : %s\n", paste(st_layers(aoi_path)$name, collapse = ", ")))
cat(sprintf("  Features : %d\n", nrow(aoi)))
cat(sprintf("  CRS      : %s\n", st_crs(aoi)$input))

# Emprise en Lambert-93
aoi_l93 <- st_transform(aoi, 2154)
bbox <- st_bbox(aoi_l93)
area_km2 <- as.numeric(st_area(st_union(aoi_l93))) / 1e6
cat(sprintf("  Emprise Lambert-93 : [%.0f, %.0f] → [%.0f, %.0f]\n",
            bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]))
cat(sprintf("  Surface  : %.2f km²\n", area_km2))

if (area_km2 > 25) {
  cat("  ATTENTION : surface > 25 km², le traitement sera long.\n")
}

# --- 3. Vérifier l'environnement Python ---
cat("\n>>> Étape 3 : Vérification de l'environnement Python\n")

library(reticulate)

# Détecter l'environnement conda (Miniforge, Miniconda, Anaconda)
conda_found <- FALSE

# 1. Si Python est déjà dans le PATH (ex: RStudio lancé depuis conda activate)
python_in_path <- tryCatch({
  cfg <- py_config()
  nzchar(cfg$python)
}, error = function(e) FALSE)

if (!python_in_path) {
  # 2. Chercher l'env conda par nom (fonctionne avec Miniforge/Miniconda/Anaconda)
  conda_envs <- tryCatch(conda_list(), error = function(e) data.frame())
  if (nrow(conda_envs) > 0) {
    match_idx <- grep("open_canopy", conda_envs$name, ignore.case = TRUE)
    if (length(match_idx) > 0) {
      use_python(conda_envs$python[match_idx[1]], required = TRUE)
      conda_found <- TRUE
      cat(sprintf("  Conda env : %s\n", conda_envs$name[match_idx[1]]))
    }
  }

  if (!conda_found) {
    # 3. Chercher via CONDA_PREFIX (variable définie quand un env est activé)
    conda_prefix <- Sys.getenv("CONDA_PREFIX", "")
    if (nzchar(conda_prefix)) {
      py_path <- if (.Platform$OS.type == "windows") {
        file.path(conda_prefix, "python.exe")
      } else {
        file.path(conda_prefix, "bin", "python")
      }
      if (file.exists(py_path)) {
        use_python(py_path, required = TRUE)
        conda_found <- TRUE
        cat(sprintf("  Conda env (CONDA_PREFIX) : %s\n", conda_prefix))
      }
    }
  }
}

if (python_in_path || conda_found) {
  cat(sprintf("  Python : %s\n", py_config()$python))
} else {
  cat("  ATTENTION : aucun environnement Python detecte.\n")
  cat("  Le pipeline necessite Python avec PyTorch pour l'inference.\n")
  cat("  Installation :\n")
  cat("    conda create -n open_canopy python=3.10\n")
  cat("    conda activate open_canopy\n")
  cat("    pip install torch torchvision numpy rasterio segmentation-models-pytorch timm\n")
  cat("  Puis lancez RStudio depuis cet environnement.\n")
}

modules_core <- c("torch", "numpy", "rasterio",
                   "segmentation_models_pytorch", "timm")
for (mod in modules_core) {
  avail <- py_module_available(mod)
  cat(sprintf("  [%s] Python %s\n", ifelse(avail, "OK", "!!"), mod))
}

# --- 4. Tester le téléchargement du modèle ---
cat("\n>>> Étape 4 : Test du téléchargement du modèle\n")
cat(sprintf("  Modèle demandé : %s\n", model_name))

source("R/04_pipeline_aoi_to_chm.R")

tryCatch({
  model_path <- download_model(model_name)
  cat(sprintf("  Modèle téléchargé : %s\n", model_path))
  cat(sprintf("  Taille : %.1f Mo\n", file.size(model_path) / 1024^2))
}, error = function(e) {
  cat(sprintf("  ERREUR téléchargement : %s\n", e$message))
  cat("  Vous pouvez fournir le chemin manuellement (voir étape 5).\n")
  model_path <<- NULL
})

# --- 5. Lancer le pipeline complet ---
cat("\n>>> Étape 5 : Lancement du pipeline AOI → CHM\n")
cat(sprintf("  AOI    : %s\n", aoi_path))
cat(sprintf("  Modèle : %s\n", model_name))
cat(sprintf("  Sortie : outputs/\n"))

t0 <- Sys.time()

tryCatch({
  result <- pipeline_aoi_to_chm(
    aoi_path   = aoi_path,
    output_dir = "outputs",
    model_name = model_name
  )

  elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)

  cat("\n==============================================================\n")
  cat("  Pipeline terminé avec succès !\n")
  cat(sprintf("  Durée : %s min\n", elapsed))
  cat("==============================================================\n\n")

  # Lister les fichiers produits
  out_files <- list.files("outputs", pattern = "\\.(tif|pdf)$", full.names = TRUE)
  if (length(out_files) > 0) {
    cat("Fichiers produits :\n")
    for (f in out_files) {
      cat(sprintf("  %s (%.1f Mo)\n", basename(f), file.size(f) / 1024^2))
    }
  }

}, error = function(e) {
  elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 0)
  cat(sprintf("\n  ERREUR après %s sec : %s\n", elapsed, e$message))
  cat("\n  Conseils :\n")
  cat("  - Vérifiez votre token HF : Sys.setenv(HF_TOKEN = 'hf_...')\n")
  cat("  - Vérifiez l'env conda : conda activate open_canopy\n")
  cat("  - Essayez avec un chemin de modèle local :\n")
  cat(sprintf("      source('R/04_pipeline_aoi_to_chm.R')\n"))
  cat(sprintf("      pipeline_aoi_to_chm('%s', model_path = 'chemin/model.ckpt')\n",
              aoi_path))
})
