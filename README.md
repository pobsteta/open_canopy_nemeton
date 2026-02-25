# Open-Canopy R — Ortho IGN 0.20m + SPOT 1.5m

Code R pour estimer la hauteur de canopée à partir des **ortho IGN** (RVB + IRC à 0.20m) en exploitant les modèles pré-entraînés **[Open-Canopy](https://huggingface.co/datasets/AI4Forest/Open-Canopy)** (SPOT 6-7 à 1.5m).

## Contexte

Le projet **Open-Canopy** fournit des modèles (UNet, PVTv2) entraînés sur images SPOT 6-7 à 1.5m pour prédire la hauteur de canopée. Ce dépôt adapte ces modèles pour fonctionner avec les **orthophotos IGN** (BD ORTHO®) à **0.20m**, qui offrent une résolution 7.5x supérieure.

### Comparaison des données sources

| | Open-Canopy (SPOT) | Ortho IGN |
|---|---|---|
| **Résolution** | 1.5 m/pixel | 0.20 m/pixel |
| **Pixels/km²** | 444 444 | 25 000 000 |
| **Source** | Satellite SPOT 6-7 | Photo aérienne IGN |
| **Bandes RVB** | R, V, B (+ PIR) | R, V, B |
| **Bandes IRC** | — | PIR, R, V |
| **Format** | GeoTIFF | JPEG2000 (.jp2) |
| **NDVI** | Non | Oui (depuis IRC) |
| **Couverture** | 87 000 km² (France) | France entière |

## Structure du projet

```
├── R/
│   ├── 01_download_open_canopy.R   # Téléchargement HF + chargement IGN
│   ├── 02_analyse_open_canopy.R    # Analyse, NDVI, visualisation
│   └── 03_prediction_open_canopy.R # Inférence modèles + post-traitement
├── data/
│   ├── open_canopy/                # Dataset HF (SPOT 1.5m)
│   └── ign/                        # Ortho IGN (RVB + IRC, 0.20m)
├── outputs/                        # Résultats et graphiques
├── .gitignore
└── README.md
```

## Prérequis

### Packages R

```r
install.packages(c("terra", "sf", "httr2", "jsonlite", "curl", "fs"))
install.packages("reticulate")  # Interface Python
```

### Environnement Python (miniforge + conda)

L'environnement conda `open_canopy` doit être installé avec PyTorch :

```bash
# Installer miniforge si nécessaire
# https://github.com/conda-forge/miniforge

conda create -n open_canopy python=3.10
conda activate open_canopy
pip install torch torchvision numpy rasterio huggingface_hub
```

### Connexion R → Python

```r
library(reticulate)
use_condaenv("open_canopy", required = TRUE)
```

## Utilisation

### 1. Charger les données

```r
source("R/01_download_open_canopy.R")

# --- Ortho IGN (placer les .jp2 dans data/ign/) ---
irc <- load_ign_ortho("data/ign/ortho_irc.jp2", type = "irc")
rvb <- load_ign_ortho("data/ign/ortho_rvb.jp2", type = "rvb")

# Télécharger via WMS IGN pour une emprise donnée (Lambert-93)
bbox <- c(843000, 6518000, 844000, 6519000)
download_ign_ortho_pair(bbox)

# --- Open-Canopy (SPOT, Hugging Face) ---
download_open_canopy_subset(split = "test", n_tiles = 5)
```

### 2. Analyser (NDVI, classification, croisement)

```r
source("R/02_analyse_open_canopy.R")

# Charger l'IRC et calculer le NDVI
irc <- load_ign_ortho("data/ign/ortho_irc.jp2", type = "irc")
ndvi <- compute_ndvi(irc)        # (PIR - R) / (PIR + R)
gndvi <- compute_gndvi(irc)      # (PIR - V) / (PIR + V)
savi <- compute_savi(irc)        # SAVI (ajusté sol)

# Masque de végétation
veg <- mask_vegetation(ndvi, threshold = 0.3)

# Visualisation multi-panneaux
plot_ign_full_comparison(rvb, irc, chm)

# Croisement NDVI x CHM
cross <- cross_ndvi_chm(ndvi, chm)

# Statistiques
irc_stats <- compute_irc_stats(irc)
chm_stats <- compute_chm_stats(chm)
```

### 3. Prédire la hauteur de canopée

```r
source("R/03_prediction_open_canopy.R")

# Configurer l'environnement Python
setup_conda_env("open_canopy")

# Télécharger un modèle pré-entraîné
model_path <- download_pretrained_model("unet")

# Pipeline complet : IGN 0.20m → agrégation 1.5m → inférence → CHM
predict_chm_from_ign("data/ign/ortho_rvb.jp2", model_path)

# Rééchantillonner manuellement
ign <- load_ign_ortho("data/ign/ortho_rvb.jp2", type = "rvb")
ign_1_5m <- resample_ign_to_spot(ign)  # 0.20m → 1.5m (facteur 8x)

# Suréchantillonner un CHM prédit vers la résolution IGN
chm_hr <- upsample_chm_to_ign(chm_predicted)  # 1.5m → 0.20m
```

## Workflow de rééchantillonnage

```
Ortho IGN (0.20m)                    Modèles Open-Canopy
5000 x 5000 px/km²                   entraînés sur SPOT 1.5m
        │                                     │
        ▼                                     │
  aggregate(fact=8)                            │
  moyenne 8x8 pixels                           │
        │                                     │
        ▼                                     ▼
  Image à 1.5m  ──────────────────►  Inférence UNet/PVTv2
  667 x 667 px/km²                            │
                                              ▼
                                    CHM prédit à 1.5m
                                              │
                                              ▼
                                    disagg(fact=8)    (optionnel)
                                              │
                                              ▼
                                    CHM à 0.20m
```

## Indices spectraux (depuis IRC IGN)

| Indice | Formule | Usage |
|---|---|---|
| NDVI | (PIR - R) / (PIR + R) | Activité végétale |
| GNDVI | (PIR - V) / (PIR + V) | Teneur en chlorophylle |
| SAVI | ((PIR - R) / (PIR + R + L)) × (1+L) | Végétation sur sol nu |

## Références

- **Open-Canopy** : Fogel et al. (2024). [arXiv:2407.09392](https://arxiv.org/abs/2407.09392)
- **Dataset HF** : [AI4Forest/Open-Canopy](https://huggingface.co/datasets/AI4Forest/Open-Canopy)
- **Code Python** : [fajwel/Open-Canopy](https://github.com/fajwel/Open-Canopy)
- **BD ORTHO® IGN** : [geoservices.ign.fr/bdortho](https://geoservices.ign.fr/bdortho)
- **Géoplateforme** : [data.geopf.fr](https://data.geopf.fr)
