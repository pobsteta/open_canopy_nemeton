# Open-Canopy R

Code R pour télécharger, analyser et visualiser le dataset **[Open-Canopy](https://huggingface.co/datasets/AI4Forest/Open-Canopy)** depuis Hugging Face.

## Description

[Open-Canopy](https://arxiv.org/abs/2407.09392) est un benchmark ouvert à l'échelle nationale pour l'estimation de la hauteur de canopée à très haute résolution (1.5 m). Le dataset couvre plus de 87 000 km² en France et combine :

- **Images SPOT 6-7** : imagerie satellite multispectrale à 1.5 m de résolution
- **CHM LiDAR** : modèles de hauteur de canopée dérivés de données LiDAR aériennes haute résolution
- **Benchmark de changement** : détection des changements de hauteur de canopée entre deux dates

## Structure du projet

```
├── R/
│   ├── 01_download_open_canopy.R   # Téléchargement depuis Hugging Face
│   ├── 02_analyse_open_canopy.R    # Analyse et visualisation
│   └── 03_prediction_open_canopy.R # Prédiction et post-traitement
├── data/                           # Données téléchargées (non versionné)
├── outputs/                        # Résultats et graphiques
├── .gitignore
└── README.md
```

## Prérequis

### Packages R requis

```r
install.packages(c("terra", "sf", "httr2", "jsonlite", "curl", "fs"))

# Optionnel (pour les modèles pré-entraînés via Python)
install.packages("reticulate")
```

### Token Hugging Face (optionnel)

Pour les datasets publics, aucun token n'est nécessaire. Si vous avez besoin d'un accès authentifié :

```bash
export HF_TOKEN="hf_votre_token_ici"
```

## Utilisation

### 1. Télécharger les données

```r
source("R/01_download_open_canopy.R")

# Télécharger un sous-ensemble de tuiles de test
download_open_canopy_subset(split = "test", n_tiles = 5)

# Télécharger des tuiles d'entraînement
download_open_canopy_subset(split = "train", n_tiles = 100, data_type = "all")

# Cloner le dataset complet (~360 Go, nécessite git-lfs)
hf_git_clone()
```

### 2. Analyser et visualiser

```r
source("R/02_analyse_open_canopy.R")

# Charger une image SPOT
spot <- load_spot_image("data/test/images/tile_001.tif")
plot_spot_rgb(spot)

# Charger et analyser un CHM
chm <- load_chm("data/test/lidar/tile_001.tif")
stats <- compute_chm_stats(chm)
plot_chm(chm)
plot_canopy_classes(chm)

# Comparer image et CHM
plot_spot_chm_comparison(spot, chm, tile_name = "Tuile 001")
```

### 3. Évaluer les prédictions

```r
source("R/03_prediction_open_canopy.R")

# Évaluer des prédictions par rapport à la référence LiDAR
prediction <- load_predictions("outputs/prediction.tif")
reference <- load_chm("data/test/lidar/tile_001.tif")
metrics <- evaluate_predictions(prediction, reference)

# Couverture forestière
forest_cover <- compute_forest_cover_grid(chm, cell_size = 100)

# Détection de perte de canopée
loss <- detect_canopy_loss(chm_t1, chm_t2, threshold = -5)
```

## Dataset Open-Canopy

| Caractéristique | Valeur |
|---|---|
| Couverture | 87 000+ km² (France) |
| Résolution | 1.5 m |
| Taille | ~360 Go |
| Train | 66 339 km² |
| Validation | 7 369 km² |
| Test | 13 675 km² |
| Buffer | 1 km entre test et train/val |

## Références

- **Paper** : Fogel et al. (2024). "Open-Canopy: Towards Very High Resolution Forest Monitoring." [arXiv:2407.09392](https://arxiv.org/abs/2407.09392)
- **Dataset** : [AI4Forest/Open-Canopy sur Hugging Face](https://huggingface.co/datasets/AI4Forest/Open-Canopy)
- **Code Python** : [fajwel/Open-Canopy sur GitHub](https://github.com/fajwel/Open-Canopy)
- **Projet** : AI4Forest, financé par l'ANR, le DLR et le BMBF
