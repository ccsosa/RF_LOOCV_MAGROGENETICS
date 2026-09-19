# Spatial Macrogenetics Prediction Workflow (Random Forest via `ranger`)

[![R](https://img.shields.io/badge/R-≥4.1.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Spatial-Analysis](https://img.shields.io/badge/Focus-Spatial_Macrogenetics-green.svg)](#)

An R pipeline for predicting continuous spatial variation in genetic diversity metrics (e.g., observed heterozygosity, $H_o$) across a species' geographic distribution using Random Forest regression (`ranger`), feature selection via Boruta, nested Leave-One-Out Cross-Validation (LOOCV), and Multivariate Environmental Similarity Surface (MESS) masking.

---

## Key Features

* **Spatial Processing & Point Aggregation**: Automatically snaps out-of-coverage sample coordinates to the nearest valid raster cell (using a 5 km buffer) and aggregates multiple samples falling in the same cell.
* **Coastal & Marine Edge Correction**: Employs focal spatial smoothing (`terra::focal`) and nearest-neighbor fallback extractions to avoid `NA` predictor values at complex coastal edges.
* **Strictly Nested LOOCV**: Variable selection (Boruta consensus) and hyperparameter grid searches (`num.trees`, `mtry`, `min.node.size`) are independently executed *within* each cross-validation fold to prevent data leakage from held-out test cells.
* **MESS Interpolation Masking**: Computes a Multivariate Environmental Similarity Surface (MESS) to restrict spatial projections strictly to regions of environmental interpolation ($MESS > 0$), masking out novel/extrapolated environmental spaces.
* **Spatial Autocorrelation Diagnostics**: Assesses spatial dependence in cross-validation residuals using Moran's $I$ with inverse-distance spatial weighting (`ape::Moran.I`).
* **Memory-Efficient Spatial Prediction**: Predicts raster outputs in chunks to manage memory overhead across large high-resolution extents.
* **Publication-Ready Exports**: Automatically generates performance comparison plots (Train vs. LOOCV), permutation-based variable importance charts, and high-resolution PDF spatial maps using `tmap`.

---

## Workflow Architecture
```mermaid
graph TD
    classDef process fill:#f9f9f9,stroke:#333,stroke-width:1px,text-align:left;

    subgraph Step1 ["1. Data Ingestion & Spatial Preprocessing"]
        A1["• Snap genetic points to environmental grid cells<br>• Aggregate duplicate cell samples (mean Ho, sample count)<br>• Extract focal-smoothed bioclimatic covariates"]
    end

    subgraph Step2 ["2. Fully Nested Leave-One-Out Cross-Validation (LOOCV)"]
        B1["For fold i in 1..N:<br> ├── Boruta Consensus Feature Selection (Training set)<br> ├── Hyperparameter Grid Search (OOB R² optimization)<br> └── Fit model on train set & predict held-out cell i"]
    end

    subgraph Step3 ["3. Production Model & Spatial Projections"]
        C1["• Re-run Boruta & grid search on complete dataset<br>• Fit final production Random Forest model<br>• Chunked prediction restricted to binary SDM footprint<br>• Mask spatial output using MESS (keep MESS > 0)<br>• Calculate Moran's I on LOOCV residuals"]
    end

    Step1 --> Step2
    Step2 --> Step3

    class A1,B1,C1 process;
```
---

## Dependencies

Required R packages:

```r
install.packages(c(
  "terra", "sf", "ranger", "patchwork", "dplyr", 
  "ggplot2", "ggpmisc", "tmap", "Boruta", "ape", 
  "RColorBrewer", "readxl", "caret", "raster", "dismo"
))
```

# # Input Requirements
- Genetic Data (Tabla_to_model.xlsx): Excel spreadsheet with a sheet named "data" containing:
- sp: Species scientific name string (e.g., "Crocodylus acutus").
- Ho: Continuous genetic response variable (Observed Heterozygosity).
- lon, lat: Georeferenced sample coordinates in WGS84 (EPSG:4326).
- Environmental Rasters (raster_dir): Directory containing continuous predictor raster layers in standard GeoTIFF format (.tif).
- Binary SDM Raster (sdm_path): A binary Species Distribution Model raster (1 = suitable/presence, 0/NA = unsuitable) to restrict spatial projections.

### Usage Example

# Load workflow function
source("RF_LOOCV.R")

# Run execution pipeline
```r
RF_LOOCV(
  outdir           = "path/to/output_dir",
  sp_name          = "Crocodylus acutus",
  raster_dir       = "path/to/raster_dir",
  data_path        = "path/to/genetic_data.xlsx",
  sdm_path         = "path/to/binary_sdm.tif",
  N_CORES          = 6L,   # Number of CPU threads for ranger
  boruta_threshold = 10L,  # Consensus threshold (e.g., 10/10 runs)
  n_seeds          = 10L   # Number of Boruta runs per iteration
)
```


###Key Methodological References
- Ranger: Wright, M. N., & Ziegler, A. (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. Journal of Statistical Software, 77(1), 1–17.
- Boruta: Kursa, M. B., & Rudnicki, W. R. (2010). Feature Selection with the Boruta Package. Journal of Statistical Software, 36(11), 1–13.
- MESS: Elith, J., Kearney, M., & Phillips, S. (2010). The art of modelling range-shifting species. Methods in Ecology and Evolution, 1(4), 330–342.
- Sosa, C.C.; Arenas, C.; García-Merchán, V.H. Human Population Density Influences Genetic Diversity of Two Rattus Species Worldwide: A Macrogenetic Approach. Genes 2023, 14, 1442. https://doi.org/10.3390/genes14071442


# Authors:
Chrystian C. Sosa, Jorge Gómez-Marulanda, Victor Hugo García-Merchán