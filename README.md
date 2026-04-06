# Pipeline_Sparse
R code for sparse spatial weight construction, benchmarking, and spatial econometric analysis on large gridded raster data in Indonesia.

# Sparse Spatial Matrix Workflow for Large Gridded Data

This repository contains two complementary R scripts for building and benchmarking a sparse spatial weights workflow for very large gridded datasets in Indonesia.

## Repository contents

### 1. Full empirical workflow
**File:** `v4 HC1 Full Spatial_Sparse_W_Full_Template_All_Processes_patched_final.R`

This script implements the full empirical pipeline, including:

- raster preparation and harmonisation
- common-support construction across GDP, population, NDVI, and MNDWI
- sparse full-lattice weight construction using Kronecker products
- projection from the full template lattice to the retained empirical support
- isolate removal and row-standardisation
- Moran’s I / Geary’s C and related diagnostics
- SAR and SDM estimation
- HC1-style inference and output export

In other words, this is the main analysis script used to move from raw gridded inputs to the final spatial econometric results.

### 2. Benchmark and stress-test workflow
**File:** `benchmark_v4_aligned_sparse_vs_cell2nb_version3_compact_metadata_space.R`

This script benchmarks the matrix-construction stage of the proposed workflow against the conventional grid-based baseline:

- **Proposed method:** sparse Kronecker construction + projection matrix + isolate removal + row-standardisation
- **Baseline:** `cell2nb + subset.nb + nb2mat`

The benchmark is aligned to the same four-layer common-support logic used in the full empirical script. Its purpose is to quantify the computational gain of the proposed working pipeline before SAR and SDM estimation.

## Main idea

The repository is built around a simple principle:

> build the spatial weights matrix in sparse form from the outset, project it onto the valid empirical support, and avoid dense or list-heavy intermediate objects wherever possible.

This design is especially important for national-scale gridded data, where a conventional neighbour-list workflow becomes slow, memory-intensive, or infeasible.

## Suggested workflow

### Option A. Run the empirical analysis
Use the full workflow script when you want to reproduce the spatial econometric analysis.

### Option B. Run the benchmark
Use the benchmark script when you want to compare:

- the proposed Table 1 sparse pipeline
- against the conventional `cell2nb + nb2mat` route

This benchmark isolates the preprocessing and spatial-weight construction stage and does **not** benchmark SAR or SDM estimation itself.

## Data requirements

The scripts assume access to the following categories of input data:

- administrative boundary shapefiles for Indonesia
- gridded GDP raster
- LandScan population raster
- NDVI raster
- MNDWI raster

The administrative boundary for Indonesia is available at GADM Maps and Data: https://gadm.org/.
We use the Gridded GDP raster based on version 2 of Kummu et al. (2018): https://doi.org/10.5281/zenodo.13943886.
For the LandScan data available at: https://landscan.ornl.gov/.
The last two gridded datasets, NDVI and MNDWI, were downloaded using Google Earth Engine, and, for reproducibility, we uploaded them to this GitHub: Indonesia_NDVI_MODIS.zip and Indonesia_MNDWI.zip.
Because file paths in the working scripts were originally written for a local research environment, you will need to update directory paths before running the code.

## Required R packages

The scripts rely on a standard spatial and matrix stack in R, including packages such as:

- `terra`
- `sf`
- `Matrix`
- `spdep`
- `spatialreg`
- `splm`
- `openxlsx` or `writexl`
- `dplyr`
- `stats`

Depending on your local configuration, additional packages may be required for exporting results or diagnostics.

## Benchmark outputs

The benchmark script is designed to report:

- elapsed time
- peak object size / memory burden
- size of intermediate and final matrix objects
- isolates removed after projection
- feasibility boundary of the conventional baseline
- compact summary tables and metadata-ready outputs

This makes it suitable for appendix tables, robustness checks, and computational reproducibility.

## Reproducibility note

The benchmark is intentionally aligned with the matrix-construction stage of the full empirical workflow. That means the benchmark reproduces the same sequence of:

1. common-support definition
2. full-grid sparse matrix construction
3. projection onto valid support
4. isolate removal
5. row-standardisation

It then stops before the formation of the regression sample and the SAR/SDM likelihood evaluation.

## Before making the repository public

Please check the following before publishing:

- Remove or replace local machine paths
- Confirm that no private data files are included
- Rename scripts into cleaner public-facing filenames if needed
- Add a license file
- Add a citation file if you want others to reference the code properly
- Optionally add a small toy example for quick testing

## Suggested public-facing filenames

For a cleaner public repository, you may wish to rename the scripts as follows:

- `full_sparse_spatial_workflow.R`
- `benchmark_sparse_vs_cell2nb.R`

## Contact/citation note

If this repository supports a manuscript submission, it is helpful to add:

- the manuscript title
- author names
- version date
- a short note on how to cite the code

---

## Short description

This repository provides a sparse spatial matrix workflow for large gridded raster data and a benchmark showing the computational advantage of sparse Kronecker-plus-projection construction relative to the conventional `cell2nb + nb2mat` approach.
