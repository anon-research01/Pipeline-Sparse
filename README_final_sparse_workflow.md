# Sparse Spatial Matrix Workflow for Large Gridded Data

This repository supports the manuscript:

**Constructing Million-Cell Spatial Weight Matrices on Irregular Support: A Kronecker-Plus-Projection Workflow**

It contains R code for constructing sparse spatial weight matrices from large regular-grid raster data observed on irregular analytical support, benchmarking the construction stage, and running the spatial econometric application for Indonesian gridded data.

## Overview

The workflow is designed for large raster systems where the full template grid is much larger than the valid empirical support. In the Indonesian application, the full national template contains about 11.3 million grid cells, while about 2.2 million cells belong to the common valid support used before isolate removal.

The main idea is to:

1. align all raster layers to a common template;
2. define the common valid support across all analytical layers;
3. construct full-lattice contiguity in sparse form using Kronecker products;
4. project the full-lattice matrix to the retained empirical support;
5. remove zero-degree retained cells;
6. row-standardise the pruned sparse matrix; and
7. use the final matrix for spatial econometric analysis.

This avoids dense or list-heavy intermediate matrix objects wherever possible.

## Main manuscript-aligned notation

The final manuscript uses the following notation consistently:

- `G`: full rectangular template grid.
- `R`: number of raster rows.
- `C`: number of raster columns.
- `N = R x C`: full template size.
- `mathcal{C}`: common valid support, defined as the intersection of valid cells across GDP, population, NDVI, and MNDWI.
- `|mathcal{C}|`: common-support size before isolate removal.
- `Q`: non-isolated retained cells after projection and isolate removal.
- `n = |Q|`: final non-isolated estimation sample and final matrix dimension.
- `W_FULL`: sparse full-lattice contiguity matrix.
- `W_COMMON`: projected matrix on the common valid support.
- `W_PRUNED`: matrix after removing zero-degree retained cells.
- `D_Q = diag(W_PRUNED 1)`: diagonal matrix of row sums after pruning.
- `W_tilde = D_Q^{-1} W_PRUNED`: final row-standardised sparse spatial weight matrix.

The final `W_tilde` matrix is an `n x n` matrix and is the matrix used in SAR, SDM, spatial lag construction, log-determinant evaluation, and direct-indirect effect calculations.

## Repository contents

### 1. Full empirical workflow

**Current working script name:** `v4 HC1 Full Spatial_Sparse_W_Full_Template_All_Processes_patched_final.R`

**Suggested public-facing name:** `full_sparse_spatial_workflow.R`

This script implements the full empirical workflow:

- raster preparation and harmonisation;
- common valid support construction across GDP, population, NDVI, and MNDWI;
- sparse full-lattice spatial weight construction using Kronecker products;
- projection from the full template lattice to the common valid support;
- isolate removal and row-standardisation;
- construction of the final `W_tilde` matrix;
- SAR and SDM estimation;
- MLE, Huber-White / HC1, and Spatial Block Subsampling standard errors;
- export of model tables and computational summaries.

### 2. Benchmark and stress-test workflow

**Current working script name:** `benchmark_v4_aligned_sparse_vs_cell2nb_version3_compact_metadata_space.R`

**Suggested public-facing name:** `benchmark_sparse_vs_spdep_vs_terra.R`

This script benchmarks the matrix-construction stage of the workflow. The benchmark compares three construction routes applied to the same retained-support scenarios and queen-contiguity rule:

1. **Proposed sparse Kronecker-projection route**
   - sparse full-lattice Kronecker construction;
   - projection to retained support;
   - isolate removal;
   - row-standardisation.

2. **Conventional dense-output `spdep` route**
   - `cell2nb`;
   - `subset.nb`;
   - isolate removal;
   - materialisation through `nb2mat`.

3. **Raster-native `terra::adjacent` route**
   - full-lattice raster adjacency pairs;
   - filtering to the retained support;
   - sparse matrix construction;
   - isolate removal;
   - row-standardisation.

The benchmark is restricted to the spatial weight matrix construction stage. It stops before regression-sample formation and before SAR or SDM likelihood estimation.

## Benchmark support scenarios

The benchmark uses the same retained-support proportions as the final manuscript:

```text
0.01%, 0.1%, 0.5%, 1%, 5%, 10%, 25%, 50%, 75%, 100%
```

In decimal form, the support set is:

```text
{0.0001, 0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00}
```

These scenarios are nested subsets of the Indonesian common valid support.

## Benchmark outputs

The benchmark records:

- support level;
- number of support cells;
- cells remaining after isolate removal;
- elapsed construction time;
- tracked peak object size;
- size of the final sparse matrix object;
- number of nonzero neighbour links;
- isolates removed;
- comparator status;
- feasibility boundary of the dense-output `spdep` route;
- pairwise time and memory ratios.

Ratios are interpreted as comparator value divided by proposed-workflow value. Values greater than one indicate that the comparator requires more elapsed time or more memory than the proposed workflow.

## Dense-reference scaling figure

The dense-reference figure in the manuscript is not a formal benchmark against an optimised dense implementation. It is an order-of-magnitude infeasibility illustration.

The final manuscript treats:

- direct dense log-determinant evaluation as scaling approximately with `O(n^3)`;
- dense storage as scaling approximately with `O(n^2)`;
- sparse storage as governed by the number of nonzero neighbour links, `nnz(W_tilde)`, plus index and pointer arrays.

The plotted sparse memory series should be interpreted as an implementation-specific sparse object size, not as an exact CSR or CSC memory formula.

## Spatial econometric estimation

The empirical application estimates SAR and SDM specifications on the final non-isolated retained support:

```text
n = |Q|
W_tilde in R^{n x n}
```

The log-determinant is written consistently as:

```text
log |I_n - rho W_tilde|
```

The Monte Carlo trace-series approximation uses:

```text
K_tr = 60 powers
m = 40 probe vectors
```

Spatial Block Subsampling is used as a dependence-aware sensitivity check for uncertainty assessment, not as a procedure that fully resolves all inferential concerns. The SBS notation used in the final manuscript is:

```text
K_SBS = 100 graph-connected blocks
m_target = 10,000 cells
m_eff = average realised block size
scale factor = (m_eff / n)^{1/2}
```

The SBS routine re-estimates the model on graph-connected blocks and rescales across-block dispersion to the full retained sample.

## Data requirements

The scripts assume access to the following input data:

- Indonesian administrative boundary shapefiles;
- gridded GDP raster;
- LandScan population raster;
- NDVI raster;
- MNDWI raster.

Indicative sources include:

- administrative boundaries: GADM;
- GDP raster: Kummu et al. downscaled gridded GDP data;
- population raster: LandScan Global;
- NDVI and MNDWI: processed through Google Earth Engine.

Local file paths in the working scripts must be updated before running the code.

## Required R packages

The workflow relies on a standard spatial and sparse-matrix stack in R, including:

- `terra`
- `raster`
- `sf`
- `spdep`
- `spatialreg`
- `Matrix`
- `igraph`
- `dplyr`
- `ggplot2`
- `openxlsx` or `writexl`
- `stats`

Additional packages may be required depending on local export, plotting, or diagnostics settings.

## Suggested workflow for users

### Option A. Reproduce the empirical analysis

Use the full empirical workflow script when the goal is to reproduce:

- the common valid support;
- the final row-standardised spatial weight matrix;
- SAR and SDM estimates;
- MLE, HC1, and SBS uncertainty measures;
- spatial econometric tables.

### Option B. Reproduce the benchmark

Use the benchmark script when the goal is to reproduce:

- Appendix A4 benchmark tables;
- comparison against `cell2nb + nb2mat`;
- comparison against `terra::adjacent`;
- support-level scalability results;
- feasibility boundary of the dense-output route.

### Option C. Use the workflow on another gridded dataset

To adapt the workflow to another raster system:

1. align all raster layers to a common grid;
2. construct the common valid support;
3. build sparse full-lattice contiguity;
4. project to the retained support;
5. remove zero-degree cells;
6. row-standardise the pruned matrix;
7. use the resulting `W_tilde` in spatial analysis.

## Reproducibility notes

The matrix-construction benchmark is intentionally aligned with the empirical workflow. Both use the same conceptual sequence:

```text
common valid support
-> sparse full-grid construction
-> projection to retained support
-> isolate removal
-> row-standardisation
```

The benchmark isolates matrix construction only. Estimation times reported in the manuscript exclude weight-matrix construction and data preparation. The reported SAR and SDM likelihood-estimation times therefore refer to the stage after the final spatial weight matrix and regression variables have already been prepared.

## Important implementation notes

Before making the repository public, check the following:

- remove or replace local machine paths;
- confirm that no private data files are included;
- rename scripts into clean public-facing filenames;
- add a license file;
- add a citation file if the code should be cited;
- add a small toy example if possible;
- confirm that package versions are documented;
- confirm that output folders are created automatically by the scripts.

## Suggested citation note

If using this repository, please cite the associated manuscript:

> Constructing Million-Cell Spatial Weight Matrices on Irregular Support: A Kronecker-Plus-Projection Workflow.

A formal citation can be added after the manuscript is accepted or assigned a DOI.

## Short description

This repository provides a sparse Kronecker-plus-projection workflow for constructing spatial weight matrices from large gridded raster data on irregular analytical support. It also provides benchmark code comparing the proposed workflow with a dense-output `spdep` route and a raster-native `terra::adjacent` route, together with scripts for spatial econometric estimation on the final non-isolated retained support.
