# ============================================================
# Benchmark v4: Proposed sparse Kronecker-projection pipeline
#               vs support-only sparse edge-list comparator
# ------------------------------------------------------------
# Purpose:
#   Compare two sparse matrix-construction strategies on the SAME
#   common support used in the v4 four-layer Indonesia process.
#
# Compared methods:
#   (1) Proposed pipeline:
#       W_full via sparse Kronecker -> P' W_full P -> drop isolates -> row-standardise
#
#   (2) Support-only sparse edge-list comparator:
#       retained indices -> row-column coordinates -> neighbour offset matching ->
#       Matrix::sparseMatrix -> drop isolates -> row-standardise
#
# Core alignment with the manuscript and previous benchmark code:
#   1) one common support across GDP, population/LandScan, NDVI, MNDWI
#   2) W_full is built by sparse Kronecker on the full template
#   3) projection matrix P has dimension N x V
#   4) W_common = t(P) %*% W_full %*% P
#   5) isolates are dropped BEFORE row-standardisation
#   6) edge-list benchmark is valid only after output-equivalence validation
#
# Important:
#   This script intentionally removes the cell2nb + nb2mat baseline.
#   It is intended to produce a stronger sparse-vs-sparse benchmark.
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(terra)
  library(sf)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ============================================================
# 1. Core sparse matrix utilities
# ============================================================

band_off1 <- function(n) {
  Matrix::bandSparse(
    n,
    k = c(-1, 1),
    diag = list(rep(1, n - 1), rep(1, n - 1)),
    symmetric = FALSE
  )
}

build_W_full <- function(r, c, type = c("queen", "rook")) {
  type <- match.arg(type)
  Ir <- Matrix::Diagonal(r)
  Ic <- Matrix::Diagonal(c)
  Br <- band_off1(r)
  Bc <- band_off1(c)

  # Row-major indexing: index = (row - 1) * c + col
  W <- kronecker(Ir, Bc) + kronecker(Br, Ic)
  if (type == "queen") W <- W + kronecker(Br, Bc)
  as(Matrix::drop0(W), "dgCMatrix")
}

row_standardize_drop0 <- function(W) {
  W <- as(W, "dgCMatrix")
  rs <- Matrix::rowSums(W)
  keep <- rs > 0
  if (!all(keep)) W <- W[keep, keep, drop = FALSE]
  rs <- Matrix::rowSums(W)
  if (length(rs) > 0 && any(rs <= 0)) {
    stop("Zero row sums remain after isolate removal.")
  }
  W_std <- Matrix::Diagonal(x = 1 / rs) %*% W
  list(W = as(Matrix::drop0(W_std), "dgCMatrix"), keep = keep)
}

cell_to_rc <- function(idx, nc) {
  idx <- as.integer(idx)
  data.frame(
    row = ((idx - 1L) %/% nc) + 1L,
    col = ((idx - 1L) %% nc) + 1L
  )
}

get_offsets <- function(contiguity = c("queen", "rook")) {
  contiguity <- match.arg(contiguity)
  rook <- rbind(
    c(-1L,  0L),
    c( 1L,  0L),
    c( 0L, -1L),
    c( 0L,  1L)
  )
  if (contiguity == "rook") return(rook)
  rbind(
    rook,
    c(-1L, -1L),
    c(-1L,  1L),
    c( 1L, -1L),
    c( 1L,  1L)
  )
}

# ============================================================
# 2. Raster and support construction utilities
#    Copied/aligned with previous v4 framework
# ============================================================

raster_to_colmajor <- function(x,
                               template = NULL,
                               method = "near",
                               var_name = deparse(substitute(x)),
                               already_aligned = FALSE) {
  if (!isTRUE(already_aligned)) {
    if (is.null(template)) {
      stop(sprintf("A template must be supplied when '%s' is not already aligned.", var_name))
    }
    x_use <- terra::resample(x, template, method = method)
  } else {
    x_use <- x
    if (is.null(template)) template <- x_use
  }

  if (terra::nlyr(x_use) != 1L) {
    stop(sprintf("'%s' must be a single-layer SpatRaster after alignment.", var_name))
  }

  if (terra::nrow(x_use) != terra::nrow(template) ||
      terra::ncol(x_use) != terra::ncol(template)) {
    stop(sprintf("Geometry mismatch for '%s' after alignment.", var_name))
  }

  vals <- terra::values(x_use, mat = FALSE)
  expected_n <- terra::ncell(template)

  if (length(vals) != expected_n) {
    stop(sprintf(
      "Value-length mismatch for '%s': expected %d cells but found %d values.",
      var_name, expected_n, length(vals)
    ))
  }

  M <- matrix(vals,
              nrow = terra::nrow(template),
              ncol = terra::ncol(template),
              byrow = TRUE)

  # This matches the previous v4 code and the row-major full-grid indexing.
  as.numeric(t(M))
}

crop_mask <- function(x, poly) {
  terra::mask(terra::crop(x, terra::vect(poly)), terra::vect(poly))
}

merge_rasters_safe <- function(r_list) {
  stopifnot(length(r_list) >= 1L)
  if (length(r_list) == 1L) return(r_list[[1L]])
  Reduce(function(a, b) terra::merge(a, b), r_list)
}

region_source_tag <- c(
  Sumatra = "sumatra",
  Java = "Java_Bali",
  Kalimantan = "kalimantan",
  Sulawesi = "sulawesi",
  Papua = "Eastern_Indonesia",
  Bali_Nusa_Tenggara = "Bali_NTB_NTT",
  Maluku = "Eastern_Indonesia"
)

build_region_polygons <- function(shp_admin2) {
  if (!("NAME_1" %in% names(shp_admin2))) {
    stop("The shapefile must contain a NAME_1 field for province names.")
  }

  province_groups <- list(
    Sumatra = c(
      "Aceh", "Sumatera Utara", "Sumatera Barat", "Jambi", "Riau",
      "Bengkulu", "Sumatera Selatan", "Lampung",
      "Bangka Belitung", "Kepulauan Riau"
    ),
    Java = c(
      "Banten", "Jakarta Raya", "Jawa Barat", "Jawa Tengah",
      "Yogyakarta", "Jawa Timur"
    ),
    Kalimantan = c(
      "Kalimantan Timur", "Kalimantan Barat", "Kalimantan Tengah",
      "Kalimantan Selatan", "Kalimantan Utara"
    ),
    Sulawesi = c(
      "Gorontalo", "Sulawesi Utara", "Sulawesi Barat",
      "Sulawesi Tengah", "Sulawesi Selatan", "Sulawesi Tenggara"
    ),
    Papua = c("Papua", "Papua Barat"),
    Bali_Nusa_Tenggara = c("Bali", "Nusa Tenggara Barat", "Nusa Tenggara Timur"),
    Maluku = c("Maluku", "Maluku Utara")
  )

  available_names <- unique(as.character(shp_admin2$NAME_1))
  missing_by_group <- lapply(province_groups, function(pnames) setdiff(pnames, available_names))
  missing_by_group <- missing_by_group[lengths(missing_by_group) > 0L]

  if (length(missing_by_group) > 0L) {
    msg <- paste(
      vapply(names(missing_by_group), function(nm) {
        paste0(nm, ": ", paste(missing_by_group[[nm]], collapse = ", "))
      }, character(1)),
      collapse = " | "
    )
    stop(paste(
      "Province name mismatch detected between the hard-coded region groups and the shapefile NAME_1 field.",
      "Please harmonize the province names before proceeding.",
      msg
    ))
  }

  polys <- lapply(names(province_groups), function(region_nm) {
    shp_admin2[shp_admin2$NAME_1 %in% province_groups[[region_nm]], ]
  })
  names(polys) <- names(province_groups)
  polys$Indonesia <- shp_admin2
  polys
}

prepare_data_sources <- function(cfg, year = 2020) {
  yr <- as.character(year)
  list(
    gdp = terra::rast(cfg$gdp_file),
    landscan = setNames(
      list(terra::rast(file.path(cfg$landscan_root, yr, paste0("landscan-global-", yr, ".tif")))),
      yr
    )
  )
}

load_single_region_inputs <- function(region_name, year, polys, sources, cfg) {
  if (identical(region_name, "Indonesia")) {
    stop("load_single_region_inputs() is for non-Indonesia island groups only.")
  }

  poly <- polys[[region_name]]
  tag  <- region_source_tag[[region_name]]
  yr   <- as.character(year)

  gdp_layer_name <- paste0("gdp_", year)
  if (!(gdp_layer_name %in% names(sources$gdp))) {
    stop(sprintf("Layer '%s' was not found in %s.", gdp_layer_name, cfg$gdp_file))
  }

  rast_gdp <- crop_mask(sources$gdp[[gdp_layer_name]], poly)
  rast_pop <- crop_mask(sources$landscan[[yr]], poly)

  ndvi_file <- file.path(cfg$data_root, "Indonesia_NDVI_MODIS", paste0("Annual_Index_", tag, "_", year, ".tif"))
  mndwi_file <- file.path(cfg$data_root, "Indonesia_MNDWI", paste0("MNDWI_", tag, "_", year, ".tif"))

  if (!file.exists(ndvi_file)) stop(sprintf("NDVI file not found: %s", ndvi_file))
  if (!file.exists(mndwi_file)) stop(sprintf("MNDWI file not found: %s", mndwi_file))

  rast_ndvi  <- crop_mask(terra::rast(ndvi_file), poly)
  rast_mndwi <- crop_mask(terra::rast(mndwi_file), poly)

  list(
    template = rast_gdp,
    gdp = rast_gdp,
    pop = rast_pop,
    ndvi = rast_ndvi,
    mndwi = rast_mndwi
  )
}

prepare_region_inputs_standalone <- function(cfg,
                                             year = 2020,
                                             cache_to_global = TRUE,
                                             global_name = "region_inputs") {
  required_cfg <- c("shapefile_admin2", "data_root", "landscan_root", "gdp_file")
  missing_cfg <- setdiff(required_cfg, names(cfg))
  if (length(missing_cfg) > 0L) {
    stop(sprintf("cfg is missing required fields: %s", paste(missing_cfg, collapse = ", ")))
  }

  if (!file.exists(cfg$shapefile_admin2)) stop(sprintf("Shapefile not found: %s", cfg$shapefile_admin2))
  if (!file.exists(cfg$gdp_file)) stop(sprintf("GDP raster not found: %s", cfg$gdp_file))

  ina_shp_lv2 <- sf::st_read(cfg$shapefile_admin2, quiet = TRUE)
  region_polys <- build_region_polygons(ina_shp_lv2)
  sources <- prepare_data_sources(cfg, year = year)

  core_regions <- c("Sumatra", "Java", "Kalimantan", "Sulawesi", "Papua", "Bali_Nusa_Tenggara", "Maluku")
  region_inputs <- setNames(vector("list", length(core_regions)), core_regions)

  for (rg in core_regions) {
    message(sprintf("Preparing region inputs for %s ...", rg))
    region_inputs[[rg]] <- load_single_region_inputs(rg, year, region_polys, sources, cfg)
  }

  region_inputs$Indonesia <- list(
    template = merge_rasters_safe(lapply(core_regions, function(rg) region_inputs[[rg]]$gdp)),
    gdp      = merge_rasters_safe(lapply(core_regions, function(rg) region_inputs[[rg]]$gdp)),
    pop      = merge_rasters_safe(lapply(core_regions, function(rg) region_inputs[[rg]]$pop)),
    ndvi     = merge_rasters_safe(lapply(core_regions, function(rg) region_inputs[[rg]]$ndvi)),
    mndwi    = merge_rasters_safe(lapply(core_regions, function(rg) region_inputs[[rg]]$mndwi))
  )

  if (isTRUE(cache_to_global)) assign(global_name, region_inputs, envir = .GlobalEnv)
  region_inputs
}

get_indonesia_inputs_v4_aligned <- function() {
  for (nm in c("region_inputs", "region_inputs_2020")) {
    if (exists(nm, inherits = TRUE)) {
      rr <- get(nm, inherits = TRUE)
      if (is.list(rr) && !is.null(rr$Indonesia) &&
          all(c("template", "gdp", "pop", "ndvi", "mndwi") %in% names(rr$Indonesia))) {
        return(rr$Indonesia)
      }
    }
  }

  needed <- c(
    "mask_indonesia_gdp_2020",
    "mask_Landscan_indonesia_2020",
    "mask_NDVI_indonesia_2020",
    "mask_MNDWI_indonesia_2020"
  )
  ok <- vapply(needed, exists, logical(1), inherits = TRUE)
  if (!all(ok)) {
    stop(paste0(
      "Could not find region_inputs$Indonesia/region_inputs_2020$Indonesia or the four mask objects: ",
      paste(needed[!ok], collapse = ", ")
    ))
  }

  list(
    template = get("mask_indonesia_gdp_2020", inherits = TRUE),
    gdp      = get("mask_indonesia_gdp_2020", inherits = TRUE),
    pop      = get("mask_Landscan_indonesia_2020", inherits = TRUE),
    ndvi     = get("mask_NDVI_indonesia_2020", inherits = TRUE),
    mndwi    = get("mask_MNDWI_indonesia_2020", inherits = TRUE)
  )
}

build_common_support_v4 <- function(template,
                                    rast_gdp, rast_pop, rast_ndvi, rast_mndwi,
                                    resample_pop = "near",
                                    resample_ndvi = "near",
                                    resample_mndwi = "near") {
  rast_gdp_a   <- terra::resample(rast_gdp,   template, method = "near")
  rast_pop_a   <- terra::resample(rast_pop,   template, method = resample_pop)
  rast_ndvi_a  <- terra::resample(rast_ndvi,  template, method = resample_ndvi)
  rast_mndwi_a <- terra::resample(rast_mndwi, template, method = resample_mndwi)

  r <- terra::nrow(template)
  c <- terra::ncol(template)
  N <- r * c

  gdp_vec   <- raster_to_colmajor(rast_gdp_a,   template, var_name = "rast_gdp",   already_aligned = TRUE)
  pop_vec   <- raster_to_colmajor(rast_pop_a,   template, var_name = "rast_pop",   already_aligned = TRUE)
  ndvi_vec  <- raster_to_colmajor(rast_ndvi_a,  template, var_name = "rast_ndvi",  already_aligned = TRUE)
  mndwi_vec <- raster_to_colmajor(rast_mndwi_a, template, var_name = "rast_mndwi", already_aligned = TRUE)

  valid_all <- !is.na(gdp_vec) & !is.na(pop_vec) & !is.na(ndvi_vec) & !is.na(mndwi_vec)
  idx_all <- which(valid_all)

  list(
    template = template,
    gdp = rast_gdp_a,
    pop = rast_pop_a,
    ndvi = rast_ndvi_a,
    mndwi = rast_mndwi_a,
    r = r,
    c = c,
    N = N,
    valid_all = valid_all,
    idx_all = idx_all,
    V = length(idx_all)
  )
}

precompute_support_v4_indonesia_standalone <- function(cfg,
                                                       year = 2020,
                                                       resample_pop = "near",
                                                       resample_ndvi = "near",
                                                       resample_mndwi = "near",
                                                       cache_to_global = TRUE,
                                                       global_name = "region_inputs",
                                                       save_region_inputs_rds = NULL,
                                                       save_support_rds = NULL) {
  region_inputs <- prepare_region_inputs_standalone(
    cfg = cfg,
    year = year,
    cache_to_global = cache_to_global,
    global_name = global_name
  )

  rr <- region_inputs$Indonesia
  support_obj <- build_common_support_v4(
    template = rr$template,
    rast_gdp = rr$gdp,
    rast_pop = rr$pop,
    rast_ndvi = rr$ndvi,
    rast_mndwi = rr$mndwi,
    resample_pop = resample_pop,
    resample_ndvi = resample_ndvi,
    resample_mndwi = resample_mndwi
  )

  out <- list(region_inputs = region_inputs, inputs = rr, support = support_obj)

  if (!is.null(save_region_inputs_rds)) saveRDS(region_inputs, save_region_inputs_rds)
  if (!is.null(save_support_rds)) saveRDS(out, save_support_rds)
  out
}

load_precomputed_support_v4 <- function(support_rds) {
  obj <- readRDS(support_rds)
  if (is.null(obj$support) || is.null(obj$inputs)) {
    stop("The precomputed support file does not contain 'support' and 'inputs'.")
  }
  obj
}

# ============================================================
# 3. Scenario generation on the common support
# ============================================================

normalize_proportions <- function(p) {
  p <- as.numeric(p)
  if (all(p > 1)) p <- p / 100
  if (any(p <= 0 | p > 1)) stop("Proportions must lie in (0, 1].")
  unique(sort(p))
}

scenario_label_pct <- function(p) {
  pct <- 100 * as.numeric(p)
  out <- ifelse(
    pct < 1,
    paste0(sub("\\.$", "", sub("0+$", "", formatC(pct, format = "f", digits = 3))), "pct"),
    paste0(sub("\\.$", "", sub("0+$", "", formatC(pct, format = "f", digits = 1))), "pct")
  )
  gsub("\\.", "p", out)
}

build_support_adjacency_for_scenarios <- function(idx_all, nr, nc, contiguity = c("queen", "rook")) {
  contiguity <- match.arg(contiguity)
  idx_all <- sort(unique(as.integer(idx_all)))
  rc <- cell_to_rc(idx_all, nc)
  V <- length(idx_all)
  N <- as.integer(nr) * as.integer(nc)

  lookup <- integer(N)
  lookup[idx_all] <- seq_len(V)
  offsets <- get_offsets(contiguity)

  adj <- vector("list", V)
  for (i in seq_len(V)) {
    rr <- rc$row[i]
    cc <- rc$col[i]
    cand_r <- rr + offsets[, 1]
    cand_c <- cc + offsets[, 2]
    ok <- cand_r >= 1L & cand_r <= nr & cand_c >= 1L & cand_c <= nc
    if (!any(ok)) {
      adj[[i]] <- integer(0)
      next
    }
    cand_idx <- (cand_r[ok] - 1L) * nc + cand_c[ok]
    pos <- lookup[cand_idx]
    adj[[i]] <- as.integer(pos[pos > 0L])
  }
  adj
}

bfs_order_adj <- function(adj, start) {
  n <- length(adj)
  seen <- rep(FALSE, n)
  queue <- integer(n)
  out <- integer(n)
  head <- 1L
  tail <- 1L
  out_len <- 0L

  queue[tail] <- start
  seen[start] <- TRUE

  while (head <= tail) {
    v <- queue[head]
    head <- head + 1L
    out_len <- out_len + 1L
    out[out_len] <- v

    nbrs <- adj[[v]]
    if (length(nbrs) > 0L) {
      new_nbrs <- nbrs[!seen[nbrs]]
      if (length(new_nbrs) > 0L) {
        seen[new_nbrs] <- TRUE
        m <- length(new_nbrs)
        queue[(tail + 1L):(tail + m)] <- new_nbrs
        tail <- tail + m
      }
    }
  }

  out[seq_len(out_len)]
}

connected_components_adj <- function(adj) {
  n <- length(adj)
  comp_id <- integer(n)
  comp <- 0L

  for (i in seq_len(n)) {
    if (comp_id[i] != 0L) next
    comp <- comp + 1L
    ord <- bfs_order_adj(adj, i)
    comp_id[ord] <- comp
  }

  list(n_comp = comp, comp_id = comp_id)
}

allocate_counts <- function(k, comp_sizes) {
  raw <- k * comp_sizes / sum(comp_sizes)
  alloc <- pmin(comp_sizes, floor(raw))
  rem <- k - sum(alloc)

  if (rem > 0L) {
    frac <- raw - floor(raw)
    ord <- order(frac, decreasing = TRUE)
    for (j in ord) {
      if (rem == 0L) break
      if (alloc[j] < comp_sizes[j]) {
        alloc[j] <- alloc[j] + 1L
        rem <- rem - 1L
      }
    }
  }

  if (rem > 0L) {
    spare <- comp_sizes - alloc
    ord2 <- order(spare, decreasing = TRUE)
    for (j in ord2) {
      if (rem == 0L) break
      add <- min(spare[j], rem)
      alloc[j] <- alloc[j] + add
      rem <- rem - add
    }
  }

  alloc
}

make_random_nested_supports <- function(idx_all, proportions, seed = 42) {
  proportions <- normalize_proportions(proportions)
  idx_all <- sort(unique(as.integer(idx_all)))
  set.seed(seed)
  ord <- sample(idx_all, length(idx_all), replace = FALSE)
  out <- vector("list", length(proportions))
  names(out) <- scenario_label_pct(proportions)
  for (i in seq_along(proportions)) {
    k <- max(1L, floor(length(ord) * proportions[i]))
    out[[i]] <- sort(ord[seq_len(k)])
  }
  out
}

make_spatial_nested_supports <- function(idx_all, nr, nc,
                                         proportions,
                                         contiguity = c("queen", "rook"),
                                         seed = 42) {
  contiguity <- match.arg(contiguity)
  proportions <- normalize_proportions(proportions)
  idx_all <- sort(unique(as.integer(idx_all)))
  adj <- build_support_adjacency_for_scenarios(idx_all, nr, nc, contiguity = contiguity)
  comp_info <- connected_components_adj(adj)
  comp_positions <- split(seq_along(idx_all), comp_info$comp_id)
  comp_sizes <- lengths(comp_positions)

  coords <- cell_to_rc(idx_all, nc)
  orders <- vector("list", length(comp_positions))
  set.seed(seed)

  for (j in seq_along(comp_positions)) {
    pos <- comp_positions[[j]]
    rc <- coords[pos, , drop = FALSE]
    ctr_row <- mean(rc$row)
    ctr_col <- mean(rc$col)
    seed_local <- pos[which.min((rc$row - ctr_row)^2 + (rc$col - ctr_col)^2)]

    local_index <- integer(length(idx_all))
    local_index[pos] <- seq_along(pos)
    adj_local <- lapply(pos, function(p) local_index[adj[[p]][adj[[p]] %in% pos]])
    ord_local <- bfs_order_adj(adj_local, start = local_index[seed_local])
    orders[[j]] <- pos[ord_local]
  }

  out <- vector("list", length(proportions))
  names(out) <- scenario_label_pct(proportions)
  N_all <- length(idx_all)

  for (i in seq_along(proportions)) {
    k <- max(1L, floor(N_all * proportions[i]))
    alloc <- allocate_counts(k, comp_sizes)
    sel_pos <- unlist(
      mapply(
        FUN = function(ord, m) if (m <= 0L) integer(0) else ord[seq_len(m)],
        orders, alloc,
        SIMPLIFY = FALSE
      ),
      use.names = FALSE
    )
    out[[i]] <- sort(idx_all[sel_pos])
  }

  out
}

prepare_scenarios_from_support_v4 <- function(support_obj,
                                              proportions = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00),
                                              contiguity = c("queen", "rook"),
                                              seed = 42,
                                              scenario_mode = c("random", "spatial")) {
  contiguity <- match.arg(contiguity)
  scenario_mode <- match.arg(scenario_mode)

  if (support_obj$V == 0L) stop("Common support is empty. Check raster alignment and masks.")

  scenarios <- switch(
    scenario_mode,
    random = make_random_nested_supports(
      idx_all = support_obj$idx_all,
      proportions = proportions,
      seed = seed
    ),
    spatial = make_spatial_nested_supports(
      idx_all = support_obj$idx_all,
      nr = support_obj$r,
      nc = support_obj$c,
      proportions = proportions,
      contiguity = contiguity,
      seed = seed
    )
  )

  list(
    scenarios = scenarios,
    scenario_sizes = data.frame(
      scenario = names(scenarios),
      n_cells = vapply(scenarios, length, integer(1)),
      pct_of_full_common_support = round(100 * vapply(scenarios, length, integer(1)) / support_obj$V, 4),
      stringsAsFactors = FALSE
    )
  )
}

prepare_benchmark_v4_indonesia_standalone <- function(cfg,
                                                      year = 2020,
                                                      proportions = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00),
                                                      contiguity = c("queen", "rook"),
                                                      seed = 42,
                                                      scenario_mode = c("random", "spatial"),
                                                      resample_pop = "near",
                                                      resample_ndvi = "near",
                                                      resample_mndwi = "near",
                                                      cache_to_global = TRUE,
                                                      global_name = "region_inputs",
                                                      support_rds = NULL,
                                                      reuse_precomputed = TRUE,
                                                      save_precomputed = TRUE,
                                                      force_rebuild = FALSE,
                                                      save_region_inputs_rds = NULL) {
  contiguity <- match.arg(contiguity)
  scenario_mode <- match.arg(scenario_mode)

  if (!is.null(support_rds) && reuse_precomputed && file.exists(support_rds) && !isTRUE(force_rebuild)) {
    prep_core <- load_precomputed_support_v4(support_rds)
  } else {
    prep_core <- precompute_support_v4_indonesia_standalone(
      cfg = cfg,
      year = year,
      resample_pop = resample_pop,
      resample_ndvi = resample_ndvi,
      resample_mndwi = resample_mndwi,
      cache_to_global = cache_to_global,
      global_name = global_name,
      save_region_inputs_rds = save_region_inputs_rds,
      save_support_rds = if (isTRUE(save_precomputed)) support_rds else NULL
    )
  }

  sc <- prepare_scenarios_from_support_v4(
    support_obj = prep_core$support,
    proportions = proportions,
    contiguity = contiguity,
    seed = seed,
    scenario_mode = scenario_mode
  )

  c(prep_core, sc)
}

prepare_benchmark_v4_indonesia <- function(proportions = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00),
                                           contiguity = c("queen", "rook"),
                                           seed = 42,
                                           scenario_mode = c("random", "spatial"),
                                           resample_pop = "near",
                                           resample_ndvi = "near",
                                           resample_mndwi = "near") {
  contiguity <- match.arg(contiguity)
  scenario_mode <- match.arg(scenario_mode)

  rr <- get_indonesia_inputs_v4_aligned()
  support_obj <- build_common_support_v4(
    template = rr$template,
    rast_gdp = rr$gdp,
    rast_pop = rr$pop,
    rast_ndvi = rr$ndvi,
    rast_mndwi = rr$mndwi,
    resample_pop = resample_pop,
    resample_ndvi = resample_ndvi,
    resample_mndwi = resample_mndwi
  )

  sc <- prepare_scenarios_from_support_v4(
    support_obj = support_obj,
    proportions = proportions,
    contiguity = contiguity,
    seed = seed,
    scenario_mode = scenario_mode
  )

  list(inputs = rr, support = support_obj, scenarios = sc$scenarios, scenario_sizes = sc$scenario_sizes)
}

# ============================================================
# 4. General benchmarking helpers
# ============================================================

bytes_num <- function(x) {
  if (is.null(x)) return(NA_real_)
  as.numeric(object.size(x))
}

mb_num <- function(x) x / 1024^2

max_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}

median_safe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  stats::median(x)
}

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(NA_character_)
  paste(x, collapse = " | ")
}

bind_rows_fill <- function(lst) {
  all_names <- unique(unlist(lapply(lst, names)))
  lst2 <- lapply(lst, function(x) {
    miss <- setdiff(all_names, names(x))
    if (length(miss) > 0L) {
      for (nm in miss) x[[nm]] <- NA
    }
    x[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, lst2)
  rownames(out) <- NULL
  out
}

critical_warning_detected <- function(warnings,
                                      patterns = c("integer overflow", "cannot allocate", "long vectors",
                                                   "negative length vectors", "memory", "overflow")) {
  if (length(warnings) == 0L) return(FALSE)
  any(vapply(patterns, function(p) any(grepl(p, warnings, ignore.case = TRUE)), logical(1)))
}

classify_failure_v4 <- function(error_message = NA_character_,
                                warning_message = NA_character_,
                                hard_patterns = c(
                                  "cannot allocate",
                                  "memory exhausted",
                                  "out of memory",
                                  "vector memory exhausted",
                                  "negative length vectors",
                                  "long vectors not supported",
                                  "integer overflow",
                                  "timeout",
                                  "reached elapsed time limit",
                                  "allocation failure",
                                  "std::bad_alloc"
                                ),
                                soft_patterns = c(
                                  "subscript out of bounds",
                                  "dimnames",
                                  "names do not match",
                                  "invalid.*argument",
                                  "missing value where true/false needed"
                                )) {
  msg <- paste(c(error_message, warning_message), collapse = " | ")
  msg_low <- tolower(msg)

  hard_hit <- any(grepl(paste(hard_patterns, collapse = "|"), msg_low, perl = TRUE))
  soft_hit <- any(grepl(paste(soft_patterns, collapse = "|"), msg_low, perl = TRUE))

  if (hard_hit) return("hard")
  if (soft_hit) return("soft")
  "soft"
}

timed_eval <- function(expr, timeout_seconds = Inf) {
  expr_sub <- substitute(expr)
  start <- proc.time()[["elapsed"]]
  warns <- character(0)

  if (is.finite(timeout_seconds)) setTimeLimit(cpu = Inf, elapsed = timeout_seconds, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)

  out <- tryCatch(
    withCallingHandlers(
      eval.parent(expr_sub),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  elapsed <- proc.time()[["elapsed"]] - start

  if (inherits(out, "error")) {
    list(ok = FALSE, value = NULL, elapsed = elapsed, error = conditionMessage(out), warnings = warns)
  } else {
    list(ok = TRUE, value = out, elapsed = elapsed, error = NA_character_, warnings = warns)
  }
}

make_fail_row <- function(method, stage, cells_input,
                          error_message = NA_character_,
                          warning_message = NA_character_,
                          elapsed_total = NA_real_) {
  data.frame(
    method = method,
    status = "failed",
    failure_stage = stage,
    failure_class = classify_failure_v4(error_message, warning_message),
    error = error_message,
    warning_message = warning_message,
    cells_input = cells_input,
    cells_after_prune = NA_real_,
    n_isolates = NA_real_,
    nnz_common = NA_real_,
    nnz_final = NA_real_,
    elapsed_total = elapsed_total,
    stringsAsFactors = FALSE
  )
}

sparse_compare_fast <- function(A, B, tol = 1e-12, compute_diff_nnz = FALSE) {
  if (is.null(A) || is.null(B)) {
    return(data.frame(
      dims_equal = FALSE,
      nnz_A = NA_real_,
      nnz_B = NA_real_,
      nnz_equal = FALSE,
      structure_equal = FALSE,
      max_abs_diff = NA_real_,
      diff_nnz = NA_real_,
      equal_within_tol = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  A <- as(A, "dgCMatrix")
  B <- as(B, "dgCMatrix")
  dims_equal <- identical(dim(A), dim(B))
  nnz_A <- length(A@x)
  nnz_B <- length(B@x)
  nnz_equal <- identical(nnz_A, nnz_B)

  structure_equal <- dims_equal && identical(A@p, B@p) && identical(A@i, B@i)
  max_abs_diff <- NA_real_
  diff_nnz <- NA_real_

  if (structure_equal) {
    max_abs_diff <- if (nnz_A == 0L) 0 else max(abs(A@x - B@x))
  }

  if (isTRUE(compute_diff_nnz) && dims_equal) {
    D <- Matrix::drop0(A - B, tol = tol)
    diff_nnz <- Matrix::nnzero(D)
    max_abs_diff <- if (length(D@x) == 0L) 0 else max(abs(D@x))
  }

  equal_within_tol <- if (isTRUE(compute_diff_nnz) && dims_equal) {
    diff_nnz == 0L
  } else {
    structure_equal && is.finite(max_abs_diff) && max_abs_diff <= tol
  }

  data.frame(
    dims_equal = dims_equal,
    nnz_A = nnz_A,
    nnz_B = nnz_B,
    nnz_equal = nnz_equal,
    structure_equal = structure_equal,
    max_abs_diff = max_abs_diff,
    diff_nnz = diff_nnz,
    equal_within_tol = equal_within_tol,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 5. Proposed method: sparse Kronecker + projection
# ============================================================

run_proposed_once_v4 <- function(nr, nc, idx_support,
                                 contiguity = c("queen", "rook"),
                                 timeout_seconds = Inf,
                                 cached_W_full = NULL,
                                 return_matrices = FALSE,
                                 critical_warning_patterns = c("integer overflow", "cannot allocate", "memory", "overflow")) {
  contiguity <- match.arg(contiguity)
  idx_support <- sort(unique(as.integer(idx_support)))
  N <- as.integer(nr) * as.integer(nc)
  V <- length(idx_support)

  if (V == 0L) stop("idx_support is empty.")
  if (any(idx_support < 1L | idx_support > N)) stop("idx_support contains out-of-range indices.")

  build_step <- if (is.null(cached_W_full)) {
    timed_eval(build_W_full(nr, nc, type = contiguity), timeout_seconds)
  } else {
    list(ok = TRUE, value = cached_W_full, elapsed = 0, error = NA_character_, warnings = character(0))
  }

  if (!build_step$ok || critical_warning_detected(build_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "proposed_kronecker_projection",
      stage = "build_W_full",
      cells_input = V,
      error_message = if (!build_step$ok) build_step$error else paste(build_step$warnings, collapse = " | "),
      warning_message = collapse_unique(build_step$warnings),
      elapsed_total = build_step$elapsed
    )
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  W_full <- as(build_step$value, "dgCMatrix")

  proj_step <- timed_eval({
    P <- Matrix::sparseMatrix(i = idx_support, j = seq_len(V), x = 1, dims = c(N, V))
    W_common <- as(Matrix::drop0(Matrix::t(P) %*% W_full %*% P), "dgCMatrix")
    list(P = P, W_common = W_common)
  }, timeout_seconds)

  if (!proj_step$ok || critical_warning_detected(proj_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "proposed_kronecker_projection",
      stage = "projection",
      cells_input = V,
      error_message = if (!proj_step$ok) proj_step$error else paste(proj_step$warnings, collapse = " | "),
      warning_message = collapse_unique(c(build_step$warnings, proj_step$warnings)),
      elapsed_total = build_step$elapsed + proj_step$elapsed
    )
    df$elapsed_build_full <- build_step$elapsed
    df$elapsed_projection <- proj_step$elapsed
    df$bytes_W_full <- bytes_num(W_full)
    df$peak_object_bytes <- bytes_num(W_full)
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  P <- proj_step$value$P
  W_common_raw <- proj_step$value$W_common

  prune_step <- timed_eval(row_standardize_drop0(W_common_raw), timeout_seconds)

  if (!prune_step$ok || critical_warning_detected(prune_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "proposed_kronecker_projection",
      stage = "prune_standardise",
      cells_input = V,
      error_message = if (!prune_step$ok) prune_step$error else paste(prune_step$warnings, collapse = " | "),
      warning_message = collapse_unique(c(build_step$warnings, proj_step$warnings, prune_step$warnings)),
      elapsed_total = build_step$elapsed + proj_step$elapsed + prune_step$elapsed
    )
    df$elapsed_build_full <- build_step$elapsed
    df$elapsed_projection <- proj_step$elapsed
    df$elapsed_prune_std <- prune_step$elapsed
    df$bytes_W_full <- bytes_num(W_full)
    df$bytes_P <- bytes_num(P)
    df$bytes_W_common <- bytes_num(W_common_raw)
    df$peak_object_bytes <- max_na(c(bytes_num(W_full), bytes_num(P), bytes_num(W_common_raw)))
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  W_final <- prune_step$value$W
  keep <- prune_step$value$keep

  df <- data.frame(
    method = "proposed_kronecker_projection",
    status = "ok",
    failure_stage = NA_character_,
    failure_class = NA_character_,
    error = NA_character_,
    warning_message = collapse_unique(c(build_step$warnings, proj_step$warnings, prune_step$warnings)),
    cells_input = V,
    cells_after_prune = nrow(W_final),
    n_isolates = sum(!keep),
    nnz_common = length(W_common_raw@x),
    nnz_final = length(W_final@x),
    elapsed_build_full = build_step$elapsed,
    elapsed_projection = proj_step$elapsed,
    elapsed_prune_std = prune_step$elapsed,
    elapsed_total = build_step$elapsed + proj_step$elapsed + prune_step$elapsed,
    bytes_W_full = bytes_num(W_full),
    bytes_P = bytes_num(P),
    bytes_W_common = bytes_num(W_common_raw),
    bytes_W_final = bytes_num(W_final),
    peak_object_bytes = max_na(c(bytes_num(W_full), bytes_num(P), bytes_num(W_common_raw), bytes_num(W_final))),
    stringsAsFactors = FALSE
  )

  list(
    result = df,
    W_common = if (isTRUE(return_matrices)) W_common_raw else NULL,
    W_final = if (isTRUE(return_matrices)) W_final else NULL,
    keep = if (isTRUE(return_matrices)) keep else NULL
  )
}

# ============================================================
# 6. Support-only sparse edge-list comparator
# ============================================================

build_W_edge_support_only <- function(nr, nc, idx_support,
                                      contiguity = c("queen", "rook"),
                                      use_integer_lookup = TRUE) {
  contiguity <- match.arg(contiguity)
  idx_support <- sort(unique(as.integer(idx_support)))
  N <- as.integer(nr) * as.integer(nc)
  V <- length(idx_support)
  if (V == 0L) stop("idx_support is empty.")
  if (any(idx_support < 1L | idx_support > N)) stop("idx_support contains out-of-range indices.")

  rc <- cell_to_rc(idx_support, nc)
  rows <- as.integer(rc$row)
  cols <- as.integer(rc$col)
  offsets <- get_offsets(contiguity)

  # Integer lookup is faster and avoids string-name lookup for large supports.
  # lookup[j] gives retained-support position of full-grid cell j, or 0 if absent.
  lookup <- integer(N)
  lookup[idx_support] <- seq_len(V)

  i_list <- vector("list", nrow(offsets))
  j_list <- vector("list", nrow(offsets))

  src_all <- seq_len(V)

  for (m in seq_len(nrow(offsets))) {
    a <- offsets[m, 1]
    b <- offsets[m, 2]

    cand_r <- rows + a
    cand_c <- cols + b

    inside <- cand_r >= 1L & cand_r <= nr & cand_c >= 1L & cand_c <= nc
    if (!any(inside)) {
      i_list[[m]] <- integer(0)
      j_list[[m]] <- integer(0)
      next
    }

    src <- src_all[inside]
    cand_idx <- as.integer((cand_r[inside] - 1L) * nc + cand_c[inside])
    dest <- lookup[cand_idx]
    ok <- dest > 0L

    i_list[[m]] <- as.integer(src[ok])
    j_list[[m]] <- as.integer(dest[ok])
  }

  ii <- unlist(i_list, use.names = FALSE)
  jj <- unlist(j_list, use.names = FALSE)

  W_edge <- Matrix::sparseMatrix(
    i = ii,
    j = jj,
    x = 1,
    dims = c(V, V),
    giveCsparse = TRUE
  )
  W_edge <- as(Matrix::drop0(W_edge), "dgCMatrix")

  list(
    W_edge = W_edge,
    lookup = lookup,
    rows = rows,
    cols = cols,
    edge_i = ii,
    edge_j = jj
  )
}

run_edge_once_v4 <- function(nr, nc, idx_support,
                             contiguity = c("queen", "rook"),
                             timeout_seconds = Inf,
                             return_matrices = FALSE,
                             critical_warning_patterns = c("integer overflow", "cannot allocate", "memory", "overflow")) {
  contiguity <- match.arg(contiguity)
  idx_support <- sort(unique(as.integer(idx_support)))
  N <- as.integer(nr) * as.integer(nc)
  V <- length(idx_support)

  if (V == 0L) stop("idx_support is empty.")
  if (any(idx_support < 1L | idx_support > N)) stop("idx_support contains out-of-range indices.")

  coord_step <- timed_eval({
    rc <- cell_to_rc(idx_support, nc)
    rows <- as.integer(rc$row)
    cols <- as.integer(rc$col)
    lookup <- integer(N)
    lookup[idx_support] <- seq_len(V)
    list(rows = rows, cols = cols, lookup = lookup)
  }, timeout_seconds)

  if (!coord_step$ok || critical_warning_detected(coord_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "edge_support_only_sparse",
      stage = "coord_lookup",
      cells_input = V,
      error_message = if (!coord_step$ok) coord_step$error else paste(coord_step$warnings, collapse = " | "),
      warning_message = collapse_unique(coord_step$warnings),
      elapsed_total = coord_step$elapsed
    )
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  rows <- coord_step$value$rows
  cols <- coord_step$value$cols
  lookup <- coord_step$value$lookup
  offsets <- get_offsets(contiguity)

  edge_step <- timed_eval({
    i_list <- vector("list", nrow(offsets))
    j_list <- vector("list", nrow(offsets))
    src_all <- seq_len(V)

    for (m in seq_len(nrow(offsets))) {
      a <- offsets[m, 1]
      b <- offsets[m, 2]
      cand_r <- rows + a
      cand_c <- cols + b
      inside <- cand_r >= 1L & cand_r <= nr & cand_c >= 1L & cand_c <= nc

      if (!any(inside)) {
        i_list[[m]] <- integer(0)
        j_list[[m]] <- integer(0)
        next
      }

      src <- src_all[inside]
      cand_idx <- as.integer((cand_r[inside] - 1L) * nc + cand_c[inside])
      dest <- lookup[cand_idx]
      ok <- dest > 0L

      i_list[[m]] <- as.integer(src[ok])
      j_list[[m]] <- as.integer(dest[ok])
    }

    list(
      edge_i = unlist(i_list, use.names = FALSE),
      edge_j = unlist(j_list, use.names = FALSE)
    )
  }, timeout_seconds)

  if (!edge_step$ok || critical_warning_detected(edge_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "edge_support_only_sparse",
      stage = "edge_matching",
      cells_input = V,
      error_message = if (!edge_step$ok) edge_step$error else paste(edge_step$warnings, collapse = " | "),
      warning_message = collapse_unique(c(coord_step$warnings, edge_step$warnings)),
      elapsed_total = coord_step$elapsed + edge_step$elapsed
    )
    df$elapsed_coord_lookup <- coord_step$elapsed
    df$elapsed_edge_match <- edge_step$elapsed
    df$bytes_lookup <- bytes_num(lookup)
    df$bytes_coord <- bytes_num(rows) + bytes_num(cols)
    df$peak_object_bytes <- max_na(c(bytes_num(lookup), bytes_num(rows) + bytes_num(cols)))
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  edge_i <- edge_step$value$edge_i
  edge_j <- edge_step$value$edge_j

  assembly_step <- timed_eval({
    W_edge <- Matrix::sparseMatrix(
      i = edge_i,
      j = edge_j,
      x = 1,
      dims = c(V, V),
      giveCsparse = TRUE
    )
    as(Matrix::drop0(W_edge), "dgCMatrix")
  }, timeout_seconds)

  if (!assembly_step$ok || critical_warning_detected(assembly_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "edge_support_only_sparse",
      stage = "sparse_assembly",
      cells_input = V,
      error_message = if (!assembly_step$ok) assembly_step$error else paste(assembly_step$warnings, collapse = " | "),
      warning_message = collapse_unique(c(coord_step$warnings, edge_step$warnings, assembly_step$warnings)),
      elapsed_total = coord_step$elapsed + edge_step$elapsed + assembly_step$elapsed
    )
    df$elapsed_coord_lookup <- coord_step$elapsed
    df$elapsed_edge_match <- edge_step$elapsed
    df$elapsed_sparse_assembly <- assembly_step$elapsed
    df$bytes_lookup <- bytes_num(lookup)
    df$bytes_coord <- bytes_num(rows) + bytes_num(cols)
    df$bytes_edge_index <- bytes_num(edge_i) + bytes_num(edge_j)
    df$peak_object_bytes <- max_na(c(bytes_num(lookup), bytes_num(rows) + bytes_num(cols), bytes_num(edge_i) + bytes_num(edge_j)))
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  W_edge_raw <- assembly_step$value

  prune_step <- timed_eval(row_standardize_drop0(W_edge_raw), timeout_seconds)

  if (!prune_step$ok || critical_warning_detected(prune_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "edge_support_only_sparse",
      stage = "prune_standardise",
      cells_input = V,
      error_message = if (!prune_step$ok) prune_step$error else paste(prune_step$warnings, collapse = " | "),
      warning_message = collapse_unique(c(coord_step$warnings, edge_step$warnings, assembly_step$warnings, prune_step$warnings)),
      elapsed_total = coord_step$elapsed + edge_step$elapsed + assembly_step$elapsed + prune_step$elapsed
    )
    df$elapsed_coord_lookup <- coord_step$elapsed
    df$elapsed_edge_match <- edge_step$elapsed
    df$elapsed_sparse_assembly <- assembly_step$elapsed
    df$elapsed_prune_std <- prune_step$elapsed
    df$bytes_lookup <- bytes_num(lookup)
    df$bytes_coord <- bytes_num(rows) + bytes_num(cols)
    df$bytes_edge_index <- bytes_num(edge_i) + bytes_num(edge_j)
    df$bytes_W_common <- bytes_num(W_edge_raw)
    df$peak_object_bytes <- max_na(c(
      bytes_num(lookup), bytes_num(rows) + bytes_num(cols), bytes_num(edge_i) + bytes_num(edge_j), bytes_num(W_edge_raw)
    ))
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  W_final <- prune_step$value$W
  keep <- prune_step$value$keep

  df <- data.frame(
    method = "edge_support_only_sparse",
    status = "ok",
    failure_stage = NA_character_,
    failure_class = NA_character_,
    error = NA_character_,
    warning_message = collapse_unique(c(coord_step$warnings, edge_step$warnings, assembly_step$warnings, prune_step$warnings)),
    cells_input = V,
    cells_after_prune = nrow(W_final),
    n_isolates = sum(!keep),
    nnz_common = length(W_edge_raw@x),
    nnz_final = length(W_final@x),
    elapsed_coord_lookup = coord_step$elapsed,
    elapsed_edge_match = edge_step$elapsed,
    elapsed_sparse_assembly = assembly_step$elapsed,
    elapsed_prune_std = prune_step$elapsed,
    elapsed_total = coord_step$elapsed + edge_step$elapsed + assembly_step$elapsed + prune_step$elapsed,
    bytes_lookup = bytes_num(lookup),
    bytes_coord = bytes_num(rows) + bytes_num(cols),
    bytes_edge_index = bytes_num(edge_i) + bytes_num(edge_j),
    bytes_W_common = bytes_num(W_edge_raw),
    bytes_W_final = bytes_num(W_final),
    peak_object_bytes = max_na(c(
      bytes_num(lookup),
      bytes_num(rows) + bytes_num(cols),
      bytes_num(edge_i) + bytes_num(edge_j),
      bytes_num(W_edge_raw),
      bytes_num(W_final)
    )),
    stringsAsFactors = FALSE
  )

  list(
    result = df,
    W_common = if (isTRUE(return_matrices)) W_edge_raw else NULL,
    W_final = if (isTRUE(return_matrices)) W_final else NULL,
    keep = if (isTRUE(return_matrices)) keep else NULL
  )
}

# ============================================================
# 7. Validation and summary helpers
# ============================================================

validate_proposed_vs_edge <- function(prop_obj, edge_obj,
                                      scenario = NA_character_,
                                      rep = NA_integer_,
                                      pct_support = NA_real_,
                                      tol = 1e-12,
                                      compute_diff_nnz = FALSE) {
  prop_ok <- is.data.frame(prop_obj$result) && identical(prop_obj$result$status[1], "ok")
  edge_ok <- is.data.frame(edge_obj$result) && identical(edge_obj$result$status[1], "ok")

  if (!prop_ok || !edge_ok) {
    return(data.frame(
      scenario = scenario,
      rep = rep,
      pct_support = pct_support,
      prop_ok = prop_ok,
      edge_ok = edge_ok,
      common_equal = NA,
      final_equal = NA,
      keep_equal = NA,
      output_equivalence_ok = FALSE,
      common_max_abs_diff = NA_real_,
      final_max_abs_diff = NA_real_,
      common_diff_nnz = NA_real_,
      final_diff_nnz = NA_real_,
      cells_after_prune_proposed = if (prop_ok) prop_obj$result$cells_after_prune[1] else NA_real_,
      cells_after_prune_edge = if (edge_ok) edge_obj$result$cells_after_prune[1] else NA_real_,
      nnz_common_proposed = if (prop_ok) prop_obj$result$nnz_common[1] else NA_real_,
      nnz_common_edge = if (edge_ok) edge_obj$result$nnz_common[1] else NA_real_,
      nnz_final_proposed = if (prop_ok) prop_obj$result$nnz_final[1] else NA_real_,
      nnz_final_edge = if (edge_ok) edge_obj$result$nnz_final[1] else NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  cmp_common <- sparse_compare_fast(
    prop_obj$W_common,
    edge_obj$W_common,
    tol = tol,
    compute_diff_nnz = compute_diff_nnz
  )
  cmp_final <- sparse_compare_fast(
    prop_obj$W_final,
    edge_obj$W_final,
    tol = tol,
    compute_diff_nnz = compute_diff_nnz
  )

  keep_equal <- identical(prop_obj$keep, edge_obj$keep)

  data.frame(
    scenario = scenario,
    rep = rep,
    pct_support = pct_support,
    prop_ok = TRUE,
    edge_ok = TRUE,
    common_equal = cmp_common$equal_within_tol[1],
    final_equal = cmp_final$equal_within_tol[1],
    keep_equal = keep_equal,
    output_equivalence_ok = isTRUE(cmp_common$equal_within_tol[1]) &&
      isTRUE(cmp_final$equal_within_tol[1]) &&
      isTRUE(keep_equal),
    common_dims_equal = cmp_common$dims_equal[1],
    common_structure_equal = cmp_common$structure_equal[1],
    common_max_abs_diff = cmp_common$max_abs_diff[1],
    common_diff_nnz = cmp_common$diff_nnz[1],
    final_dims_equal = cmp_final$dims_equal[1],
    final_structure_equal = cmp_final$structure_equal[1],
    final_max_abs_diff = cmp_final$max_abs_diff[1],
    final_diff_nnz = cmp_final$diff_nnz[1],
    cells_after_prune_proposed = prop_obj$result$cells_after_prune[1],
    cells_after_prune_edge = edge_obj$result$cells_after_prune[1],
    nnz_common_proposed = prop_obj$result$nnz_common[1],
    nnz_common_edge = edge_obj$result$nnz_common[1],
    nnz_final_proposed = prop_obj$result$nnz_final[1],
    nnz_final_edge = edge_obj$result$nnz_final[1],
    stringsAsFactors = FALSE
  )
}

summarize_validation <- function(validation_df) {
  if (is.null(validation_df) || nrow(validation_df) == 0L) return(data.frame())
  out <- lapply(split(validation_df, validation_df$scenario, drop = TRUE), function(d) {
    data.frame(
      scenario = d$scenario[1],
      pct_support = d$pct_support[1],
      runs = nrow(d),
      prop_ok_runs = sum(d$prop_ok, na.rm = TRUE),
      edge_ok_runs = sum(d$edge_ok, na.rm = TRUE),
      equivalence_ok_runs = sum(d$output_equivalence_ok, na.rm = TRUE),
      all_equivalent = all(d$output_equivalence_ok, na.rm = TRUE),
      median_common_max_abs_diff = median_safe(d$common_max_abs_diff),
      median_final_max_abs_diff = median_safe(d$final_max_abs_diff),
      median_cells_after_prune_proposed = median_safe(d$cells_after_prune_proposed),
      median_cells_after_prune_edge = median_safe(d$cells_after_prune_edge),
      median_nnz_common_proposed = median_safe(d$nnz_common_proposed),
      median_nnz_common_edge = median_safe(d$nnz_common_edge),
      median_nnz_final_proposed = median_safe(d$nnz_final_proposed),
      median_nnz_final_edge = median_safe(d$nnz_final_edge),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$pct_support), , drop = FALSE]
}

benchmark_metadata_proposed_edge <- function() {
  data.frame(
    variable_name = c(
      "scenario", "method", "pct_support", "runs", "ok_runs", "failed_runs",
      "cells_input", "cells_after_prune", "n_isolates", "nnz_common", "nnz_final",
      "elapsed_total", "elapsed_build_full", "elapsed_projection", "elapsed_coord_lookup",
      "elapsed_edge_match", "elapsed_sparse_assembly", "elapsed_prune_std",
      "bytes_W_full_mb", "bytes_P_mb", "bytes_lookup_mb", "bytes_coord_mb",
      "bytes_edge_index_mb", "bytes_W_common_mb", "bytes_W_final_mb", "peak_object_mb",
      "output_equivalence_ok", "common_max_abs_diff", "final_max_abs_diff"
    ),
    description = c(
      "Scenario label for retained-support benchmark level.",
      "Benchmark method identifier.",
      "Retained support as share of full common support.",
      "Number of runs combined in summary row.",
      "Number of successful runs.",
      "Number of failed runs.",
      "Number of cells entering matrix-construction stage before isolate removal.",
      "Number of cells remaining after isolate removal.",
      "Number of isolates removed after projection or support-only construction.",
      "Number of nonzero entries in unstandardised common-support matrix.",
      "Number of nonzero entries in final row-standardised matrix.",
      "Total elapsed time for method-specific matrix-construction workflow.",
      "Elapsed time for W_FULL construction; proposed method only.",
      "Elapsed time for projection t(P) %*% W_full %*% P; proposed method only.",
      "Elapsed time for row-column coordinate and lookup construction; edge method only.",
      "Elapsed time for neighbour-offset matching; edge method only.",
      "Elapsed time for Matrix::sparseMatrix assembly; edge method only.",
      "Elapsed time for isolate removal and row-standardisation.",
      "Object size of W_FULL in MB; proposed method only.",
      "Object size of projection matrix P in MB; proposed method only.",
      "Object size of integer lookup vector in MB; edge method only.",
      "Object size of retained row-column coordinate vectors in MB; edge method only.",
      "Object size of edge-index vectors in MB; edge method only.",
      "Object size of unstandardised common-support sparse matrix in MB.",
      "Object size of final row-standardised sparse matrix in MB.",
      "Largest tracked object size in MB.",
      "Validation flag: final benchmark comparison is valid only when output equivalence is TRUE.",
      "Maximum absolute difference between proposed and edge common-support matrices.",
      "Maximum absolute difference between proposed and edge final row-standardised matrices."
    ),
    stringsAsFactors = FALSE
  )
}

summarize_benchmark_raw_v4 <- function(raw) {
  summary_list <- lapply(split(raw, list(raw$scenario, raw$method), drop = TRUE), function(d) {
    data.frame(
      scenario = d$scenario[1],
      method = d$method[1],
      pct_support = d$pct_support[1],
      runs = nrow(d),
      ok_runs = sum(d$status == "ok", na.rm = TRUE),
      failed_runs = sum(d$status == "failed", na.rm = TRUE),
      cells_input = median_safe(d$cells_input),
      cells_after_prune = median_safe(d$cells_after_prune),
      n_isolates = median_safe(d$n_isolates),
      nnz_common = median_safe(d$nnz_common),
      nnz_final = median_safe(d$nnz_final),
      elapsed_total = median_safe(d$elapsed_total),
      elapsed_build_full = median_safe(d$elapsed_build_full),
      elapsed_projection = median_safe(d$elapsed_projection),
      elapsed_coord_lookup = median_safe(d$elapsed_coord_lookup),
      elapsed_edge_match = median_safe(d$elapsed_edge_match),
      elapsed_sparse_assembly = median_safe(d$elapsed_sparse_assembly),
      elapsed_prune_std = median_safe(d$elapsed_prune_std),
      peak_object_mb = median_safe(d$peak_object_mb),
      bytes_W_full_mb = median_safe(d$bytes_W_full_mb),
      bytes_P_mb = median_safe(d$bytes_P_mb),
      bytes_lookup_mb = median_safe(d$bytes_lookup_mb),
      bytes_coord_mb = median_safe(d$bytes_coord_mb),
      bytes_edge_index_mb = median_safe(d$bytes_edge_index_mb),
      bytes_W_common_mb = median_safe(d$bytes_W_common_mb),
      bytes_W_final_mb = median_safe(d$bytes_W_final_mb),
      failure_message = collapse_unique(d$error[d$status == "failed"]),
      warning_message = collapse_unique(d$warning_message),
      failure_class = collapse_unique(d$failure_class[d$status == "failed"]),
      failure_stage = collapse_unique(d$failure_stage[d$status == "failed"]),
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  summary_df[order(summary_df$pct_support, summary_df$method), , drop = FALSE]
}

make_compact_summary_v4 <- function(summary_df) {
  keep <- c(
    "scenario", "method", "pct_support",
    "runs", "ok_runs", "failed_runs",
    "cells_input", "cells_after_prune", "n_isolates",
    "nnz_common", "nnz_final", "elapsed_total", "peak_object_mb",
    "bytes_W_full_mb", "bytes_P_mb", "bytes_lookup_mb", "bytes_edge_index_mb",
    "bytes_W_common_mb", "bytes_W_final_mb", "failure_message"
  )
  keep <- keep[keep %in% names(summary_df)]
  summary_df[, keep, drop = FALSE]
}

capture_benchmark_env <- function() {
  list(
    sys_info = Sys.info(),
    r_version = R.version.string,
    session_info = utils::capture.output(sessionInfo())
  )
}

write_benchmark_outputs <- function(raw, summary_df, compact_summary,
                                    validation_df, validation_summary,
                                    metadata_df, out_dir,
                                    prefix = "benchmark_v4_proposed_vs_edge") {
  if (is.null(out_dir) || !nzchar(out_dir)) return(invisible(NULL))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  write.csv(raw, file.path(out_dir, paste0(prefix, "_raw.csv")), row.names = FALSE)
  write.csv(summary_df, file.path(out_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
  write.csv(compact_summary, file.path(out_dir, paste0(prefix, "_summary_compact.csv")), row.names = FALSE)
  write.csv(validation_df, file.path(out_dir, paste0(prefix, "_validation.csv")), row.names = FALSE)
  write.csv(validation_summary, file.path(out_dir, paste0(prefix, "_validation_summary.csv")), row.names = FALSE)
  write.csv(metadata_df, file.path(out_dir, paste0(prefix, "_metadata.csv")), row.names = FALSE)

  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(
      list(
        raw = raw,
        summary = summary_df,
        compact_summary = compact_summary,
        validation = validation_df,
        validation_summary = validation_summary,
        metadata = metadata_df
      ),
      path = file.path(out_dir, paste0(prefix, "_results.xlsx"))
    )
  }

  if (requireNamespace("xlsx", quietly = TRUE)) {
    xlsx::write.xlsx(raw, file.path(out_dir, paste0(prefix, "_raw.xlsx")), row.names = FALSE)
    xlsx::write.xlsx(summary_df, file.path(out_dir, paste0(prefix, "_summary.xlsx")), row.names = FALSE)
    xlsx::write.xlsx(compact_summary, file.path(out_dir, paste0(prefix, "_summary_compact.xlsx")), row.names = FALSE)
    xlsx::write.xlsx(validation_df, file.path(out_dir, paste0(prefix, "_validation.xlsx")), row.names = FALSE)
    xlsx::write.xlsx(validation_summary, file.path(out_dir, paste0(prefix, "_validation_summary.xlsx")), row.names = FALSE)
    xlsx::write.xlsx(metadata_df, file.path(out_dir, paste0(prefix, "_metadata.xlsx")), row.names = FALSE)
  }

  invisible(NULL)
}

# ============================================================
# 8. Main benchmark runner: Proposed vs edge-list comparator
# ============================================================

benchmark_proposed_vs_edge_v4 <- function(nr, nc,
                                          scenario_list,
                                          full_common_V = NULL,
                                          contiguity = c("queen", "rook"),
                                          reps = 1L,
                                          timeout_seconds = Inf,
                                          cache_W_full = FALSE,
                                          validate_equivalence = TRUE,
                                          validation_tol = 1e-12,
                                          compute_diff_nnz = FALSE,
                                          verbose = TRUE,
                                          out_dir = NULL,
                                          write_each_scenario = TRUE,
                                          log_prefix = "benchmark_v4_proposed_vs_edge") {
  contiguity <- match.arg(contiguity)
  if (is.null(names(scenario_list))) names(scenario_list) <- paste0("scenario_", seq_along(scenario_list))

  full_common_V <- full_common_V %||% max(vapply(scenario_list, length, integer(1)))

  cached_W_full <- NULL
  if (isTRUE(cache_W_full)) {
    message("Building cached W_FULL once ...")
    cached_W_full <- build_W_full(nr, nc, type = contiguity)
  }

  rows_out <- list()
  val_out <- list()
  k <- 1L
  v <- 1L

  for (sc_name in names(scenario_list)) {
    idx_support <- sort(unique(as.integer(scenario_list[[sc_name]])))
    pct_support <- length(idx_support) / full_common_V

    for (rep_id in seq_len(reps)) {
      if (verbose) {
        cat(sprintf(
          "Scenario %s | rep %d/%d | cells = %d | pct = %.4f%%\n",
          sc_name, rep_id, reps, length(idx_support), 100 * pct_support
        ))
      }

      return_mats <- isTRUE(validate_equivalence)

      prop_obj <- run_proposed_once_v4(
        nr = nr,
        nc = nc,
        idx_support = idx_support,
        contiguity = contiguity,
        timeout_seconds = timeout_seconds,
        cached_W_full = cached_W_full,
        return_matrices = return_mats
      )

      edge_obj <- run_edge_once_v4(
        nr = nr,
        nc = nc,
        idx_support = idx_support,
        contiguity = contiguity,
        timeout_seconds = timeout_seconds,
        return_matrices = return_mats
      )

      prop_res <- prop_obj$result
      edge_res <- edge_obj$result

      prop_res$scenario <- sc_name
      prop_res$rep <- rep_id
      prop_res$pct_support <- pct_support
      edge_res$scenario <- sc_name
      edge_res$rep <- rep_id
      edge_res$pct_support <- pct_support

      # MB columns
      for (nm in names(prop_res)) {
        if (startsWith(nm, "bytes_")) prop_res[[paste0(nm, "_mb")]] <- mb_num(prop_res[[nm]])
      }
      prop_res$peak_object_mb <- mb_num(prop_res$peak_object_bytes)

      for (nm in names(edge_res)) {
        if (startsWith(nm, "bytes_")) edge_res[[paste0(nm, "_mb")]] <- mb_num(edge_res[[nm]])
      }
      edge_res$peak_object_mb <- mb_num(edge_res$peak_object_bytes)

      rows_out[[k]] <- prop_res
      k <- k + 1L
      rows_out[[k]] <- edge_res
      k <- k + 1L

      if (isTRUE(validate_equivalence)) {
        val_out[[v]] <- validate_proposed_vs_edge(
          prop_obj = prop_obj,
          edge_obj = edge_obj,
          scenario = sc_name,
          rep = rep_id,
          pct_support = pct_support,
          tol = validation_tol,
          compute_diff_nnz = compute_diff_nnz
        )
        v <- v + 1L
      }

      rm(prop_obj, edge_obj)
      gc()
    }

    if (isTRUE(write_each_scenario) && !is.null(out_dir)) {
      raw_now <- bind_rows_fill(rows_out)
      num_cols <- setdiff(names(raw_now), c("method", "status", "error", "warning_message", "failure_class", "failure_stage", "scenario"))
      for (nm in num_cols) raw_now[[nm]] <- suppressWarnings(as.numeric(raw_now[[nm]]))
      summary_now <- summarize_benchmark_raw_v4(raw_now)
      compact_now <- make_compact_summary_v4(summary_now)

      val_now <- if (length(val_out)) bind_rows_fill(val_out) else data.frame()
      val_sum_now <- summarize_validation(val_now)
      meta_now <- benchmark_metadata_proposed_edge()

      write_benchmark_outputs(
        raw = raw_now,
        summary_df = summary_now,
        compact_summary = compact_now,
        validation_df = val_now,
        validation_summary = val_sum_now,
        metadata_df = meta_now,
        out_dir = out_dir,
        prefix = paste0(log_prefix, "_progress_", sc_name)
      )
    }
  }

  raw <- bind_rows_fill(rows_out)
  raw <- raw[order(raw$scenario, raw$method, raw$rep), , drop = FALSE]
  num_cols <- setdiff(names(raw), c("method", "status", "error", "warning_message", "failure_class", "failure_stage", "scenario"))
  for (nm in num_cols) raw[[nm]] <- suppressWarnings(as.numeric(raw[[nm]]))

  summary_df <- summarize_benchmark_raw_v4(raw)
  compact_summary <- make_compact_summary_v4(summary_df)
  validation_df <- if (length(val_out)) bind_rows_fill(val_out) else data.frame()
  validation_summary <- summarize_validation(validation_df)
  metadata_df <- benchmark_metadata_proposed_edge()

  if (!is.null(out_dir)) {
    write_benchmark_outputs(
      raw = raw,
      summary_df = summary_df,
      compact_summary = compact_summary,
      validation_df = validation_df,
      validation_summary = validation_summary,
      metadata_df = metadata_df,
      out_dir = out_dir,
      prefix = log_prefix
    )
    saveRDS(
      list(
        raw = raw,
        summary = summary_df,
        compact_summary = compact_summary,
        validation = validation_df,
        validation_summary = validation_summary,
        metadata = metadata_df,
        scenarios = scenario_list,
        env = capture_benchmark_env()
      ),
      file.path(out_dir, paste0(log_prefix, "_results.rds"))
    )
  }

  list(
    raw = raw,
    summary = summary_df,
    compact_summary = compact_summary,
    validation = validation_df,
    validation_summary = validation_summary,
    metadata = metadata_df,
    scenarios = scenario_list,
    env = capture_benchmark_env()
  )
}

# ============================================================
# 9. Standalone example
#    Edit paths as needed.
# ============================================================

run_standalone_example <- function() {
  cfg <- list(
    shapefile_admin2 = "D:/SHAPEFILE/gadm41_IDN_2.shp",
    data_root = "E:/Spatial Sparse W Matrix Research/DATA",
    landscan_root = "F:/LANDSCAN DATA",
    gdp_file = "E:/Spatial Sparse W Matrix Research/DATA/rast_gdp_tot_1990_2020_30arcsec.tif"
  )

  root_dir <- "E:/Spatial Sparse W Matrix Research/v3 Outputs/benchmark_v4_proposed_vs_edge"
  cache_dir <- file.path(root_dir, "cache")
  out_dir <- file.path(root_dir, "outputs")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  support_rds <- file.path(cache_dir, "support_obj_indonesia_2020.rds")
  region_inputs_rds <- file.path(cache_dir, "region_inputs_indonesia_2020.rds")

  prep <- prepare_benchmark_v4_indonesia_standalone(
    cfg = cfg,
    year = 2020,
    proportions = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00),
    contiguity = "queen",
    seed = 42,
    scenario_mode = "random",
    cache_to_global = TRUE,
    global_name = "region_inputs",
    support_rds = support_rds,
    reuse_precomputed = TRUE,
    save_precomputed = TRUE,
    force_rebuild = FALSE,
    save_region_inputs_rds = region_inputs_rds
  )

  print(prep$scenario_sizes)
  cat("Rows x Cols         :", prep$support$r, "x", prep$support$c, "\n")
  cat("Template cells      :", prep$support$N, "\n")
  cat("Full common support :", prep$support$V, "\n")

  bm <- benchmark_proposed_vs_edge_v4(
    nr = prep$support$r,
    nc = prep$support$c,
    scenario_list = prep$scenarios,
    full_common_V = prep$support$V,
    contiguity = "queen",
    reps = 2,
    timeout_seconds = 7200,
    cache_W_full = FALSE,
    validate_equivalence = TRUE,
    validation_tol = 1e-12,
    compute_diff_nnz = FALSE,
    verbose = TRUE,
    out_dir = out_dir,
    write_each_scenario = TRUE,
    log_prefix = "benchmark_v4_proposed_vs_edge"
  )

  print(bm$summary)
  print(bm$validation_summary)

  saveRDS(prep, file.path(root_dir, "benchmark_inputs_v4_proposed_vs_edge.rds"))
  saveRDS(bm, file.path(root_dir, "benchmark_results_v4_proposed_vs_edge.rds"))

  invisible(list(prep = prep, bm = bm))
}

# Uncomment to run:
# res <- run_standalone_example()



# ============================================================
# 10. NEW: Raster-native terra::adjacent full-lattice comparator
#     and three-method benchmark runner
# ------------------------------------------------------------
# Compared methods in this new runner:
#   (1) proposed_kronecker_projection
#   (2) edge_support_only_sparse
#   (3) terra_adjacent_full_lattice
#
# The terra comparator is intentionally designed as a middle-ground
# outside-spdep benchmark. It scans/constructs adjacency from the
# full raster lattice using terra::adjacent and only afterwards filters
# to the same retained support. This keeps it comparable with the
# full-lattice-to-support logic of the proposed Kronecker-projection route,
# while avoiding spdep neighbour-list objects.
# ============================================================

make_empty_sparse <- function(n) {
  Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    x = numeric(0),
    dims = c(n, n),
    giveCsparse = TRUE
  )
}

make_terra_template <- function(nr, nc) {
  terra::rast(
    nrows = as.integer(nr),
    ncols = as.integer(nc),
    xmin = 0,
    xmax = as.numeric(nc),
    ymin = 0,
    ymax = as.numeric(nr)
  )
}

terra_directions_from_contiguity <- function(contiguity = c("queen", "rook")) {
  contiguity <- match.arg(contiguity)
  if (contiguity == "queen") 8L else 4L
}

# Build W_COMMON by scanning the full raster lattice with terra::adjacent,
# then filtering the resulting full-lattice adjacency pairs to the retained
# common support. If chunk_size is not NULL, all full cells are still scanned,
# but terra::adjacent is called in chunks to reduce memory pressure.
run_terra_adjacent_full_once_v4 <- function(nr, nc, idx_support,
                                            contiguity = c("queen", "rook"),
                                            timeout_seconds = Inf,
                                            return_matrices = FALSE,
                                            terra_chunk_size = 500000L,
                                            critical_warning_patterns = c("integer overflow", "cannot allocate", "memory", "overflow")) {
  contiguity <- match.arg(contiguity)
  idx_support <- sort(unique(as.integer(idx_support)))
  N <- as.integer(nr) * as.integer(nc)
  V <- length(idx_support)

  if (V == 0L) stop("idx_support is empty.")
  if (any(idx_support < 1L | idx_support > N)) stop("idx_support contains out-of-range indices.")

  lookup_step <- timed_eval({
    template <- make_terra_template(nr, nc)
    lookup <- integer(N)
    lookup[idx_support] <- seq_len(V)
    list(template = template, lookup = lookup)
  }, timeout_seconds)

  if (!lookup_step$ok || critical_warning_detected(lookup_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "terra_adjacent_full_lattice",
      stage = "template_lookup",
      cells_input = V,
      error_message = if (!lookup_step$ok) lookup_step$error else paste(lookup_step$warnings, collapse = " | "),
      warning_message = collapse_unique(lookup_step$warnings),
      elapsed_total = lookup_step$elapsed
    )
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  template <- lookup_step$value$template
  lookup <- lookup_step$value$lookup
  directions <- terra_directions_from_contiguity(contiguity)

  adj_step <- timed_eval({
    if (is.null(terra_chunk_size) || length(terra_chunk_size) == 0L || is.na(terra_chunk_size) || !is.finite(terra_chunk_size)) {
      # Full materialisation of terra adjacent pairs in one call.
      adj <- terra::adjacent(
        x = template,
        cells = seq_len(N),
        directions = directions,
        pairs = TRUE,
        include = FALSE,
        symmetrical = FALSE
      )
      adj <- as.matrix(adj)
      if (ncol(adj) < 2L) stop("terra::adjacent did not return a two-column pair matrix.")
      from <- as.integer(adj[, 1])
      to   <- as.integer(adj[, 2])
      ii <- lookup[from]
      jj <- lookup[to]
      ok <- ii > 0L & jj > 0L
      list(
        edge_i = as.integer(ii[ok]),
        edge_j = as.integer(jj[ok]),
        n_full_pairs_scanned = nrow(adj),
        max_pair_chunk_bytes = bytes_num(adj),
        terra_chunk_size = NA_real_,
        terra_mode = "one_call_full_pairs"
      )
    } else {
      # Full-lattice scan in chunks. This still scans every full-template cell,
      # but it avoids storing the full pair matrix at once.
      chunk_size <- as.integer(terra_chunk_size)
      if (chunk_size <= 0L) stop("terra_chunk_size must be positive or NULL.")

      starts <- seq.int(1L, N, by = chunk_size)
      i_list <- vector("list", length(starts))
      j_list <- vector("list", length(starts))
      n_full_pairs_total <- 0
      max_pair_chunk_bytes <- 0

      for (g in seq_along(starts)) {
        s <- starts[g]
        e <- min(N, s + chunk_size - 1L)
        cells_g <- s:e

        adj_g <- terra::adjacent(
          x = template,
          cells = cells_g,
          directions = directions,
          pairs = TRUE,
          include = FALSE,
            symmetrical = FALSE
        )
        adj_g <- as.matrix(adj_g)
        if (length(adj_g) == 0L) {
          i_list[[g]] <- integer(0)
          j_list[[g]] <- integer(0)
          next
        }
        if (ncol(adj_g) < 2L) stop("terra::adjacent did not return a two-column pair matrix.")

        n_full_pairs_total <- n_full_pairs_total + nrow(adj_g)
        max_pair_chunk_bytes <- max(max_pair_chunk_bytes, bytes_num(adj_g), na.rm = TRUE)

        from <- as.integer(adj_g[, 1])
        to   <- as.integer(adj_g[, 2])
        ii <- lookup[from]
        jj <- lookup[to]
        ok <- ii > 0L & jj > 0L
        i_list[[g]] <- as.integer(ii[ok])
        j_list[[g]] <- as.integer(jj[ok])
      }

      list(
        edge_i = unlist(i_list, use.names = FALSE),
        edge_j = unlist(j_list, use.names = FALSE),
        n_full_pairs_scanned = n_full_pairs_total,
        max_pair_chunk_bytes = max_pair_chunk_bytes,
        terra_chunk_size = chunk_size,
        terra_mode = "chunked_full_lattice_scan"
      )
    }
  }, timeout_seconds)

  if (!adj_step$ok || critical_warning_detected(adj_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "terra_adjacent_full_lattice",
      stage = "terra_adjacent_filter",
      cells_input = V,
      error_message = if (!adj_step$ok) adj_step$error else paste(adj_step$warnings, collapse = " | "),
      warning_message = collapse_unique(c(lookup_step$warnings, adj_step$warnings)),
      elapsed_total = lookup_step$elapsed + adj_step$elapsed
    )
    df$elapsed_template_lookup <- lookup_step$elapsed
    df$elapsed_terra_adjacent_filter <- adj_step$elapsed
    df$bytes_lookup <- bytes_num(lookup)
    df$peak_object_bytes <- bytes_num(lookup)
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  edge_i <- adj_step$value$edge_i
  edge_j <- adj_step$value$edge_j
  n_full_pairs_scanned <- adj_step$value$n_full_pairs_scanned
  max_pair_chunk_bytes <- adj_step$value$max_pair_chunk_bytes
  terra_chunk_size_used <- adj_step$value$terra_chunk_size
  terra_mode <- adj_step$value$terra_mode

  assembly_step <- timed_eval({
    if (length(edge_i) == 0L) {
      W_terra <- make_empty_sparse(V)
    } else {
      W_terra <- Matrix::sparseMatrix(
        i = edge_i,
        j = edge_j,
        x = 1,
        dims = c(V, V),
        giveCsparse = TRUE
      )
    }
    as(Matrix::drop0(W_terra), "dgCMatrix")
  }, timeout_seconds)

  if (!assembly_step$ok || critical_warning_detected(assembly_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "terra_adjacent_full_lattice",
      stage = "sparse_assembly",
      cells_input = V,
      error_message = if (!assembly_step$ok) assembly_step$error else paste(assembly_step$warnings, collapse = " | "),
      warning_message = collapse_unique(c(lookup_step$warnings, adj_step$warnings, assembly_step$warnings)),
      elapsed_total = lookup_step$elapsed + adj_step$elapsed + assembly_step$elapsed
    )
    df$elapsed_template_lookup <- lookup_step$elapsed
    df$elapsed_terra_adjacent_filter <- adj_step$elapsed
    df$elapsed_sparse_assembly <- assembly_step$elapsed
    df$n_full_pairs_scanned <- n_full_pairs_scanned
    df$terra_chunk_size <- terra_chunk_size_used
    df$terra_mode <- terra_mode
    df$bytes_lookup <- bytes_num(lookup)
    df$bytes_terra_max_pair_chunk <- max_pair_chunk_bytes
    df$bytes_edge_index <- bytes_num(edge_i) + bytes_num(edge_j)
    df$peak_object_bytes <- max_na(c(bytes_num(lookup), max_pair_chunk_bytes, bytes_num(edge_i) + bytes_num(edge_j)))
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  W_terra_raw <- assembly_step$value

  prune_step <- timed_eval(row_standardize_drop0(W_terra_raw), timeout_seconds)

  if (!prune_step$ok || critical_warning_detected(prune_step$warnings, critical_warning_patterns)) {
    df <- make_fail_row(
      method = "terra_adjacent_full_lattice",
      stage = "prune_standardise",
      cells_input = V,
      error_message = if (!prune_step$ok) prune_step$error else paste(prune_step$warnings, collapse = " | "),
      warning_message = collapse_unique(c(lookup_step$warnings, adj_step$warnings, assembly_step$warnings, prune_step$warnings)),
      elapsed_total = lookup_step$elapsed + adj_step$elapsed + assembly_step$elapsed + prune_step$elapsed
    )
    df$elapsed_template_lookup <- lookup_step$elapsed
    df$elapsed_terra_adjacent_filter <- adj_step$elapsed
    df$elapsed_sparse_assembly <- assembly_step$elapsed
    df$elapsed_prune_std <- prune_step$elapsed
    df$n_full_pairs_scanned <- n_full_pairs_scanned
    df$terra_chunk_size <- terra_chunk_size_used
    df$terra_mode <- terra_mode
    df$bytes_lookup <- bytes_num(lookup)
    df$bytes_terra_max_pair_chunk <- max_pair_chunk_bytes
    df$bytes_edge_index <- bytes_num(edge_i) + bytes_num(edge_j)
    df$bytes_W_common <- bytes_num(W_terra_raw)
    df$peak_object_bytes <- max_na(c(bytes_num(lookup), max_pair_chunk_bytes, bytes_num(edge_i) + bytes_num(edge_j), bytes_num(W_terra_raw)))
    return(list(result = df, W_common = NULL, W_final = NULL, keep = NULL))
  }

  W_final <- prune_step$value$W
  keep <- prune_step$value$keep

  df <- data.frame(
    method = "terra_adjacent_full_lattice",
    status = "ok",
    failure_stage = NA_character_,
    failure_class = NA_character_,
    error = NA_character_,
    warning_message = collapse_unique(c(lookup_step$warnings, adj_step$warnings, assembly_step$warnings, prune_step$warnings)),
    cells_input = V,
    cells_after_prune = nrow(W_final),
    n_isolates = sum(!keep),
    nnz_common = length(W_terra_raw@x),
    nnz_final = length(W_final@x),
    elapsed_template_lookup = lookup_step$elapsed,
    elapsed_terra_adjacent_filter = adj_step$elapsed,
    elapsed_sparse_assembly = assembly_step$elapsed,
    elapsed_prune_std = prune_step$elapsed,
    elapsed_total = lookup_step$elapsed + adj_step$elapsed + assembly_step$elapsed + prune_step$elapsed,
    n_full_pairs_scanned = n_full_pairs_scanned,
    terra_chunk_size = terra_chunk_size_used,
    terra_mode = terra_mode,
    bytes_lookup = bytes_num(lookup),
    bytes_terra_max_pair_chunk = max_pair_chunk_bytes,
    bytes_edge_index = bytes_num(edge_i) + bytes_num(edge_j),
    bytes_W_common = bytes_num(W_terra_raw),
    bytes_W_final = bytes_num(W_final),
    peak_object_bytes = max_na(c(
      bytes_num(lookup),
      max_pair_chunk_bytes,
      bytes_num(edge_i) + bytes_num(edge_j),
      bytes_num(W_terra_raw),
      bytes_num(W_final)
    )),
    stringsAsFactors = FALSE
  )

  list(
    result = df,
    W_common = if (isTRUE(return_matrices)) W_terra_raw else NULL,
    W_final = if (isTRUE(return_matrices)) W_final else NULL,
    keep = if (isTRUE(return_matrices)) keep else NULL
  )
}

# Generic validation: compare any comparator output against the proposed output.
validate_against_proposed <- function(prop_obj, comp_obj,
                                      comparator_method,
                                      scenario = NA_character_,
                                      rep = NA_integer_,
                                      pct_support = NA_real_,
                                      tol = 1e-12,
                                      compute_diff_nnz = FALSE) {
  prop_ok <- is.data.frame(prop_obj$result) && identical(prop_obj$result$status[1], "ok")
  comp_ok <- is.data.frame(comp_obj$result) && identical(comp_obj$result$status[1], "ok")

  if (!prop_ok || !comp_ok) {
    return(data.frame(
      scenario = scenario,
      rep = rep,
      pct_support = pct_support,
      comparator_method = comparator_method,
      prop_ok = prop_ok,
      comparator_ok = comp_ok,
      common_equal = NA,
      final_equal = NA,
      keep_equal = NA,
      output_equivalence_ok = FALSE,
      common_max_abs_diff = NA_real_,
      final_max_abs_diff = NA_real_,
      common_diff_nnz = NA_real_,
      final_diff_nnz = NA_real_,
      cells_after_prune_proposed = if (prop_ok) prop_obj$result$cells_after_prune[1] else NA_real_,
      cells_after_prune_comparator = if (comp_ok) comp_obj$result$cells_after_prune[1] else NA_real_,
      nnz_common_proposed = if (prop_ok) prop_obj$result$nnz_common[1] else NA_real_,
      nnz_common_comparator = if (comp_ok) comp_obj$result$nnz_common[1] else NA_real_,
      nnz_final_proposed = if (prop_ok) prop_obj$result$nnz_final[1] else NA_real_,
      nnz_final_comparator = if (comp_ok) comp_obj$result$nnz_final[1] else NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  cmp_common <- sparse_compare_fast(
    prop_obj$W_common,
    comp_obj$W_common,
    tol = tol,
    compute_diff_nnz = compute_diff_nnz
  )
  cmp_final <- sparse_compare_fast(
    prop_obj$W_final,
    comp_obj$W_final,
    tol = tol,
    compute_diff_nnz = compute_diff_nnz
  )

  keep_equal <- identical(prop_obj$keep, comp_obj$keep)

  data.frame(
    scenario = scenario,
    rep = rep,
    pct_support = pct_support,
    comparator_method = comparator_method,
    prop_ok = TRUE,
    comparator_ok = TRUE,
    common_equal = cmp_common$equal_within_tol[1],
    final_equal = cmp_final$equal_within_tol[1],
    keep_equal = keep_equal,
    output_equivalence_ok = isTRUE(cmp_common$equal_within_tol[1]) &&
      isTRUE(cmp_final$equal_within_tol[1]) &&
      isTRUE(keep_equal),
    common_dims_equal = cmp_common$dims_equal[1],
    common_structure_equal = cmp_common$structure_equal[1],
    common_max_abs_diff = cmp_common$max_abs_diff[1],
    common_diff_nnz = cmp_common$diff_nnz[1],
    final_dims_equal = cmp_final$dims_equal[1],
    final_structure_equal = cmp_final$structure_equal[1],
    final_max_abs_diff = cmp_final$max_abs_diff[1],
    final_diff_nnz = cmp_final$diff_nnz[1],
    cells_after_prune_proposed = prop_obj$result$cells_after_prune[1],
    cells_after_prune_comparator = comp_obj$result$cells_after_prune[1],
    nnz_common_proposed = prop_obj$result$nnz_common[1],
    nnz_common_comparator = comp_obj$result$nnz_common[1],
    nnz_final_proposed = prop_obj$result$nnz_final[1],
    nnz_final_comparator = comp_obj$result$nnz_final[1],
    stringsAsFactors = FALSE
  )
}

summarize_validation_generic <- function(validation_df) {
  if (is.null(validation_df) || nrow(validation_df) == 0L) return(data.frame())
  split_key <- interaction(validation_df$scenario, validation_df$comparator_method, drop = TRUE)
  out <- lapply(split(validation_df, split_key, drop = TRUE), function(d) {
    data.frame(
      scenario = d$scenario[1],
      comparator_method = d$comparator_method[1],
      pct_support = d$pct_support[1],
      runs = nrow(d),
      prop_ok_runs = sum(d$prop_ok, na.rm = TRUE),
      comparator_ok_runs = sum(d$comparator_ok, na.rm = TRUE),
      equivalence_ok_runs = sum(d$output_equivalence_ok, na.rm = TRUE),
      all_equivalent = all(d$output_equivalence_ok, na.rm = TRUE),
      median_common_max_abs_diff = median_safe(d$common_max_abs_diff),
      median_final_max_abs_diff = median_safe(d$final_max_abs_diff),
      median_cells_after_prune_proposed = median_safe(d$cells_after_prune_proposed),
      median_cells_after_prune_comparator = median_safe(d$cells_after_prune_comparator),
      median_nnz_common_proposed = median_safe(d$nnz_common_proposed),
      median_nnz_common_comparator = median_safe(d$nnz_common_comparator),
      median_nnz_final_proposed = median_safe(d$nnz_final_proposed),
      median_nnz_final_comparator = median_safe(d$nnz_final_comparator),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$pct_support, out$comparator_method), , drop = FALSE]
}

# General summary that keeps all method-specific timing and object-size fields.
summarize_benchmark_raw_general <- function(raw) {
  if (is.null(raw) || nrow(raw) == 0L) return(data.frame())

  groups <- split(raw, list(raw$scenario, raw$method), drop = TRUE)
  out <- lapply(groups, function(d) {
    base <- data.frame(
      scenario = d$scenario[1],
      method = d$method[1],
      pct_support = d$pct_support[1],
      runs = nrow(d),
      ok_runs = sum(d$status == "ok", na.rm = TRUE),
      failed_runs = sum(d$status == "failed", na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    numeric_names <- names(d)[vapply(d, is.numeric, logical(1))]
    numeric_names <- setdiff(numeric_names, c("rep", "pct_support"))
    for (nm in numeric_names) {
      base[[nm]] <- median_safe(d[[nm]])
    }

    base$failure_message <- collapse_unique(d$error[d$status == "failed"])
    base$warning_message <- collapse_unique(d$warning_message)
    base$failure_class <- collapse_unique(d$failure_class[d$status == "failed"])
    base$failure_stage <- collapse_unique(d$failure_stage[d$status == "failed"])
    base
  })

  out <- bind_rows_fill(out)
  out[order(out$pct_support, out$method), , drop = FALSE]
}

make_compact_summary_general <- function(summary_df) {
  keep <- c(
    "scenario", "method", "pct_support",
    "runs", "ok_runs", "failed_runs",
    "cells_input", "cells_after_prune", "n_isolates",
    "nnz_common", "nnz_final",
    "elapsed_total", "peak_object_mb",
    "elapsed_build_full", "elapsed_projection",
    "elapsed_coord_lookup", "elapsed_edge_match", "elapsed_template_lookup",
    "elapsed_terra_adjacent_filter", "elapsed_sparse_assembly", "elapsed_prune_std",
    "n_full_pairs_scanned", "terra_chunk_size",
    "bytes_W_full_mb", "bytes_P_mb", "bytes_lookup_mb",
    "bytes_coord_mb", "bytes_edge_index_mb", "bytes_terra_max_pair_chunk_mb",
    "bytes_W_common_mb", "bytes_W_final_mb",
    "failure_message", "failure_stage", "failure_class"
  )
  keep <- keep[keep %in% names(summary_df)]
  summary_df[, keep, drop = FALSE]
}

make_public_comparison_general <- function(summary_df) {
  keep <- c(
    "scenario", "method", "pct_support", "runs", "ok_runs", "failed_runs",
    "cells_input", "cells_after_prune", "n_isolates", "nnz_common", "nnz_final",
    "elapsed_total", "peak_object_mb", "bytes_W_common_mb", "bytes_W_final_mb",
    "failure_message"
  )
  keep <- keep[keep %in% names(summary_df)]
  out <- summary_df[, keep, drop = FALSE]
  out
}

make_gain_vs_proposed_general <- function(summary_df) {
  prop <- summary_df[summary_df$method == "proposed_kronecker_projection", , drop = FALSE]
  alt <- summary_df[summary_df$method != "proposed_kronecker_projection", , drop = FALSE]
  if (nrow(prop) == 0L || nrow(alt) == 0L) return(data.frame())

  merged <- merge(
    alt,
    prop,
    by = "scenario",
    suffixes = c("_alt", "_prop"),
    all.x = TRUE
  )

  merged$time_ratio_alt_over_prop <- merged$elapsed_total_alt / merged$elapsed_total_prop
  merged$memory_ratio_alt_over_prop <- merged$peak_object_mb_alt / merged$peak_object_mb_prop
  merged$faster_method <- ifelse(
    is.finite(merged$time_ratio_alt_over_prop) & merged$time_ratio_alt_over_prop < 1,
    merged$method_alt,
    "proposed_kronecker_projection"
  )
  merged$lower_memory_method <- ifelse(
    is.finite(merged$memory_ratio_alt_over_prop) & merged$memory_ratio_alt_over_prop < 1,
    merged$method_alt,
    "proposed_kronecker_projection"
  )

  keep <- c(
    "scenario", "pct_support_alt", "method_alt",
    "cells_input_alt", "cells_after_prune_alt", "nnz_final_alt",
    "elapsed_total_prop", "elapsed_total_alt", "time_ratio_alt_over_prop", "faster_method",
    "peak_object_mb_prop", "peak_object_mb_alt", "memory_ratio_alt_over_prop", "lower_memory_method",
    "ok_runs_alt", "failed_runs_alt", "failure_message_alt"
  )
  keep <- keep[keep %in% names(merged)]
  out <- merged[, keep, drop = FALSE]
  names(out) <- sub("_alt$", "", names(out))
  names(out) <- sub("pct_support", "pct_support", names(out))
  out[order(out$pct_support_alt %||% out$pct_support, out$method), , drop = FALSE]
}

benchmark_metadata_proposed_edge_terra <- function() {
  data.frame(
    variable_name = c(
      "method", "proposed_kronecker_projection", "edge_support_only_sparse", "terra_adjacent_full_lattice",
      "elapsed_terra_adjacent_filter", "n_full_pairs_scanned", "terra_chunk_size", "output_equivalence_ok"
    ),
    description = c(
      "Benchmark method identifier.",
      "Proposed full-lattice sparse Kronecker construction followed by projection, isolate removal, and row-standardisation.",
      "Direct support-only sparse edge-list comparator using retained-cell coordinate offsets.",
      "Raster-native outside-spdep comparator that scans full-lattice adjacency using terra::adjacent, filters to retained support, constructs a sparse Matrix, removes isolates, and row-standardises.",
      "Elapsed time for the terra::adjacent full-lattice scan and support filtering stage.",
      "Number of full-lattice adjacency pairs scanned by terra::adjacent before support filtering.",
      "Chunk size used for terra::adjacent full-lattice scan; NA indicates one-call full pair materialisation.",
      "Validation flag indicating whether comparator output equals the proposed output under identical support and contiguity rule."
    ),
    stringsAsFactors = FALSE
  )
}

write_benchmark_outputs_general <- function(raw, summary_df, compact_summary,
                                            public_comparison, gain_vs_proposed,
                                            validation_df, validation_summary,
                                            metadata_df, out_dir,
                                            prefix = "benchmark_v4_proposed_edge_terra") {
  if (is.null(out_dir) || !nzchar(out_dir)) return(invisible(NULL))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  write.csv(raw, file.path(out_dir, paste0(prefix, "_raw.csv")), row.names = FALSE)
  write.csv(summary_df, file.path(out_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
  write.csv(compact_summary, file.path(out_dir, paste0(prefix, "_summary_compact.csv")), row.names = FALSE)
  write.csv(public_comparison, file.path(out_dir, paste0(prefix, "_public_comparison.csv")), row.names = FALSE)
  write.csv(gain_vs_proposed, file.path(out_dir, paste0(prefix, "_gain_vs_proposed.csv")), row.names = FALSE)
  write.csv(validation_df, file.path(out_dir, paste0(prefix, "_validation.csv")), row.names = FALSE)
  write.csv(validation_summary, file.path(out_dir, paste0(prefix, "_validation_summary.csv")), row.names = FALSE)
  write.csv(metadata_df, file.path(out_dir, paste0(prefix, "_metadata.csv")), row.names = FALSE)

  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(
      list(
        raw = raw,
        summary = summary_df,
        compact_summary = compact_summary,
        public_comparison = public_comparison,
        gain_vs_proposed = gain_vs_proposed,
        validation = validation_df,
        validation_summary = validation_summary,
        metadata = metadata_df
      ),
      path = file.path(out_dir, paste0(prefix, "_results.xlsx"))
    )
  }

  invisible(NULL)
}


# Create a stable output structure. Final outputs are written to
# <out_dir>/final_outputs and scenario-level progress CSV/XLSX files are
# written to <out_dir>/progress_csv. Set use_output_subfolders = FALSE in
# benchmark_proposed_edge_terra_v4() to recover the old flat output layout.
setup_benchmark_output_dirs <- function(out_dir,
                                        use_output_subfolders = TRUE,
                                        final_subdir = "final_outputs",
                                        progress_subdir = "progress_csv") {
  if (is.null(out_dir) || !nzchar(out_dir)) {
    return(list(
      base_out_dir = NULL,
      final_out_dir = NULL,
      progress_out_dir = NULL
    ))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(use_output_subfolders)) {
    final_out_dir <- file.path(out_dir, final_subdir)
    progress_out_dir <- file.path(out_dir, progress_subdir)
    dir.create(final_out_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(progress_out_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    final_out_dir <- out_dir
    progress_out_dir <- out_dir
  }

  list(
    base_out_dir = out_dir,
    final_out_dir = final_out_dir,
    progress_out_dir = progress_out_dir
  )
}

benchmark_proposed_edge_terra_v4 <- function(nr, nc,
                                             scenario_list,
                                             full_common_V = NULL,
                                             contiguity = c("queen", "rook"),
                                             reps = 1L,
                                             timeout_seconds = Inf,
                                             cache_W_full = FALSE,
                                             run_proposed = TRUE,
                                             run_edge = TRUE,
                                             run_terra = TRUE,
                                             terra_chunk_size = 500000L,
                                             validate_equivalence = TRUE,
                                             validation_tol = 1e-12,
                                             compute_diff_nnz = FALSE,
                                             verbose = TRUE,
                                             out_dir = NULL,
                                             write_each_scenario = TRUE,
                                             log_prefix = "benchmark_v4_proposed_edge_terra",
                                             use_output_subfolders = TRUE,
                                             final_subdir = "final_outputs",
                                             progress_subdir = "progress_csv") {
  contiguity <- match.arg(contiguity)
  if (is.null(names(scenario_list))) names(scenario_list) <- paste0("scenario_", seq_along(scenario_list))

  full_common_V <- full_common_V %||% max(vapply(scenario_list, length, integer(1)))

  out_dirs <- setup_benchmark_output_dirs(
    out_dir = out_dir,
    use_output_subfolders = use_output_subfolders,
    final_subdir = final_subdir,
    progress_subdir = progress_subdir
  )
  final_out_dir <- out_dirs$final_out_dir
  progress_out_dir <- out_dirs$progress_out_dir

  cached_W_full <- NULL
  if (isTRUE(cache_W_full)) {
    message("Building cached W_FULL once ...")
    cached_W_full <- build_W_full(nr, nc, type = contiguity)
  }

  rows_out <- list()
  val_out <- list()
  k <- 1L
  v <- 1L

  for (sc_name in names(scenario_list)) {
    idx_support <- sort(unique(as.integer(scenario_list[[sc_name]])))
    pct_support <- length(idx_support) / full_common_V

    for (rep_id in seq_len(reps)) {
      if (verbose) {
        cat(sprintf(
          "Scenario %s | rep %d/%d | cells = %d | pct = %.4f%%\n",
          sc_name, rep_id, reps, length(idx_support), 100 * pct_support
        ))
      }

      return_mats <- isTRUE(validate_equivalence)

      prop_obj <- NULL
      if (isTRUE(run_proposed)) {
        prop_obj <- run_proposed_once_v4(
          nr = nr,
          nc = nc,
          idx_support = idx_support,
          contiguity = contiguity,
          timeout_seconds = timeout_seconds,
          cached_W_full = cached_W_full,
          return_matrices = return_mats
        )
        prop_res <- prop_obj$result
        prop_res$scenario <- sc_name
        prop_res$rep <- rep_id
        prop_res$pct_support <- pct_support
        for (nm in names(prop_res)) {
          if (startsWith(nm, "bytes_")) prop_res[[paste0(nm, "_mb")]] <- mb_num(prop_res[[nm]])
        }
        prop_res$peak_object_mb <- mb_num(prop_res$peak_object_bytes)
        rows_out[[k]] <- prop_res
        k <- k + 1L
      }

      if (isTRUE(run_edge)) {
        edge_obj <- run_edge_once_v4(
          nr = nr,
          nc = nc,
          idx_support = idx_support,
          contiguity = contiguity,
          timeout_seconds = timeout_seconds,
          return_matrices = return_mats
        )
        edge_res <- edge_obj$result
        edge_res$scenario <- sc_name
        edge_res$rep <- rep_id
        edge_res$pct_support <- pct_support
        for (nm in names(edge_res)) {
          if (startsWith(nm, "bytes_")) edge_res[[paste0(nm, "_mb")]] <- mb_num(edge_res[[nm]])
        }
        edge_res$peak_object_mb <- mb_num(edge_res$peak_object_bytes)
        rows_out[[k]] <- edge_res
        k <- k + 1L

        if (isTRUE(validate_equivalence) && !is.null(prop_obj)) {
          val_out[[v]] <- validate_against_proposed(
            prop_obj = prop_obj,
            comp_obj = edge_obj,
            comparator_method = "edge_support_only_sparse",
            scenario = sc_name,
            rep = rep_id,
            pct_support = pct_support,
            tol = validation_tol,
            compute_diff_nnz = compute_diff_nnz
          )
          v <- v + 1L
        }
        rm(edge_obj)
      }

      if (isTRUE(run_terra)) {
        terra_obj <- run_terra_adjacent_full_once_v4(
          nr = nr,
          nc = nc,
          idx_support = idx_support,
          contiguity = contiguity,
          timeout_seconds = timeout_seconds,
          return_matrices = return_mats,
          terra_chunk_size = terra_chunk_size
        )
        terra_res <- terra_obj$result
        terra_res$scenario <- sc_name
        terra_res$rep <- rep_id
        terra_res$pct_support <- pct_support
        for (nm in names(terra_res)) {
          if (startsWith(nm, "bytes_")) terra_res[[paste0(nm, "_mb")]] <- mb_num(terra_res[[nm]])
        }
        terra_res$peak_object_mb <- mb_num(terra_res$peak_object_bytes)
        rows_out[[k]] <- terra_res
        k <- k + 1L

        if (isTRUE(validate_equivalence) && !is.null(prop_obj)) {
          val_out[[v]] <- validate_against_proposed(
            prop_obj = prop_obj,
            comp_obj = terra_obj,
            comparator_method = "terra_adjacent_full_lattice",
            scenario = sc_name,
            rep = rep_id,
            pct_support = pct_support,
            tol = validation_tol,
            compute_diff_nnz = compute_diff_nnz
          )
          v <- v + 1L
        }
        rm(terra_obj)
      }

      rm(prop_obj)
      gc()
    }

    if (isTRUE(write_each_scenario) && !is.null(progress_out_dir)) {
      raw_now <- bind_rows_fill(rows_out)
      raw_now <- raw_now[order(raw_now$scenario, raw_now$method, raw_now$rep), , drop = FALSE]
      num_cols <- setdiff(names(raw_now), c("method", "status", "error", "warning_message", "failure_class", "failure_stage", "scenario", "terra_mode"))
      for (nm in num_cols) raw_now[[nm]] <- suppressWarnings(as.numeric(raw_now[[nm]]))
      summary_now <- summarize_benchmark_raw_general(raw_now)
      compact_now <- make_compact_summary_general(summary_now)
      public_now <- make_public_comparison_general(summary_now)
      gain_now <- make_gain_vs_proposed_general(summary_now)
      val_now <- if (length(val_out)) bind_rows_fill(val_out) else data.frame()
      val_sum_now <- summarize_validation_generic(val_now)
      meta_now <- benchmark_metadata_proposed_edge_terra()

      write_benchmark_outputs_general(
        raw = raw_now,
        summary_df = summary_now,
        compact_summary = compact_now,
        public_comparison = public_now,
        gain_vs_proposed = gain_now,
        validation_df = val_now,
        validation_summary = val_sum_now,
        metadata_df = meta_now,
        out_dir = progress_out_dir,
        prefix = paste0(log_prefix, "_progress_", sc_name)
      )
    }
  }

  raw <- bind_rows_fill(rows_out)
  raw <- raw[order(raw$scenario, raw$method, raw$rep), , drop = FALSE]
  num_cols <- setdiff(names(raw), c("method", "status", "error", "warning_message", "failure_class", "failure_stage", "scenario", "terra_mode"))
  for (nm in num_cols) raw[[nm]] <- suppressWarnings(as.numeric(raw[[nm]]))

  summary_df <- summarize_benchmark_raw_general(raw)
  compact_summary <- make_compact_summary_general(summary_df)
  public_comparison <- make_public_comparison_general(summary_df)
  gain_vs_proposed <- make_gain_vs_proposed_general(summary_df)
  validation_df <- if (length(val_out)) bind_rows_fill(val_out) else data.frame()
  validation_summary <- summarize_validation_generic(validation_df)
  metadata_df <- benchmark_metadata_proposed_edge_terra()

  if (!is.null(final_out_dir)) {
    write_benchmark_outputs_general(
      raw = raw,
      summary_df = summary_df,
      compact_summary = compact_summary,
      public_comparison = public_comparison,
      gain_vs_proposed = gain_vs_proposed,
      validation_df = validation_df,
      validation_summary = validation_summary,
      metadata_df = metadata_df,
      out_dir = final_out_dir,
      prefix = log_prefix
    )
    saveRDS(
      list(
        raw = raw,
        summary = summary_df,
        compact_summary = compact_summary,
        public_comparison = public_comparison,
        gain_vs_proposed = gain_vs_proposed,
        validation = validation_df,
        validation_summary = validation_summary,
        metadata = metadata_df,
        scenarios = scenario_list,
        env = capture_benchmark_env()
      ),
      file.path(final_out_dir, paste0(log_prefix, "_results.rds"))
    )
  }

  list(
    raw = raw,
    summary = summary_df,
    compact_summary = compact_summary,
    public_comparison = public_comparison,
    gain_vs_proposed = gain_vs_proposed,
    validation = validation_df,
    validation_summary = validation_summary,
    metadata = metadata_df,
    scenarios = scenario_list,
    env = capture_benchmark_env()
  )
}

# ============================================================
# 11. Standalone example for the three-method benchmark
#     Edit paths as needed before running.
# ============================================================

run_standalone_example_proposed_edge_terra <- function() {
  cfg <- list(
    shapefile_admin2 = "D:/SHAPEFILE/gadm41_IDN_2.shp",
    data_root = "E:/Spatial Sparse W Matrix Research/DATA",
    landscan_root = "F:/LANDSCAN DATA",
    gdp_file = "E:/Spatial Sparse W Matrix Research/DATA/rast_gdp_tot_1990_2020_30arcsec.tif"
  )

  root_dir <- "E:/Spatial Sparse W Matrix Research/Latest v2 Outputs/benchmark_v4_proposed_edge_terra"
  cache_dir <- file.path(root_dir, "cache")
  out_dir <- file.path(root_dir, "outputs")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  # With use_output_subfolders = TRUE, outputs are written to:
  #   file.path(out_dir, "final_outputs")   for final CSV/XLSX/RDS files
  #   file.path(out_dir, "progress_csv")    for scenario-level progress files

  support_rds <- file.path(cache_dir, "support_obj_indonesia_2020.rds")
  region_inputs_rds <- file.path(cache_dir, "region_inputs_indonesia_2020.rds")

  prep <- prepare_benchmark_v4_indonesia_standalone(
    cfg = cfg,
    year = 2020,
    proportions = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00),
    contiguity = "queen",
    seed = 42,
    scenario_mode = "random",
    cache_to_global = TRUE,
    global_name = "region_inputs",
    support_rds = support_rds,
    reuse_precomputed = TRUE,
    save_precomputed = TRUE,
    force_rebuild = FALSE,
    save_region_inputs_rds = region_inputs_rds
  )

  print(prep$scenario_sizes)
  cat("Rows x Cols         :", prep$support$r, "x", prep$support$c, "\n")
  cat("Template cells      :", prep$support$N, "\n")
  cat("Full common support :", prep$support$V, "\n")

  bm <- benchmark_proposed_edge_terra_v4(
    nr = prep$support$r,
    nc = prep$support$c,
    scenario_list = prep$scenarios,
    full_common_V = prep$support$V,
    contiguity = "queen",
    reps = 2,
    timeout_seconds = 7200,
    cache_W_full = FALSE,
    run_proposed = TRUE,
    run_edge = TRUE,
    run_terra = TRUE,
    terra_chunk_size = 500000L,
    validate_equivalence = TRUE,
    validation_tol = 1e-12,
    compute_diff_nnz = FALSE,
    verbose = TRUE,
    out_dir = out_dir,
    write_each_scenario = TRUE,
    log_prefix = "benchmark_v4_proposed_edge_terra",
    use_output_subfolders = TRUE,
    final_subdir = "final_outputs",
    progress_subdir = "progress_csv"
  )

  print(bm$summary)
  print(bm$validation_summary)
  print(bm$gain_vs_proposed)

  saveRDS(prep, file.path(root_dir, "benchmark_inputs_v4_proposed_edge_terra.rds"))
  saveRDS(bm, file.path(root_dir, "benchmark_results_v4_proposed_edge_terra.rds"))

  invisible(list(prep = prep, bm = bm))
}

# Uncomment to run the new three-method benchmark:
res <- run_standalone_example_proposed_edge_terra()
