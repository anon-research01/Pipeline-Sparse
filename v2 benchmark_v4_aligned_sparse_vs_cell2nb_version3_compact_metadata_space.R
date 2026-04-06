suppressPackageStartupMessages({
  library(Matrix)
  library(terra)
  library(spdep)
  library(sf)
  library(xlsx)
})

# ============================================================
# Benchmark script aligned to the v4 full process
# ------------------------------------------------------------
# Goal:
#   Compare the proposed Table-1 sparse pipeline against the
#   conventional cell2nb + nb2mat baseline on the SAME common
#   support used in the v4 four-layer Indonesia process.
#
# Core alignment with v4:
#   1) one common support across GDP, population, NDVI, MNDWI
#   2) W_full built by sparse Kronecker on the full template
#   3) projection matrix P has dimension N x V
#   4) W_common = t(P) %*% W_full %*% P
#   5) isolates dropped BEFORE row-standardisation
#
# Patch highlights in this version:
#   A) preparation is split into a one-time support precompute step
#   B) precomputed support can be saved/loaded from .rds
#   C) scenario generation is separated from support building
#   D) spatial nested scenarios no longer call full-grid cell2nb
#      during preparation; adjacency is built only on the support
# ============================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

# =========================
# Utilities copied/aligned from v4
# =========================

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

  W <- kronecker(Ir, Bc) + kronecker(Br, Ic)
  if (type == "queen") W <- W + kronecker(Br, Bc)
  as(W, "dgCMatrix")
}

row_standardize_drop0 <- function(W) {
  rs <- Matrix::rowSums(W)
  keep <- rs > 0
  if (!all(keep)) W <- W[keep, keep, drop = FALSE]
  rs <- Matrix::rowSums(W)
  W <- Matrix::Diagonal(x = 1 / rs) %*% W
  list(W = as(Matrix::drop0(W), "dgCMatrix"), keep = keep)
}

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

  as.numeric(t(M))
}

# =========================
# Standalone region-input preparation aligned to v4
# =========================

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

prepare_region_inputs_standalone <- function(cfg, year = 2020, cache_to_global = TRUE, global_name = "region_inputs") {
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

# =========================
# Small helpers
# =========================

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
  out <- gsub("\\.", "p", out)
  out
}


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

cell_to_rc <- function(idx, nc) {
  idx <- as.integer(idx)
  data.frame(
    row = ((idx - 1L) %/% nc) + 1L,
    col = ((idx - 1L) %% nc) + 1L
  )
}

# =========================
# Build common support exactly as in v4
# =========================

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

  needed <- c("mask_indonesia_gdp_2020", "mask_Landscan_indonesia_2020", "mask_NDVI_indonesia_2020", "mask_MNDWI_indonesia_2020")
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

# =========================
# Support precompute / caching layer
# =========================

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

prepare_scenarios_from_support_v4 <- function(support_obj,
                                              proportions = c(0.05, 0.10, 0.25, 0.35, 0.50, 0.75, 1.00),
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
                                                      proportions = c(0.05, 0.10, 0.25, 0.35, 0.50, 0.75, 1.00),
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

# =========================
# Nested scenario generation on the v4 common support
# Fast spatial mode without full-grid cell2nb
# =========================

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

build_support_adjacency <- function(idx_all, nr, nc, contiguity = c("queen", "rook")) {
  contiguity <- match.arg(contiguity)
  idx_all <- sort(unique(as.integer(idx_all)))
  rc <- cell_to_rc(idx_all, nc)
  keys <- as.character(idx_all)
  pos_lookup <- seq_along(idx_all)
  names(pos_lookup) <- keys

  offsets <- if (contiguity == "queen") {
    rbind(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L),
          c(-1L, -1L), c(-1L, 1L), c(1L, -1L), c(1L, 1L))
  } else {
    rbind(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L))
  }

  adj <- vector("list", length(idx_all))
  for (i in seq_along(idx_all)) {
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
    pos <- unname(pos_lookup[as.character(cand_idx)])
    adj[[i]] <- as.integer(pos[!is.na(pos)])
  }
  adj
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

make_spatial_nested_supports <- function(idx_all, nr, nc,
                                         proportions,
                                         contiguity = c("queen", "rook"),
                                         seed = 42) {
  contiguity <- match.arg(contiguity)
  proportions <- normalize_proportions(proportions)
  idx_all <- sort(unique(as.integer(idx_all)))
  adj <- build_support_adjacency(idx_all, nr, nc, contiguity = contiguity)
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

    # remap local positions for the component
    local_index <- integer(length(idx_all))
    local_index[pos] <- seq_along(pos)
    adj_local <- lapply(pos, function(p) local_index[adj[[p]][adj[[p]] %in% pos]])
    ord_local <- bfs_order_adj(adj_local, start = local_index[seed_local])
    orders[[j]] <- pos[ord_local]
  }

  out <- vector("list", length(proportions))
  names(out) <- scenario_label_pct(proportions)
  N <- length(idx_all)

  for (i in seq_along(proportions)) {
    k <- max(1L, floor(N * proportions[i]))
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

# =========================
# Proposed method: EXACT Table-1/v4 matrix stage
# =========================

run_proposed_once_v4 <- function(nr, nc, idx_support,
                                 contiguity = c("queen", "rook"),
                                 timeout_seconds = Inf,
                                 cached_W_full = NULL,
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
    return(data.frame(
      method = "proposed_table1_v4",
      status = "failed",
      error = if (!build_step$ok) build_step$error else paste(build_step$warnings, collapse = " | "),
      warning_message = collapse_unique(build_step$warnings),
      cells_input = V,
      cells_after_prune = NA_real_,
      n_isolates = NA_real_,
      nnz_common = NA_real_,
      nnz_final = NA_real_,
      elapsed_build_full = build_step$elapsed,
      elapsed_projection = NA_real_,
      elapsed_prune_std = NA_real_,
      elapsed_total = build_step$elapsed,
      bytes_W_full = NA_real_,
      bytes_P = NA_real_,
      bytes_W_common = NA_real_,
      bytes_W_final = NA_real_,
      peak_object_bytes = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  W_full <- as(build_step$value, "dgCMatrix")

  proj_step <- timed_eval({
    P <- Matrix::sparseMatrix(i = idx_support, j = seq_len(V), x = 1, dims = c(N, V))
    W_common <- as(Matrix::drop0(Matrix::t(P) %*% W_full %*% P), "dgCMatrix")
    list(P = P, W_common = W_common)
  }, timeout_seconds)

  if (!proj_step$ok || critical_warning_detected(proj_step$warnings, critical_warning_patterns)) {
    return(data.frame(
      method = "proposed_table1_v4",
      status = "failed",
      error = if (!proj_step$ok) proj_step$error else paste(proj_step$warnings, collapse = " | "),
      warning_message = collapse_unique(c(build_step$warnings, proj_step$warnings)),
      cells_input = V,
      cells_after_prune = NA_real_,
      n_isolates = NA_real_,
      nnz_common = NA_real_,
      nnz_final = NA_real_,
      elapsed_build_full = build_step$elapsed,
      elapsed_projection = proj_step$elapsed,
      elapsed_prune_std = NA_real_,
      elapsed_total = build_step$elapsed + proj_step$elapsed,
      bytes_W_full = bytes_num(W_full),
      bytes_P = NA_real_,
      bytes_W_common = NA_real_,
      bytes_W_final = NA_real_,
      peak_object_bytes = bytes_num(W_full),
      stringsAsFactors = FALSE
    ))
  }

  P <- proj_step$value$P
  W_common_raw <- proj_step$value$W_common
  prune_step <- timed_eval(row_standardize_drop0(W_common_raw), timeout_seconds)

  if (!prune_step$ok || critical_warning_detected(prune_step$warnings, critical_warning_patterns)) {
    return(data.frame(
      method = "proposed_table1_v4",
      status = "failed",
      error = if (!prune_step$ok) prune_step$error else paste(prune_step$warnings, collapse = " | "),
      warning_message = collapse_unique(c(build_step$warnings, proj_step$warnings, prune_step$warnings)),
      cells_input = V,
      cells_after_prune = NA_real_,
      n_isolates = NA_real_,
      nnz_common = length(W_common_raw@x),
      nnz_final = NA_real_,
      elapsed_build_full = build_step$elapsed,
      elapsed_projection = proj_step$elapsed,
      elapsed_prune_std = prune_step$elapsed,
      elapsed_total = build_step$elapsed + proj_step$elapsed + prune_step$elapsed,
      bytes_W_full = bytes_num(W_full),
      bytes_P = bytes_num(P),
      bytes_W_common = bytes_num(W_common_raw),
      bytes_W_final = NA_real_,
      peak_object_bytes = max_na(c(bytes_num(W_full), bytes_num(P), bytes_num(W_common_raw))),
      stringsAsFactors = FALSE
    ))
  }

  W_final <- prune_step$value$W
  keep <- prune_step$value$keep

  data.frame(
    method = "proposed_table1_v4",
    status = "ok",
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
}

# =========================
# Baseline method: cell2nb + subset.nb + nb2mat
# =========================

prune_nb <- function(nb_obj) {
  deg <- spdep::card(nb_obj)
  keep <- deg > 0L
  nb_pruned <- spdep::subset.nb(nb_obj, subset = keep)
  list(nb_pruned = nb_pruned, keep = keep, n_isolates = sum(!keep))
}

rowstd_dense <- function(W) {
  if (!length(W)) return(W)
  rs <- rowSums(W)
  if (any(rs <= 0)) stop("Zero row sums encountered during dense row-standardisation.")
  W / rs
}

run_baseline_once_v4 <- function(nr, nc, idx_support,
                                 contiguity = c("queen", "rook"),
                                 nb2mat_style = c("B", "W"),
                                 timeout_seconds = Inf,
                                 cached_nb_full = NULL,
                                 critical_warning_patterns = c("integer overflow", "cannot allocate", "memory", "overflow")) {
  contiguity <- match.arg(contiguity)
  nb2mat_style <- match.arg(nb2mat_style)
  idx_support <- sort(unique(as.integer(idx_support)))
  N <- as.integer(nr) * as.integer(nc)
  V <- length(idx_support)

  if (V == 0L) stop("idx_support is empty.")
  if (any(idx_support < 1L | idx_support > N)) stop("idx_support contains out-of-range indices.")

  keep_full <- rep(FALSE, N)
  keep_full[idx_support] <- TRUE

  fail_row <- function(stage, error_message = NA_character_, warning_message = NA_character_,
                       elapsed_cell2nb = NA_real_, elapsed_subset_nb = NA_real_,
                       elapsed_prune_nb = NA_real_, elapsed_nb2mat = NA_real_,
                       bytes_nb_full = NA_real_, bytes_nb_common = NA_real_,
                       bytes_nb_pruned = NA_real_, bytes_W_final = NA_real_,
                       peak_object_bytes = NA_real_, n_isolates = NA_real_,
                       nnz_common = NA_real_) {
    data.frame(
      method = "baseline_cell2nb_nb2mat_v4",
      status = "failed",
      error = error_message,
      warning_message = warning_message,
      failure_class = classify_failure_v4(
        error_message = error_message,
        warning_message = warning_message
      ),
      failure_stage = stage,
      cells_input = V,
      cells_after_prune = NA_real_,
      n_isolates = n_isolates,
      nnz_common = nnz_common,
      nnz_final = NA_real_,
      elapsed_cell2nb = elapsed_cell2nb,
      elapsed_subset_nb = elapsed_subset_nb,
      elapsed_prune_nb = elapsed_prune_nb,
      elapsed_nb2mat = elapsed_nb2mat,
      elapsed_total = sum(c(elapsed_cell2nb, elapsed_subset_nb, elapsed_prune_nb, elapsed_nb2mat), na.rm = TRUE),
      bytes_nb_full = bytes_nb_full,
      bytes_nb_common = bytes_nb_common,
      bytes_nb_pruned = bytes_nb_pruned,
      bytes_W_final = bytes_W_final,
      peak_object_bytes = peak_object_bytes,
      stringsAsFactors = FALSE
    )
  }

  step_cell2nb <- if (is.null(cached_nb_full)) {
    timed_eval(spdep::cell2nb(nrow = nr, ncol = nc, type = contiguity, torus = FALSE), timeout_seconds)
  } else {
    list(ok = TRUE, value = cached_nb_full, elapsed = 0, error = NA_character_, warnings = character(0))
  }

  warn_cell2nb <- collapse_unique(step_cell2nb$warnings)
  if (!step_cell2nb$ok || critical_warning_detected(step_cell2nb$warnings, critical_warning_patterns)) {
    return(fail_row(
      stage = "cell2nb",
      error_message = if (!step_cell2nb$ok) step_cell2nb$error else warn_cell2nb,
      warning_message = warn_cell2nb,
      elapsed_cell2nb = step_cell2nb$elapsed
    ))
  }

  nb_full <- step_cell2nb$value
  step_subset <- timed_eval(spdep::subset.nb(nb_full, subset = keep_full), timeout_seconds)
  warn_subset <- collapse_unique(c(step_cell2nb$warnings, step_subset$warnings))

  if (!step_subset$ok || critical_warning_detected(step_subset$warnings, critical_warning_patterns)) {
    return(fail_row(
      stage = "subset_nb",
      error_message = if (!step_subset$ok) step_subset$error else collapse_unique(step_subset$warnings),
      warning_message = warn_subset,
      elapsed_cell2nb = step_cell2nb$elapsed,
      elapsed_subset_nb = step_subset$elapsed,
      bytes_nb_full = bytes_num(nb_full),
      peak_object_bytes = bytes_num(nb_full)
    ))
  }

  nb_common <- step_subset$value
  nnz_common_approx <- sum(spdep::card(nb_common))
  step_prune <- timed_eval(prune_nb(nb_common), timeout_seconds)
  warn_prune <- collapse_unique(c(step_cell2nb$warnings, step_subset$warnings, step_prune$warnings))

  if (!step_prune$ok || critical_warning_detected(step_prune$warnings, critical_warning_patterns)) {
    return(fail_row(
      stage = "prune_nb",
      error_message = if (!step_prune$ok) step_prune$error else collapse_unique(step_prune$warnings),
      warning_message = warn_prune,
      elapsed_cell2nb = step_cell2nb$elapsed,
      elapsed_subset_nb = step_subset$elapsed,
      elapsed_prune_nb = step_prune$elapsed,
      bytes_nb_full = bytes_num(nb_full),
      bytes_nb_common = bytes_num(nb_common),
      peak_object_bytes = max_na(c(bytes_num(nb_full), bytes_num(nb_common))),
      nnz_common = nnz_common_approx
    ))
  }

  prune_obj <- step_prune$value
  nb_pruned <- prune_obj$nb_pruned
  n_isolates <- prune_obj$n_isolates

  step_nb2mat <- timed_eval({
    if (length(nb_pruned) == 0L) {
      matrix(0, 0, 0)
    } else if (nb2mat_style == "W") {
      spdep::nb2mat(nb_pruned, style = "W", zero.policy = TRUE)
    } else {
      W_bin <- spdep::nb2mat(nb_pruned, style = "B", zero.policy = TRUE)
      rowstd_dense(W_bin)
    }
  }, timeout_seconds)
  warn_nb2mat <- collapse_unique(c(step_cell2nb$warnings, step_subset$warnings, step_prune$warnings, step_nb2mat$warnings))

  if (!step_nb2mat$ok || critical_warning_detected(step_nb2mat$warnings, critical_warning_patterns)) {
    return(fail_row(
      stage = "nb2mat",
      error_message = if (!step_nb2mat$ok) step_nb2mat$error else collapse_unique(step_nb2mat$warnings),
      warning_message = warn_nb2mat,
      elapsed_cell2nb = step_cell2nb$elapsed,
      elapsed_subset_nb = step_subset$elapsed,
      elapsed_prune_nb = step_prune$elapsed,
      elapsed_nb2mat = step_nb2mat$elapsed,
      bytes_nb_full = bytes_num(nb_full),
      bytes_nb_common = bytes_num(nb_common),
      bytes_nb_pruned = bytes_num(nb_pruned),
      peak_object_bytes = max_na(c(bytes_num(nb_full), bytes_num(nb_common), bytes_num(nb_pruned))),
      n_isolates = n_isolates,
      nnz_common = nnz_common_approx
    ))
  }

  W_final <- step_nb2mat$value
  nnz_final <- if (!length(W_final)) 0 else sum(W_final != 0)

  data.frame(
    method = "baseline_cell2nb_nb2mat_v4",
    status = "ok",
    error = NA_character_,
    warning_message = warn_nb2mat,
    failure_class = NA_character_,
    failure_stage = NA_character_,
    cells_input = V,
    cells_after_prune = nrow(W_final),
    n_isolates = n_isolates,
    nnz_common = nnz_common_approx,
    nnz_final = nnz_final,
    elapsed_cell2nb = step_cell2nb$elapsed,
    elapsed_subset_nb = step_subset$elapsed,
    elapsed_prune_nb = step_prune$elapsed,
    elapsed_nb2mat = step_nb2mat$elapsed,
    elapsed_total = step_cell2nb$elapsed + step_subset$elapsed + step_prune$elapsed + step_nb2mat$elapsed,
    bytes_nb_full = bytes_num(nb_full),
    bytes_nb_common = bytes_num(nb_common),
    bytes_nb_pruned = bytes_num(nb_pruned),
    bytes_W_final = bytes_num(W_final),
    peak_object_bytes = max_na(c(bytes_num(nb_full), bytes_num(nb_common), bytes_num(nb_pruned), bytes_num(W_final))),
    stringsAsFactors = FALSE
  )
}

# =========================
# Benchmark logging helpers
# =========================

benchmark_metadata_v4 <- function() {
  data.frame(
    variable_name = c(
      "scenario", "method", "pct_support", "runs", "ok_runs", "failed_runs", "skipped_runs",
      "hard_failed_runs", "soft_failed_runs", "cells_input", "cells_after_prune", "n_isolates",
      "nnz_final", "elapsed_total", "peak_object_mb",
      "bytes_W_full_mb", "bytes_W_common_mb", "bytes_W_final_mb",
      "bytes_nb_full_mb", "bytes_nb_common_mb", "bytes_nb_pruned_mb",
      "failure_message", "failure_class", "failure_stage"
    ),
    description = c(
      "Scenario label for the retained-support benchmark level.",
      "Benchmark method identifier.",
      "Retained support as a share of the full common support.",
      "Number of runs combined in the summary row.",
      "Number of successful runs.",
      "Number of failed runs.",
      "Number of baseline runs skipped after a prior hard failure.",
      "Number of failed runs classified as hard failures.",
      "Number of failed runs classified as soft failures.",
      "Number of cells entering the matrix-construction stage before isolate removal.",
      "Number of cells remaining after isolate removal.",
      "Number of isolates removed after projection or support restriction.",
      "Number of non-zero entries in the final matrix.",
      "Total elapsed time for the matrix-construction workflow.",
      "Largest object size observed among tracked objects, in MB.",
      "Object size of W_FULL in MB; proposed method only.",
      "Object size of W_COMMON in MB; proposed method only.",
      "Object size of final row-standardised W in MB; proposed for final sparse W, baseline for dense/matrix output only when applicable.",
      "Object size of the full neighbour-list object produced by cell2nb, in MB; baseline only.",
      "Object size of the support-restricted neighbour-list object, in MB; baseline only.",
      "Object size of the pruned neighbour-list object after isolate removal, in MB; baseline only.",
      "Error message recorded for failed runs.",
      "Failure classification used by the skip rule: hard or soft.",
      "Stage at which a failure occurred."
    ),
    method_applicability = c(
      "both", "both", "both", "both", "both", "both", "baseline only",
      "baseline only", "baseline only", "both", "both", "both",
      "both", "both", "both",
      "proposed only", "proposed only", "method-specific",
      "baseline only", "baseline only", "baseline only",
      "failed runs only", "failed baseline runs only", "failed baseline runs only"
    ),
    notes = c(
      "Generated from the retained-support proportions.",
      "proposed_table1_v4 or baseline_cell2nb_nb2mat_v4.",
      "Computed relative to the full common support used in the benchmark.",
      "Typically equals the number of repetitions for the scenario-method pair.",
      "Summary count.",
      "Summary count.",
      "Applies when skip_larger_baseline_after_fail is enabled and triggered by a hard failure.",
      "Hard failures are typically size-related, for example memory exhaustion or timeout.",
      "Soft failures are non-size-related or ambiguous failures.",
      "Median across runs in the summary table.",
      "Median across runs in the summary table.",
      "Median across runs in the summary table.",
      "Median across runs in the summary table.",
      "Median across runs in the summary table.",
      "Median across runs in the summary table.",
      "NA for the baseline rows.",
      "NA for the baseline rows.",
      "For proposed rows this is the final sparse matrix; for baseline rows this is the final matrix object produced after nb2mat and standardisation.",
      "NA for the proposed rows.",
      "NA for the proposed rows.",
      "NA for the proposed rows.",
      "Collapsed across failed runs in the summary table.",
      "Only reported for failed baseline runs.",
      "Only reported for failed baseline runs."
    ),
    stringsAsFactors = FALSE
  )
}

summarize_benchmark_raw_v4 <- function(raw) {
  summary_list <- lapply(split(raw, list(raw$scenario, raw$method), drop = TRUE), function(d) {
    data.frame(
      scenario = d$scenario[1], method = d$method[1], pct_support = d$pct_support[1],
      runs = nrow(d), ok_runs = sum(d$status == "ok", na.rm = TRUE), failed_runs = sum(d$status == "failed", na.rm = TRUE),
      skipped_runs = sum(d$status == "skipped_after_hard_failure", na.rm = TRUE),
      hard_failed_runs = sum(d$status == "failed" & d$failure_class == "hard", na.rm = TRUE),
      soft_failed_runs = sum(d$status == "failed" & d$failure_class == "soft", na.rm = TRUE),
      cells_input = median_safe(d$cells_input), cells_after_prune = median_safe(d$cells_after_prune),
      n_isolates = median_safe(d$n_isolates), nnz_final = median_safe(d$nnz_final),
      elapsed_total = median_safe(d$elapsed_total),
      elapsed_build_full = median_safe(d$elapsed_build_full), elapsed_projection = median_safe(d$elapsed_projection),
      elapsed_prune_std = median_safe(d$elapsed_prune_std), elapsed_cell2nb = median_safe(d$elapsed_cell2nb),
      elapsed_subset_nb = median_safe(d$elapsed_subset_nb), elapsed_prune_nb = median_safe(d$elapsed_prune_nb),
      elapsed_nb2mat = median_safe(d$elapsed_nb2mat),
      peak_object_mb = median_safe(d$peak_object_mb),
      bytes_W_full_mb = median_safe(d$bytes_W_full_mb),
      bytes_P_mb = median_safe(d$bytes_P_mb),
      bytes_W_common_mb = median_safe(d$bytes_W_common_mb),
      bytes_W_final_mb = median_safe(d$bytes_W_final_mb),
      bytes_nb_full_mb = median_safe(d$bytes_nb_full_mb),
      bytes_nb_common_mb = median_safe(d$bytes_nb_common_mb),
      bytes_nb_pruned_mb = median_safe(d$bytes_nb_pruned_mb),
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
    "runs", "ok_runs", "failed_runs", "skipped_runs",
    "cells_input", "cells_after_prune", "n_isolates", "nnz_final",
    "elapsed_total", "peak_object_mb",
    "bytes_W_full_mb", "bytes_W_common_mb", "bytes_W_final_mb",
    "bytes_nb_full_mb", "bytes_nb_common_mb", "bytes_nb_pruned_mb",
    "failure_message"
  )
  keep <- keep[keep %in% names(summary_df)]
  summary_df[, keep, drop = FALSE]
}

write_benchmark_progress_v4 <- function(raw, summary_df, out_dir,
                                        prefix = "benchmark_v4_aligned",
                                        scenario_name = NULL,
                                        write_excel = TRUE,
                                        write_csv = TRUE) {
  if (is.null(out_dir) || !nzchar(out_dir)) return(invisible(NULL))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  stamp <- if (is.null(scenario_name)) "latest" else gsub("[^A-Za-z0-9_-]", "_", scenario_name)
  compact_summary <- make_compact_summary_v4(summary_df)
  metadata_df <- benchmark_metadata_v4()

  if (isTRUE(write_csv)) {
    write.csv(raw, file.path(out_dir, paste0(prefix, "_raw_", stamp, ".csv")), row.names = FALSE)
    write.csv(summary_df, file.path(out_dir, paste0(prefix, "_summary_", stamp, ".csv")), row.names = FALSE)
    write.csv(compact_summary, file.path(out_dir, paste0(prefix, "_summary_compact_", stamp, ".csv")), row.names = FALSE)
    write.csv(metadata_df, file.path(out_dir, paste0(prefix, "_metadata_", stamp, ".csv")), row.names = FALSE)

    write.csv(raw, file.path(out_dir, paste0(prefix, "_raw_latest.csv")), row.names = FALSE)
    write.csv(summary_df, file.path(out_dir, paste0(prefix, "_summary_latest.csv")), row.names = FALSE)
    write.csv(compact_summary, file.path(out_dir, paste0(prefix, "_summary_compact_latest.csv")), row.names = FALSE)
    write.csv(metadata_df, file.path(out_dir, paste0(prefix, "_metadata_latest.csv")), row.names = FALSE)
  }

  if (isTRUE(write_excel) && requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(
      list(raw = raw, summary = summary_df, compact_summary = compact_summary, metadata = metadata_df),
      path = file.path(out_dir, paste0(prefix, "_progress_", stamp, ".xlsx"))
    )
    writexl::write_xlsx(
      list(raw = raw, summary = summary_df, compact_summary = compact_summary, metadata = metadata_df),
      path = file.path(out_dir, paste0(prefix, "_progress_latest.xlsx"))
    )
  }

  invisible(NULL)
}

# =========================
# Benchmark runner
# =========================

capture_benchmark_env <- function() {
  list(sys_info = Sys.info(), r_version = R.version.string, session_info = utils::capture.output(sessionInfo()))
}

benchmark_irregular_support_v4 <- function(nr, nc,
                                           scenario_list,
                                           contiguity = c("queen", "rook"),
                                           reps = 1L,
                                           timeout_seconds = Inf,
                                           baseline_nb2mat_style = c("B", "W"),
                                           cache_full = FALSE,
                                           cache_baseline = FALSE,
                                           verbose = TRUE,
                                           skip_larger_baseline_after_fail = TRUE,
                                           critical_warning_patterns = c("integer overflow", "cannot allocate", "memory", "overflow"),
                                           out_dir = NULL,
                                           log_progress = TRUE,
                                           write_excel = TRUE,
                                           write_csv = TRUE,
                                           log_prefix = "benchmark_v4_aligned") {
  contiguity <- match.arg(contiguity)
  baseline_nb2mat_style <- match.arg(baseline_nb2mat_style)

  if (is.null(names(scenario_list))) names(scenario_list) <- paste0("scenario_", seq_along(scenario_list))
  cached_W_full <- NULL
  cached_nb_full <- NULL
  if (cache_full) cached_W_full <- build_W_full(nr, nc, type = contiguity)
  if (cache_baseline) cached_nb_full <- spdep::cell2nb(nrow = nr, ncol = nc, type = contiguity, torus = FALSE)

  full_n <- max(vapply(scenario_list, length, integer(1)))
  rows_out <- list()
  k <- 1L
  baseline_hard_failed_from_here <- FALSE

  for (sc_name in names(scenario_list)) {
    idx_support <- sort(unique(as.integer(scenario_list[[sc_name]])))
    pct_support <- length(idx_support) / full_n

    for (rep_id in seq_len(reps)) {
      if (verbose) cat(sprintf("Scenario %s | rep %d/%d | cells = %d
", sc_name, rep_id, reps, length(idx_support)))

      prop_res <- run_proposed_once_v4(
        nr = nr, nc = nc, idx_support = idx_support,
        contiguity = contiguity, timeout_seconds = timeout_seconds,
        cached_W_full = cached_W_full,
        critical_warning_patterns = critical_warning_patterns
      )
      prop_res$scenario <- sc_name
      prop_res$rep <- rep_id
      prop_res$pct_support <- pct_support
      prop_res$bytes_W_full_mb <- mb_num(prop_res$bytes_W_full)
      prop_res$bytes_P_mb <- mb_num(prop_res$bytes_P)
      prop_res$bytes_W_common_mb <- mb_num(prop_res$bytes_W_common)
      prop_res$bytes_W_final_mb <- mb_num(prop_res$bytes_W_final)
      prop_res$peak_object_mb <- mb_num(prop_res$peak_object_bytes)
      rows_out[[k]] <- prop_res
      k <- k + 1L

      if (skip_larger_baseline_after_fail && baseline_hard_failed_from_here) {
        base_res <- data.frame(
          method = "baseline_cell2nb_nb2mat_v4",
          status = "skipped_after_hard_failure",
          error = NA_character_, warning_message = NA_character_,
          failure_class = "hard", failure_stage = NA_character_,
          cells_input = length(idx_support), cells_after_prune = NA_real_, n_isolates = NA_real_,
          nnz_common = NA_real_, nnz_final = NA_real_,
          elapsed_cell2nb = NA_real_, elapsed_subset_nb = NA_real_, elapsed_prune_nb = NA_real_, elapsed_nb2mat = NA_real_,
          elapsed_total = NA_real_, bytes_nb_full = NA_real_, bytes_nb_common = NA_real_, bytes_nb_pruned = NA_real_,
          bytes_W_final = NA_real_, peak_object_bytes = NA_real_,
          scenario = sc_name, rep = rep_id, pct_support = pct_support,
          bytes_nb_full_mb = NA_real_, bytes_nb_common_mb = NA_real_,
          bytes_nb_pruned_mb = NA_real_, bytes_W_final_mb = NA_real_,
          peak_object_mb = NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        base_res <- run_baseline_once_v4(
          nr = nr, nc = nc, idx_support = idx_support,
          contiguity = contiguity, nb2mat_style = baseline_nb2mat_style,
          timeout_seconds = timeout_seconds,
          cached_nb_full = cached_nb_full,
          critical_warning_patterns = critical_warning_patterns
        )
        base_res$scenario <- sc_name
        base_res$rep <- rep_id
        base_res$pct_support <- pct_support
        base_res$bytes_nb_full_mb <- mb_num(base_res$bytes_nb_full)
        base_res$bytes_nb_common_mb <- mb_num(base_res$bytes_nb_common)
        base_res$bytes_nb_pruned_mb <- mb_num(base_res$bytes_nb_pruned)
        base_res$bytes_W_final_mb <- mb_num(base_res$bytes_W_final)
        base_res$peak_object_mb <- mb_num(base_res$peak_object_bytes)
        if (skip_larger_baseline_after_fail &&
            identical(base_res$status[1], "failed") &&
            identical(base_res$failure_class[1], "hard")) {
          baseline_hard_failed_from_here <- TRUE
        }
      }

      rows_out[[k]] <- base_res
      k <- k + 1L
      gc()
    }

    if (isTRUE(log_progress) && !is.null(out_dir)) {
      raw_now <- bind_rows_fill(rows_out)
      raw_now <- raw_now[order(raw_now$scenario, raw_now$method, raw_now$rep), , drop = FALSE]
      num_cols <- setdiff(names(raw_now), c("method", "status", "error", "warning_message", "failure_class", "failure_stage", "scenario"))
      for (nm in num_cols) raw_now[[nm]] <- suppressWarnings(as.numeric(raw_now[[nm]]))
      summary_now <- summarize_benchmark_raw_v4(raw_now)
      write_benchmark_progress_v4(
        raw = raw_now,
        summary_df = summary_now,
        out_dir = out_dir,
        prefix = log_prefix,
        scenario_name = sc_name,
        write_excel = write_excel,
        write_csv = write_csv
      )
    }
  }

  raw <- bind_rows_fill(rows_out)
  raw <- raw[order(raw$scenario, raw$method, raw$rep), , drop = FALSE]
  num_cols <- setdiff(names(raw), c("method", "status", "error", "warning_message", "failure_class", "failure_stage", "scenario"))
  for (nm in num_cols) raw[[nm]] <- suppressWarnings(as.numeric(raw[[nm]]))

  summary_df <- summarize_benchmark_raw_v4(raw)

  if (!is.null(out_dir)) {
    write_benchmark_progress_v4(
      raw = raw,
      summary_df = summary_df,
      out_dir = out_dir,
      prefix = log_prefix,
      scenario_name = "final",
      write_excel = write_excel,
      write_csv = write_csv
    )
    saveRDS(list(raw = raw, summary = summary_df, scenarios = scenario_list),
            file.path(out_dir, paste0(log_prefix, "_results_latest.rds")))
  }

  list(raw = raw, summary = summary_df, env = capture_benchmark_env(), scenarios = scenario_list)
}

# =========================
# Convenience wrapper to build v4-aligned inputs and scenarios
# =========================

prepare_benchmark_v4_indonesia <- function(proportions = c(0.05, 0.10, 0.25, 0.35, 0.50, 0.75, 1.00),
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

# =========================
# Standalone example without loading the v4 script
# =========================

run_standalone_example <- function() {
  cfg <- list(
    shapefile_admin2 = "D:/SHAPEFILE/gadm41_IDN_2.shp",
    data_root = "E:/Spatial Sparse W Matrix Research/DATA",
    landscan_root = "F:/LANDSCAN DATA",
    gdp_file = "E:/Spatial Sparse W Matrix Research/DATA/rast_gdp_tot_1990_2020_30arcsec.tif"
  )

  root_dir <- "E:/Spatial Sparse W Matrix Research/v4 Outputs/benchmark_v4_aligned"
  cache_dir <- file.path(root_dir, "cache")
  log_dir <- file.path(root_dir, "progress_logs")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

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
  cat("Rows x Cols         :", prep$support$r, "x", prep$support$c, "
")
  cat("Template cells      :", prep$support$N, "
")
  cat("Full common support :", prep$support$V, "
")

  bm <- benchmark_irregular_support_v4(
    nr = prep$support$r,
    nc = prep$support$c,
    scenario_list = prep$scenarios,
    contiguity = "queen",
    reps = 2,
    timeout_seconds = 7200,
    baseline_nb2mat_style = "B",
    cache_full = FALSE,
    cache_baseline = TRUE,
    verbose = TRUE,
    skip_larger_baseline_after_fail = TRUE,
    out_dir = log_dir,
    log_progress = TRUE,
    write_excel = TRUE,
    write_csv = TRUE,
    log_prefix = "benchmark_v4_aligned"
  )

  print(bm$summary)

  compact_summary <- make_compact_summary_v4(bm$summary)
  metadata_df <- benchmark_metadata_v4()

  write.csv(bm$raw, file.path(root_dir, "benchmark_raw_v4_aligned.csv"), row.names = FALSE)
  write.csv(bm$summary, file.path(root_dir, "benchmark_summary_v4_aligned.csv"), row.names = FALSE)
  write.csv(compact_summary, file.path(root_dir, "benchmark_summary_compact_v4_aligned.csv"), row.names = FALSE)
  write.csv(metadata_df, file.path(root_dir, "benchmark_metadata_v4_aligned.csv"), row.names = FALSE)
  
  write.xlsx(bm$raw, file.path(root_dir, "benchmark_raw_v4_aligned.xlsx"))
  write.xlsx(bm$summary, file.path(root_dir, "benchmark_summary_v4_aligned.xlsx"))
  write.xlsx(compact_summary, file.path(root_dir, "benchmark_summary_compact_v4_aligned.xlsx"))
  write.xlsx(metadata_df, file.path(root_dir, "benchmark_metadata_v4_aligned.xlsx"))
  
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(
      list(raw = bm$raw, summary = bm$summary, compact_summary = compact_summary, metadata = metadata_df),
      file.path(root_dir, "benchmark_v4_aligned_results.xlsx")
    )
  }
  saveRDS(bm, file.path(root_dir, "benchmark_results_v4_aligned.rds"))
  saveRDS(prep, file.path(root_dir, "benchmark_inputs_v4_aligned.rds"))

  invisible(list(prep = prep, bm = bm, compact_summary = compact_summary, metadata = metadata_df))
}

res <- run_standalone_example()