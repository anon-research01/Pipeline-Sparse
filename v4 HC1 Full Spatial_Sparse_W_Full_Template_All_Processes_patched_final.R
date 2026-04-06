# ============================================================
# Full reviewer-ready template: preprocessing + sparse common-W
# SAR and SDM estimation for the 2020 four-layer application
# ------------------------------------------------------------
# This version integrates the shapefile/raster preprocessing
# steps with the cleaned sparse SAR/SDM estimation workflow.
# It is written so that the full pipeline can be run from the
# source rasters to the exported workbooks in `out_path`.
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(terra)
  library(spdep)
  library(spatialreg)
  library(openxlsx)
})

out_path <- "E:/Spatial Sparse W Matrix Research/v4 Full Elapsed HC1 Outputs/sbs_m = 500 SBS_K = 30"

# =========================
# Basic utilities
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
  Ir <- Diagonal(r)
  Ic <- Diagonal(c)
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
  W <- Diagonal(x = 1 / rs) %*% W
  list(W = W, keep = keep)
}

dgC_to_listw <- function(W_dgC) {
  stopifnot(inherits(W_dgC, "dgCMatrix"))
  n <- nrow(W_dgC)
  p <- W_dgC@p
  i <- W_dgC@i
  x <- W_dgC@x
  
  nb_row <- vector("list", n)
  gl_row <- vector("list", n)
  
  for (col in seq_len(n)) {
    idx <- (p[col] + 1L):p[col + 1L]
    if (!length(idx)) next
    rows <- i[idx] + 1L
    vals <- x[idx]
    
    for (k in seq_along(rows)) {
      r <- rows[k]
      nb_row[[r]] <- c(nb_row[[r]], col)
      gl_row[[r]] <- c(gl_row[[r]], vals[k])
    }
  }
  
  # Convert empty neighbour sets into explicit no-neighbour entries
  for (r in seq_len(n)) {
    if (length(nb_row[[r]]) == 0L) {
      nb_row[[r]] <- 0L
      gl_row[[r]] <- numeric(0)
    }
  }
  
  attr(nb_row, "region.id") <- as.character(seq_len(n))
  attr(nb_row, "call") <- match.call()
  attr(nb_row, "type") <- "user"
  attr(nb_row, "sym") <- FALSE
  class(nb_row) <- "nb"
  
  nb2listw(nb_row, glist = gl_row, style = "W", zero.policy = TRUE)
}

.round_numeric_df <- function(df, digits = 6) {
  out <- df
  is_num <- vapply(out, is.numeric, logical(1))
  out[is_num] <- lapply(out[is_num], round, digits = digits)
  out
}

.format_p <- function(p, digits = 4, eps = 1e-16) {
  format.pval(p, digits = digits, eps = eps)
}

.p_from_z <- function(z) {
  2 * pnorm(abs(z), lower.tail = FALSE)
}

# Convert a single-layer raster to the same column-major vector convention
# used in the sparse projection step. By default the raster is first aligned
# to the template; if the input has already been resampled upstream, set
# already_aligned = TRUE to avoid a second resampling pass.
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
# Trace / log-determinant helpers
# =========================

trace_powers_mc <- function(W, kmax = 50, m = 30, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  traces <- numeric(kmax)

  for (b in seq_len(m)) {
    v <- sample(c(-1, 1), nrow(W), replace = TRUE)
    z <- v
    for (k in seq_len(kmax)) {
      z <- W %*% z
      traces[k] <- traces[k] + sum(v * z)
    }
  }

  traces / m
}

logdet_I_minus_rhoW <- function(W, rho, method = c("MC", "LU"), trW = NULL) {
  method <- match.arg(method)

  if (method == "MC") {
    if (is.null(trW)) stop("trW must be supplied for method = 'MC'.")
    k <- seq_along(trW)
    return(-sum((rho^k / k) * trW))
  }

  S <- Diagonal(nrow(W)) - rho * W
  dd <- try(determinant(S, logarithm = TRUE), silent = TRUE)
  if (inherits(dd, "try-error")) return(-Inf)
  as.numeric(dd$modulus)
}

solve_S <- function(W, rho, b) {
  S <- Diagonal(nrow(W)) - rho * W
  drop(solve(S, b, sparse = TRUE))
}

avg_diag_Sinv_series <- function(rho, trW, n) {
  (n + sum((rho^(seq_along(trW))) * trW)) / n
}

avg_rowsum_Sinv <- function(W, rho, row_standardized = TRUE) {
  if (row_standardized) return(1 / (1 - rho))
  mean(solve_S(W, rho, rep(1, nrow(W))))
}

# =========================
# Manual SAR MLE
# =========================

sar_mle_manual <- function(y, X, W,
                           rho_bounds = c(-0.99, 0.99),
                           logdet_method = c("MC", "LU"),
                           trW = NULL,
                           kmax = 50,
                           m_mc = 30,
                           seed = 123,
                           row_standardized = TRUE,
                           compute_impacts = TRUE,
                           R_impacts = 0) {
  logdet_method <- match.arg(logdet_method)
  n <- length(y)
  p <- ncol(X)

  if (logdet_method == "MC" && is.null(trW)) {
    trW <- trace_powers_mc(W, kmax = kmax, m = m_mc, seed = seed)
  }

  XtX <- crossprod(X)
  Xt  <- t(X)
  
  .solve_beta_sar <- function(Sy) {
    rhs <- Xt %*% Sy
    
    beta <- tryCatch(
      solve(XtX, rhs),
      error = function(e) {
        tryCatch(
          qr.solve(X, Sy, tol = 1e-10),
          error = function(e2) {
            if (!requireNamespace("MASS", quietly = TRUE)) stop(e2)
            MASS::ginv(X) %*% Sy
          }
        )
      }
    )
    
    beta <- as.numeric(beta)
    
    if (length(beta) != ncol(X)) {
      tmp <- rep(NA_real_, ncol(X))
      tmp[seq_len(min(length(beta), ncol(X)))] <- beta[seq_len(min(length(beta), ncol(X)))]
      beta <- tmp
    }
    
    beta[!is.finite(beta)] <- 0
    beta
  }

  ll_prof <- function(rho) {
    Sy <- y - rho * drop(W %*% y)
    beta <- .solve_beta_sar(Sy)
    e <- Sy - drop(X %*% beta)
    s2 <- sum(e * e) / n
    if (!is.finite(s2) || s2 <= 0) return(-Inf)
    ldS <- logdet_I_minus_rhoW(W, rho, method = logdet_method, trW = trW)
    if (!is.finite(ldS)) return(-Inf)
    ldS - (n / 2) * log(s2)
  }

  obj_prof <- function(r) {
    val <- ll_prof(r)
    if (is.finite(val)) -val else 1e100
  }
  
  opt <- optimize(obj_prof, interval = rho_bounds, maximum = FALSE)
  rho_hat <- opt$minimum
  ll_hat  <- ll_prof(rho_hat)

  Sy <- y - rho_hat * drop(W %*% y)
  beta_hat <- .solve_beta_sar(Sy)
  e <- Sy - drop(X %*% beta_hat)
  sigma2 <- sum(e * e) / n

  Vbeta <- tryCatch(
    sigma2 * solve(XtX),
    error = function(e) {
      if (!requireNamespace("MASS", quietly = TRUE)) stop(e)
      sigma2 * MASS::ginv(as.matrix(XtX))
    }
  )
  se_beta <- sqrt(diag(Vbeta))

  se_rho <- NA_real_
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    H <- try(numDeriv::hessian(ll_prof, x = rho_hat), silent = TRUE)
    if (!inherits(H, "try-error") && is.finite(H)) {
      info <- -H
      if (info > 0) se_rho <- sqrt(1 / info)
    }
  }

  k_par <- p + 2
  logLik <- ll_hat - (n / 2) * log(2 * pi) - n / 2
  AIC <- -2 * logLik + 2 * k_par
  BIC <- -2 * logLik + log(n) * k_par
  fitted <- drop(X %*% beta_hat)

  impacts <- NULL
  if (compute_impacts) {
    direct_factor <- avg_diag_Sinv_series(rho_hat, trW, n)
    total_factor  <- avg_rowsum_Sinv(W, rho_hat, row_standardized = row_standardized)
    indirect_factor <- total_factor - direct_factor

    names_x <- colnames(X)
    if (is.null(names_x)) names_x <- paste0("x", seq_len(p))
    idx_coef <- which(names_x != "Intercept")

    impacts <- data.frame(
      variable = names_x[idx_coef],
      direct   = direct_factor * beta_hat[idx_coef],
      indirect = indirect_factor * beta_hat[idx_coef],
      total    = total_factor * beta_hat[idx_coef],
      row.names = NULL
    )

    if (R_impacts > 0) {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' is required when R_impacts > 0.")
      }
      set.seed(seed)
      draws <- MASS::mvrnorm(
        R_impacts,
        mu = c(beta_hat[idx_coef], rho_hat),
        Sigma = diag(c(diag(Vbeta)[idx_coef], ifelse(is.na(se_rho), 0, se_rho^2))),
        empirical = FALSE
      )
      bd <- ncol(draws) - 1
      dir_draw <- ind_draw <- tot_draw <- matrix(NA_real_, R_impacts, bd)
      for (r in seq_len(R_impacts)) {
        br <- draws[r, seq_len(bd)]
        rr <- draws[r, bd + 1]
        df <- avg_diag_Sinv_series(rr, trW, n)
        tf <- avg_rowsum_Sinv(W, rr, row_standardized)
        idf <- tf - df
        dir_draw[r, ] <- df * br
        ind_draw[r, ] <- idf * br
        tot_draw[r, ] <- tf * br
      }
      impacts$se_direct   <- apply(dir_draw, 2, sd, na.rm = TRUE)
      impacts$se_indirect <- apply(ind_draw, 2, sd, na.rm = TRUE)
      impacts$se_total    <- apply(tot_draw, 2, sd, na.rm = TRUE)
      impacts$z_direct    <- impacts$direct / impacts$se_direct
      impacts$z_indirect  <- impacts$indirect / impacts$se_indirect
      impacts$z_total     <- impacts$total / impacts$se_total
      impacts$p_direct    <- .p_from_z(impacts$z_direct)
      impacts$p_indirect  <- .p_from_z(impacts$z_indirect)
      impacts$p_total     <- .p_from_z(impacts$z_total)
    }
  }

  out <- list(
    model = "SAR",
    n = n,
    k = p,
    y = y,
    X = X,
    W = W,
    rho = rho_hat,
    beta = beta_hat,
    sigma2 = sigma2,
    se_beta = se_beta,
    se_rho = se_rho,
    vcov_beta = Vbeta,
    logLik = logLik,
    AIC = AIC,
    BIC = BIC,
    residuals = e,
    fitted = fitted,
    logdet_method = logdet_method,
    trW = trW,
    impacts = impacts
  )
  class(out) <- "sar_manual"
  out
}

# =========================
# Manual SDM MLE
# =========================

sdm_mle_manual <- function(y, X, W,
                           rho_bounds = c(-0.99, 0.99),
                           logdet_method = c("MC", "LU"),
                           trW = NULL,
                           kmax = 50,
                           m_mc = 30,
                           seed = 123,
                           row_standardized = TRUE,
                           compute_impacts = TRUE,
                           R_impacts = 0) {
  
  logdet_method <- match.arg(logdet_method)
  n <- length(y)
  
  if (logdet_method == "MC" && is.null(trW)) {
    trW <- trace_powers_mc(W, kmax = kmax, m = m_mc, seed = seed)
  }
  
  cnX <- colnames(X)
  if (is.null(cnX)) cnX <- paste0("x", seq_len(ncol(X)))
  
  int_idx <- which(tolower(cnX) %in% c("(intercept)", "intercept"))
  if (!length(int_idx)) int_idx <- 1L
  
  X0 <- X[, -int_idx, drop = FALSE]
  WX0 <- as.matrix(W %*% X0)
  Z <- cbind(X, WX0)
  
  cnZ <- c(cnX, paste0("W_", colnames(X0)))
  if (length(cnZ) != ncol(Z)) {
    cnZ <- paste0("z", seq_len(ncol(Z)))
  }
  colnames(Z) <- cnZ
  
  ZtZ <- crossprod(Z)
  Zt  <- t(Z)
  
  .solve_theta_sdm <- function(Sy) {
    rhs <- Zt %*% Sy
    
    theta <- tryCatch(
      solve(ZtZ, rhs),
      error = function(e) {
        tryCatch(
          qr.solve(Z, Sy, tol = 1e-10),
          error = function(e2) {
            if (!requireNamespace("MASS", quietly = TRUE)) stop(e2)
            MASS::ginv(Z) %*% Sy
          }
        )
      }
    )
    
    theta <- as.numeric(theta)
    
    if (length(theta) != ncol(Z)) {
      tmp <- rep(NA_real_, ncol(Z))
      tmp[seq_len(min(length(theta), ncol(Z)))] <- theta[seq_len(min(length(theta), ncol(Z)))]
      theta <- tmp
    }
    
    theta[!is.finite(theta)] <- 0
    theta
  }
  
  ll_prof <- function(rho) {
    Sy <- y - rho * drop(W %*% y)
    theta <- .solve_theta_sdm(Sy)
    e <- Sy - drop(Z %*% theta)
    s2 <- sum(e * e) / n
    
    if (!is.finite(s2) || s2 <= 0) return(-Inf)
    
    ldS <- logdet_I_minus_rhoW(W, rho, method = logdet_method, trW = trW)
    if (!is.finite(ldS)) return(-Inf)
    
    ldS - (n / 2) * log(s2)
  }
  
  obj_prof <- function(r) {
    val <- ll_prof(r)
    if (is.finite(val)) -val else 1e100
  }
  
  opt <- optimize(obj_prof, interval = rho_bounds, maximum = FALSE)
  rho_hat <- opt$minimum
  
  Sy <- y - rho_hat * drop(W %*% y)
  theta_hat <- .solve_theta_sdm(Sy)
  e <- Sy - drop(Z %*% theta_hat)
  sigma2 <- sum(e * e) / n
  
  Vtheta <- tryCatch(
    sigma2 * solve(ZtZ),
    error = function(e) {
      if (!requireNamespace("MASS", quietly = TRUE)) stop(e)
      sigma2 * MASS::ginv(as.matrix(ZtZ))
    }
  )
  se_theta <- sqrt(pmax(diag(Vtheta), 0))
  
  se_rho <- NA_real_
  if (requireNamespace("numDeriv", quietly = TRUE)) {
    H <- try(numDeriv::hessian(ll_prof, x = rho_hat), silent = TRUE)
    if (!inherits(H, "try-error") && length(H) == 1L && is.finite(H)) {
      info <- -as.numeric(H)
      if (is.finite(info) && info > 0) {
        se_rho <- sqrt(1 / info)
      }
    }
  }
  
  k_par <- length(theta_hat) + 2
  ll_hat <- ll_prof(rho_hat)
  logLik <- ll_hat - (n / 2) * log(2 * pi) - n / 2
  AIC <- -2 * logLik + 2 * k_par
  BIC <- -2 * logLik + log(n) * k_par
  fitted <- drop(Z %*% theta_hat)
  
  impacts <- NULL
  if (compute_impacts) {
    direct_Sinv <- avg_diag_Sinv_series(rho_hat, trW, n)
    
    tr_SinvW <- sum((rho_hat^(seq_along(trW))) * c(trW[-length(trW)], 0))
    direct_SinvW <- tr_SinvW / n
    
    total_Sinv <- avg_rowsum_Sinv(W, rho_hat, row_standardized)
    total_SinvW <- total_Sinv
    
    pX <- ncol(X)
    pWX <- ncol(X0)
    idx_X <- seq_len(pX)
    idx_WX <- pX + seq_len(pWX)
    
    vars <- colnames(X0)
    if (is.null(vars)) vars <- paste0("x", seq_len(pWX))
    
    x_no_intercept_idx <- setdiff(seq_len(pX), int_idx)
    
    get_row <- function(j) {
      beta_j  <- theta_hat[idx_X[x_no_intercept_idx[j]]]
      theta_j <- theta_hat[idx_WX[j]]
      
      direct   <- direct_Sinv  * beta_j + direct_SinvW * theta_j
      total    <- total_Sinv   * beta_j + total_SinvW  * theta_j
      indirect <- total - direct
      
      c(direct = direct, indirect = indirect, total = total)
    }
    
    M <- t(vapply(seq_len(pWX), get_row, numeric(3)))
    impacts <- data.frame(
      variable = vars,
      direct = M[, 1],
      indirect = M[, 2],
      total = M[, 3]
    )
  }
  
  out <- list(
    model = "SDM",
    n = n,
    y = y,
    X = X,
    W = W,
    rho = rho_hat,
    theta = theta_hat,
    sigma2 = sigma2,
    se_theta = se_theta,
    se_rho = se_rho,
    vcov_theta = Vtheta,
    logLik = logLik,
    AIC = AIC,
    BIC = BIC,
    residuals = e,
    fitted = fitted,
    logdet_method = logdet_method,
    trW = trW,
    impacts = impacts,
    coef_names = cnZ
  )
  class(out) <- "sdm_manual"
  out
}

# =========================
# Print methods
# =========================

print.sar_manual <- function(x, ...) {
  cat("Manual SAR (asymmetric W allowed)\n")
  cat(sprintf("Observations: %d, Parameters: %d + rho + sigma2\n", x$n, x$k))
  cat(sprintf("Method logdet: %s\n", x$logdet_method))
  cat(sprintf("LogLik: %.3f  AIC: %.3f  BIC: %.3f\n\n", x$logLik, x$AIC, x$BIC))

  cn <- colnames(x$X)
  if (is.null(cn)) cn <- paste0("x", seq_len(x$k))
  tb <- cbind(
    Estimate = x$beta,
    `Std. Error` = x$se_beta,
    `z value` = x$beta / x$se_beta,
    `Pr(>|z|)` = .p_from_z(x$beta / x$se_beta)
  )
  rownames(tb) <- cn
  print(.round_numeric_df(as.data.frame(tb), 4))
  cat("\n")

  cat(sprintf("rho = %.6f", x$rho))
  if (!is.na(x$se_rho)) {
    z <- x$rho / x$se_rho
    cat(sprintf("  (SE = %.6f, z = %.3f, p = %.3g)\n", x$se_rho, z, .p_from_z(z)))
  } else {
    cat("\n")
  }
  cat(sprintf("sigma^2 = %.6f\n", x$sigma2))

  if (!is.null(x$impacts)) {
    cat("\nAverage Impacts (LeSage-Pace):\n")
    print(.round_numeric_df(x$impacts, 6))
  }
  invisible(x)
}

print.sdm_manual <- function(x, ...) {
  cat("Manual SDM (asymmetric W allowed)\n")
  cat(sprintf("Observations: %d\n", x$n))
  cat(sprintf("Method logdet: %s\n", x$logdet_method))
  cat(sprintf("LogLik: %.3f  AIC: %.3f  BIC: %.3f\n\n", x$logLik, x$AIC, x$BIC))

  names_theta <- x$coef_names
  if (is.null(names_theta)) {
    cnX <- colnames(x$X)
    if (is.null(cnX)) cnX <- paste0("x", seq_len(ncol(x$X)))
    int_idx <- which(tolower(cnX) %in% c("intercept", "(intercept)"))
    if (!length(int_idx)) int_idx <- 1L
    names_theta <- c(cnX, paste0("W_", cnX[-int_idx]))
  }

  tb <- cbind(
    Estimate = x$theta,
    `Std. Error` = x$se_theta,
    `z value` = x$theta / x$se_theta,
    `Pr(>|z|)` = .p_from_z(x$theta / x$se_theta)
  )
  rownames(tb) <- names_theta
  print(.round_numeric_df(as.data.frame(tb), 4))
  cat("\n")

  cat(sprintf("rho = %.6f", x$rho))
  if (!is.na(x$se_rho)) {
    z <- x$rho / x$se_rho
    cat(sprintf("  (SE = %.6f, z = %.3f, p = %.3g)\n", x$se_rho, z, .p_from_z(z)))
  } else {
    cat("\n")
  }
  cat(sprintf("sigma^2 = %.6f\n", x$sigma2))

  if (!is.null(x$impacts)) {
    cat("\nAverage Impacts (LeSage-Pace):\n")
    print(.round_numeric_df(x$impacts, 6))
  }
  invisible(x)
}

# =========================
# Generic region wrapper
# =========================
# The functions below apply the same sparse-matrix workflow to any study area
# defined on a raster template. This design is intentional: Indonesia-wide and
# island-level applications differ only in the raster inputs, not in the
# construction of the common support, the projection step, the row-standardized
# sparse weights matrix, or the SAR/SDM estimation logic.
#
# Key arguments:
# - template: reference raster defining the full rectangular lattice on which
#   the initial sparse neighbourhood matrix is constructed.
# - rast_gdp: dependent-variable raster. In the manuscript's illustration, GDP
#   is transformed using log1p() after the retained common-support cells are set.
# - rast_pop, rast_ndvi: explanatory rasters aligned to the same template.
# - logdet_method: method used for the log-determinant term in the manual sparse
#   likelihood. "MC" denotes the Monte Carlo trace approximation.
# - type_W: neighbourhood rule used to construct the initial lattice-based
#   weights matrix on the full template ("queen" or "rook").
# - engine: estimation backend. "manual" uses the custom sparse likelihood
#   implementation, while "spatialreg" uses the package backend where feasible.
# - seed: random seed used for reproducibility of stochastic steps.
#
# The outputs report three sample sizes that should be distinguished in the
# manuscript and reviewer materials:
# 1. total_template_cells: full rectangular grid size;
# 2. common_support_cells: cells jointly observed for all regression variables;
# 3. retained_cells: common-support cells remaining after isolates are removed.
# This distinction helps avoid conflating the full raster footprint with the
# final econometric estimation support.

run_sar_sdm_region <- function(
    template,
    rast_gdp, rast_pop, rast_ndvi, rast_mndwi,
    region_name = "Region",
    resample_pop = "near",
    resample_ndvi = "near",
    resample_mndwi = "near",
    logdet_method = c("MC", "LU"),
    type_W = c("queen", "rook"),
    transform_y = function(z) log1p(z),
    scale_X = TRUE,
    engine = c("manual", "spatialreg"),
    use_MC_trace = TRUE,
    impacts_R = 1000,
    seed = 123) {
  
  logdet_method <- match.arg(logdet_method)
  type_W <- match.arg(type_W)
  engine <- match.arg(engine)
  
  # All variables are aligned to the same template so that the empirical
  # support is defined as one common intersection across GDP, population,
  # NDVI, and MNDWI.
  rast_gdp_a   <- terra::resample(rast_gdp,   template, method = "near")
  rast_pop_a   <- terra::resample(rast_pop,   template, method = resample_pop)
  rast_ndvi_a  <- terra::resample(rast_ndvi,  template, method = resample_ndvi)
  rast_mndwi_a <- terra::resample(rast_mndwi, template, method = resample_mndwi)
  
  r <- nrow(template)
  c <- ncol(template)
  N <- r * c
  
  # Sparse lattice matrix on the full rectangular template
  W_full <- build_W_full(r, c, type_W)
  
  # Extract aligned values in the same column-major order used by W
  gdp_vec   <- raster_to_colmajor(rast_gdp_a,   template, var_name = "rast_gdp",   already_aligned = TRUE)
  pop_vec   <- raster_to_colmajor(rast_pop_a,   template, var_name = "rast_pop",   already_aligned = TRUE)
  ndvi_vec  <- raster_to_colmajor(rast_ndvi_a,  template, var_name = "rast_ndvi",  already_aligned = TRUE)
  mndwi_vec <- raster_to_colmajor(rast_mndwi_a, template, var_name = "rast_mndwi", already_aligned = TRUE)
  
  valid_gdp   <- !is.na(gdp_vec)
  valid_pop   <- !is.na(pop_vec)
  valid_ndvi  <- !is.na(ndvi_vec)
  valid_mndwi <- !is.na(mndwi_vec)
  
  # One common support used throughout
  valid_all <- valid_gdp & valid_pop & valid_ndvi & valid_mndwi
  
  idx_all <- which(valid_all)
  V <- length(idx_all)
  if (V == 0L) stop("Common support is empty. Check masks and raster alignment.")
  
  # Projection from full template to observed common support
  P <- sparseMatrix(i = idx_all, j = seq_len(V), x = 1, dims = c(N, V))
  W_common <- t(P) %*% W_full %*% P
  
  # Drop isolates before row-standardization
  rs_obj <- row_standardize_drop0(W_common)
  W_common <- rs_obj$W
  idx_kept <- idx_all[rs_obj$keep]
  
  y <- transform_y(gdp_vec[idx_kept])
  X_core <- cbind(
    pop   = pop_vec[idx_kept],
    ndvi  = ndvi_vec[idx_kept],
    mndwi = mndwi_vec[idx_kept]
  )
  if (scale_X) X_core <- scale(X_core)
  X <- cbind(Intercept = 1, X_core)
  
  set.seed(seed)
  
  if (engine == "spatialreg") {
    lw <- dgC_to_listw(W_common)
    trs <- NULL
    
    if (use_MC_trace) {
      tmp <- try(spatialreg::trW(lw, type = "MC"), silent = TRUE)
      if (!inherits(tmp, "try-error")) trs <- tmp
    }
    
    t0_sar <- Sys.time()
    fit_sar <- spatialreg::lagsarlm(
      y ~ X[, -1, drop = FALSE],
      listw = lw,
      method = "LU",
      trs = trs,
      control = list(super = TRUE)
    )
    t1_sar <- Sys.time()
    
    t0_sdm <- Sys.time()
    fit_sdm <- spatialreg::lagsarlm(
      y ~ X[, -1, drop = FALSE],
      listw = lw,
      type = "mixed",
      Durbin = TRUE,
      method = "LU",
      trs = trs,
      control = list(super = TRUE)
    )
    t1_sdm <- Sys.time()
    
    t0_imp_sar <- Sys.time()
    impacts_sar <- impacts(fit_sar, listw = lw, R = impacts_R)
    t1_imp_sar <- Sys.time()
    
    t0_imp_sdm <- Sys.time()
    impacts_sdm <- impacts(fit_sdm, listw = lw, R = impacts_R)
    t1_imp_sdm <- Sys.time()
    
    elapsed_sar_sec <- as.numeric(difftime(t1_sar, t0_sar, units = "secs"))
    elapsed_sdm_sec <- as.numeric(difftime(t1_sdm, t0_sdm, units = "secs"))
    elapsed_impacts_sar_sec <- as.numeric(difftime(t1_imp_sar, t0_imp_sar, units = "secs"))
    elapsed_impacts_sdm_sec <- as.numeric(difftime(t1_imp_sdm, t0_imp_sdm, units = "secs"))
    
  } else {
    # Under the custom engine, Monte Carlo traces are computed once and reused
    trW <- trace_powers_mc(W_common, kmax = 60, m = 40, seed = seed)
    
    t0_sar <- Sys.time()
    fit_sar <- sar_mle_manual(
      y, X, W_common,
      logdet_method = logdet_method,
      trW = trW,
      rho_bounds = c(-0.99, 0.99),
      row_standardized = TRUE,
      compute_impacts = TRUE,
      R_impacts = 0,
      seed = seed
    )
    t1_sar <- Sys.time()
    
    t0_sdm <- Sys.time()
    fit_sdm <- sdm_mle_manual(
      y, X, W_common,
      logdet_method = logdet_method,
      trW = trW,
      rho_bounds = c(-0.99, 0.99),
      row_standardized = TRUE,
      compute_impacts = TRUE,
      R_impacts = 0,
      seed = seed
    )
    t1_sdm <- Sys.time()
    
    impacts_sar <- fit_sar$impacts
    impacts_sdm <- fit_sdm$impacts
    
    elapsed_sar_sec <- as.numeric(difftime(t1_sar, t0_sar, units = "secs"))
    elapsed_sdm_sec <- as.numeric(difftime(t1_sdm, t0_sdm, units = "secs"))
    
    # For the manual engine, impacts are computed inside the estimator
    elapsed_impacts_sar_sec <- NA_real_
    elapsed_impacts_sdm_sec <- NA_real_
  }
  
  list(
    region_name = region_name,
    engine = engine,
    W_common = W_common,
    y = as.numeric(y),
    X = as.matrix(X),
    idx_kept = idx_kept,
    valid_all = valid_all,
    dims = c(
      r = r, c = c, N = N,
      V_common = V,
      V_after_keep = nrow(W_common)
    ),
    fit_sar = fit_sar,
    fit_sdm = fit_sdm,
    impacts_sar = impacts_sar,
    impacts_sdm = impacts_sdm,
    seed = seed,
    type_W = type_W,
    scale_X = scale_X,
    mndwi_included = TRUE,
    total_template_cells = N,
    common_support_cells = V,
    retained_cells = nrow(W_common),
    dropped_isolates = V - nrow(W_common),
    share_common_support = V / N,
    share_retained_cells = nrow(W_common) / N,
    elapsed_sar_sec = elapsed_sar_sec,
    elapsed_sdm_sec = elapsed_sdm_sec,
    elapsed_impacts_sar_sec = elapsed_impacts_sar_sec,
    elapsed_impacts_sdm_sec = elapsed_impacts_sdm_sec,
    elapsed_models_total_sec = elapsed_sar_sec + elapsed_sdm_sec,
    elapsed_full_modeling_sec = sum(
      c(elapsed_sar_sec, elapsed_sdm_sec, elapsed_impacts_sar_sec, elapsed_impacts_sdm_sec),
      na.rm = TRUE
    )
  )
}

# Backward-compatible aliases are kept so that legacy region-specific calls can
# still be used in reviewer materials without altering the estimation logic.
run_sar_sdm_indonesia <- function(...) {
  run_sar_sdm_region(..., region_name = "Indonesia")
}

run_sar_sdm_java <- function(...) {
  run_sar_sdm_region(..., region_name = "Java")
}

# Convenience wrapper
# -------------------
# This helper performs estimation and workbook export in one call. It is used
# only to reduce repetition across regions. The underlying sample construction,
# estimator choice, and inference workflow remain unchanged.
#
# Additional arguments for reviewer interpretation:
# - dataset_name: label used in the exported workbook and worksheet metadata.
# - rast_mndwi: MNDWI raster treated as a regular explanatory variable, aligned
#   to the same template as population and NDVI.
# - sbs_K: number of graph-connected spatial subsamples used in the SBS step.
# - sbs_m: target number of retained cells within each connected subsample.
# - run_sbs: if TRUE, the workbook includes SBS-based inference alongside the
#   model-based MLE and HC1 standard errors.
run_and_export_region <- function(
    dataset_name,
    template,
    rast_gdp, rast_pop, rast_ndvi, rast_mndwi,
    out_dir,
    sbs_K = 30,
    sbs_m = 500,
    logdet_method = "MC",
    type_W = "queen",
    engine = "manual",
    seed = 42,
    run_sbs = TRUE,
    export_block_draws = TRUE,
    hc1_mode = c("conditional", "full", "both"),
    ...) {
  
  hc1_mode <- match.arg(hc1_mode)
  
  res <- run_sar_sdm_region(
    template = template,
    rast_gdp = rast_gdp,
    rast_pop = rast_pop,
    rast_ndvi = rast_ndvi,
    rast_mndwi = rast_mndwi,
    region_name = dataset_name,
    logdet_method = logdet_method,
    type_W = type_W,
    engine = engine,
    seed = seed,
    ...
  )
  
  export_obj <- export_sparse_outputs(
    res = res,
    dataset_name = dataset_name,
    out_dir = out_dir,
    sbs_K = sbs_K,
    sbs_m = sbs_m,
    seed = seed,
    run_sbs = run_sbs,
    export_block_draws = export_block_draws,
    hc1_mode = hc1_mode
  )
  
  invisible(list(
    res = res,
    export_obj = export_obj
  ))
}

# =========================
# HC1 robust SEs
# =========================

sar_add_sandwich_se <- function(fit_sar, y, X, W, type = c("HC0", "HC1")) {
  type <- match.arg(type)
  n <- length(y)
  k <- ncol(X)
  rho <- fit_sar$rho
  Sy <- as.numeric(y - rho * (W %*% y))
  e  <- as.numeric(Sy - X %*% fit_sar$beta)

  XtX <- crossprod(X)
  bread <- tryCatch(
    solve(XtX),
    error = function(e) {
      if (!requireNamespace("MASS", quietly = TRUE)) stop(e)
      MASS::ginv(XtX)
    }
  )
  meat <- crossprod(X * e, X * e)
  if (type == "HC1") meat <- meat * (n / (n - k))
  Vrob <- bread %*% meat %*% bread

  fit_sar$se_beta_hc <- sqrt(diag(Vrob))
  fit_sar$vcov_beta_hc <- Vrob
  fit_sar$hc_type <- type
  fit_sar
}

sdm_add_sandwich_se <- function(fit_sdm, y, X, W, type = c("HC0", "HC1")) {
  type <- match.arg(type)
  n <- length(y)
  rho <- fit_sdm$rho

  cnX <- colnames(X)
  if (is.null(cnX)) cnX <- paste0("x", seq_len(ncol(X)))
  int_idx <- which(tolower(cnX) %in% c("intercept", "(intercept)"))
  if (!length(int_idx)) int_idx <- 1L

  X0 <- X[, -int_idx, drop = FALSE]
  WX0 <- as.matrix(W %*% X0)
  Z <- cbind(X, WX0)

  Sy <- as.numeric(y - rho * (W %*% y))
  e  <- as.numeric(Sy - Z %*% fit_sdm$theta)

  ZtZ <- crossprod(Z)
  bread <- tryCatch(
    solve(ZtZ),
    error = function(e) {
      if (!requireNamespace("MASS", quietly = TRUE)) stop(e)
      MASS::ginv(ZtZ)
    }
  )
  meat <- crossprod(Z * e, Z * e)
  if (type == "HC1") meat <- meat * (n / (n - ncol(Z)))
  Vrob <- bread %*% meat %*% bread

  fit_sdm$se_theta_hc <- sqrt(diag(Vrob))
  fit_sdm$vcov_theta_hc <- Vrob
  fit_sdm$hc_type <- type
  fit_sdm
}

# Optional full-parameter HC1 extension
# -------------------------------------
# The original HC1 adjustment above is conditional on rho-hat and therefore
# applies only to the slope coefficients. The functions below implement an
# optional full-parameter sandwich approximation for c(rho, slopes), so that
# rho can also receive an HC1-type standard error when the user explicitly
# requests it. This extension is intended for robustness reporting and is not
# required for reproducing the main manuscript tables.

.safe_solve <- function(A) {
  tryCatch(
    solve(A),
    error = function(e) {
      if (!requireNamespace("MASS", quietly = TRUE)) stop(e)
      MASS::ginv(A)
    }
  )
}

.tr_SinvW_from_traces <- function(rho, trW) {
  if (is.null(trW) || !length(trW)) return(NA_real_)
  k <- seq_along(trW)
  sum((rho^(k - 1)) * trW)
}

sar_add_sandwich_se_full <- function(fit_sar, y, X, W, type = c("HC0", "HC1")) {
  type <- match.arg(type)
  
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for full HC1 on rho.")
  }
  
  n <- length(y)
  p <- ncol(X)
  q <- 1 + p
  
  rho_hat  <- as.numeric(fit_sar$rho)
  beta_hat <- as.numeric(fit_sar$beta)
  par_hat  <- c(rho = rho_hat, beta_hat)
  par_names <- c("rho", colnames(X))
  if (is.null(par_names[-1])) par_names <- c("rho", paste0("x", seq_len(p)))
  names(par_hat) <- par_names
  
  Wy <- as.numeric(W %*% y)
  
  ll_full <- function(par) {
    rho  <- as.numeric(par[1])
    beta <- as.numeric(par[-1])
    
    if (!is.finite(rho) || rho <= -0.99 || rho >= 0.99) return(-1e12)
    
    e  <- as.numeric(y - rho * Wy - X %*% beta)
    s2 <- mean(e^2)
    if (!is.finite(s2) || s2 <= 0) return(-1e12)
    
    ldS <- logdet_I_minus_rhoW(
      W, rho,
      method = fit_sar$logdet_method,
      trW = fit_sar$trW
    )
    if (!is.finite(ldS)) return(-1e12)
    
    ldS - (n / 2) * log(s2)
  }
  
  H <- try(numDeriv::hessian(ll_full, x = par_hat), silent = TRUE)
  if (inherits(H, "try-error")) {
    warning("Full HC1 Hessian failed for SAR. Returning fit unchanged.")
    return(fit_sar)
  }
  
  A <- -H
  Ainv <- .safe_solve(A)
  
  e  <- as.numeric(y - rho_hat * Wy - X %*% beta_hat)
  s2 <- mean(e^2)
  
  tr_SinvW <- .tr_SinvW_from_traces(rho_hat, fit_sar$trW)
  if (!is.finite(tr_SinvW)) {
    ld_fun <- function(r) logdet_I_minus_rhoW(
      W, r,
      method = fit_sar$logdet_method,
      trW = fit_sar$trW
    )
    d_ld <- try(numDeriv::grad(ld_fun, x = rho_hat), silent = TRUE)
    if (inherits(d_ld, "try-error") || !is.finite(d_ld)) {
      warning("Could not evaluate full HC1 rho score for SAR. Returning fit unchanged.")
      return(fit_sar)
    }
    tr_SinvW <- -as.numeric(d_ld)
  }
  
  S <- matrix(0, nrow = n, ncol = q)
  colnames(S) <- par_names
  
  S[, 1] <- (-tr_SinvW / n) + (Wy * e) / s2
  S[, -1] <- sweep(X, 1, e / s2, FUN = "*")
  
  B <- crossprod(S)
  if (type == "HC1") B <- B * (n / (n - q))
  
  Vrob <- Ainv %*% B %*% Ainv
  Vrob <- (Vrob + t(Vrob)) / 2
  
  se_full <- sqrt(pmax(diag(Vrob), 0))
  names(se_full) <- par_names
  
  fit_sar$hc_type_full <- type
  fit_sar$vcov_par_hc_full <- Vrob
  fit_sar$se_par_hc_full <- se_full
  fit_sar$se_rho_hc_full <- as.numeric(se_full["rho"])
  fit_sar$se_beta_hc_full <- as.numeric(se_full[-1])
  
  fit_sar
}

sdm_add_sandwich_se_full <- function(fit_sdm, y, X, W, type = c("HC0", "HC1")) {
  type <- match.arg(type)
  
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for full HC1 on rho.")
  }
  
  n <- length(y)
  rho_hat <- as.numeric(fit_sdm$rho)
  
  cnX <- colnames(X)
  if (is.null(cnX)) cnX <- paste0("x", seq_len(ncol(X)))
  int_idx <- which(tolower(cnX) %in% c("intercept", "(intercept)"))
  if (!length(int_idx)) int_idx <- 1L
  
  X0  <- X[, -int_idx, drop = FALSE]
  WX0 <- as.matrix(W %*% X0)
  Z   <- cbind(X, WX0)
  
  cnZ <- c(cnX, paste0("W_", colnames(X0)))
  if (any(is.na(cnZ)) || !length(cnZ)) cnZ <- paste0("z", seq_len(ncol(Z)))
  colnames(Z) <- cnZ
  
  theta_hat <- as.numeric(fit_sdm$theta)
  q <- 1 + length(theta_hat)
  
  par_hat <- c(rho = rho_hat, theta_hat)
  par_names <- c("rho", cnZ)
  names(par_hat) <- par_names
  
  Wy <- as.numeric(W %*% y)
  
  ll_full <- function(par) {
    rho   <- as.numeric(par[1])
    theta <- as.numeric(par[-1])
    
    if (!is.finite(rho) || rho <= -0.99 || rho >= 0.99) return(-1e12)
    
    e  <- as.numeric(y - rho * Wy - Z %*% theta)
    s2 <- mean(e^2)
    if (!is.finite(s2) || s2 <= 0) return(-1e12)
    
    ldS <- logdet_I_minus_rhoW(
      W, rho,
      method = fit_sdm$logdet_method,
      trW = fit_sdm$trW
    )
    if (!is.finite(ldS)) return(-1e12)
    
    ldS - (n / 2) * log(s2)
  }
  
  H <- try(numDeriv::hessian(ll_full, x = par_hat), silent = TRUE)
  if (inherits(H, "try-error")) {
    warning("Full HC1 Hessian failed for SDM. Returning fit unchanged.")
    return(fit_sdm)
  }
  
  A <- -H
  Ainv <- .safe_solve(A)
  
  e  <- as.numeric(y - rho_hat * Wy - Z %*% theta_hat)
  s2 <- mean(e^2)
  
  tr_SinvW <- .tr_SinvW_from_traces(rho_hat, fit_sdm$trW)
  if (!is.finite(tr_SinvW)) {
    ld_fun <- function(r) logdet_I_minus_rhoW(
      W, r,
      method = fit_sdm$logdet_method,
      trW = fit_sdm$trW
    )
    d_ld <- try(numDeriv::grad(ld_fun, x = rho_hat), silent = TRUE)
    if (inherits(d_ld, "try-error") || !is.finite(d_ld)) {
      warning("Could not evaluate full HC1 rho score for SDM. Returning fit unchanged.")
      return(fit_sdm)
    }
    tr_SinvW <- -as.numeric(d_ld)
  }
  
  S <- matrix(0, nrow = n, ncol = q)
  colnames(S) <- par_names
  
  S[, 1] <- (-tr_SinvW / n) + (Wy * e) / s2
  S[, -1] <- sweep(Z, 1, e / s2, FUN = "*")
  
  B <- crossprod(S)
  if (type == "HC1") B <- B * (n / (n - q))
  
  Vrob <- Ainv %*% B %*% Ainv
  Vrob <- (Vrob + t(Vrob)) / 2
  
  se_full <- sqrt(pmax(diag(Vrob), 0))
  names(se_full) <- par_names
  
  fit_sdm$hc_type_full <- type
  fit_sdm$vcov_par_hc_full <- Vrob
  fit_sdm$se_par_hc_full <- se_full
  fit_sdm$se_rho_hc_full <- as.numeric(se_full["rho"])
  fit_sdm$se_theta_hc_full <- as.numeric(se_full[-1])
  
  fit_sdm
}

# =========================
# Spatial block subsampling
# =========================
# The subsampling procedure is implemented on graph-connected spatial blocks.
# This is not an i.i.d. bootstrap. Instead, each subsample is grown over the
# symmetrized spatial graph to preserve local neighbourhood connectivity.
#
# Interpretation of the main tuning parameters:
# - K: number of graph-connected subsamples drawn from the retained support;
# - m: target number of retained cells in each subsample.
#
# The across-block standard deviation is later scaled by sqrt(m / n) before it
# is reported in the coefficient tables. This is the convention used in the
# present script when SBS is compared with model-based and HC1 standard errors.

.grow_block <- function(W_sym, start, m) {
  visited <- rep(FALSE, nrow(W_sym))
  q <- start
  visited[start] <- TRUE
  block <- integer(0)

  get_neighbors <- function(v) {
    p <- W_sym@p
    i <- W_sym@i
    if (p[v] == p[v + 1L]) integer(0) else (i[(p[v] + 1L):p[v + 1L]] + 1L)
  }

  while (length(block) < m && length(q) > 0L) {
    v <- q[1]
    q <- q[-1]
    block <- c(block, v)
    if (length(block) >= m) break
    nb <- get_neighbors(v)
    nb <- nb[!visited[nb]]
    visited[nb] <- TRUE
    q <- c(q, nb)
  }

  block
}

subsample_se <- function(y, X, W,
                         model = c("SAR", "SDM"),
                         K = 30,
                         m = 10000,
                         seed = 123,
                         logdet_method = "MC") {
  model <- match.arg(model)
  set.seed(seed)
  W_sym <- sign(W + t(W))
  est_list <- vector("list", K)
  block_sizes <- integer(K)

  for (k in seq_len(K)) {
    start <- sample.int(length(y), 1L)
    idx <- .grow_block(W_sym, start, m)
    block_sizes[k] <- length(idx)
    yb <- y[idx]
    Xb <- X[idx, , drop = FALSE]
    Wb <- W[idx, idx, drop = FALSE]
    trWb <- trace_powers_mc(Wb, kmax = 40, m = 20, seed = seed + k)

    if (model == "SAR") {
      fit <- sar_mle_manual(yb, Xb, Wb, logdet_method = logdet_method, trW = trWb, compute_impacts = FALSE)
      est_list[[k]] <- c(rho = fit$rho, coef = fit$beta)
    } else {
      fit <- sdm_mle_manual(yb, Xb, Wb, logdet_method = logdet_method, trW = trWb, compute_impacts = FALSE)
      est_list[[k]] <- c(rho = fit$rho, coef = fit$theta)
    }
  }

  E <- do.call(rbind, est_list)
  list(
    estimates = E,
    mean = colMeans(E),
    se_sub = apply(E, 2, sd),
    K = K,
    m_target = m,
    block_sizes = block_sizes,
    effective_m = mean(block_sizes),
    note = "SEs are across-block standard deviations from graph-connected spatial subsamples; scaling uses the average realized block size."
  )
}

.map_subse <- function(se_sub, coef_names, has_rho = TRUE) {
  nm <- names(se_sub)
  if (is.null(nm)) {
    idx <- seq_along(coef_names)
    return(setNames(se_sub[if (has_rho) 1 + idx else idx], coef_names))
  }
  idx <- match(paste0("coef", seq_along(coef_names)), nm)
  if (all(!is.na(idx))) return(setNames(se_sub[idx], coef_names))
  vec <- se_sub[setdiff(nm, "rho")]
  vec <- vec[seq_along(coef_names)]
  setNames(vec, coef_names)
}

.map_block_means <- function(M_est, coef_names, include_rho = TRUE) {
  nm <- colnames(M_est)
  out <- numeric(length(coef_names))
  names(out) <- coef_names

  if (is.null(nm)) {
    idx <- seq_along(coef_names)
    out[] <- colMeans(M_est[, if (include_rho) 1 + idx else idx, drop = FALSE], na.rm = TRUE)
    return(out)
  }

  idx <- match(paste0("coef", seq_along(coef_names)), nm)
  if (all(!is.na(idx))) {
    out[] <- colMeans(M_est[, idx, drop = FALSE], na.rm = TRUE)
    return(out)
  }

  cols <- setdiff(seq_len(ncol(M_est)), if (include_rho && "rho" %in% nm) match("rho", nm) else integer(0))
  cols <- cols[seq_along(coef_names)]
  out[] <- colMeans(M_est[, cols, drop = FALSE], na.rm = TRUE)
  out
}

scale_sub_stats <- function(ss_obj, n, m = NULL, coef_names, include_rho = TRUE) {
  m_eff <- if (!is.null(ss_obj$block_sizes)) mean(ss_obj$block_sizes, na.rm = TRUE) else m
  if (is.null(m_eff) || !is.finite(m_eff) || m_eff <= 0) {
    stop("A valid effective subsample size is required for SBS scaling.")
  }
  fac <- sqrt(m_eff / n)
  list(
    factor = fac,
    effective_m = m_eff,
    se = list(
      rho  = if (include_rho && "rho" %in% colnames(ss_obj$estimates)) as.numeric(ss_obj$se_sub["rho"] * fac) else NA_real_,
      coef = .map_subse(ss_obj$se_sub * fac, coef_names, has_rho = include_rho)
    ),
    mean = list(
      rho  = if (include_rho && "rho" %in% colnames(ss_obj$estimates)) mean(ss_obj$estimates[, "rho"], na.rm = TRUE) else NA_real_,
      coef = .map_block_means(ss_obj$estimates, coef_names, include_rho = include_rho)
    )
  )
}

# =========================
# Support and table helpers
# =========================

build_support_summary <- function(res, dataset_name = "Indonesia") {
  data.frame(
    dataset = dataset_name,
    mndwi_included = if (!is.null(res$mndwi_included)) res$mndwi_included else TRUE,
    total_template_cells = res$total_template_cells,
    common_support_cells = res$common_support_cells,
    retained_cells = res$retained_cells,
    dropped_isolates = res$dropped_isolates,
    share_common_support = res$share_common_support,
    share_retained_cells = res$share_retained_cells,
    stringsAsFactors = FALSE
  )
}

make_sar_table <- function(res, ss_sar = NULL, m_block = NULL, digits_p = 4) {
  beta <- res$fit_sar$beta
  se_mle <- res$fit_sar$se_beta
  se_hc1 <- if (!is.null(res$fit_sar$se_beta_hc)) res$fit_sar$se_beta_hc else rep(NA_real_, length(beta))
  cn <- colnames(res$X)

  out <- data.frame(
    term = cn,
    estimate = as.numeric(beta),
    se_mle = as.numeric(se_mle),
    z_mle = as.numeric(beta / se_mle),
    p_mle = .format_p(.p_from_z(beta / se_mle), digits = digits_p),
    se_hc1 = as.numeric(se_hc1),
    z_hc1 = as.numeric(beta / se_hc1),
    p_hc1 = .format_p(.p_from_z(beta / se_hc1), digits = digits_p),
    stringsAsFactors = FALSE
  )

  rho_row <- data.frame(
    term = "rho",
    estimate = as.numeric(res$fit_sar$rho),
    se_mle = as.numeric(res$fit_sar$se_rho),
    z_mle = as.numeric(res$fit_sar$rho / res$fit_sar$se_rho),
    p_mle = .format_p(.p_from_z(res$fit_sar$rho / res$fit_sar$se_rho), digits = digits_p),
    se_hc1 = NA_real_,
    z_hc1 = NA_real_,
    p_hc1 = NA_character_,
    stringsAsFactors = FALSE
  )

  if (!is.null(ss_sar) && !is.null(m_block)) {
    ss <- scale_sub_stats(ss_sar, n = length(res$y), coef_names = cn, include_rho = TRUE)
    out$mean_block <- as.numeric(ss$mean$coef[cn])
    out$delta_block <- out$estimate - out$mean_block
    out$se_sbs <- as.numeric(ss$se$coef[cn])
    out$z_sbs <- out$estimate / out$se_sbs
    out$p_sbs <- .format_p(.p_from_z(out$z_sbs), digits = digits_p)

    rho_row$mean_block <- as.numeric(ss$mean$rho)
    rho_row$delta_block <- rho_row$estimate - rho_row$mean_block
    rho_row$se_sbs <- as.numeric(ss$se$rho)
    rho_row$z_sbs <- rho_row$estimate / rho_row$se_sbs
    rho_row$p_sbs <- .format_p(.p_from_z(rho_row$z_sbs), digits = digits_p)
  }

  rbind(out, rho_row)
}

make_sdm_table <- function(res, ss_sdm = NULL, m_block = NULL, digits_p = 4) {
  theta <- res$fit_sdm$theta
  se_mle <- res$fit_sdm$se_theta
  se_hc1 <- if (!is.null(res$fit_sdm$se_theta_hc)) res$fit_sdm$se_theta_hc else rep(NA_real_, length(theta))
  names_theta <- res$fit_sdm$coef_names

  out <- data.frame(
    term = names_theta,
    estimate = as.numeric(theta),
    se_mle = as.numeric(se_mle),
    z_mle = as.numeric(theta / se_mle),
    p_mle = .format_p(.p_from_z(theta / se_mle), digits = digits_p),
    se_hc1 = as.numeric(se_hc1),
    z_hc1 = as.numeric(theta / se_hc1),
    p_hc1 = .format_p(.p_from_z(theta / se_hc1), digits = digits_p),
    stringsAsFactors = FALSE
  )

  rho_row <- data.frame(
    term = "rho",
    estimate = as.numeric(res$fit_sdm$rho),
    se_mle = as.numeric(res$fit_sdm$se_rho),
    z_mle = as.numeric(res$fit_sdm$rho / res$fit_sdm$se_rho),
    p_mle = .format_p(.p_from_z(res$fit_sdm$rho / res$fit_sdm$se_rho), digits = digits_p),
    se_hc1 = NA_real_,
    z_hc1 = NA_real_,
    p_hc1 = NA_character_,
    stringsAsFactors = FALSE
  )

  if (!is.null(ss_sdm) && !is.null(m_block)) {
    ss <- scale_sub_stats(ss_sdm, n = length(res$y), coef_names = names_theta, include_rho = TRUE)
    out$mean_block <- as.numeric(ss$mean$coef[names_theta])
    out$delta_block <- out$estimate - out$mean_block
    out$se_sbs <- as.numeric(ss$se$coef[names_theta])
    out$z_sbs <- out$estimate / out$se_sbs
    out$p_sbs <- .format_p(.p_from_z(out$z_sbs), digits = digits_p)

    rho_row$mean_block <- as.numeric(ss$mean$rho)
    rho_row$delta_block <- rho_row$estimate - rho_row$mean_block
    rho_row$se_sbs <- as.numeric(ss$se$rho)
    rho_row$z_sbs <- rho_row$estimate / rho_row$se_sbs
    rho_row$p_sbs <- .format_p(.p_from_z(rho_row$z_sbs), digits = digits_p)
  }

  rbind(out, rho_row)
}

extract_impacts_table <- function(impacts_obj) {
  if (inherits(impacts_obj, "LagImpact") || inherits(impacts_obj, "WXImpact")) {
    sm <- summary(impacts_obj, zstats = TRUE, short = TRUE)
    out <- as.data.frame(sm$res)
    out$variable <- rownames(out)
    rownames(out) <- NULL
    return(out[, c("variable", setdiff(names(out), "variable"))])
  }
  if (is.data.frame(impacts_obj)) return(impacts_obj)
  stop("Unsupported impacts object.")
}

# Full-HC1 table helpers
# ----------------------
# These helpers are used only when the optional full-parameter HC1 patch is
# requested. They preserve the existing workbook layout but replace the HC1
# columns with the full-parameter version so that rho also receives HC1-based
# inference.

make_sar_table_full_hc <- function(res, ss_sar = NULL, m_block = NULL, digits_p = 4) {
  beta <- res$fit_sar$beta
  se_mle <- res$fit_sar$se_beta
  se_hc1 <- if (!is.null(res$fit_sar$se_beta_hc_full)) res$fit_sar$se_beta_hc_full else rep(NA_real_, length(beta))
  cn <- colnames(res$X)
  
  out <- data.frame(
    term = cn,
    estimate = as.numeric(beta),
    se_mle = as.numeric(se_mle),
    z_mle = as.numeric(beta / se_mle),
    p_mle = .format_p(.p_from_z(beta / se_mle), digits = digits_p),
    se_hc1 = as.numeric(se_hc1),
    z_hc1 = as.numeric(beta / se_hc1),
    p_hc1 = .format_p(.p_from_z(beta / se_hc1), digits = digits_p),
    stringsAsFactors = FALSE
  )
  
  rho_hc <- if (!is.null(res$fit_sar$se_rho_hc_full)) as.numeric(res$fit_sar$se_rho_hc_full) else NA_real_
  rho_row <- data.frame(
    term = "rho",
    estimate = as.numeric(res$fit_sar$rho),
    se_mle = as.numeric(res$fit_sar$se_rho),
    z_mle = as.numeric(res$fit_sar$rho / res$fit_sar$se_rho),
    p_mle = .format_p(.p_from_z(res$fit_sar$rho / res$fit_sar$se_rho), digits = digits_p),
    se_hc1 = rho_hc,
    z_hc1 = as.numeric(res$fit_sar$rho / rho_hc),
    p_hc1 = .format_p(.p_from_z(res$fit_sar$rho / rho_hc), digits = digits_p),
    stringsAsFactors = FALSE
  )
  
  if (is.na(rho_hc)) {
    rho_row$z_hc1 <- NA_real_
    rho_row$p_hc1 <- NA_character_
  }
  
  if (!is.null(ss_sar) && !is.null(m_block)) {
    ss <- scale_sub_stats(ss_sar, n = length(res$y), m = m_block, coef_names = cn, include_rho = TRUE)
    out$mean_block <- as.numeric(ss$mean$coef[cn])
    out$delta_block <- out$estimate - out$mean_block
    out$se_sbs <- as.numeric(ss$se$coef[cn])
    out$z_sbs <- out$estimate / out$se_sbs
    out$p_sbs <- .format_p(.p_from_z(out$z_sbs), digits = digits_p)
    
    rho_row$mean_block <- as.numeric(ss$mean$rho)
    rho_row$delta_block <- rho_row$estimate - rho_row$mean_block
    rho_row$se_sbs <- as.numeric(ss$se$rho)
    rho_row$z_sbs <- rho_row$estimate / rho_row$se_sbs
    rho_row$p_sbs <- .format_p(.p_from_z(rho_row$z_sbs), digits = digits_p)
  }
  
  rbind(out, rho_row)
}

make_sdm_table_full_hc <- function(res, ss_sdm = NULL, m_block = NULL, digits_p = 4) {
  theta <- res$fit_sdm$theta
  se_mle <- res$fit_sdm$se_theta
  se_hc1 <- if (!is.null(res$fit_sdm$se_theta_hc_full)) res$fit_sdm$se_theta_hc_full else rep(NA_real_, length(theta))
  names_theta <- res$fit_sdm$coef_names
  
  out <- data.frame(
    term = names_theta,
    estimate = as.numeric(theta),
    se_mle = as.numeric(se_mle),
    z_mle = as.numeric(theta / se_mle),
    p_mle = .format_p(.p_from_z(theta / se_mle), digits = digits_p),
    se_hc1 = as.numeric(se_hc1),
    z_hc1 = as.numeric(theta / se_hc1),
    p_hc1 = .format_p(.p_from_z(theta / se_hc1), digits = digits_p),
    stringsAsFactors = FALSE
  )
  
  rho_hc <- if (!is.null(res$fit_sdm$se_rho_hc_full)) as.numeric(res$fit_sdm$se_rho_hc_full) else NA_real_
  rho_row <- data.frame(
    term = "rho",
    estimate = as.numeric(res$fit_sdm$rho),
    se_mle = as.numeric(res$fit_sdm$se_rho),
    z_mle = as.numeric(res$fit_sdm$rho / res$fit_sdm$se_rho),
    p_mle = .format_p(.p_from_z(res$fit_sdm$rho / res$fit_sdm$se_rho), digits = digits_p),
    se_hc1 = rho_hc,
    z_hc1 = as.numeric(res$fit_sdm$rho / rho_hc),
    p_hc1 = .format_p(.p_from_z(res$fit_sdm$rho / rho_hc), digits = digits_p),
    stringsAsFactors = FALSE
  )
  
  if (is.na(rho_hc)) {
    rho_row$z_hc1 <- NA_real_
    rho_row$p_hc1 <- NA_character_
  }
  
  if (!is.null(ss_sdm) && !is.null(m_block)) {
    ss <- scale_sub_stats(ss_sdm, n = length(res$y), m = m_block, coef_names = names_theta, include_rho = TRUE)
    out$mean_block <- as.numeric(ss$mean$coef[names_theta])
    out$delta_block <- out$estimate - out$mean_block
    out$se_sbs <- as.numeric(ss$se$coef[names_theta])
    out$z_sbs <- out$estimate / out$se_sbs
    out$p_sbs <- .format_p(.p_from_z(out$z_sbs), digits = digits_p)
    
    rho_row$mean_block <- as.numeric(ss$mean$rho)
    rho_row$delta_block <- rho_row$estimate - rho_row$mean_block
    rho_row$se_sbs <- as.numeric(ss$se$rho)
    rho_row$z_sbs <- rho_row$estimate / rho_row$se_sbs
    rho_row$p_sbs <- .format_p(.p_from_z(rho_row$z_sbs), digits = digits_p)
  }
  
  rbind(out, rho_row)
}

extract_impacts_table <- function(impacts_obj) {
  if (inherits(impacts_obj, "LagImpact") || inherits(impacts_obj, "WXImpact")) {
    sm <- summary(impacts_obj, zstats = TRUE, short = TRUE)
    out <- as.data.frame(sm$res)
    out$variable <- rownames(out)
    rownames(out) <- NULL
    return(out[, c("variable", setdiff(names(out), "variable"))])
  }
  if (is.data.frame(impacts_obj)) return(impacts_obj)
  stop("Unsupported impacts object.")
}

# =========================
# Export workbook
# =========================
# The workbook export is intended as a reviewer-facing summary of the full
# sparse-matrix workflow. It writes out:
# - support_size: template size, common support, retained cells, isolates;
# - sar_coefficients / sdm_coefficients: MLE, HC1, and, if requested, SBS;
# - sar_impacts / sdm_impacts: average direct, indirect, and total effects;
# - model_fit: sample size and information criteria;
# - sbs_draws_*: optional block-level estimates from the SBS procedure.
#
# Key arguments:
# - sbs_K: number of graph-connected subsamples used in SBS inference;
# - sbs_m: target connected block size used in each subsample;
# - run_sbs: whether SBS-based uncertainty measures are computed and exported.
# - hc1_mode: "conditional" reproduces the manuscript-oriented HC1 treatment,
#   "full" reports the optional full-parameter HC1 including rho, and "both"
#   exports both versions side by side in separate sheets.

export_sparse_outputs <- function(res,
                                  dataset_name = "Indonesia",
                                  out_dir = out_path,
                                  sbs_K = 30,
                                  sbs_m = 500,
                                  seed = 42,
                                  run_sbs = TRUE,
                                  export_block_draws = TRUE,
                                  hc1_mode = c("conditional", "full", "both")) {
  
  hc1_mode <- match.arg(hc1_mode)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # HC1 inference
  if (hc1_mode %in% c("conditional", "both")) {
    res$fit_sar <- sar_add_sandwich_se(res$fit_sar, y = res$y, X = res$X, W = res$W_common, type = "HC1")
    res$fit_sdm <- sdm_add_sandwich_se(res$fit_sdm, y = res$y, X = res$X, W = res$W_common, type = "HC1")
  }
  
  if (hc1_mode %in% c("full", "both")) {
    res$fit_sar <- sar_add_sandwich_se_full(res$fit_sar, y = res$y, X = res$X, W = res$W_common, type = "HC1")
    res$fit_sdm <- sdm_add_sandwich_se_full(res$fit_sdm, y = res$y, X = res$X, W = res$W_common, type = "HC1")
  }
  
  # SBS inference + elapsed time
  sbs_sar <- NULL
  sbs_sdm <- NULL
  elapsed_sbs_sar_sec <- NA_real_
  elapsed_sbs_sdm_sec <- NA_real_
  
  if (isTRUE(run_sbs)) {
    t0_sbs_sar <- Sys.time()
    sbs_sar <- subsample_se(
      res$y, res$X, res$W_common,
      model = "SAR", K = sbs_K, m = sbs_m,
      seed = seed, logdet_method = "MC"
    )
    t1_sbs_sar <- Sys.time()
    elapsed_sbs_sar_sec <- as.numeric(difftime(t1_sbs_sar, t0_sbs_sar, units = "secs"))
    
    t0_sbs_sdm <- Sys.time()
    sbs_sdm <- subsample_se(
      res$y, res$X, res$W_common,
      model = "SDM", K = sbs_K, m = sbs_m,
      seed = seed, logdet_method = "MC"
    )
    t1_sbs_sdm <- Sys.time()
    elapsed_sbs_sdm_sec <- as.numeric(difftime(t1_sbs_sdm, t0_sbs_sdm, units = "secs"))
  }
  
  # Tables
  support_df <- build_support_summary(res, dataset_name = dataset_name)
  
  if (hc1_mode == "full") {
    sar_coef_df <- make_sar_table_full_hc(res, ss_sar = sbs_sar, m_block = if (run_sbs) sbs_m else NULL)
    sdm_coef_df <- make_sdm_table_full_hc(res, ss_sdm = sbs_sdm, m_block = if (run_sbs) sbs_m else NULL)
    sar_coef_df_cond <- NULL
    sdm_coef_df_cond <- NULL
  } else {
    sar_coef_df <- make_sar_table(res, ss_sar = sbs_sar, m_block = if (run_sbs) sbs_m else NULL)
    sdm_coef_df <- make_sdm_table(res, ss_sdm = sbs_sdm, m_block = if (run_sbs) sbs_m else NULL)
    sar_coef_df_cond <- if (hc1_mode == "both") sar_coef_df else NULL
    sdm_coef_df_cond <- if (hc1_mode == "both") sdm_coef_df else NULL
    
    if (hc1_mode == "both") {
      sar_coef_df <- make_sar_table_full_hc(res, ss_sar = sbs_sar, m_block = if (run_sbs) sbs_m else NULL)
      sdm_coef_df <- make_sdm_table_full_hc(res, ss_sdm = sbs_sdm, m_block = if (run_sbs) sbs_m else NULL)
    }
  }
  
  sar_imp_df <- extract_impacts_table(res$impacts_sar)
  sdm_imp_df <- extract_impacts_table(res$impacts_sdm)
  
  # Effective SBS block size for metadata, if available
  effective_m_sar <- if (!is.null(sbs_sar) && !is.null(sbs_sar$block_sizes)) mean(sbs_sar$block_sizes, na.rm = TRUE) else NA_real_
  effective_m_sdm <- if (!is.null(sbs_sdm) && !is.null(sbs_sdm$block_sizes)) mean(sbs_sdm$block_sizes, na.rm = TRUE) else NA_real_
  effective_m <- mean(c(effective_m_sar, effective_m_sdm), na.rm = TRUE)
  if (!is.finite(effective_m)) effective_m <- NA_real_
  
  # Model fit + elapsed time
  gof_df <- data.frame(
    model = c("SAR", "SDM"),
    n = c(res$fit_sar$n, res$fit_sdm$n),
    logLik = c(res$fit_sar$logLik, res$fit_sdm$logLik),
    AIC = c(res$fit_sar$AIC, res$fit_sdm$AIC),
    BIC = c(res$fit_sar$BIC, res$fit_sdm$BIC),
    elapsed_estimation_sec = c(res$elapsed_sar_sec, res$elapsed_sdm_sec),
    elapsed_impacts_sec = c(res$elapsed_impacts_sar_sec, res$elapsed_impacts_sdm_sec),
    elapsed_total_modeling_sec = c(
      sum(c(res$elapsed_sar_sec, res$elapsed_impacts_sar_sec), na.rm = TRUE),
      sum(c(res$elapsed_sdm_sec, res$elapsed_impacts_sdm_sec), na.rm = TRUE)
    ),
    elapsed_sbs_sec = c(elapsed_sbs_sar_sec, elapsed_sbs_sdm_sec),
    stringsAsFactors = FALSE
  )
  
  meta_df <- data.frame(
    item = c(
      "dataset_name",
      "engine",
      "W_type",
      "seed",
      "HC1_mode",
      "SBS_blocks_K",
      "SBS_block_size_m_nominal",
      "SBS_block_size_m_effective",
      "SBS_scale_factor_effective",
      "HC1_full_note",
      "Timing_note_1",
      "Timing_note_2",
      "SBS_timing_note"
    ),
    value = c(
      dataset_name,
      res$engine,
      res$type_W,
      res$seed,
      hc1_mode,
      if (run_sbs) sbs_K else NA,
      if (run_sbs) sbs_m else NA,
      if (run_sbs) effective_m else NA,
      if (run_sbs && is.finite(effective_m)) sqrt(effective_m / length(res$y)) else NA,
      if (hc1_mode %in% c("full", "both")) {
        "Full HC1 includes rho via optional numerical sandwich approximation."
      } else {
        "HC1 is conditional on rho-hat and applies to slopes only."
      },
      "elapsed_estimation_sec records model-fitting time for SAR and SDM separately.",
      "Under the manual engine, impacts are computed inside the estimator, so elapsed_impacts_sec may be NA and estimation time may already include impact calculation.",
      "elapsed_sbs_sec records the wall-clock time required to run graph-connected spatial block subsampling for each model."
    ),
    stringsAsFactors = FALSE
  )
  
  # Round for cleaner reviewer-facing workbook
  support_df <- .round_numeric_df(support_df, 6)
  if (!is.null(sar_coef_df_cond)) sar_coef_df_cond <- .round_numeric_df(sar_coef_df_cond, 6)
  if (!is.null(sdm_coef_df_cond)) sdm_coef_df_cond <- .round_numeric_df(sdm_coef_df_cond, 6)
  sar_coef_df <- .round_numeric_df(sar_coef_df, 6)
  sdm_coef_df <- .round_numeric_df(sdm_coef_df, 6)
  sar_imp_df <- .round_numeric_df(sar_imp_df, 6)
  sdm_imp_df <- .round_numeric_df(sdm_imp_df, 6)
  gof_df <- .round_numeric_df(gof_df, 6)
  
  wb <- createWorkbook()
  addWorksheet(wb, "meta")
  addWorksheet(wb, "support_size")
  
  if (hc1_mode == "both") {
    addWorksheet(wb, "sar_coeff_cond")
    addWorksheet(wb, "sdm_coeff_cond")
    addWorksheet(wb, "sar_coeff_full")
    addWorksheet(wb, "sdm_coeff_full")
  } else {
    addWorksheet(wb, "sar_coefficients")
    addWorksheet(wb, "sdm_coefficients")
  }
  
  addWorksheet(wb, "sar_impacts")
  addWorksheet(wb, "sdm_impacts")
  addWorksheet(wb, "model_fit")
  
  writeData(wb, "meta", meta_df)
  writeData(wb, "support_size", support_df)
  
  if (hc1_mode == "both") {
    writeData(wb, "sar_coeff_cond", sar_coef_df_cond)
    writeData(wb, "sdm_coeff_cond", sdm_coef_df_cond)
    writeData(wb, "sar_coeff_full", sar_coef_df)
    writeData(wb, "sdm_coeff_full", sdm_coef_df)
  } else {
    writeData(wb, "sar_coefficients", sar_coef_df)
    writeData(wb, "sdm_coefficients", sdm_coef_df)
  }
  
  writeData(wb, "sar_impacts", sar_imp_df)
  writeData(wb, "sdm_impacts", sdm_imp_df)
  writeData(wb, "model_fit", gof_df)
  
  if (isTRUE(run_sbs) && isTRUE(export_block_draws)) {
    if (!is.null(sbs_sar) && !is.null(sbs_sar$estimates)) {
      addWorksheet(wb, "sbs_draws_sar")
      writeData(wb, "sbs_draws_sar", .round_numeric_df(as.data.frame(sbs_sar$estimates), 6))
    }
    if (!is.null(sbs_sdm) && !is.null(sbs_sdm$estimates)) {
      addWorksheet(wb, "sbs_draws_sdm")
      writeData(wb, "sbs_draws_sdm", .round_numeric_df(as.data.frame(sbs_sdm$estimates), 6))
    }
  }
  
  out_file <- file.path(out_dir, paste0(dataset_name, "_SAR_SDM_outputs.xlsx"))
  saveWorkbook(wb, out_file, overwrite = TRUE)
  
  invisible(list(
    workbook = out_file,
    support_size = support_df,
    model_fit = gof_df,
    sar_coefficients = sar_coef_df,
    sdm_coefficients = sdm_coef_df,
    sar_impacts = sar_imp_df,
    sdm_impacts = sdm_imp_df,
    sbs_sar = sbs_sar,
    sbs_sdm = sbs_sdm,
    elapsed_sbs_sar_sec = elapsed_sbs_sar_sec,
    elapsed_sbs_sdm_sec = elapsed_sbs_sdm_sec,
    elapsed_sbs_total_sec = sum(c(elapsed_sbs_sar_sec, elapsed_sbs_sdm_sec), na.rm = TRUE)
  ))
}

# =========================
# Full preprocessing pipeline for the 2020 application
# =========================
# The sections below integrate the preprocessing steps that were previously
# stored in separate scripts. They are kept in the same file so that the full
# workflow can be run from raw shapefile/raster inputs to the exported
# tabulations in `out_path`.
#
# This consolidated version is aligned with the current four-layer design used
# in the manuscript revision:
#   - GDP (dependent variable)
#   - LandScan population
#   - NDVI
#   - MNDWI
#
# The older SAVI and CO preprocessing blocks are not included in the default
# run here because the current reviewer-facing regression specification and the
# revised data description focus on the four-layer workflow above.

# -------------------------
# User configuration
# -------------------------
cfg <- list(
  shapefile_admin2 = "D:/SHAPEFILE/gadm41_IDN_2.shp",
  data_root = "E:/Spatial Sparse W Matrix Research/DATA",
  landscan_root = "F:/LANDSCAN DATA",
  gdp_file = "E:/Spatial Sparse W Matrix Research/DATA/rast_gdp_tot_1990_2020_30arcsec.tif",
  out_dir = out_path
)

analysis_year <- 2020

# Control flags for the full pipeline
RUN_AUTOCORR_SUMMARY <- TRUE
RUN_REGRESSION_SUITE <- TRUE
RUN_SINGLE_REGION_EXAMPLE <- FALSE

# Regions to process in the full run
regions_to_run <- c(
  "Sumatra", "Java", "Kalimantan", "Sulawesi",
  "Papua", "Bali_Nusa_Tenggara", "Maluku", "Indonesia"
)

# -------------------------
# Polygon helpers
# -------------------------
build_region_polygons <- function(shp_admin2) {
  if (!("NAME_1" %in% names(shp_admin2))) {
    stop("The shapefile must contain a 'NAME_1' field for province names.")
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
  missing_by_group <- missing_by_group[lengths(missing_by_group) > 0]

  if (length(missing_by_group) > 0) {
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
    pnames <- province_groups[[region_nm]]
    shp_admin2[shp_admin2$NAME_1 %in% pnames, ]
  })
  names(polys) <- names(province_groups)

  empty_regions <- names(polys)[vapply(polys, nrow, integer(1)) == 0L]
  if (length(empty_regions) > 0) {
    stop(sprintf("No polygons were found for the following regions: %s", paste(empty_regions, collapse = ", ")))
  }

  # Indonesia is the combined national polygon used only for convenience in
  # descriptive summaries; island-level merges are still used where source
  # rasters are stored by island group.
  polys$Indonesia <- shp_admin2

  polys
}

# Source tags used in the existing raster filenames
region_source_tag <- c(
  Sumatra = "sumatra",
  Java = "Java_Bali",
  Kalimantan = "kalimantan",
  Sulawesi = "sulawesi",
  Papua = "Eastern_Indonesia",
  Bali_Nusa_Tenggara = "Bali_NTB_NTT",
  Maluku = "Eastern_Indonesia"
)

`%||%` <- function(x, y) if (is.null(x)) y else x

crop_mask <- function(x, poly) {
  terra::mask(terra::crop(x, vect(poly)), vect(poly))
}

merge_rasters_safe <- function(r_list) {
  stopifnot(length(r_list) >= 1)
  if (length(r_list) == 1) return(r_list[[1]])
  Reduce(function(a, b) terra::merge(a, b), r_list)
}

# -------------------------
# Raster loading helpers
# -------------------------
prepare_data_sources <- function(cfg, year = 2020) {
  yr <- as.character(year)
  list(
    gdp = terra::rast(cfg$gdp_file),
    landscan = setNames(list(terra::rast(file.path(cfg$landscan_root, yr, paste0("landscan-global-", yr, ".tif")))), yr)
  )
}

load_single_region_inputs <- function(region_name, year, polys, sources, cfg) {
  if (identical(region_name, "Indonesia")) {
    stop("load_single_region_inputs() is for non-Indonesia island groups only.")
  }

  poly <- polys[[region_name]]
  tag  <- region_source_tag[[region_name]]
  yr   <- as.character(year)

  # GDP and LandScan come from national rasters, then are cropped and masked to
  # the island polygon.
  gdp_layer_name <- paste0("gdp_", year)
  if (!(gdp_layer_name %in% names(sources$gdp))) {
    stop(sprintf("Layer '%s' was not found in %s.", gdp_layer_name, cfg$gdp_file))
  }

  rast_gdp <- crop_mask(sources$gdp[[gdp_layer_name]], poly)
  rast_pop <- crop_mask(sources$landscan[[yr]], poly)

  # NDVI and MNDWI are stored as island-group-specific files in the user's
  # current directory structure. Java uses the Java-Bali source file and is
  # then masked back to the Java polygon, which matches the earlier scripts.
  ndvi_file <- file.path(
    cfg$data_root, "Indonesia_NDVI_MODIS",
    paste0("Annual_Index_", tag, "_", year, ".tif")
  )
  mndwi_file <- file.path(
    cfg$data_root, "Indonesia_MNDWI",
    paste0("MNDWI_", tag, "_", year, ".tif")
  )

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

prepare_analysis_rasters <- function(year, polys, cfg) {
  sources <- prepare_data_sources(cfg, year = year)

  core_regions <- c(
    "Sumatra", "Java", "Kalimantan", "Sulawesi",
    "Papua", "Bali_Nusa_Tenggara", "Maluku"
  )

  region_inputs <- setNames(vector("list", length(core_regions)), core_regions)
  for (rg in core_regions) {
    region_inputs[[rg]] <- load_single_region_inputs(rg, year, polys, sources, cfg)
  }

  # Indonesia is constructed by merging the already-masked island-group rasters.
  region_inputs$Indonesia <- list(
    template = merge_rasters_safe(lapply(core_regions, function(rg) region_inputs[[rg]]$gdp)),
    gdp      = merge_rasters_safe(lapply(core_regions, function(rg) region_inputs[[rg]]$gdp)),
    pop      = merge_rasters_safe(lapply(core_regions, function(rg) region_inputs[[rg]]$pop)),
    ndvi     = merge_rasters_safe(lapply(core_regions, function(rg) region_inputs[[rg]]$ndvi)),
    mndwi    = merge_rasters_safe(lapply(core_regions, function(rg) region_inputs[[rg]]$mndwi))
  )

  region_inputs
}

# -------------------------
# Descriptive support summary
# -------------------------
build_region_dimension_table <- function(region_inputs) {
  out <- lapply(names(region_inputs), function(region_name) {
    rr <- region_inputs[[region_name]]
    
    template <- rr$template
    r <- terra::nrow(template)
    c <- terra::ncol(template)
    N <- terra::ncell(template)
    
    gdp_vec   <- raster_to_colmajor(rr$gdp,   template, method = "near", var_name = paste0(region_name, "_gdp"))
    pop_vec   <- raster_to_colmajor(rr$pop,   template, method = "near", var_name = paste0(region_name, "_pop"))
    ndvi_vec  <- raster_to_colmajor(rr$ndvi,  template, method = "near", var_name = paste0(region_name, "_ndvi"))
    mndwi_vec <- raster_to_colmajor(rr$mndwi, template, method = "near", var_name = paste0(region_name, "_mndwi"))
    
    valid_gdp   <- !is.na(gdp_vec)
    valid_pop   <- !is.na(pop_vec)
    valid_ndvi  <- !is.na(ndvi_vec)
    valid_mndwi <- !is.na(mndwi_vec)
    
    common_support <- valid_gdp & valid_pop & valid_ndvi & valid_mndwi
    
    data.frame(
      region = region_name,
      rows = r,
      cols = c,
      total_grids = N,
      non_na_gdp = sum(valid_gdp),
      non_na_pop = sum(valid_pop),
      non_na_ndvi = sum(valid_ndvi),
      non_na_mndwi = sum(valid_mndwi),
      common_support_cells = sum(common_support)
    )
  })
  
  do.call(rbind, out)
}

# -------------------------
# Sparse autocorrelation summary
# -------------------------
# This section reproduces, in a generic form, the sparse W-tilde construction
# used earlier for the variable-specific spatial autocorrelation checks. It is
# written as a reusable helper so that the same procedure can be applied to each
# region and each variable without duplicating large code blocks.
compute_sparse_autocorr <- function(rast_obj,
                                    region_name,
                                    variable_name,
                                    out_dir,
                                    type_W = "queen",
                                    save_w_tilde = TRUE) {
  r <- nrow(rast_obj)
  c <- ncol(rast_obj)
  N <- r * c
  
  start_sparse <- Sys.time()
  W_full <- build_W_full(r, c, type = type_W)
  end_sparse <- Sys.time()
  
  if (terra::nlyr(rast_obj) != 1L) {
    stop(sprintf("'%s' must be a single-layer SpatRaster.", variable_name))
  }
  
  vals <- terra::values(rast_obj, mat = FALSE)
  
  if (length(vals) != N) {
    stop(sprintf(
      "Value-length mismatch for %s in %s: expected %d cells but found %d values.",
      variable_name, region_name, N, length(vals)
    ))
  }
  
  M <- matrix(vals, nrow = r, ncol = c, byrow = TRUE)
  
  # Flatten using the same column-major convention used in the sparse workflow
  vals_col <- as.numeric(t(M))
  
  idx_obs <- which(!is.na(vals_col))
  V <- length(idx_obs)
  
  if (V == 0L) {
    return(data.frame(
      region = region_name,
      variable = variable_name,
      rows = r,
      cols = c,
      total_grids = N,
      non_na = 0L,
      retained_non_isolates = 0L,
      moran_I = NA_real_,
      moran_p = NA_real_,
      geary_C = NA_real_,
      geary_p = NA_real_,
      sparse_build_sec = as.numeric(difftime(end_sparse, start_sparse, units = "secs")),
      adjusted_build_sec = NA_real_,
      w_tilde_mb = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  start_adjust <- Sys.time()
  P <- sparseMatrix(i = idx_obs, j = seq_len(V), x = 1, dims = c(N, V))
  W_tilde <- t(P) %*% W_full %*% P
  
  # Drop isolates before autocorrelation testing
  deg <- Matrix::rowSums(W_tilde != 0)
  keep <- deg > 0
  
  if (!any(keep)) {
    end_adjust <- Sys.time()
    return(data.frame(
      region = region_name,
      variable = variable_name,
      rows = r,
      cols = c,
      total_grids = N,
      non_na = V,
      retained_non_isolates = 0L,
      moran_I = NA_real_,
      moran_p = NA_real_,
      geary_C = NA_real_,
      geary_p = NA_real_,
      sparse_build_sec = as.numeric(difftime(end_sparse, start_sparse, units = "secs")),
      adjusted_build_sec = as.numeric(difftime(end_adjust, start_adjust, units = "secs")),
      w_tilde_mb = as.numeric(object.size(W_tilde)) / 1024^2,
      stringsAsFactors = FALSE
    ))
  }
  
  W_tilde <- W_tilde[keep, keep, drop = FALSE]
  z <- vals_col[idx_obs][keep]
  end_adjust <- Sys.time()
  
  lw <- dgC_to_listw(as(W_tilde, "dgCMatrix"))
  
  mor <- try(moran.test(z, lw, zero.policy = TRUE), silent = TRUE)
  gea <- try(geary.test(z, lw, zero.policy = TRUE), silent = TRUE)
  
  moran_I <- if (inherits(mor, "try-error")) NA_real_ else unname(mor$estimate[["Moran I statistic"]])
  moran_p <- if (inherits(mor, "try-error")) NA_real_ else mor$p.value
  geary_C <- if (inherits(gea, "try-error")) NA_real_ else unname(gea$estimate[["Geary C statistic"]])
  geary_p <- if (inherits(gea, "try-error")) NA_real_ else gea$p.value
  
  if (isTRUE(save_w_tilde)) {
    rds_dir <- file.path(out_dir, "rds_w_tilde")
    dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      W_tilde,
      file = file.path(
        rds_dir,
        paste0(tolower(region_name), "_", tolower(variable_name), "_W_tilde.rds")
      )
    )
  }
  
  data.frame(
    region = region_name,
    variable = variable_name,
    rows = r,
    cols = c,
    total_grids = N,
    non_na = V,
    retained_non_isolates = sum(keep),
    moran_I = moran_I,
    moran_p = moran_p,
    geary_C = geary_C,
    geary_p = geary_p,
    sparse_build_sec = as.numeric(difftime(end_sparse, start_sparse, units = "secs")),
    adjusted_build_sec = as.numeric(difftime(end_adjust, start_adjust, units = "secs")),
    w_tilde_mb = as.numeric(object.size(W_tilde)) / 1024^2,
    stringsAsFactors = FALSE
  )
}

run_autocorr_suite <- function(region_inputs,
                               regions = names(region_inputs),
                               out_dir,
                               type_W = "queen",
                               save_w_tilde = TRUE) {
  all_rows <- list()
  k <- 1L

  var_map <- list(
    Grid_GDP = "gdp",
    Pop_dis = "pop",
    Veg_quality = "ndvi",
    Water_quality = "mndwi"
  )

  for (rg in regions) {
    obj <- region_inputs[[rg]]
    for (vn in names(var_map)) {
      all_rows[[k]] <- compute_sparse_autocorr(
        rast_obj = obj[[var_map[[vn]]]],
        region_name = rg,
        variable_name = vn,
        out_dir = out_dir,
        type_W = type_W,
        save_w_tilde = save_w_tilde
      )
      k <- k + 1L
    }
  }

  summary_df <- do.call(rbind, all_rows)
  summary_df <- .round_numeric_df(summary_df, 6)

  wb <- createWorkbook()
  addWorksheet(wb, "autocorr_summary")
  writeData(wb, "autocorr_summary", summary_df)

  out_file <- file.path(out_dir, paste0("Autocorrelation_Summary_", analysis_year, ".xlsx"))
  saveWorkbook(wb, out_file, overwrite = TRUE)

  invisible(list(summary = summary_df, workbook = out_file))
}

# -------------------------
# Regression suite wrapper
# -------------------------
run_regression_suite <- function(region_inputs,
                                 regions = names(region_inputs),
                                 out_dir,
                                 sbs_K = 30,
                                 sbs_m = 500,
                                 logdet_method = "MC",
                                 type_W = "queen",
                                 engine = "manual",
                                 seed = 42,
                                 run_sbs = TRUE,
                                 hc1_mode = c("conditional", "full", "both"),
                                 analysis_year = 2020) {
  
  hc1_mode <- match.arg(hc1_mode)
  
  out_list <- setNames(vector("list", length(regions)), regions)
  support_rows <- list()
  fit_rows <- list()
  timing_rows <- list()
  
  for (rg in regions) {
    obj <- region_inputs[[rg]]
    
    out_list[[rg]] <- run_and_export_region(
      dataset_name = rg,
      template = obj$template,
      rast_gdp = obj$gdp,
      rast_pop = obj$pop,
      rast_ndvi = obj$ndvi,
      rast_mndwi = obj$mndwi,
      out_dir = out_dir,
      sbs_K = sbs_K,
      sbs_m = sbs_m,
      logdet_method = logdet_method,
      type_W = type_W,
      engine = engine,
      seed = seed,
      run_sbs = run_sbs,
      hc1_mode = hc1_mode
    )
    
    res <- out_list[[rg]]$res
    exp_obj <- out_list[[rg]]$export_obj
    
    support_rows[[rg]] <- data.frame(
      region = rg,
      total_template_cells = res$total_template_cells,
      common_support_cells = res$common_support_cells,
      retained_cells = res$retained_cells,
      dropped_isolates = res$dropped_isolates,
      share_common_support = res$share_common_support,
      share_retained_cells = res$share_retained_cells,
      stringsAsFactors = FALSE
    )
    
    fit_rows[[paste0(rg, "_SAR")]] <- data.frame(
      region = rg,
      model = "SAR",
      logLik = res$fit_sar$logLik,
      AIC = res$fit_sar$AIC,
      BIC = res$fit_sar$BIC,
      rho = res$fit_sar$rho,
      sigma2 = res$fit_sar$sigma2,
      elapsed_estimation_sec = res$elapsed_sar_sec,
      elapsed_impacts_sec = res$elapsed_impacts_sar_sec,
      elapsed_total_modeling_sec = sum(
        c(res$elapsed_sar_sec, res$elapsed_impacts_sar_sec),
        na.rm = TRUE
      ),
      elapsed_sbs_sec = exp_obj$elapsed_sbs_sar_sec,
      stringsAsFactors = FALSE
    )
    
    fit_rows[[paste0(rg, "_SDM")]] <- data.frame(
      region = rg,
      model = "SDM",
      logLik = res$fit_sdm$logLik,
      AIC = res$fit_sdm$AIC,
      BIC = res$fit_sdm$BIC,
      rho = res$fit_sdm$rho,
      sigma2 = res$fit_sdm$sigma2,
      elapsed_estimation_sec = res$elapsed_sdm_sec,
      elapsed_impacts_sec = res$elapsed_impacts_sdm_sec,
      elapsed_total_modeling_sec = sum(
        c(res$elapsed_sdm_sec, res$elapsed_impacts_sdm_sec),
        na.rm = TRUE
      ),
      elapsed_sbs_sec = exp_obj$elapsed_sbs_sdm_sec,
      stringsAsFactors = FALSE
    )
    
    timing_rows[[paste0(rg, "_SAR")]] <- data.frame(
      region = rg,
      model = "SAR",
      total_template_cells = res$total_template_cells,
      common_support_cells = res$common_support_cells,
      retained_cells = res$retained_cells,
      elapsed_estimation_sec = res$elapsed_sar_sec,
      elapsed_impacts_sec = res$elapsed_impacts_sar_sec,
      elapsed_total_modeling_sec = sum(
        c(res$elapsed_sar_sec, res$elapsed_impacts_sar_sec),
        na.rm = TRUE
      ),
      retained_cells_per_sec = if (is.finite(res$elapsed_sar_sec) && res$elapsed_sar_sec > 0) {
        res$retained_cells / res$elapsed_sar_sec
      } else {
        NA_real_
      },
      elapsed_sbs_sec = exp_obj$elapsed_sbs_sar_sec,
      stringsAsFactors = FALSE
    )
    
    timing_rows[[paste0(rg, "_SDM")]] <- data.frame(
      region = rg,
      model = "SDM",
      total_template_cells = res$total_template_cells,
      common_support_cells = res$common_support_cells,
      retained_cells = res$retained_cells,
      elapsed_estimation_sec = res$elapsed_sdm_sec,
      elapsed_impacts_sec = res$elapsed_impacts_sdm_sec,
      elapsed_total_modeling_sec = sum(
        c(res$elapsed_sdm_sec, res$elapsed_impacts_sdm_sec),
        na.rm = TRUE
      ),
      retained_cells_per_sec = if (is.finite(res$elapsed_sdm_sec) && res$elapsed_sdm_sec > 0) {
        res$retained_cells / res$elapsed_sdm_sec
      } else {
        NA_real_
      },
      elapsed_sbs_sec = exp_obj$elapsed_sbs_sdm_sec,
      stringsAsFactors = FALSE
    )
  }
  
  support_df <- .round_numeric_df(do.call(rbind, support_rows), 6)
  fit_df <- .round_numeric_df(do.call(rbind, fit_rows), 6)
  timing_df <- .round_numeric_df(do.call(rbind, timing_rows), 6)
  
  wb <- createWorkbook()
  addWorksheet(wb, "support_summary")
  addWorksheet(wb, "model_fit_summary")
  addWorksheet(wb, "timing_summary")
  addWorksheet(wb, "timing_note")
  
  writeData(wb, "support_summary", support_df)
  writeData(wb, "model_fit_summary", fit_df)
  writeData(wb, "timing_summary", timing_df)
  
  timing_note_df <- data.frame(
    item = c(
      "elapsed_estimation_sec",
      "elapsed_impacts_sec",
      "elapsed_total_modeling_sec",
      "SBS_timing_note",
      "retained_cells_per_sec",
      "manual_engine_note"
    ),
    description = c(
      "Elapsed wall-clock time for fitting the SAR or SDM model.",
      "Elapsed wall-clock time for post-estimation impact calculations. This may be NA under the manual engine.",
      "Sum of estimation and impact-calculation time for the reported model.",
      "elapsed_sbs_sec records the wall-clock time required to run graph-connected spatial block subsampling for the reported model.",
      "Retained common-support cells divided by elapsed_estimation_sec. This is a descriptive throughput indicator.",
      "Under the manual engine, impact calculations are embedded in the estimation routine. Therefore, elapsed_estimation_sec may already include impact computation, while elapsed_impacts_sec can remain NA."
    ),
    stringsAsFactors = FALSE
  )
  writeData(wb, "timing_note", timing_note_df)
  
  out_file <- file.path(out_dir, paste0("Regression_Suite_Summary_", analysis_year, ".xlsx"))
  saveWorkbook(wb, out_file, overwrite = TRUE)
  
  invisible(list(
    results = out_list,
    support_summary = support_df,
    model_fit_summary = fit_df,
    timing_summary = timing_df,
    workbook = out_file
  ))
}

# =========================
# Main execution block
# =========================
# The block below can be run as-is once the paths above match the local machine.
# It will:
# 1. read the shapefile,
# 2. build island-group polygons,
# 3. prepare the four analysis rasters for each region,
# 4. export a table of grid dimensions/common support,
# 5. optionally export sparse autocorrelation summaries, and
# 6. run SAR/SDM estimation plus workbook export for each region.

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

ina_shp_lv2 <- sf::st_read(cfg$shapefile_admin2, quiet = TRUE)
region_polys <- build_region_polygons(ina_shp_lv2)
region_inputs_2020 <- prepare_analysis_rasters(
  year = analysis_year,
  polys = region_polys,
  cfg = cfg
)

# Table 2-style dimensional summary based on the final four-layer common support
dimension_table_2020 <- build_region_dimension_table(region_inputs_2020)
wb_dims <- createWorkbook()
addWorksheet(wb_dims, "dimensions_common_support")
writeData(wb_dims, "dimensions_common_support", .round_numeric_df(dimension_table_2020, 6))
saveWorkbook(
  wb_dims,
  file.path(cfg$out_dir, paste0("Region_Dimensions_CommonSupport_", analysis_year, ".xlsx")),
  overwrite = TRUE
)

if (isTRUE(RUN_AUTOCORR_SUMMARY)) {
  autocorr_outputs <- run_autocorr_suite(
    region_inputs = region_inputs_2020,
    regions = regions_to_run,
    out_dir = cfg$out_dir,
    type_W = "queen",
    save_w_tilde = TRUE
  )
}

if (isTRUE(RUN_REGRESSION_SUITE)) {
  regression_outputs <- run_regression_suite(
    region_inputs = region_inputs_2020,
    regions = regions_to_run,
    out_dir = cfg$out_dir,
    sbs_K = 30,
    sbs_m = 500,
    logdet_method = "MC",
    type_W = "queen",
    engine = "manual",
    seed = 42,
    run_sbs = TRUE,
    hc1_mode = "both",
    analysis_year = 2020
  )
}

# -------------------------
# Single-region example calls
# -------------------------
# These examples are optional and are disabled by default to avoid rerunning
# Indonesia after the full regression suite has already written the same output.
if (isTRUE(RUN_SINGLE_REGION_EXAMPLE)) {
out_idn <- run_and_export_region(
  dataset_name = "Indonesia",
  
  # The template defines the full rectangular grid on which the initial
  # lattice-based spatial weights matrix is constructed.
  template  = region_inputs_2020$Indonesia$template,
  
  # Gridded GDP is used as the dependent variable in the regression.
  rast_gdp  = region_inputs_2020$Indonesia$gdp,
  
  # Population, NDVI, and MNDWI are the explanatory raster layers.
  rast_pop  = region_inputs_2020$Indonesia$pop,
  rast_ndvi = region_inputs_2020$Indonesia$ndvi,
  rast_mndwi = region_inputs_2020$Indonesia$mndwi,
  
  # Output folder for the exported workbook.
  out_dir   = cfg$out_dir,
  
  # Spatial block subsampling settings:
  # sbs_K = number of graph-connected subsamples
  # sbs_m = target number of retained cells in each subsample
  sbs_K = 30,
  sbs_m = 500,
  
  # Estimation settings:
  # "MC" uses Monte Carlo trace approximation for the log-determinant.
  # "queen" uses queen contiguity on the full template grid.
  # "manual" uses the custom sparse likelihood implementation.
  logdet_method = "MC",
  type_W = "queen",
  engine = "manual",
  
  # Random seed for reproducibility of stochastic steps.
  seed = 42,
  
  # If TRUE, SBS standard errors are computed and exported
  # alongside model-based MLE and HC1 standard errors.
  run_sbs = TRUE,
  export_block_draws = TRUE,
  # HC1 reporting mode:
  # "conditional" = HC1 for slope coefficients only; rho is evaluated using MLE and SBS
  # "full"        = approximate full-parameter HC1, including rho
  # "both"        = export both versions for transparency and robustness comparison
  hc1_mode = "both"
)


#
# out_java <- run_and_export_region(
#   dataset_name = "Java",
#   template  = region_inputs_2020$Java$template,
#   rast_gdp  = region_inputs_2020$Java$gdp,
#   rast_pop  = region_inputs_2020$Java$pop,
#   rast_ndvi = region_inputs_2020$Java$ndvi,
#   rast_mndwi = region_inputs_2020$Java$mndwi,
#   out_dir   = cfg$out_dir,
#   sbs_K = 30,
#   sbs_m = 500,
#   logdet_method = "MC",
#   type_W = "queen",
#   engine = "manual",
#   seed = 42,
#   run_sbs = TRUE,
#   export_block_draws = TRUE,
#   hc1_mode = "both"
# )
}
