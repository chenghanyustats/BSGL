library(MCMCpack)
library(splines)
library(GIGrvg)

# =========================================================
# Basis
# =========================================================

make_basis <- function(coords, L = 50, saved_knots = NULL) {
  n <- nrow(coords)

  if (!is.null(saved_knots)) {
    u_knots <- saved_knots$u_knots
    v_knots <- saved_knots$v_knots

    Bx <- bs(coords[, 1],
             knots = u_knots,
             Boundary.knots = saved_knots$u_boundary,
             intercept = TRUE,
             degree = saved_knots$degree)
    By <- bs(coords[, 2],
             knots = v_knots,
             Boundary.knots = saved_knots$v_boundary,
             intercept = TRUE,
             degree = saved_knots$degree)
  } else {
    L_sqrt <- ceiling(sqrt(L))
    Bx <- bs(coords[, 1], df = L_sqrt, intercept = TRUE)
    By <- bs(coords[, 2], df = L_sqrt, intercept = TRUE)
  }

  Psi <- matrix(0, nrow = n, ncol = L)
  idx <- 1
  for (i in 1:ncol(Bx)) {
    for (j in 1:ncol(By)) {
      if (idx <= L) {
        Psi[, idx] <- Bx[, i] * By[, j]
        idx <- idx + 1
      }
    }
  }
  while (idx <= L) {
    Psi[, idx] <- Bx[, ncol(Bx)] * By[, ncol(By)]
    idx <- idx + 1
  }

  if (is.null(saved_knots)) {
    attr(Psi, "saved_knots") <- list(
      u_knots = attr(Bx, "knots"),
      v_knots = attr(By, "knots"),
      u_boundary = attr(Bx, "Boundary.knots"),
      v_boundary = attr(By, "Boundary.knots"),
      degree = attr(Bx, "degree")
    )
  }

  Psi
}

# =========================================================
# BSGL fit
# =========================================================

#' Fit a Bayesian Spatial Group Lasso model
#'
#' Fits the proposed BSGL spatially varying coefficient model using
#' tensor-product spline basis functions and Bayesian group-lasso shrinkage.
#'
#' @param y Numeric response vector.
#' @param X Numeric covariate matrix with one column per covariate.
#' @param coords Numeric coordinate matrix with two columns.
#' @param L Number of tensor-product basis functions.
#' @param niter Total number of MCMC iterations.
#' @param nburn Number of burn-in iterations.
#' @param a_lam Shape hyperparameter for the global shrinkage prior.
#' @param b_lam Rate hyperparameter for the global shrinkage prior.
#' @param verbose Logical; if `TRUE`, print progress.
#'
#' @return An object of class `"bsgl_fit"` containing posterior samples and
#'   fitted basis information.
fit_bsgl <- function(y, X, coords, L = 36,
                     niter = 5000, nburn = 500,
                     a_lam = 15, b_lam = 1,
                     verbose = TRUE) {

  n <- length(y)
  m <- ncol(X)
  nsave <- niter - nburn

  if (verbose) {
    cat("BSGL\n")
    cat("n =", n, " p =", m, " L =", L, "\n")
    cat("lambda2 prior: Gamma(", a_lam, ", ", b_lam, ")\n", sep = "")
  }

  # basis
  Psi <- make_basis(coords, L)
  saved_knots <- attr(Psi, "saved_knots")

  # precompute X_j * Psi and XtX_j
  XP <- vector("list", m)
  XtX <- vector("list", m)
  for (j in 1:m) {
    XP[[j]] <- X[, j] * Psi
    XtX[[j]] <- crossprod(XP[[j]])
  }

  # ===== initialization: keep original =====
  B <- matrix(0, L, m)
  tau2 <- rep(0.01, m)
  s2 <- 0.1
  lam2 <- 1

  B_smp <- array(0, c(L, m, nsave))
  tau2_smp <- matrix(0, nsave, m)
  s2_smp <- numeric(nsave)
  lam2_smp <- numeric(nsave)

  # fitted/residual
  fit <- rep(0, n)
  res <- y - fit

  t0 <- Sys.time()

  for (it in 1:niter) {

    # ---- update beta_j with residual cycling
    for (j in 1:m) {
      b_old <- B[, j]
      Xb_old <- XP[[j]] %*% b_old

      # partial residual: add old contribution back
      rj <- res + Xb_old

      # ===== keep original posterior formula =====
      Prec <- (1 / s2) * XtX[[j]] + (1 / (s2 * tau2[j])) * diag(L)
      R <- chol(Prec)
      Vj <- chol2inv(R)
      mj <- (1 / s2) * Vj %*% crossprod(XP[[j]], rj)

      b_new <- as.numeric(rmvnorm(1, mj, Vj))
      B[, j] <- b_new

      Xb_new <- XP[[j]] %*% b_new
      res <- rj - Xb_new
    }

    # ---- update tau2_j: EXACTLY your original parameterization
    for (j in 1:m) {
      b2 <- sum(B[, j]^2)
      chi_j <- max(lam2 / 2, 1e-6)
      psi_j <- max(b2 / (2 * s2), 1e-6)
      gam_j <- rgig(1, lambda = -0.5, chi = chi_j, psi = psi_j)
      tau2[j] <- 1 / max(gam_j, 0.01)
    }

    # ---- update sigma2: same as original
    pen <- sum(sapply(1:m, function(j) sum(B[, j]^2) / tau2[j]))
    shp_s <- 1 + (n + m * L) / 2
    rate_s <- 1 + 0.5 * (sum(res^2) + pen)
    s2 <- rinvgamma(1, shape = shp_s, scale = rate_s)

    # ---- update lambda2: same as original
    shp_l <- a_lam + m * (L + 1) / 2
    rate_l <- b_lam + 0.5 * sum(tau2)
    lam2 <- rgamma(1, shape = shp_l, rate = rate_l)

    # ---- store
    if (it > nburn) {
      idx <- it - nburn
      B_smp[, , idx] <- B
      tau2_smp[idx, ] <- tau2
      s2_smp[idx] <- s2
      lam2_smp[idx] <- lam2
    }

    if (verbose && it %% 100 == 0) {
      elap <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      rate <- it / max(elap, 1e-8)
      eta <- (niter - it) / rate / 60
      bnorm <- apply(B, 2, function(z) sqrt(sum(z^2)))
      cat("Iter", it,
          "| s2:", round(s2, 4),
          "| lam2:", round(lam2, 4),
          "| ||B_j||:", paste(round(bnorm, 2), collapse = " "),
          "| ETA:", round(eta, 1), "min\n")
    }
  }

  out <- list(
    beta = B_smp,
    tau2 = tau2_smp,
    s2 = s2_smp,
    lam2 = lam2_smp,
    Psi = Psi,
    saved_knots = saved_knots,
    n = n,
    m = m,
    L = L
  )
  class(out) <- "bsgl_fit"
  out
}

# =========================================================
# Prediction
# =========================================================

#' Predict from a fitted BSGL model
#'
#' Computes posterior predictive summaries at new covariate and coordinate
#' values.
#'
#' @param fit A fitted object returned by [fit_bsgl()] or
#'   `fit_bsgl_multichain()`.
#' @param X_new Numeric matrix of new covariates.
#' @param coords_new Numeric matrix of new coordinates.
#' @param return_samples Logical; if `TRUE`, include posterior sample matrices.
#'
#' @return A list with posterior predictive means, credible limits, and
#'   optionally posterior samples.
pred_bsgl <- function(fit, X_new, coords_new, return_samples = FALSE) {
  Psi_new <- make_basis(
    coords_new,
    fit$L,
    saved_knots = fit$saved_knots
  )

  X_new <- as.matrix(X_new)
  n_new <- nrow(X_new)
  ns <- dim(fit$beta)[3]

  mu_smp <- matrix(0, nrow = n_new, ncol = ns)

  for (s in 1:ns) {
    mu <- rep(0, n_new)
    for (j in 1:fit$m) {
      mu <- mu + X_new[, j] * as.vector(Psi_new %*% fit$beta[, j, s])
    }
    mu_smp[, s] <- mu
  }

  y_smp <- mu_smp
  for (s in 1:ns) {
    y_smp[, s] <- mu_smp[, s] + rnorm(n_new, mean = 0, sd = sqrt(fit$s2[s]))
  }

  out <- list(
    mean = rowMeans(mu_smp),
    lower = apply(y_smp, 1, quantile, 0.025),
    upper = apply(y_smp, 1, quantile, 0.975),
    latent_lower = apply(mu_smp, 1, quantile, 0.025),
    latent_upper = apply(mu_smp, 1, quantile, 0.975)
  )

  if (return_samples) {
    out$mu_samples <- mu_smp
    out$y_samples <- y_smp
  }

  return(out)
}

# =========================================================
# Reconstruct coefficient surfaces
# =========================================================

betas_bsgl <- function(fit, coords_or_grid) {
  Bbar <- apply(fit$beta, c(1, 2), mean)

  if (is.list(coords_or_grid) && "grid_points" %in% names(coords_or_grid)) {
    coords <- as.matrix(coords_or_grid$grid_points)
  } else {
    coords <- as.matrix(coords_or_grid)
  }

  Psi_new <- make_basis(coords, fit$L, saved_knots = fit$saved_knots)

  out <- matrix(0, nrow(Psi_new), fit$m)
  for (j in 1:fit$m) {
    out[, j] <- Psi_new %*% Bbar[, j]
  }
  out
}

# =========================================================
# Summary
# =========================================================

sum_bsgl <- function(fit, thresh = 0.1) {
  Bbar <- apply(fit$beta, c(1, 2), mean)
  norms <- matrix(0, dim(fit$beta)[3], fit$m)

  for (j in 1:fit$m) {
    norms[, j] <- apply(fit$beta[, j, ], 2, function(z) sqrt(sum(z^2)))
  }

  pip <- colMeans(norms > thresh)

  list(
    Bbar = Bbar,
    mean_norm = colMeans(norms),
    pip = pip,
    s2_mean = mean(fit$s2),
    lam2_mean = mean(fit$lam2)
  )
}

# =========================================================
# Evaluation
# =========================================================
eval_bsgl <- function(fit, test_data) {
  pr <- pred_bsgl(fit, test_data$X, test_data$coords)

  cover <- mean(test_data$y >= pr$lower & test_data$y <= pr$upper)
  mspe  <- mean((test_data$y - pr$mean)^2)

  list(
    mspe = mspe,
    cover = cover
  )
}


fit_bsgl_multichain <- function(y, X, coords, L = 36,
                                niter = 5000, nburn = 500,
                                a_lam = 15, b_lam = 1,
                                nchains = 4,
                                parallel = TRUE,
                                seeds = NULL,
                                verbose = TRUE) {

  if (is.null(seeds)) {
    seeds <- 1000 + seq_len(nchains)
  }
  stopifnot(length(seeds) >= nchains)

  if (verbose) {
    cat("=== Multi-chain BSGL Fitting ===\n")
    cat("Chains:", nchains,
        "| Iterations per chain:", niter,
        "| Burn-in:", nburn, "\n")
  }

  run_one_chain <- function(chain_id, seed) {
    set.seed(seed)
    if (verbose) cat("Starting chain", chain_id, "with seed", seed, "\n")

    fit_bsgl(
      y = y,
      X = X,
      coords = coords,
      L = L,
      niter = niter,
      nburn = nburn,
      a_lam = a_lam,
      b_lam = b_lam,
      verbose = FALSE
    )
  }

  if (parallel) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(nchains)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      parallel::clusterExport(
        cl,
        varlist = c(
          "make_basis", "fit_bsgl",
          "y", "X", "coords", "L", "niter", "nburn",
          "a_lam", "b_lam", "verbose"
        ),
        envir = environment()
      )

      parallel::clusterEvalQ(cl, {
        library(MCMCpack)
        library(splines)
        library(GIGrvg)
        NULL
      })

      fits <- parallel::parLapply(
        cl,
        X = seq_len(nchains),
        fun = function(i) {
          set.seed(seeds[i])
          fit_bsgl(
            y = y,
            X = X,
            coords = coords,
            L = L,
            niter = niter,
            nburn = nburn,
            a_lam = a_lam,
            b_lam = b_lam,
            verbose = FALSE
          )
        }
      )
    } else {
      fits <- parallel::mclapply(
        seq_len(nchains),
        function(i) run_one_chain(i, seeds[i]),
        mc.cores = nchains
      )
    }
  } else {
    fits <- lapply(seq_len(nchains), function(i) run_one_chain(i, seeds[i]))
  }

  if (verbose) cat("Combining chains...\n")

  beta_all <- abind::abind(lapply(fits, `[[`, "beta"), along = 3)
  tau2_all <- do.call(rbind, lapply(fits, `[[`, "tau2"))
  s2_all   <- unlist(lapply(fits, `[[`, "s2"))
  lam2_all <- unlist(lapply(fits, `[[`, "lam2"))

  out <- list(
    beta = beta_all,
    tau2 = tau2_all,
    s2 = s2_all,
    lam2 = lam2_all,
    Psi = fits[[1]]$Psi,
    saved_knots = fits[[1]]$saved_knots,
    n = fits[[1]]$n,
    m = fits[[1]]$m,
    L = fits[[1]]$L,
    nchains = nchains,
    chain_seeds = seeds[seq_len(nchains)]
  )

  class(out) <- "bsgl_fit"

  if (verbose) {
    cat("Done. Total posterior draws kept:", dim(out$beta)[3], "\n")
  }

  out
}
