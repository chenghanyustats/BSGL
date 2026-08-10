library(mvtnorm)
library(MCMCpack)
library(splines)
library(ggplot2)
library(cowplot)
library(viridis)

# =========================================================
# Gaussian-prior basis SVC
# =========================================================

#' Fit a Gaussian-prior spatially varying coefficient model
#'
#' Fits the GSVC comparison model using the same spline basis representation as
#' BSGL but without group-lasso shrinkage.
#'
#' @param y Numeric response vector.
#' @param X Numeric covariate matrix with one column per covariate.
#' @param coords Numeric coordinate matrix with two columns.
#' @param L Number of tensor-product basis functions.
#' @param niter Total number of MCMC iterations.
#' @param nburn Number of burn-in iterations.
#' @param a_k,b_k Hyperparameters for the coefficient-scale prior.
#' @param a_s,b_s Hyperparameters for the error-variance prior.
#' @param verbose Logical; if `TRUE`, print progress.
#'
#' @return An object of class `"gs_svc"` containing posterior samples and fitted
#'   basis information.
fit_gs <- function(y, X, coords, L = 36,
                   niter = 5000, nburn = 500,
                   a_k = 2, b_k = 1,
                   a_s = 2, b_s = 1,
                   verbose = TRUE) {
  
  n <- length(y)
  m <- ncol(X)
  nsave <- niter - nburn
  
  if (verbose) {
    cat("Gaussian basis SVC\n")
    cat("n =", n, " p =", m, " L =", L, "\n")
  }
  
  # basis
  Psi <- make_basis(coords, L)
  saved_knots <- attr(Psi, "saved_knots")
  
  # precompute X_j * Psi and crossprod
  XP <- vector("list", m)
  XtX <- vector("list", m)
  Xty <- vector("list", m)
  
  for (j in 1:m) {
    XP[[j]] <- X[, j] * Psi
    XtX[[j]] <- crossprod(XP[[j]])
    Xty[[j]] <- crossprod(XP[[j]], y)
  }
  
  # init
  A <- matrix(0, L, m)
  s2 <- var(y)
  k2 <- 1
  
  A_smp <- array(0, c(L, m, nsave))
  s2_smp <- numeric(nsave)
  k2_smp <- numeric(nsave)
  
  # fitted and residual
  fit <- rep(0, n)
  res <- y - fit
  
  t0 <- Sys.time()
  
  for (it in 1:niter) {
    
    # ---- update A[,j]
    for (j in 1:m) {
      a_old <- A[, j]
      Xa_old <- XP[[j]] %*% a_old
      
      # add back old contribution
      rj <- res + Xa_old
      
      # posterior
      Prec <- XtX[[j]] + diag(1 / k2, L)
      R <- chol(Prec)
      V0 <- chol2inv(R)
      mj <- V0 %*% crossprod(XP[[j]], rj)
      Vj <- s2 * V0
      
      a_new <- as.numeric(rmvnorm(1, mj, Vj))
      A[, j] <- a_new
      
      Xa_new <- XP[[j]] %*% a_new
      
      # update residual efficiently
      res <- rj - Xa_new
    }
    
    # ---- update k2
    A2 <- sum(A^2)
    shp_k <- a_k + m * L / 2
    scl_k <- b_k + A2 / (2 * s2)
    k2 <- rinvgamma(1, shape = shp_k, scale = scl_k)
    
    # ---- update s2
    SSE <- sum(res^2)
    S <- 0.5 * SSE + 0.5 * A2 / k2
    shp_s <- a_s + (n + m * L) / 2
    scl_s <- b_s + S
    s2 <- rinvgamma(1, shape = shp_s, scale = scl_s)
    
    # ---- store
    if (it > nburn) {
      idx <- it - nburn
      A_smp[, , idx] <- A
      s2_smp[idx] <- s2
      k2_smp[idx] <- k2
    }
    
    if (verbose && it %% 100 == 0) {
      elap <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      rate <- it / max(elap, 1e-8)
      eta <- (niter - it) / rate / 60
      anorm <- apply(A, 2, function(z) sqrt(sum(z^2)))
      cat("Iter", it,
          "| s2:", round(s2, 4),
          "| k2:", round(k2, 4),
          "| ||A_j||:", paste(round(anorm, 2), collapse = " "),
          "| ETA:", round(eta, 1), "min\n")
    }
  }
  
  out <- list(
    A = A_smp,
    s2 = s2_smp,
    k2 = k2_smp,
    Psi = Psi,
    saved_knots = saved_knots,
    n = n,
    m = m,
    L = L
  )
  class(out) <- "gs_svc"
  out
}

# =========================================================
# prediction
# =========================================================

#' Predict from a fitted GSVC model
#'
#' Computes posterior predictive summaries at new covariate and coordinate
#' values.
#'
#' @param fit A fitted object returned by [fit_gs()].
#' @param X_new Numeric matrix of new covariates.
#' @param coords_new Numeric matrix of new coordinates.
#' @param return_samples Logical; if `TRUE`, include posterior sample matrices.
#'
#' @return A list with posterior predictive means, credible limits, and
#'   optionally posterior samples.
pred_gs <- function(fit, X_new, coords_new, return_samples = FALSE) {
  Psi_new <- make_basis(
    coords_new,
    fit$L,
    saved_knots = fit$saved_knots
  )
  
  X_new <- as.matrix(X_new)
  n_new <- nrow(X_new)
  ns <- dim(fit$A)[3]
  
  mu_smp <- matrix(0, nrow = n_new, ncol = ns)
  
  for (s in 1:ns) {
    mu <- rep(0, n_new)
    for (j in 1:fit$m) {
      mu <- mu + X_new[, j] * as.vector(Psi_new %*% fit$A[, j, s])
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
# coefficient surfaces
# =========================================================

betas_gs <- function(fit, coords_or_grid) {
  Abar <- apply(fit$A, c(1, 2), mean)
  
  if (is.list(coords_or_grid) && "grid_points" %in% names(coords_or_grid)) {
    coords <- as.matrix(coords_or_grid$grid_points)
  } else {
    coords <- as.matrix(coords_or_grid)
  }
  
  Psi_new <- make_basis(coords, fit$L, saved_knots = fit$saved_knots)
  
  B <- matrix(0, nrow(Psi_new), fit$m)
  for (j in 1:fit$m) {
    B[, j] <- Psi_new %*% Abar[, j]
  }
  B
}

# =========================================================
# SCP and credible-interval maps
# =========================================================

gs_surface_ci <- function(fit, coords_or_grid, ci_level = 0.95) {
  if (is.list(coords_or_grid) && "grid_points" %in% names(coords_or_grid)) {
    coords <- as.matrix(coords_or_grid$grid_points)
  } else {
    coords <- as.matrix(coords_or_grid)
  }

  Psi_new <- make_basis(coords, fit$L, saved_knots = fit$saved_knots)
  n_grid <- nrow(Psi_new)
  n_samples <- dim(fit$A)[3]

  alpha <- 1 - ci_level
  lower_prob <- alpha / 2
  upper_prob <- 1 - alpha / 2

  mean_beta <- matrix(0, n_grid, fit$m)
  lower_beta <- matrix(0, n_grid, fit$m)
  upper_beta <- matrix(0, n_grid, fit$m)
  excludes0 <- matrix(FALSE, n_grid, fit$m)

  for (j in seq_len(fit$m)) {
    beta_samples <- matrix(0, n_grid, n_samples)
    for (s in seq_len(n_samples)) {
      beta_samples[, s] <- Psi_new %*% fit$A[, j, s]
    }

    mean_beta[, j] <- rowMeans(beta_samples)
    lower_beta[, j] <- apply(beta_samples, 1, quantile, lower_prob)
    upper_beta[, j] <- apply(beta_samples, 1, quantile, upper_prob)
    excludes0[, j] <- (lower_beta[, j] > 0) | (upper_beta[, j] < 0)
  }

  colnames(mean_beta) <- paste0("V", seq_len(fit$m))
  colnames(lower_beta) <- colnames(mean_beta)
  colnames(upper_beta) <- colnames(mean_beta)
  colnames(excludes0) <- colnames(mean_beta)

  list(
    mean = mean_beta,
    lower = lower_beta,
    upper = upper_beta,
    excludes0 = excludes0
  )
}

#' Compute GSVC spatial coverage probabilities
#'
#' Computes the spatial coverage probability (SCP) for each covariate surface as
#' the proportion of grid locations where the pointwise credible interval
#' excludes zero.
#'
#' @param fit A fitted object returned by [fit_gs()].
#' @param coords_or_grid A coordinate matrix or a list with element
#'   `grid_points`.
#' @param ci_level Credible interval level.
#'
#' @return A named numeric vector of SCP values.
calc_scp_gs <- function(fit, coords_or_grid, ci_level = 0.95) {
  ci <- gs_surface_ci(fit, coords_or_grid, ci_level = ci_level)
  scp <- colMeans(ci$excludes0)
  names(scp) <- paste0("β", seq_len(fit$m))
  scp
}

# =========================================================
# summary
# =========================================================

sum_gs <- function(fit) {
  Abar <- apply(fit$A, c(1, 2), mean)
  norms <- matrix(0, dim(fit$A)[3], fit$m)
  
  for (j in 1:fit$m) {
    norms[, j] <- apply(fit$A[, j, ], 2, function(z) sqrt(sum(z^2)))
  }
  
  list(
    Abar = Abar,
    mean_norm = colMeans(norms),
    s2_mean = mean(fit$s2),
    k2_mean = mean(fit$k2)
  )
}

# =========================================================
# one-shot evaluation
# =========================================================

eval_gs <- function(fit, test_data) {
  pr <- pred_gs(fit, test_data$X, test_data$coords)
  
  cover <- mean(test_data$y >= pr$lower & test_data$y <= pr$upper)
  mspe  <- mean((test_data$y - pr$mean)^2)
  
  list(
    mspe = mspe,
    cover = cover
  )
}

plot_comp <- function(bsgl_fit, gs_fit, gam_model, mgwr_model,
                      grid_data, meta_info, j,
                      n_val, p_val,
                      grayscale = FALSE, show_title = TRUE,
                      save_plot = TRUE, save_dir = "plots",
                      width = 13, height = 2.6, dpi = 300) {
  library(ggplot2)
  library(cowplot)
  library(viridis)
  library(grid)
  
  # -----------------------------
  # create save dir
  # -----------------------------
  if (save_plot) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # -----------------------------
  # grid points
  # -----------------------------
  grid_points <- as.data.frame(grid_data$grid_points)
  colnames(grid_points) <- c("x", "y")
  
  # -----------------------------
  # beta surfaces
  # -----------------------------
  beta_true <- calc_true_beta(as.matrix(grid_points), meta_info$p)[, j]
  beta_bsgl <- betas_bsgl(bsgl_fit, grid_data)[, j]
  beta_gam  <- get_gam_betas(gam_model, as.matrix(grid_points))[, j]
  beta_mgwr <- get_mgwr_betas(mgwr_model, as.matrix(grid_points))[, j]
  beta_gs   <- betas_gs(gs_fit, grid_data)[, j]
  
  plot_data <- list(
    "True"         = data.frame(x = grid_points$x, y = grid_points$y, beta = beta_true),
    "BSGL"         = data.frame(x = grid_points$x, y = grid_points$y, beta = beta_bsgl),
    "GGP-GAM"      = data.frame(x = grid_points$x, y = grid_points$y, beta = beta_gam),
    "GSVC"         = data.frame(x = grid_points$x, y = grid_points$y, beta = beta_gs),
    "MGWR"         = data.frame(x = grid_points$x, y = grid_points$y, beta = beta_mgwr)
  )
  
  lims <- range(unlist(lapply(plot_data, function(df) df$beta)), na.rm = TRUE)
  
  make_panel <- function(df, ttl, keep_legend = FALSE) {
    p <- ggplot(df, aes(x = x, y = y, fill = beta)) +
      geom_raster() +
      coord_fixed() +
      theme_void(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.margin = margin(1, 1, 1, 1),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0, "pt")
      )
    
    if (show_title) p <- p + ggtitle(ttl)
    
    if (grayscale) {
      p <- p + scale_fill_gradient(
        low = "white", high = "black",
        limits = lims, name = NULL,
        guide = guide_colorbar(
          barheight = unit(4.2, "cm"),   # 比之前更高
          barwidth  = unit(0.22, "cm"),
          ticks = TRUE
        )
      )
    } else {
      p <- p + scale_fill_viridis_c(
        limits = lims, name = NULL,
        guide = guide_colorbar(
          barheight = unit(4.2, "cm"),   # 比之前更高
          barwidth  = unit(0.22, "cm"),
          ticks = TRUE
        )
      )
    }
    
    if (!keep_legend) {
      p <- p + theme(legend.position = "none")
    } else {
      p <- p + theme(
        legend.position = "right",
        legend.text = element_text(size = 10, face = "bold")
      )
    }
    
    p
  }
  
  p1 <- make_panel(plot_data[["True"]], "True", FALSE)
  p2 <- make_panel(plot_data[["BSGL"]], "BSGL", FALSE)
  p3 <- make_panel(plot_data[["GGP-GAM"]], "GGP-GAM", FALSE)
  p4 <- make_panel(plot_data[["GSVC"]], "GSVC", FALSE)
  p5_leg <- make_panel(plot_data[["MGWR"]], "MGWR", TRUE)
  
  lgd <- cowplot::get_legend(p5_leg)
  p5  <- p5_leg + theme(legend.position = "none")
  
  main_row <- cowplot::plot_grid(
    p1, p2, p3, p4, p5,
    nrow = 1,
    rel_widths = c(1, 1, 1, 1, 1),
    align = "h",
    axis = "tb"
  )
  
  final_plot <- cowplot::plot_grid(
    main_row, lgd,
    nrow = 1,
    rel_widths = c(1, 0.04)
  )
  
  if (save_plot) {
    file_name <- paste0("beta", j, "_compare_n", n_val, "_p", p_val, ".png")
    ggsave(
      filename = file.path(save_dir, file_name),
      plot = final_plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }
  
  return(final_plot)
}


# 
# plot_comp <- function(bsgl_fit, gs_fit, ssgl_model,
#                               grid_data, meta_info, j,
#                               n_basis, global_settings,
#                               grayscale = FALSE, show_title = TRUE) {
#   library(ggplot2)
#   library(cowplot)
#   library(viridis)
#   
#   grid_points <- grid_data$grid_points
#   
#   beta_true <- calc_true_beta(grid_points, meta_info$p)[, j]
#   beta_bsgl <- betas_bsgl(bsgl_fit, grid_data)[, j]
#   beta_gs   <- betas_gs(gs_fit, grid_data)[, j]
#   beta_ssgl <- betas_ssgl(ssgl_model, grid_data, n_basis, global_settings)[, j]
#   
#   plot_data <- list(
#     "True" = data.frame(x = grid_points$x, y = grid_points$y, beta = beta_true),
#     "BSGL" = data.frame(x = grid_points$x, y = grid_points$y, beta = beta_bsgl),
#     "Gaussian-SVC" = data.frame(x = grid_points$x, y = grid_points$y, beta = beta_gs),
#     "SSGL" = data.frame(x = grid_points$x, y = grid_points$y, beta = beta_ssgl)
#   )
#   
#   lims <- range(unlist(lapply(plot_data, function(df) df$beta)), na.rm = TRUE)
#   
#   make_panel <- function(df, ttl, keep_legend = FALSE) {
#     p <- ggplot(df, aes(x = x, y = y, fill = beta)) +
#       geom_raster() +
#       coord_fixed() +
#       theme_bw() +
#       theme(
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
#         plot.margin = margin(2, 2, 2, 2)
#       )
#     
#     if (show_title) p <- p + ggtitle(ttl)
#     
#     if (grayscale) {
#       p <- p + scale_fill_gradient(
#         low = "white", high = "black",
#         limits = lims, name = NULL
#       )
#     } else {
#       p <- p + scale_fill_viridis_c(
#         limits = lims, name = NULL
#       )
#     }
#     
#     if (!keep_legend) {
#       p <- p + theme(legend.position = "none")
#     } else {
#       p <- p + theme(
#         legend.position = "right",
#         legend.key.width = unit(0.25, "cm"),
#         legend.key.height = unit(1.0, "cm"),
#         legend.text = element_text(size = 10, face = "bold")
#       )
#     }
#     
#     p
#   }
#   
#   p1 <- make_panel(plot_data[["True"]], "True", FALSE)
#   p2 <- make_panel(plot_data[["BSGL"]], "BSGL", FALSE)
#   p3 <- make_panel(plot_data[["Gaussian-SVC"]], "Gaussian-SVC", FALSE)
#   p4_leg <- make_panel(plot_data[["SSGL"]], "SSGL", TRUE)
#   
#   lgd <- cowplot::get_legend(p4_leg)
#   p4 <- p4_leg + theme(legend.position = "none")
#   
#   main_row <- cowplot::plot_grid(
#     p1, p2, p3, p4,
#     nrow = 1,
#     rel_widths = c(1, 1, 1, 1)
#   )
#   
#   final_plot <- cowplot::plot_grid(
#     main_row, lgd,
#     nrow = 1,
#     rel_widths = c(1, 0.06)
#   )
#   
#   final_plot
# }
