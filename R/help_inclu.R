# 4.28 update

#' Compute BSGL spatial coverage probabilities
#'
#' Computes the spatial coverage probability (SCP) for each BSGL coefficient
#' surface as the proportion of grid locations where the pointwise credible
#' interval excludes zero.
#'
#' @param bsgl_model A fitted object returned by [fit_bsgl()] or
#'   `fit_bsgl_multichain()`.
#' @param grid_data A list with element `grid_points` containing grid
#'   coordinates.
#' @param ci_level Credible interval level.
#'
#' @return A named numeric vector of SCP values.
calc_scp_bsgl <- function(bsgl_model, grid_data, ci_level = 0.95) {
  
  m <- bsgl_model$m
  n_samples <- dim(bsgl_model$beta)[3]
  
  grid_coords <- as.matrix(grid_data$grid_points)
  n_grid <- nrow(grid_coords)
  
  Psi_grid <- make_basis(
    grid_coords,
    bsgl_model$L,
    saved_knots = bsgl_model$saved_knots
  )
  
  scp <- numeric(m)
  names(scp) <- paste0("β", 1:m)
  
  alpha <- 1 - ci_level
  lower_prob <- alpha / 2
  upper_prob <- 1 - alpha / 2
  
  for (j in 1:m) {
    beta_grid_samples <- matrix(0, n_grid, n_samples)
    
    for (s in 1:n_samples) {
      beta_grid_samples[, s] <- Psi_grid %*% bsgl_model$beta[, j, s]
    }
    
    lower <- apply(beta_grid_samples, 1, quantile, lower_prob)
    upper <- apply(beta_grid_samples, 1, quantile, upper_prob)
    
    excludes0 <- (lower > 0) | (upper < 0)
    scp[j] <- mean(excludes0)
  }
  
  return(scp)
}

# Function to calculate scp (Posterior Inclusion Probability) for Lasso
calc_scp <- function(lasso_model, grid_data, ci_level = 0.95) {  # 加ci_level参数
  
  m <- lasso_model$m
  n_samples <- dim(lasso_model$beta_samples)[3]
  grid_coords <- as.matrix(grid_data$grid_points)
  n_grid <- nrow(grid_coords)
  
  # Generate basis functions using saved knots
  Psi_grid <- make_basis(grid_coords, lasso_model$L, saved_knots = lasso_model$saved_knots)
  
  # Calculate scp for each variable
  scp <- numeric(m)
  names(scp) <- paste0("β", 1:m)
  
  # 计算CI的上下界
  alpha <- 1 - ci_level
  lower_prob <- alpha / 2
  upper_prob <- 1 - alpha / 2
  
  for (j in 1:m) {
    # Get beta samples for variable j at all grid points
    beta_grid_samples <- matrix(0, n_grid, n_samples)
    
    for (s in 1:n_samples) {
      beta_grid_samples[, s] <- Psi_grid %*% lasso_model$beta_samples[, j, s]
    }
    
    # Calculate CI for each grid point
    lower <- apply(beta_grid_samples, 1, quantile, lower_prob)
    upper <- apply(beta_grid_samples, 1, quantile, upper_prob)
    
    # Check if CI EXCLUDES 0
    excludes0 <- (lower > 0) | (upper < 0)
    scp[j] <- mean(excludes0)
  }
  
  return(scp)
}

# Function to calculate exclusion map (for visualization)
calc_exclusion_map <- function(lasso_model, grid_data) {
  
  m <- lasso_model$m
  n_samples <- dim(lasso_model$beta_samples)[3]
  grid_coords <- as.matrix(grid_data$grid_points)
  n_grid <- nrow(grid_coords)
  
  # Generate basis functions using saved knots
  Psi_grid <- make_basis(grid_coords, lasso_model$L, saved_knots = lasso_model$saved_knots)
  
  # Store exclusion maps for all variables
  exclusion_maps <- array(NA, dim = c(n_grid, m))
  colnames(exclusion_maps) <- paste0("beta", 1:m)
  
  for (j in 1:m) {
    # Get beta samples for variable j at all grid points
    beta_grid_samples <- matrix(0, n_grid, n_samples)
    
    for (s in 1:n_samples) {
      beta_grid_samples[, s] <- Psi_grid %*% lasso_model$beta_samples[, j, s]
    }
    
    # Calculate 95% CI for each grid point
    lower <- apply(beta_grid_samples, 1, quantile, 0.025)
    upper <- apply(beta_grid_samples, 1, quantile, 0.975)
    
    # Check if CI EXCLUDES 0 (TRUE = significant, FALSE = includes 0)
    exclusion_maps[, j] <- (lower > 0) | (upper < 0)  # 改这里！
  }
  
  return(exclusion_maps)
}

# Function to plot beta credible intervals for Lasso
plot_beta_ci <- function(lasso_model, grid_data) {
  
  library(ggplot2)
  library(viridis)
  
  m <- lasso_model$m
  n_samples <- dim(lasso_model$beta_samples)[3]
  grid_coords <- as.matrix(grid_data$grid_points)
  
  # Generate basis functions using saved knots
  Psi_grid <- make_basis(grid_coords, lasso_model$L, saved_knots = lasso_model$saved_knots)
  
  plots <- list()
  
  for (j in 1:m) {
    # Get beta samples for variable j at all grid points
    beta_grid_samples <- matrix(0, nrow(grid_coords), n_samples)
    
    for (s in 1:n_samples) {
      beta_grid_samples[, s] <- Psi_grid %*% lasso_model$beta_samples[, j, s]
    }
    
    # Calculate statistics
    lower <- apply(beta_grid_samples, 1, quantile, 0.025)
    upper <- apply(beta_grid_samples, 1, quantile, 0.975)
    mean_beta <- apply(beta_grid_samples, 1, mean)
    incl0 <- (lower <= 0) & (upper >= 0)
    incl0_rate <- mean(incl0)
    
    # Create plot data
    df <- data.frame(
      x = grid_coords[, 1],
      y = grid_coords[, 2],
      mean_beta = mean_beta,
      incl0 = incl0
    )
    
    # Create plot
    p <- ggplot(df, aes(x, y, color = mean_beta, shape = incl0)) +
      geom_point(size = 1.5, alpha = 0.8) +
      scale_color_viridis_c(name = "Mean β") +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), 
                         name = "95% CI\nincludes 0",
                         labels = c("No", "Yes")) +
      labs(title = sprintf("β%d (%.1f%% include 0)", j, incl0_rate * 100),
           x = "X", y = "Y") +
      theme_minimal() +
      coord_fixed()
    
    plots[[j]] <- p
  }
  
  # Print plots
  for (p in plots) {
    if (!is.null(p)) print(p)
  }
  
  return(plots)
}

# Function to create proportion maps for Lasso replicates
prop_maps <- function(incl0_maps_list, grid_data) {
  
  library(ggplot2)
  library(viridis)
  
  n_reps <- length(incl0_maps_list)
  grid_coords <- as.matrix(grid_data$grid_points)
  n_grid <- nrow(grid_coords)
  m <- ncol(incl0_maps_list[[1]])
  
  # Calculate proportion of replicates where each point includes 0
  prop_maps <- array(0, dim = c(n_grid, m))
  colnames(prop_maps) <- paste0("beta", 1:m)
  
  # Sum across replicates
  for (rep in 1:n_reps) {
    prop_maps <- prop_maps + incl0_maps_list[[rep]]
  }
  
  # Convert to proportions
  prop_maps <- prop_maps / n_reps
  
  # Create plots
  plots <- list()
  
  for (j in 1:m) {
    df <- data.frame(
      x = grid_coords[, 1],
      y = grid_coords[, 2],
      prop_incl0 = prop_maps[, j]
    )
    
    p <- ggplot(df, aes(x, y, color = prop_incl0)) +
      geom_point(size = 1.5, alpha = 0.8) +
      scale_color_viridis_c(name = "Proportion\nincluding 0", 
                            limits = c(0, 1)) +
      labs(title = sprintf("β%d: Proportion of CI including 0 (n=%d reps)", j, n_reps),
           x = "X", y = "Y") +
      theme_minimal() +
      coord_fixed()
    
    plots[[j]] <- p
  }
  
  # Print plots
  for (p in plots) {
    print(p)
  }
  
  return(list(
    prop_maps = prop_maps,
    plots = plots
  ))
}

# Function to save include-zero maps for Lasso
save_incl0_maps <- function(incl0_maps_list, filename) {
  # Convert list of matrices to 3D array for efficient storage
  n_grid <- nrow(incl0_maps_list[[1]])
  m <- ncol(incl0_maps_list[[1]])
  n_reps <- length(incl0_maps_list)
  
  incl0_array <- array(NA, dim = c(n_grid, m, n_reps))
  
  for (rep in 1:n_reps) {
    incl0_array[, , rep] <- incl0_maps_list[[rep]]
  }
  
  save(incl0_array, file = filename)
  cat("Saved include-zero maps to", filename, "\n")
}


plot_beta_ci_facet <- function(lasso_model, grid_data, beta_indices = 1:3, 
                               n = NULL, p = NULL, ncol = 3, 
                               width = NULL, height = NULL, show_legend = TRUE,
                               grayscale = FALSE) {  # 新增参数
  library(ggplot2)
  library(viridis)
  library(dplyr)
  
  dir.create("plots", showWarnings = FALSE)
  
  n_samples <- dim(lasso_model$beta_samples)[3]
  grid_coords <- as.matrix(grid_data$grid_points)
  Psi_grid <- make_basis(grid_coords, lasso_model$L, saved_knots = lasso_model$saved_knots)
  
  all_data <- data.frame()
  labels <- c()
  
  for (i in seq_along(beta_indices)) {
    j <- beta_indices[i]
    
    beta_grid_samples <- matrix(0, nrow(grid_coords), n_samples)
    for (s in 1:n_samples) {
      beta_grid_samples[, s] <- Psi_grid %*% lasso_model$beta_samples[, j, s]
    }
    
    lower <- apply(beta_grid_samples, 1, quantile, 0.025)
    upper <- apply(beta_grid_samples, 1, quantile, 0.975)
    mean_beta <- apply(beta_grid_samples, 1, mean)
    
    excludes0 <- (lower > 0) | (upper < 0)
    scp <- mean(excludes0)
    
    df_j <- data.frame(
      x = grid_coords[, 1],
      y = grid_coords[, 2],
      mean_beta = mean_beta,
      excludes0 = excludes0,
      beta_var = paste0("beta[", j, "]")
    )
    
    all_data <- rbind(all_data, df_j)
    labels[df_j$beta_var[1]] <- paste0("β", j, " (SCP = ", sprintf("%.1f%%", scp*100), ")")
  }
  
  p <- ggplot(all_data, aes(x, y)) +
    geom_point(
      aes(color = mean_beta, shape = as.factor(excludes0)),
      size = 1.3
    ) +
    scale_shape_manual(
      values = c(`FALSE` = 1, `TRUE` = 16),
      name = "95% CI excludes 0",
      labels = c("No", "Yes")
    ) +
    facet_wrap(~ beta_var, ncol = ncol, labeller = labeller(beta_var = labels)) +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = if(show_legend) "bottom" else "none",
      legend.key.size = unit(0.5, "cm"),
      legend.title = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, size = 0.4),
      panel.grid = element_blank(),
      plot.margin = margin(0.2, 0.2, 0.1, 0.2, "cm")
    )
  
  # 根据参数选择颜色方案
  if (grayscale) {
    p <- p + scale_color_gradient(
      low = "white", 
      high = "black",
      guide = "none"
    )
  } else {
    p <- p + scale_color_viridis_c(guide = "none")
  }
  
  beta_str <- paste(beta_indices, collapse = "")
  if (!is.null(n) && !is.null(p)) {
    filename <- paste0("beta_ci_facet_", beta_str, "_n", n, "_p", p, 
                       if(grayscale) "_gray" else "", ".png")
  } else {
    filename <- paste0("beta_ci_facet_", beta_str, 
                       if(grayscale) "_gray" else "", ".png")
  }
  save_path <- file.path("plots", filename)
  
  if (is.null(width) || is.null(height)) {
    if (length(beta_indices) == 1) {
      width <- 800 / 300
      height <- 340 / 300
    } else if (length(beta_indices) <= 3) {
      width <- 800 / 300
      height <- 340 / 300
    } else {
      width <- 800 / 300
      height <- 400 / 300
    }
  }
  ggsave(save_path, p, width = width, height = height, dpi = 300, bg = "white")
  
  return(p)
}

# Function to calculate F1 score with detailed breakdown including FPR
calc_f1_beta <- function(lasso_model, grid_data, true_beta, ci_level = 0.95, verbose = FALSE) {
  
  m <- lasso_model$m
  n_samples <- dim(lasso_model$beta_samples)[3]
  grid_coords <- as.matrix(grid_data$grid_points)
  n_grid <- nrow(grid_coords)
  
  # Generate basis functions
  Psi_grid <- make_basis(grid_coords, lasso_model$L, saved_knots = lasso_model$saved_knots)
  
  # Calculate CI bounds
  alpha <- 1 - ci_level
  lower_prob <- alpha / 2
  upper_prob <- 1 - alpha / 2
  
  # Results
  f1_scores <- numeric(m)
  names(f1_scores) <- paste0("β", 1:m)
  
  # Store detailed results
  details <- list()
  
  for (j in 1:m) {
    # Get beta samples at grid points
    beta_grid_samples <- matrix(0, n_grid, n_samples)
    for (s in 1:n_samples) {
      beta_grid_samples[, s] <- Psi_grid %*% lasso_model$beta_samples[, j, s]
    }
    
    # Calculate CI
    lower <- apply(beta_grid_samples, 1, quantile, lower_prob)
    upper <- apply(beta_grid_samples, 1, quantile, upper_prob)
    
    # Predicted: CI excludes 0
    predicted_signal <- (lower > 0) | (upper < 0)
    
    # True: true beta is non-zero
    true_signal <- abs(true_beta[, j]) > 1e-8
    
    # Confusion matrix
    TP <- sum(predicted_signal & true_signal)
    FP <- sum(predicted_signal & !true_signal)
    FN <- sum(!predicted_signal & true_signal)
    TN <- sum(!predicted_signal & !true_signal)
    
    # Metrics
    precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
    recall <- ifelse(TP + FN > 0, TP / (TP + FN), 0)
    specificity <- ifelse(TN + FP > 0, TN / (TN + FP), 0)
    fpr <- ifelse(FP + TN > 0, FP / (FP + TN), 0)
    
    f1_scores[j] <- ifelse(precision + recall > 0, 
                           2 * precision * recall / (precision + recall), 0)
    
    # Store details
    details[[paste0("β", j)]] <- list(
      TP = TP, FP = FP, FN = FN, TN = TN,
      precision = precision,
      recall = recall,
      specificity = specificity,
      fpr = fpr,
      f1 = f1_scores[j],
      n_true_signal = sum(true_signal),
      n_pred_signal = sum(predicted_signal),
      n_grid = n_grid
    )
    
    if (verbose) {
      cat("\n=== β", j, " ===\n", sep = "")
      cat("True Positive (correct signal):", TP, "\n")
      cat("False Positive (false alarm):", FP, "\n")
      cat("False Negative (missed signal):", FN, "\n")
      cat("True Negative (correct null):", TN, "\n")
      cat("Precision:", round(precision, 4), "\n")
      cat("Recall:", round(recall, 4), "\n")
      cat("Specificity:", round(specificity, 4), "\n")
      cat("FPR:", round(fpr, 4), "\n")
      cat("F1:", round(f1_scores[j], 4), "\n")
    }
  }
  
  if (verbose) {
    return(list(f1_scores = f1_scores, details = details))
  } else {
    return(f1_scores)
  }
}
