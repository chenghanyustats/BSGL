cv_bsgl <- function(y, X, coords,
                    L_range = c(16, 25, 36, 49, 81),
                    a_lambda_range = c(1, 2, 5),
                    b_lambda_range = c(0.01, 0.1, 1),
                    niter = 1000, nburn = 200,
                    k_folds = 5, verbose = TRUE) {

  n <- length(y)
  fold_indices <- sample(rep(1:k_folds, length.out = n))

  param_grid <- expand.grid(
    L = L_range,
    a_lambda = a_lambda_range,
    b_lambda = b_lambda_range
  )

  n_params <- nrow(param_grid)
  cv_results <- data.frame(
    L = param_grid$L,
    a_lambda = param_grid$a_lambda,
    b_lambda = param_grid$b_lambda,
    mean_mspe = NA,
    sd_mspe = NA
  )

  if (verbose) {
    cat("Starting", k_folds, "-fold CV with", n_params, "parameter combinations...\n")
  }

  for (i in 1:n_params) {
    L <- param_grid$L[i]
    a_lambda <- param_grid$a_lambda[i]
    b_lambda <- param_grid$b_lambda[i]

    if (verbose) {
      cat(sprintf("Testing params %d/%d: L=%d, a_lambda=%.2f, b_lambda=%.2f\n",
                  i, n_params, L, a_lambda, b_lambda))
    }

    fold_mspe <- numeric(k_folds)

    for (fold in 1:k_folds) {
      test_idx <- which(fold_indices == fold)
      train_idx <- which(fold_indices != fold)

      y_train <- y[train_idx]
      X_train <- X[train_idx, , drop = FALSE]
      coords_train <- coords[train_idx, , drop = FALSE]

      y_test <- y[test_idx]
      X_test <- X[test_idx, , drop = FALSE]
      coords_test <- coords[test_idx, , drop = FALSE]

      tryCatch({
        model <- fit_bsgl(
          y = y_train,
          X = X_train,
          coords = coords_train,
          L = L,
          niter = niter,
          nburn = nburn,
          a_lam = a_lambda,
          b_lam = b_lambda,
          verbose = FALSE
        )

        pred <- pred_bsgl(model, X_test, coords_test)
        fold_mspe[fold] <- mean((y_test - pred$mean)^2)

      }, error = function(e) {
        if (verbose) cat("Error in fold", fold, ":", e$message, "\n")
        fold_mspe[fold] <- Inf
      })
    }

    cv_results$mean_mspe[i] <- mean(fold_mspe)
    cv_results$sd_mspe[i] <- sd(fold_mspe)

    if (verbose) {
      cat(sprintf("  Mean MSPE: %.4f (+/-%.4f)\n",
                  cv_results$mean_mspe[i], cv_results$sd_mspe[i]))
    }
  }

  best_idx <- which.min(cv_results$mean_mspe)
  best_params <- cv_results[best_idx, ]

  if (verbose) {
    cat("\n=== CV Results ===\n")
    print(cv_results)
    cat("\nBest parameters:\n")
    print(best_params)
  }

  list(
    cv_results = cv_results,
    best_params = best_params,
    best_L = best_params$L,
    best_a_lambda = best_params$a_lambda,
    best_b_lambda = best_params$b_lambda
  )
}

fit_final_bsgl <- function(y, X, coords, cv_results, niter = 3000, nburn = 500) {
  best_params <- cv_results$best_params

  cat("Fitting final model with best parameters:\n")
  cat("L =", best_params$L, "\n")
  cat("a_lambda =", best_params$a_lambda, "\n")
  cat("b_lambda =", best_params$b_lambda, "\n")

  final_model <- fit_bsgl(
    y = y,
    X = X,
    coords = coords,
    L = best_params$L,
    niter = niter,
    nburn = nburn,
    a_lam = best_params$a_lambda,
    b_lam = best_params$b_lambda,
    verbose = TRUE
  )

  final_model
}
