# Real EVI interaction validation.
#
# This keeps the observed response y = log(EVI + 1) and the original ten main
# predictors. It only augments the real-data design matrix with one hypothesized
# useful interaction and one hypothesized null interaction.

source(file.path("R", "helpers_gplasso.R"))
source(file.path("R", "help_inclu.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

get_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "") return(default)
  as.integer(value)
}

get_env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "") return(default)
  as.numeric(value)
}

get_env_bool <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "") return(default)
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

scale01_train_test <- function(train_x, test_x) {
  train_scaled <- train_x
  test_scaled <- test_x

  for (j in seq_len(ncol(train_x))) {
    lo <- min(train_x[, j], na.rm = TRUE)
    hi <- max(train_x[, j], na.rm = TRUE)
    denom <- hi - lo

    if (!is.finite(denom) || denom == 0) {
      train_scaled[, j] <- 0
      test_scaled[, j] <- 0
    } else {
      train_scaled[, j] <- (train_x[, j] - lo) / denom
      test_scaled[, j] <- (test_x[, j] - lo) / denom
    }
  }

  colnames(train_scaled) <- colnames(train_x)
  colnames(test_scaled) <- colnames(test_x)
  list(train = train_scaled, test = test_scaled)
}

scale_one_by_train <- function(train_x, test_x) {
  lo <- min(train_x, na.rm = TRUE)
  hi <- max(train_x, na.rm = TRUE)
  denom <- hi - lo

  if (!is.finite(denom) || denom == 0) {
    return(list(train = rep(0, length(train_x)), test = rep(0, length(test_x))))
  }

  list(train = (train_x - lo) / denom, test = (test_x - lo) / denom)
}

make_grid_data <- function(train_coords, train_X, grid_n) {
  grid_coords <- expand.grid(
    x = seq(min(train_coords[, 1]), max(train_coords[, 1]), length.out = grid_n),
    y = seq(min(train_coords[, 2]), max(train_coords[, 2]), length.out = grid_n)
  )
  grid_X <- matrix(
    rep(apply(train_X, 2, median, na.rm = TRUE), each = nrow(grid_coords)),
    nrow = nrow(grid_coords),
    ncol = ncol(train_X)
  )
  colnames(grid_X) <- colnames(train_X)

  list(grid_points = as.matrix(grid_coords), X = grid_X)
}

surface_ci_bsgl <- function(fit, grid_data, ci_level = 0.95) {
  grid_coords <- as.matrix(grid_data$grid_points)
  Psi_grid <- make_basis(grid_coords, fit$L, saved_knots = fit$saved_knots)
  n_grid <- nrow(grid_coords)
  n_samples <- dim(fit$beta)[3]

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
      beta_samples[, s] <- Psi_grid %*% fit$beta[, j, s]
    }
    mean_beta[, j] <- rowMeans(beta_samples)
    lower_beta[, j] <- apply(beta_samples, 1, quantile, lower_prob)
    upper_beta[, j] <- apply(beta_samples, 1, quantile, upper_prob)
    excludes0[, j] <- (lower_beta[, j] > 0) | (upper_beta[, j] < 0)
  }

  list(mean = mean_beta, lower = lower_beta, upper = upper_beta, excludes0 = excludes0)
}

plot_scp_bar <- function(scp_df, out_file) {
  p <- ggplot(scp_df, aes(x = reorder(label, SCP), y = SCP, fill = term_type)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(values = c(Main = "#3B6EA8", Interaction = "#2F7D4F")) +
    labs(x = NULL, y = "SCP", fill = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank())

  ggsave(out_file, p, width = 8.5, height = 6.2, dpi = 300, bg = "white")
  invisible(p)
}

plot_interaction_surfaces <- function(surface_df, out_file) {
  p <- ggplot(surface_df, aes(x, y, fill = beta_mean)) +
    geom_raster() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    coord_fixed() +
    facet_wrap(~ label, nrow = 1) +
    labs(x = NULL, y = NULL, fill = "Mean beta") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank())

  ggsave(out_file, p, width = 9.5, height = 4.2, dpi = 300, bg = "white")
  invisible(p)
}

fit_eval_bsgl <- function(y_train, X_train, coords_train, y_test, X_test, coords_test,
                          L, niter, nburn, a_lam, b_lam) {
  fit <- fit_bsgl(
    y = y_train,
    X = X_train,
    coords = coords_train,
    L = L,
    niter = niter,
    nburn = nburn,
    a_lam = a_lam,
    b_lam = b_lam,
    verbose = TRUE
  )
  pred <- pred_bsgl(fit, X_test, coords_test)
  list(
    fit = fit,
    pred = pred,
    mspe = mean((y_test - pred$mean)^2),
    coverage = mean(y_test >= pred$lower & y_test <= pred$upper)
  )
}

make_lc_dummy_design <- function(train_lc, test_lc) {
  train_levels <- sort(unique(train_lc))
  ref_level <- train_levels[1]
  dummy_levels <- train_levels[train_levels != ref_level]

  train_dummy <- sapply(dummy_levels, function(level) as.numeric(train_lc == level))
  test_dummy <- sapply(dummy_levels, function(level) as.numeric(test_lc == level))

  if (length(dummy_levels) == 1) {
    train_dummy <- matrix(train_dummy, ncol = 1)
    test_dummy <- matrix(test_dummy, ncol = 1)
  }

  colnames(train_dummy) <- paste0("LC_Type4_", dummy_levels)
  colnames(test_dummy) <- paste0("LC_Type4_", dummy_levels)

  list(
    train = train_dummy,
    test = test_dummy,
    reference = ref_level,
    levels = train_levels
  )
}

out_dir <- Sys.getenv("REAL_EVI_INT_OUT_DIR", unset = "results/interaction/evi_interaction_validation_rerun")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_size <- get_env_int("REAL_EVI_INT_SAMPLE_SIZE", 1000)
sample_seed <- get_env_int("REAL_EVI_INT_SAMPLE_SEED", 20260531)
split_seed <- get_env_int("REAL_EVI_INT_SPLIT_SEED", 456)
niter <- get_env_int("REAL_EVI_INT_NITER", 1000)
nburn <- get_env_int("REAL_EVI_INT_NBURN", 250)
fit_L <- get_env_int("REAL_EVI_INT_L", 25)
a_lam <- get_env_num("REAL_EVI_INT_A_LAMBDA", 20)
b_lam <- get_env_num("REAL_EVI_INT_B_LAMBDA", 1)
grid_n <- get_env_int("REAL_EVI_INT_GRID_N", 35)
run_baseline <- get_env_bool("REAL_EVI_INT_RUN_BASELINE", TRUE)
save_fit <- get_env_bool("REAL_EVI_INT_SAVE_FIT", TRUE)
lc_mode <- Sys.getenv("REAL_EVI_INT_LC_MODE", unset = "numeric")
if (!(lc_mode %in% c("numeric", "drop", "dummy"))) {
  stop("REAL_EVI_INT_LC_MODE must be one of: numeric, drop, dummy.")
}

continuous_main_vars <- c(
  "red_reflectance",
  "NIR_reflectance",
  "blue_reflectance",
  "MIR_reflectance",
  "GPP",
  "LE",
  "view_zenith_angle",
  "sun_zenith_angle",
  "relative_azimuth_angle"
)
main_vars <- continuous_main_vars
if (lc_mode == "numeric") {
  main_vars <- c(main_vars, "LC_Type4")
}

useful_anchor_var <- Sys.getenv("REAL_EVI_INT_USEFUL_ANCHOR_VAR", unset = "NIR_reflectance")
useful_base_var <- Sys.getenv("REAL_EVI_INT_USEFUL_VAR", unset = "GPP")
null_anchor_var <- Sys.getenv("REAL_EVI_INT_NULL_ANCHOR_VAR", unset = "NIR_reflectance")
null_base_var <- Sys.getenv("REAL_EVI_INT_NULL_VAR", unset = "relative_azimuth_angle")
allowed_null_vars <- c(continuous_main_vars, if (lc_mode == "numeric") "LC_Type4")
if (!(useful_anchor_var %in% allowed_null_vars)) {
  stop("REAL_EVI_INT_USEFUL_ANCHOR_VAR must be one of the main predictor names.")
}
if (!(useful_base_var %in% allowed_null_vars)) {
  stop("REAL_EVI_INT_USEFUL_VAR must be one of the main predictor names.")
}
if (!(null_base_var %in% allowed_null_vars)) {
  stop("REAL_EVI_INT_NULL_VAR must be one of the main predictor names.")
}
if (useful_anchor_var == useful_base_var) {
  stop("REAL_EVI_INT_USEFUL_ANCHOR_VAR and REAL_EVI_INT_USEFUL_VAR must be different.")
}
if (null_anchor_var == null_base_var) {
  stop("REAL_EVI_INT_NULL_ANCHOR_VAR and REAL_EVI_INT_NULL_VAR must be different.")
}
if (setequal(c(null_anchor_var, null_base_var), c(useful_anchor_var, useful_base_var))) {
  stop("The null interaction cannot duplicate the useful interaction.")
}
useful_interaction <- paste0(useful_anchor_var, "_x_", useful_base_var)
null_interaction <- paste0(null_anchor_var, "_x_", null_base_var)
augmented_vars <- c(main_vars, useful_interaction, null_interaction)

labels <- c(
  red_reflectance = "Red Reflectance",
  NIR_reflectance = "NIR Reflectance",
  blue_reflectance = "Blue Reflectance",
  MIR_reflectance = "MIR Reflectance",
  GPP = "GPP",
  LE = "LE",
  view_zenith_angle = "View Zenith Angle",
  sun_zenith_angle = "Sun Zenith Angle",
  relative_azimuth_angle = "Relative Azimuth Angle",
  LC_Type4 = "LC Type4",
  NIR_reflectance_x_GPP = "NIR x GPP",
  NIR_reflectance_x_red_reflectance = "NIR x Red Reflectance",
  NIR_reflectance_x_blue_reflectance = "NIR x Blue Reflectance",
  NIR_reflectance_x_MIR_reflectance = "NIR x MIR Reflectance",
  NIR_reflectance_x_LE = "NIR x LE",
  NIR_reflectance_x_view_zenith_angle = "NIR x View Zenith",
  NIR_reflectance_x_sun_zenith_angle = "NIR x Sun Zenith",
  NIR_reflectance_x_relative_azimuth_angle = "NIR x Relative Azimuth",
  NIR_reflectance_x_LC_Type4 = "NIR x LC Type4",
  blue_reflectance_x_relative_azimuth_angle = "Blue x Relative Azimuth",
  blue_reflectance_x_view_zenith_angle = "Blue x View Zenith",
  MIR_reflectance_x_relative_azimuth_angle = "MIR x Relative Azimuth",
  MIR_reflectance_x_view_zenith_angle = "MIR x View Zenith",
  view_zenith_angle_x_relative_azimuth_angle = "View Zenith x Relative Azimuth",
  sun_zenith_angle_x_relative_azimuth_angle = "Sun Zenith x Relative Azimuth",
  GPP_x_relative_azimuth_angle = "GPP x Relative Azimuth",
  LE_x_relative_azimuth_angle = "LE x Relative Azimuth"
)

load(file.path("data", "data_cleaned_small_expanded.RData"))
real_data <- data_cleaned_small

set.seed(sample_seed)
idx <- if (nrow(real_data) > sample_size) {
  sample(nrow(real_data), sample_size)
} else {
  seq_len(nrow(real_data))
}
sample_data <- real_data[idx, ]

set.seed(split_seed)
test_indices <- sample(seq_len(nrow(sample_data)), size = floor(0.2 * nrow(sample_data)))
train_indices <- setdiff(seq_len(nrow(sample_data)), test_indices)

boundary_indices <- unique(c(
  which(sample_data$scaled_x == min(sample_data$scaled_x)),
  which(sample_data$scaled_x == max(sample_data$scaled_x)),
  which(sample_data$scaled_y == min(sample_data$scaled_y)),
  which(sample_data$scaled_y == max(sample_data$scaled_y))
))
test_indices <- setdiff(test_indices, boundary_indices)
train_indices <- sort(unique(c(train_indices, boundary_indices)))

y <- log(sample_data$EVI + 1)
train_y <- y[train_indices]
test_y <- y[test_indices]

train_X_main_raw <- as.matrix(sample_data[train_indices, main_vars])
test_X_main_raw <- as.matrix(sample_data[test_indices, main_vars])
scaled_main <- scale01_train_test(train_X_main_raw, test_X_main_raw)
train_X_main <- scaled_main$train
test_X_main <- scaled_main$test

lc_meta <- NULL
if (lc_mode == "dummy") {
  lc_meta <- make_lc_dummy_design(
    train_lc = sample_data$LC_Type4[train_indices],
    test_lc = sample_data$LC_Type4[test_indices]
  )
  train_X_main <- cbind(train_X_main, lc_meta$train)
  test_X_main <- cbind(test_X_main, lc_meta$test)
  main_vars <- colnames(train_X_main)

  dummy_labels <- setNames(
    paste0("LC Type4 = ", sub("^LC_Type4_", "", colnames(lc_meta$train))),
    colnames(lc_meta$train)
  )
  labels <- c(labels, dummy_labels)
}
augmented_vars <- c(main_vars, useful_interaction, null_interaction)

train_center <- colMeans(
  train_X_main[, unique(c(useful_anchor_var, useful_base_var, null_anchor_var, null_base_var)), drop = FALSE]
)
train_useful_raw <- (train_X_main[, useful_anchor_var] - train_center[useful_anchor_var]) *
  (train_X_main[, useful_base_var] - train_center[useful_base_var])
test_useful_raw <- (test_X_main[, useful_anchor_var] - train_center[useful_anchor_var]) *
  (test_X_main[, useful_base_var] - train_center[useful_base_var])
train_null_raw <- (train_X_main[, null_anchor_var] - train_center[null_anchor_var]) *
  (train_X_main[, null_base_var] - train_center[null_base_var])
test_null_raw <- (test_X_main[, null_anchor_var] - train_center[null_anchor_var]) *
  (test_X_main[, null_base_var] - train_center[null_base_var])

useful_scaled <- scale_one_by_train(train_useful_raw, test_useful_raw)
null_scaled <- scale_one_by_train(train_null_raw, test_null_raw)

train_X_aug <- cbind(
  train_X_main,
  useful_scaled$train,
  null_scaled$train
)
colnames(train_X_aug)[ncol(train_X_aug) - 1] <- useful_interaction
colnames(train_X_aug)[ncol(train_X_aug)] <- null_interaction
test_X_aug <- cbind(
  test_X_main,
  useful_scaled$test,
  null_scaled$test
)
colnames(test_X_aug)[ncol(test_X_aug) - 1] <- useful_interaction
colnames(test_X_aug)[ncol(test_X_aug)] <- null_interaction

train_coords <- as.matrix(sample_data[train_indices, c("scaled_x", "scaled_y")])
test_coords <- as.matrix(sample_data[test_indices, c("scaled_x", "scaled_y")])
colnames(train_coords) <- c("coords_x", "coords_y")
colnames(test_coords) <- c("coords_x", "coords_y")

cat("Observed-EVI real interaction validation\n")
cat("n sample:", nrow(sample_data), "\n")
cat("train/test:", nrow(train_X_aug), "/", nrow(test_X_aug), "\n")
cat("LC handling:", lc_mode, "\n")
if (!is.null(lc_meta)) {
  cat("LC reference level:", lc_meta$reference, "\n")
}
cat("hypothesized useful interaction:", useful_interaction, "\n")
cat("hypothesized null interaction:", null_interaction, "\n")
cat("useful anchor variable:", useful_anchor_var, "\n")
cat("useful base variable:", useful_base_var, "\n")
cat("cor(useful anchor, useful interaction):", round(cor(train_X_aug[, useful_anchor_var], train_X_aug[, useful_interaction]), 3), "\n")
cat("cor(useful base, useful interaction):", round(cor(train_X_aug[, useful_base_var], train_X_aug[, useful_interaction]), 3), "\n")
cat("null anchor variable:", null_anchor_var, "\n")
cat("null base variable:", null_base_var, "\n")
cat("cor(null anchor, null interaction):", round(cor(train_X_aug[, null_anchor_var], train_X_aug[, null_interaction]), 3), "\n")
cat("cor(null base, null interaction):", round(cor(train_X_aug[, null_base_var], train_X_aug[, null_interaction]), 3), "\n\n")

aug_res <- fit_eval_bsgl(
  y_train = train_y,
  X_train = train_X_aug,
  coords_train = train_coords,
  y_test = test_y,
  X_test = test_X_aug,
  coords_test = test_coords,
  L = fit_L,
  niter = niter,
  nburn = nburn,
  a_lam = a_lam,
  b_lam = b_lam
)

pred_summary <- data.frame(
  Method = "BSGL + interactions",
  MSPE = aug_res$mspe,
  Coverage = aug_res$coverage,
  n_sample = nrow(sample_data),
  n_train = nrow(train_X_aug),
  n_test = nrow(test_X_aug),
  p = ncol(train_X_aug),
  L = fit_L,
  niter = niter,
  nburn = nburn
)

baseline_res <- NULL
if (run_baseline) {
  cat("\nFitting baseline BSGL without interactions...\n")
  baseline_res <- fit_eval_bsgl(
    y_train = train_y,
    X_train = train_X_main,
    coords_train = train_coords,
    y_test = test_y,
    X_test = test_X_main,
    coords_test = test_coords,
    L = fit_L,
    niter = niter,
    nburn = nburn,
    a_lam = a_lam,
    b_lam = b_lam
  )
  pred_summary <- rbind(
    data.frame(
      Method = "BSGL baseline",
      MSPE = baseline_res$mspe,
      Coverage = baseline_res$coverage,
      n_sample = nrow(sample_data),
      n_train = nrow(train_X_aug),
      n_test = nrow(test_X_aug),
      p = ncol(train_X_main),
      L = fit_L,
      niter = niter,
      nburn = nburn
    ),
    pred_summary
  )
}

grid_data <- make_grid_data(train_coords, train_X_aug, grid_n)
surface_ci <- surface_ci_bsgl(aug_res$fit, grid_data, ci_level = 0.95)
colnames(surface_ci$mean) <- augmented_vars
colnames(surface_ci$excludes0) <- augmented_vars
scp <- colMeans(surface_ci$excludes0)
label_values <- unname(labels[augmented_vars])
missing_labels <- is.na(label_values)
label_values[missing_labels] <- gsub("_", " ", augmented_vars[missing_labels])

scp_df <- data.frame(
  variable = augmented_vars,
  label = label_values,
  term_type = ifelse(augmented_vars %in% main_vars, "Main", "Interaction"),
  SCP = as.numeric(scp[augmented_vars]),
  mean_norm = sum_bsgl(aug_res$fit)$mean_norm,
  row.names = NULL
)
scp_df <- scp_df[order(scp_df$SCP, decreasing = TRUE), ]

surface_df <- data.frame(
  x = grid_data$grid_points[, 1],
  y = grid_data$grid_points[, 2],
  useful = surface_ci$mean[, useful_interaction],
  null = surface_ci$mean[, null_interaction]
)
surface_long <- rbind(
  data.frame(x = surface_df$x, y = surface_df$y, beta_mean = surface_df$useful, label = unname(labels[useful_interaction])),
  data.frame(x = surface_df$x, y = surface_df$y, beta_mean = surface_df$null, label = unname(labels[null_interaction]))
)

write.csv(pred_summary, file.path(out_dir, "evi_interaction_prediction_summary.csv"), row.names = FALSE)
write.csv(scp_df, file.path(out_dir, "evi_interaction_scp_summary.csv"), row.names = FALSE)
write.csv(surface_df, file.path(out_dir, "evi_interaction_surfaces.csv"), row.names = FALSE)

plot_scp_bar(scp_df, file.path(out_dir, "evi_interaction_scp_bar.png"))
plot_interaction_surfaces(surface_long, file.path(out_dir, "evi_interaction_surfaces.png"))

if (save_fit) {
  saveRDS(
    list(
      augmented_fit = aug_res$fit,
      baseline_fit = if (!is.null(baseline_res)) baseline_res$fit else NULL,
      prediction_summary = pred_summary,
      scp_summary = scp_df,
      grid_data = grid_data,
      surface_ci = surface_ci,
      variables = augmented_vars
    ),
    file.path(out_dir, "evi_interaction_validation_fit.rds")
  )
}

cat("\nPrediction summary\n")
print(pred_summary)
cat("\nSCP summary\n")
print(scp_df)
cat("\nOutputs saved to:", normalizePath(out_dir), "\n")
