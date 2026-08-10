# Real-data no-interaction SCP comparison: BSGL vs Gaussian SVC.
# LC_Type4 is retained as a categorical predictor using dummy variables.

source(file.path("R", "helpers_gplasso.R"))
source(file.path("R", "helpers_gpprior.R"))
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

  list(train = train_dummy, test = test_dummy, reference = ref_level, levels = train_levels)
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

plot_scp_compare <- function(scp_df, out_file) {
  p <- ggplot(scp_df, aes(x = reorder(label, SCP), y = SCP, fill = Method)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(values = c(BSGL = "#2F7D4F", `Gaussian SVC` = "#3B6EA8")) +
    labs(x = NULL, y = "SCP", fill = NULL) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank())
  ggsave(out_file, p, width = 9.2, height = 7.2, dpi = 300, bg = "white")
  invisible(p)
}

out_dir <- Sys.getenv("REAL_NOINT_OUT_DIR", unset = "results/real_data/nointeraction_scp_compare_n10000_rerun")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_size <- get_env_int("REAL_NOINT_SAMPLE_SIZE", 10000)
sample_seed <- get_env_int("REAL_NOINT_SAMPLE_SEED", 20260531)
split_seed <- get_env_int("REAL_NOINT_SPLIT_SEED", 456)
niter <- get_env_int("REAL_NOINT_NITER", 1000)
nburn <- get_env_int("REAL_NOINT_NBURN", 250)
fit_L <- get_env_int("REAL_NOINT_L", 25)
a_lam <- get_env_num("REAL_NOINT_A_LAMBDA", 20)
b_lam <- get_env_num("REAL_NOINT_B_LAMBDA", 1)
grid_n <- get_env_int("REAL_NOINT_GRID_N", 40)
save_fit <- get_env_bool("REAL_NOINT_SAVE_FIT", TRUE)

continuous_vars <- c(
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

load(file.path("data", "data_cleaned_small_expanded.RData"))
real_data <- data_cleaned_small

set.seed(sample_seed)
idx <- if (nrow(real_data) > sample_size) sample(nrow(real_data), sample_size) else seq_len(nrow(real_data))
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

train_y <- log(sample_data$EVI[train_indices] + 1)
test_y <- log(sample_data$EVI[test_indices] + 1)

scaled_cont <- scale01_train_test(
  as.matrix(sample_data[train_indices, continuous_vars]),
  as.matrix(sample_data[test_indices, continuous_vars])
)

lc <- make_lc_dummy_design(
  train_lc = sample_data$LC_Type4[train_indices],
  test_lc = sample_data$LC_Type4[test_indices]
)

train_X <- cbind(scaled_cont$train, lc$train)
test_X <- cbind(scaled_cont$test, lc$test)
var_names <- colnames(train_X)

train_coords <- as.matrix(sample_data[train_indices, c("scaled_x", "scaled_y")])
test_coords <- as.matrix(sample_data[test_indices, c("scaled_x", "scaled_y")])
colnames(train_coords) <- c("coords_x", "coords_y")
colnames(test_coords) <- c("coords_x", "coords_y")

grid_data <- make_grid_data(train_coords, train_X, grid_n)

cat("Real no-interaction SCP comparison\n")
cat("n sample:", nrow(sample_data), "\n")
cat("train/test:", nrow(train_X), "/", nrow(test_X), "\n")
cat("p:", ncol(train_X), "\n")
cat("LC reference level:", lc$reference, "\n\n")

cat("Fitting BSGL...\n")
bsgl_fit <- fit_bsgl(
  y = train_y,
  X = train_X,
  coords = train_coords,
  L = fit_L,
  niter = niter,
  nburn = nburn,
  a_lam = a_lam,
  b_lam = b_lam,
  verbose = TRUE
)

cat("Fitting Gaussian SVC...\n")
gs_fit <- fit_gs(
  y = train_y,
  X = train_X,
  coords = train_coords,
  L = fit_L,
  niter = niter,
  nburn = nburn,
  a_k = 2,
  b_k = 1,
  a_s = 2,
  b_s = 1,
  verbose = TRUE
)

bsgl_pred <- pred_bsgl(bsgl_fit, test_X, test_coords)
gs_pred <- pred_gs(gs_fit, test_X, test_coords)

pred_summary <- data.frame(
  Method = c("BSGL", "Gaussian SVC"),
  MSPE = c(
    mean((test_y - bsgl_pred$mean)^2),
    mean((test_y - gs_pred$mean)^2)
  ),
  Coverage = c(
    mean(test_y >= bsgl_pred$lower & test_y <= bsgl_pred$upper),
    mean(test_y >= gs_pred$lower & test_y <= gs_pred$upper)
  ),
  n_sample = nrow(sample_data),
  n_train = nrow(train_X),
  n_test = nrow(test_X),
  p = ncol(train_X),
  L = fit_L,
  niter = niter,
  nburn = nburn
)

bsgl_scp <- calc_scp_bsgl(bsgl_fit, grid_data, ci_level = 0.95)
gs_scp <- calc_scp_gs(gs_fit, grid_data, ci_level = 0.95)
names(bsgl_scp) <- var_names
names(gs_scp) <- var_names

labels <- setNames(gsub("_", " ", var_names), var_names)
labels[continuous_vars] <- c(
  "Red Reflectance", "NIR Reflectance", "Blue Reflectance", "MIR Reflectance",
  "GPP", "LE", "View Zenith Angle", "Sun Zenith Angle", "Relative Azimuth Angle"
)
labels[colnames(lc$train)] <- paste0("LC Type4 = ", sub("^LC_Type4_", "", colnames(lc$train)))

scp_df <- rbind(
  data.frame(Method = "BSGL", variable = var_names, label = unname(labels[var_names]), SCP = as.numeric(bsgl_scp[var_names])),
  data.frame(Method = "Gaussian SVC", variable = var_names, label = unname(labels[var_names]), SCP = as.numeric(gs_scp[var_names]))
)

lc_vars <- colnames(lc$train)
group_summary <- data.frame(
  Method = rep(c("BSGL", "Gaussian SVC"), each = 3),
  group = rep(c("LC_Type4_max", "LC_Type4_mean", "Non_LC_top_mean"), times = 2),
  SCP = c(
    max(bsgl_scp[lc_vars]), mean(bsgl_scp[lc_vars]), mean(sort(bsgl_scp[continuous_vars], decreasing = TRUE)[1:4]),
    max(gs_scp[lc_vars]), mean(gs_scp[lc_vars]), mean(sort(gs_scp[continuous_vars], decreasing = TRUE)[1:4])
  )
)

write.csv(pred_summary, file.path(out_dir, "real_nointeraction_prediction_summary.csv"), row.names = FALSE)
write.csv(scp_df, file.path(out_dir, "real_nointeraction_scp_by_column.csv"), row.names = FALSE)
write.csv(group_summary, file.path(out_dir, "real_nointeraction_scp_group_summary.csv"), row.names = FALSE)
plot_scp_compare(scp_df, file.path(out_dir, "real_nointeraction_scp_compare.png"))

if (save_fit) {
  saveRDS(
    list(
      bsgl_fit = bsgl_fit,
      gs_fit = gs_fit,
      prediction_summary = pred_summary,
      scp_by_column = scp_df,
      group_summary = group_summary,
      lc_reference = lc$reference,
      variables = var_names
    ),
    file.path(out_dir, "real_nointeraction_scp_compare_fit.rds")
  )
}

cat("\nPrediction summary\n")
print(pred_summary)
cat("\nSCP by column\n")
print(scp_df)
cat("\nGroup summary\n")
print(group_summary)
cat("\nOutputs saved to:", normalizePath(out_dir), "\n")
