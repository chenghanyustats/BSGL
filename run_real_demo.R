# run_real_demo.R - Package-based MODIS EVI reviewer example.
# The script installs the local package if needed, loads the prepared MODIS EVI
# subset, fits the four comparison methods, and writes one comparison table.

cat("===============================================\n")
cat("BSGL MODIS EVI Example\n")
cat("===============================================\n\n")

cat("1. PACKAGE: Installing/loading BSGL\n")
cat("   -------------------------------\n")

force_install <- tolower(Sys.getenv("BSGL_DEMO_INSTALL", unset = "true")) %in%
  c("1", "true", "yes")

if (force_install || !requireNamespace("BSGL", quietly = TRUE)) {
  cat("   Installing local package from repository root...\n")
  install.packages(".", repos = NULL, type = "source")
}

suppressPackageStartupMessages(library(BSGL))
cat("   Loaded package: BSGL\n")

metric_row <- function(method, observed, predicted) {
  data.frame(
    Method = method,
    Test_MSE = mean((observed - predicted)^2),
    Test_MAE = mean(abs(observed - predicted)),
    row.names = NULL
  )
}

na_metric_row <- function(method) {
  data.frame(
    Method = method,
    Test_MSE = NA_real_,
    Test_MAE = NA_real_,
    row.names = NULL
  )
}

quiet_eval <- function(expr) {
  out_file <- tempfile()
  msg_file <- tempfile()
  out_con <- file(out_file, open = "wt")
  msg_con <- file(msg_file, open = "wt")
  sink(out_con)
  sink(msg_con, type = "message")
  on.exit({
    while (sink.number(type = "message") > 0) sink(type = "message")
    while (sink.number() > 0) sink()
    close(out_con)
    close(msg_con)
    unlink(c(out_file, msg_file))
  }, add = TRUE)
  force(expr)
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
  list(train = train_dummy, test = test_dummy, reference = ref_level,
       levels = train_levels)
}

cat("\n2. DATA: Loading prepared MODIS EVI subset\n")
cat("   --------------------------------------\n")

data_file <- file.path("data", "data_cleaned_small_expanded.RData")
if (!file.exists(data_file)) {
  stop("Prepared MODIS data file is required at: ", data_file)
}

load(data_file)
real_data <- data_cleaned_small

sample_size <- as.integer(Sys.getenv("REAL_DEMO_SAMPLE_SIZE", "1000"))
sample_seed <- as.integer(Sys.getenv("REAL_DEMO_SAMPLE_SEED", "20260531"))
split_seed <- as.integer(Sys.getenv("REAL_DEMO_SPLIT_SEED", "456"))

set.seed(sample_seed)
idx <- if (nrow(real_data) > sample_size) {
  sample(nrow(real_data), sample_size)
} else {
  seq_len(nrow(real_data))
}
sample_data <- real_data[idx, ]

set.seed(split_seed)
test_indices <- sample(seq_len(nrow(sample_data)),
                       size = floor(0.2 * nrow(sample_data)))
train_indices <- setdiff(seq_len(nrow(sample_data)), test_indices)

boundary_indices <- unique(c(
  which(sample_data$scaled_x == min(sample_data$scaled_x)),
  which(sample_data$scaled_x == max(sample_data$scaled_x)),
  which(sample_data$scaled_y == min(sample_data$scaled_y)),
  which(sample_data$scaled_y == max(sample_data$scaled_y))
))
test_indices <- setdiff(test_indices, boundary_indices)
train_indices <- sort(unique(c(train_indices, boundary_indices)))

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

train_y <- log(sample_data$EVI[train_indices] + 1)
test_y <- log(sample_data$EVI[test_indices] + 1)

scaled_cont <- scale01_train_test(
  as.matrix(sample_data[train_indices, continuous_vars]),
  as.matrix(sample_data[test_indices, continuous_vars])
)
lc <- make_lc_dummy_design(
  sample_data$LC_Type4[train_indices],
  sample_data$LC_Type4[test_indices]
)

train_X <- cbind(scaled_cont$train, lc$train)
test_X <- cbind(scaled_cont$test, lc$test)
train_coords <- as.matrix(sample_data[train_indices, c("scaled_x", "scaled_y")])
test_coords <- as.matrix(sample_data[test_indices, c("scaled_x", "scaled_y")])

cat("   Training data: n =", length(train_y),
    ", conceptual m = 10, design columns =", ncol(train_X), "\n")
cat("   Test data: n =", length(test_y), "\n")

cat("\n3. MODELS: Fitting package methods\n")
cat("   -------------------------------\n")

L <- as.integer(Sys.getenv("REAL_DEMO_L", "25"))
niter <- as.integer(Sys.getenv("REAL_DEMO_NITER", "600"))
nburn <- as.integer(Sys.getenv("REAL_DEMO_NBURN", "100"))
a_lam <- as.numeric(Sys.getenv("REAL_DEMO_A_LAM", "20"))
b_lam <- as.numeric(Sys.getenv("REAL_DEMO_B_LAM", "1"))
run_mgwr <- tolower(Sys.getenv("REAL_DEMO_RUN_MGWR", unset = "false")) %in%
  c("1", "true", "yes")
mgwr_max_train <- as.integer(Sys.getenv("REAL_DEMO_MGWR_MAX_TRAIN", "1000"))
mgwr_skip_message <- paste(
  "MGWR is skipped by default for the real-data smoke test because",
  "GWmodel::gwr.multiscale can be slow on real-data subsets.",
  "Set REAL_DEMO_RUN_MGWR=true to attempt the optional local MGWR fit."
)

cat("   L =", L, "| niter =", niter, "| nburn =", nburn,
    "| a_lam =", a_lam, "| b_lam =", b_lam, "\n\n")

set.seed(20260531)

cat("   Fitting BSGL...\n")
bsgl_fit <- fit_bsgl(
  y = train_y,
  X = train_X,
  coords = train_coords,
  L = L,
  niter = niter,
  nburn = nburn,
  a_lam = a_lam,
  b_lam = b_lam,
  verbose = FALSE
)
bsgl_pred <- pred_bsgl(bsgl_fit, test_X, test_coords)$mean

cat("   Fitting GSVC...\n")
gsvc_fit <- fit_gs(
  y = train_y,
  X = train_X,
  coords = train_coords,
  L = L,
  niter = niter,
  nburn = nburn,
  verbose = FALSE
)
gsvc_pred <- pred_gs(gsvc_fit, test_X, test_coords)$mean

cat("   Fitting GGP-GAM...\n")
gam_fit <- fit_gam(train_X, train_y, train_coords)
gam_pred <- pred_gam(gam_fit, test_X, test_coords)

if (run_mgwr && length(train_y) > mgwr_max_train) {
  cat("   MGWR note: optional local fit skipped because training n = ",
      length(train_y), " exceeds REAL_DEMO_MGWR_MAX_TRAIN = ",
      mgwr_max_train, ". Use a smaller REAL_DEMO_SAMPLE_SIZE or increase ",
      "REAL_DEMO_MGWR_MAX_TRAIN to attempt it.\n", sep = "")
  mgwr_row <- na_metric_row("MGWR")
} else if (run_mgwr) {
  cat("   Fitting MGWR...\n")
  mgwr_row <- tryCatch(
    {
      mgwr_var_n <- ncol(train_X) + 1L
      mgwr_bws0 <- rep(max(length(train_y) - 1L, 2L), mgwr_var_n)
      mgwr_fit <- fit_mgwr(
        train_X,
        train_y,
        train_coords,
        bws0 = mgwr_bws0,
        bw.seled = rep(TRUE, mgwr_var_n),
        hatmatrix = FALSE
      )
      mgwr_pred <- pred_mgwr(mgwr_fit, test_X, test_coords)
      metric_row("MGWR", test_y, mgwr_pred)
    },
    error = function(e) {
      cat("   MGWR note: optional local comparison recorded as NA.\n")
      na_metric_row("MGWR")
    }
  )
} else {
  cat("   ", mgwr_skip_message, "\n", sep = "")
  mgwr_row <- na_metric_row("MGWR")
}

cat("\n4. RESULTS: Test-set comparison\n")
cat("   ----------------------------\n")

comparison <- rbind(
  metric_row("BSGL", test_y, bsgl_pred),
  metric_row("GSVC", test_y, gsvc_pred),
  metric_row("GGP-GAM", test_y, gam_pred),
  mgwr_row
)

comparison$Test_MSE <- round(comparison$Test_MSE, 8)
comparison$Test_MAE <- round(comparison$Test_MAE, 8)

print(comparison, row.names = FALSE)

dir.create("results", recursive = TRUE, showWarnings = FALSE)
out_file <- file.path("results", "real_demo_method_comparison.csv")
write.csv(comparison, out_file, row.names = FALSE)
cat("\n   Saved comparison table to:", out_file, "\n")

cat("\n===============================================\n")
cat("RESULTS WRITTEN SUCCESSFULLY\n")
cat("===============================================\n")
