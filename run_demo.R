# run_demo.R - Package-based reviewer example for the BSGL paper code.
# The script installs the local package if needed, loads BSGL, fits the four
# methods on the prepared simulation subset, and writes one comparison table.

cat("===============================================\n")
cat("BSGL Package Example\n")
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

cat("\n2. DATA: Loading prepared simulation subset\n")
cat("   ---------------------------------------\n")

data_path <- file.path("data", "simulation_demo", "n1000_p5")
train_file <- file.path(data_path, "train_data.RData")
test_file <- file.path(data_path, "test_data.RData")

if (!file.exists(train_file) || !file.exists(test_file)) {
  stop("Prepared data not found. Expected train/test files under: ", data_path)
}

load(train_file)
load(test_file)

cat("   Training data: n =", length(train_data$y),
    ", m =", ncol(train_data$X), "\n")
cat("   Test data: n =", length(test_data$y), "\n")

cat("\n3. MODELS: Fitting package methods\n")
cat("   -------------------------------\n")

L <- as.integer(Sys.getenv("DEMO_L", "25"))
niter <- as.integer(Sys.getenv("DEMO_NITER", "1000"))
nburn <- as.integer(Sys.getenv("DEMO_NBURN", "100"))
a_lam <- as.numeric(Sys.getenv("DEMO_A_LAM", "15"))
b_lam <- as.numeric(Sys.getenv("DEMO_B_LAM", "1"))

cat("   L =", L, "| niter =", niter, "| nburn =", nburn,
    "| a_lam =", a_lam, "| b_lam =", b_lam, "\n\n")

set.seed(123)

cat("   Fitting BSGL...\n")
bsgl_fit <- fit_bsgl(
  y = train_data$y,
  X = train_data$X,
  coords = train_data$coords,
  L = L,
  niter = niter,
  nburn = nburn,
  a_lam = a_lam,
  b_lam = b_lam,
  verbose = FALSE
)
bsgl_pred <- pred_bsgl(bsgl_fit, test_data$X, test_data$coords)$mean

cat("   Fitting GSVC...\n")
gsvc_fit <- fit_gs(
  y = train_data$y,
  X = train_data$X,
  coords = train_data$coords,
  L = L,
  niter = niter,
  nburn = nburn,
  verbose = FALSE
)
gsvc_pred <- pred_gs(gsvc_fit, test_data$X, test_data$coords)$mean

cat("   Fitting GGP-GAM...\n")
gam_fit <- fit_gam(train_data$X, train_data$y, train_data$coords)
gam_pred <- pred_gam(gam_fit, test_data$X, test_data$coords)

cat("   Fitting MGWR...\n")
mgwr_row <- tryCatch(
  {
    mgwr_fit <- fit_mgwr(train_data$X, train_data$y, train_data$coords)
    mgwr_pred <- pred_mgwr(mgwr_fit, test_data$X, test_data$coords)
    metric_row("MGWR", test_data$y, mgwr_pred)
  },
  error = function(e) {
    cat("   MGWR note: optional local comparison recorded as NA.\n")
    na_metric_row("MGWR")
  }
)

cat("\n4. RESULTS: Test-set comparison\n")
cat("   ----------------------------\n")

comparison <- rbind(
  metric_row("BSGL", test_data$y, bsgl_pred),
  metric_row("GSVC", test_data$y, gsvc_pred),
  metric_row("GGP-GAM", test_data$y, gam_pred),
  mgwr_row
)

comparison$Test_MSE <- round(comparison$Test_MSE, 6)
comparison$Test_MAE <- round(comparison$Test_MAE, 6)

print(comparison, row.names = FALSE)

dir.create("results", recursive = TRUE, showWarnings = FALSE)
out_file <- file.path("results", "demo_method_comparison.csv")
write.csv(comparison, out_file, row.names = FALSE)
cat("\n   Saved comparison table to:", out_file, "\n")

cat("\n===============================================\n")
cat("RESULTS WRITTEN SUCCESSFULLY\n")
cat("===============================================\n")
