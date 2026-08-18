# ============================================================================
# BSGL core runtime smoke tests
#
# PURPOSE
#   Execute shortened copies of the main full-reproduction producer workflows
#   to detect runtime-only path/function/object mismatches without overwriting
#   retained manuscript results.
#
# PLACEMENT
#   Save as:
#       scripts/05_runtime_smoke_tests.R
#
# RUN FROM REPOSITORY ROOT
#       source("scripts/05_runtime_smoke_tests.R")
#
# IMPORTANT
#   - Original scripts are NOT edited.
#   - All smoke-test output is written under tempdir().
#   - Retained results/ and figures/ are NOT overwritten.
#   - Short MCMC settings are for software validation only, not inference.
# ============================================================================

cat("============================================================\n")
cat("BSGL core runtime smoke tests\n")
cat("============================================================\n\n")

if (!file.exists("DESCRIPTION") ||
    !dir.exists("R") ||
    !dir.exists(file.path("scripts", "full_reproduction"))) {
  stop("Run this script from the BSGL repository root.")
}

project_root <- normalizePath(".", mustWork = TRUE)
oldwd <- getwd()
on.exit(setwd(oldwd), add = TRUE)
setwd(project_root)

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
smoke_root <- file.path(
  tempdir(),
  paste0("BSGL_runtime_smoke_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(smoke_root, recursive = TRUE, showWarnings = FALSE)

cat("Repository root:", project_root, "\n")
cat("Smoke output root:", smoke_root, "\n\n")

rscript <- file.path(R.home("bin"), "Rscript")
if (!file.exists(rscript)) stop("Could not locate Rscript.")

# Keep a snapshot of git status; the smoke runner should not alter the repo.
git_status <- function() {
  out <- tryCatch(
    system2("git", c("status", "--short"), stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  )
  sort(out)
}
git_before <- git_status()

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
read_script <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

count_fixed <- function(text, pattern) {
  hits <- gregexpr(pattern, text, fixed = TRUE)[[1]]
  if (length(hits) == 1L && hits[1] == -1L) return(0L)
  length(hits)
}

replace_exact <- function(text, old, new, expected = 1L, label = old) {
  n <- count_fixed(text, old)
  if (!identical(n, as.integer(expected))) {
    stop(
      "Patch mismatch for ", label, ": expected ", expected,
      " occurrence(s), found ", n, "."
    )
  }
  gsub(old, new, text, fixed = TRUE)
}

r_path <- function(path) {
  # macOS/Linux paths are used here; normalize slashes defensively.
  gsub("\\\\", "/", normalizePath(path, mustWork = FALSE))
}

r_string <- function(path) {
  paste0('"', gsub('"', '\\"', r_path(path), fixed = TRUE), '"')
}

run_modified_script <- function(
    test_name,
    source_rel,
    patch_fun,
    expected_files = character(),
    env = character()) {

  cat("\n------------------------------------------------------------\n")
  cat("TEST:", test_name, "\n")
  cat("SOURCE:", source_rel, "\n")
  cat("------------------------------------------------------------\n")

  source_path <- file.path(project_root, source_rel)
  if (!file.exists(source_path)) {
    stop("Missing source script: ", source_rel)
  }

  text <- read_script(source_path)
  text <- patch_fun(text)

  work_dir <- file.path(smoke_root, "patched_scripts")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  temp_script <- file.path(
    work_dir,
    paste0(gsub("[^A-Za-z0-9_]+", "_", test_name), ".R")
  )
  writeLines(
    strsplit(text, "\n", fixed = TRUE)[[1]],
    temp_script,
    useBytes = TRUE
  )

  # Parse the patched copy before executing it.
  parse(file = temp_script)

  log_dir <- file.path(smoke_root, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(
    log_dir,
    paste0(gsub("[^A-Za-z0-9_]+", "_", test_name), ".log")
  )

  args <- c("--vanilla", shQuote(temp_script))
  status <- system2(
    rscript,
    args = args,
    stdout = log_file,
    stderr = log_file,
    env = env
  )
  if (is.null(status)) status <- 0L

  expected_ok <- if (length(expected_files)) {
    file.exists(expected_files)
  } else {
    logical()
  }

  if (status == 0L && (length(expected_ok) == 0L || all(expected_ok))) {
    cat("PASS:", test_name, "\n")
    if (length(expected_files)) {
      cat("Expected outputs found:\n")
      for (f in expected_files) cat("  ", f, "\n", sep = "")
    }
    return(data.frame(
      test = test_name,
      status = "PASS",
      exit_status = status,
      log = log_file,
      stringsAsFactors = FALSE
    ))
  }

  cat("FAIL:", test_name, "\n")
  cat("Exit status:", status, "\n")
  if (length(expected_files)) {
    missing <- expected_files[!expected_ok]
    if (length(missing)) {
      cat("Missing expected output(s):\n")
      for (f in missing) cat("  ", f, "\n", sep = "")
    }
  }

  if (file.exists(log_file)) {
    log_lines <- readLines(log_file, warn = FALSE, encoding = "UTF-8")
    cat("\nLast log lines:\n")
    cat(tail(log_lines, 60), sep = "\n")
    cat("\n")
  }

  data.frame(
    test = test_name,
    status = "FAIL",
    exit_status = status,
    log = log_file,
    stringsAsFactors = FALSE
  )
}

run_original_script <- function(
    test_name,
    source_rel,
    expected_files = character(),
    env = character()) {

  run_modified_script(
    test_name = test_name,
    source_rel = source_rel,
    patch_fun = identity,
    expected_files = expected_files,
    env = env
  )
}

results <- list()

# ===========================================================================
# 1. sim_10.R
#    One retained data combination, one replicate, one CV combination.
#    MGWR is skipped here and tested separately by rerun_mgwr.R.
# ===========================================================================
sim10_out <- file.path(smoke_root, "sim_10")
sim10_summary <- file.path(
  sim10_out, "summary", "rep10_summary_mean_sd.csv"
)
sim10_result <- file.path(
  sim10_out, "single_runs", "result_n1000_p5_rep01.rds"
)

results[[length(results) + 1L]] <- run_modified_script(
  test_name = "sim_10 core simulation",
  source_rel = file.path("scripts", "full_reproduction", "sim_10.R"),
  patch_fun = function(x) {
    x <- replace_exact(
      x,
      'out_dir <- "results/simulation/full_rep10_run"',
      paste0("out_dir <- ", r_string(sim10_out)),
      label = "sim_10 out_dir"
    )
    x <- replace_exact(x, "n_vec <- c(1000, 5000, 10000)",
                       "n_vec <- c(1000)", label = "sim_10 n_vec")
    x <- replace_exact(x, "p_vec <- c(5, 10, 20)",
                       "p_vec <- c(5)", label = "sim_10 p_vec")
    x <- replace_exact(x, "rep_vec <- 1:10",
                       "rep_vec <- 1L", label = "sim_10 rep_vec")
    x <- replace_exact(x, "cv_L_range <- c(25, 36)",
                       "cv_L_range <- c(25)", label = "sim_10 CV L")
    x <- replace_exact(x, "cv_a_lambda_range <- c(15, 30, 35)",
                       "cv_a_lambda_range <- c(15)", label = "sim_10 CV a")
    x <- replace_exact(x, "cv_b_lambda_range <- c(0.1, 1, 2)",
                       "cv_b_lambda_range <- c(1)", label = "sim_10 CV b")
    x <- replace_exact(x, "cv_niter <- 800",
                       "cv_niter <- 80", label = "sim_10 CV iterations")
    x <- replace_exact(x, "cv_nburn <- 100",
                       "cv_nburn <- 20", label = "sim_10 CV burn")
    x <- replace_exact(x, "cv_k_folds <- 5",
                       "cv_k_folds <- 2", label = "sim_10 CV folds")
    x <- replace_exact(x, "mcmc_niter <- 5000",
                       "mcmc_niter <- 120", label = "sim_10 MCMC iterations")
    x <- replace_exact(x, "mcmc_nburn <- 500",
                       "mcmc_nburn <- 20", label = "sim_10 MCMC burn")
    x <- replace_exact(x, "skip_mgwr <- FALSE",
                       "skip_mgwr <- TRUE", label = "sim_10 skip MGWR")
    x
  },
  expected_files = c(sim10_summary, sim10_result)
)

# ===========================================================================
# 2. rerun_mgwr.R
#    Uses the partial/core result generated above and exercises the MGWR
#    completion branch on n=1000, p=5, replicate 1.
# ===========================================================================
mgwr_result <- file.path(
  sim10_out, "single_runs", "result_full_n1000_p5_rep01.rds"
)

if (identical(results[[length(results)]]$status, "PASS")) {
  results[[length(results) + 1L]] <- run_modified_script(
    test_name = "rerun_mgwr completion",
    source_rel = file.path("scripts", "full_reproduction", "rerun_mgwr.R"),
    patch_fun = function(x) {
      x <- replace_exact(
        x,
        'out_dir <- "results/simulation/full_rep10_run"',
        paste0("out_dir <- ", r_string(sim10_out)),
        label = "rerun_mgwr out_dir"
      )
      x <- replace_exact(x, "n_vec <- c(1000, 5000, 10000)",
                         "n_vec <- c(1000)", label = "rerun_mgwr n_vec")
      x <- replace_exact(x, "p_vec <- c(5, 10, 20)",
                         "p_vec <- c(5)", label = "rerun_mgwr p_vec")
      x <- replace_exact(x, "rep_vec <- 1:10",
                         "rep_vec <- 1L", label = "rerun_mgwr rep_vec")
      x
    },
    expected_files = mgwr_result
  )
} else {
  cat("\nSKIP: rerun_mgwr completion because sim_10 failed.\n")
  results[[length(results) + 1L]] <- data.frame(
    test = "rerun_mgwr completion",
    status = "SKIP",
    exit_status = NA_integer_,
    log = NA_character_,
    stringsAsFactors = FALSE
  )
}

# ===========================================================================
# 3. sim_gpsvc.R
#    One retained independent-data combination.
# ===========================================================================
gpsvc_out <- file.path(smoke_root, "sim_gpsvc")
gpsvc_csv <- file.path(gpsvc_out, "gaussian_svc_all_results.csv")

results[[length(results) + 1L]] <- run_modified_script(
  test_name = "sim_gpsvc GSVC-only",
  source_rel = file.path("scripts", "full_reproduction", "sim_gpsvc.R"),
  patch_fun = function(x) {
    # This assumes the repository-audit patch has already aligned the paper grid.
    x <- replace_exact(x, "n_list <- c(1000, 5000, 10000)",
                       "n_list <- c(1000)", label = "sim_gpsvc n_list")
    x <- replace_exact(x, "p_list <- c(5, 10, 20)",
                       "p_list <- c(5)", label = "sim_gpsvc p_list")
    x <- replace_exact(
      x,
      'gs_out_dir <- file.path("results", "simulation", "gs_independent_runs")',
      paste0("gs_out_dir <- ", r_string(gpsvc_out)),
      label = "sim_gpsvc output directory"
    )
    x <- replace_exact(x, "niter = 5000,",
                       "niter = 120,", label = "sim_gpsvc iterations")
    x <- replace_exact(x, "nburn = 500,",
                       "nburn = 20,", label = "sim_gpsvc burn")
    x
  },
  expected_files = gpsvc_csv
)

# ===========================================================================
# 4. run_abrupt_boundary_mainsettings.R
#    Smaller generated data and short BSGL/GSVC fits.
# ===========================================================================
abrupt_out <- file.path(smoke_root, "abrupt_boundary")
abrupt_fit <- file.path(abrupt_out, "abrupt_boundary_mainsettings_fit.rds")

results[[length(results) + 1L]] <- run_modified_script(
  test_name = "abrupt-boundary producer",
  source_rel = file.path(
    "scripts", "full_reproduction", "run_abrupt_boundary_mainsettings.R"
  ),
  patch_fun = function(x) {
    x <- replace_exact(
      x,
      'out_dir <- "results/abrupt_boundary/main_simulation_settings_rerun"',
      paste0("out_dir <- ", r_string(abrupt_out)),
      label = "abrupt output directory"
    )
    x <- replace_exact(
      x,
      "data <- make_abrupt_main_data(n = 1000, p = 5, sigma2 = 0.1, seed = 123)",
      "data <- make_abrupt_main_data(n = 250, p = 5, sigma2 = 0.1, seed = 123)",
      label = "abrupt sample size"
    )
    x <- replace_exact(x, "L_range = c(25, 36),",
                       "L_range = c(25),", label = "abrupt CV L")
    x <- replace_exact(x, "a_lambda_range = c(15, 30, 35),",
                       "a_lambda_range = c(15),", label = "abrupt CV a")
    x <- replace_exact(x, "b_lambda_range = c(0.1, 1, 2),",
                       "b_lambda_range = c(1),", label = "abrupt CV b")
    x <- replace_exact(x, "niter = 800,",
                       "niter = 80,", label = "abrupt CV iterations")
    x <- replace_exact(x, "nburn = 100,",
                       "nburn = 20,", label = "abrupt CV burn")
    x <- replace_exact(x, "k_folds = 5",
                       "k_folds = 2", label = "abrupt CV folds")
    x <- replace_exact(x, "niter = 5000,",
                       "niter = 120,", expected = 2L,
                       label = "abrupt BSGL/GSVC iterations")
    x <- replace_exact(x, "nburn = 500,",
                       "nburn = 20,", expected = 2L,
                       label = "abrupt BSGL/GSVC burn")
    x
  },
  expected_files = abrupt_fit
)

# ===========================================================================
# 5. interaction_pilot_main_effects.R
#    Smaller generated data; exercises BSGL, GSVC, GGP-GAM, and MGWR.
# ===========================================================================
interaction_out <- file.path(smoke_root, "interaction")
interaction_fit <- file.path(interaction_out, "main_interaction_pilot_fit.rds")

results[[length(results) + 1L]] <- run_modified_script(
  test_name = "simulated-interaction producer",
  source_rel = file.path(
    "scripts", "full_reproduction", "interaction_pilot_main_effects.R"
  ),
  patch_fun = function(x) {
    x <- replace_exact(
      x,
      'out_dir <- "results/interaction/main_effects_rerun"',
      paste0("out_dir <- ", r_string(interaction_out)),
      label = "interaction output directory"
    )
    x <- replace_exact(
      x,
      "dat <- make_interaction_data_main(n = 1000, p_main = 5, sigma2 = 0.1, seed = 20260531)",
      "dat <- make_interaction_data_main(n = 250, p_main = 5, sigma2 = 0.1, seed = 20260531)",
      label = "interaction sample size"
    )
    x <- replace_exact(x, "niter <- 1200",
                       "niter <- 120", label = "interaction iterations")
    x <- replace_exact(x, "nburn <- 300",
                       "nburn <- 20", label = "interaction burn")
    x
  },
  expected_files = interaction_fit
)

# ===========================================================================
# 6. real_nointeraction_scp_compare.R
#    This script already exposes smoke-friendly environment variables, so the
#    original file is executed directly.
# ===========================================================================
real_no_out <- file.path(smoke_root, "real_nointeraction")
real_no_fit <- file.path(real_no_out, "real_nointeraction_scp_compare_fit.rds")

real_no_env <- c(
  paste0("REAL_NOINT_OUT_DIR=", r_path(real_no_out)),
  "REAL_NOINT_SAMPLE_SIZE=250",
  "REAL_NOINT_SAMPLE_SEED=20260531",
  "REAL_NOINT_SPLIT_SEED=456",
  "REAL_NOINT_NITER=120",
  "REAL_NOINT_NBURN=20",
  "REAL_NOINT_L=25",
  "REAL_NOINT_A_LAMBDA=20",
  "REAL_NOINT_B_LAMBDA=1",
  "REAL_NOINT_GRID_N=12",
  "REAL_NOINT_SAVE_FIT=true"
)

results[[length(results) + 1L]] <- run_original_script(
  test_name = "MODIS no-interaction producer",
  source_rel = file.path(
    "scripts", "full_reproduction", "real_nointeraction_scp_compare.R"
  ),
  expected_files = real_no_fit,
  env = real_no_env
)

# ===========================================================================
# 7. real_interaction_evi_validation.R
#    Short augmented BSGL fit; baseline disabled to avoid duplicating work.
#    Interaction choices match the retained downstream map workflow.
# ===========================================================================
real_int_out <- file.path(smoke_root, "real_interaction")
real_int_fit <- file.path(real_int_out, "evi_interaction_validation_fit.rds")

real_int_env <- c(
  paste0("REAL_EVI_INT_OUT_DIR=", r_path(real_int_out)),
  "REAL_EVI_INT_SAMPLE_SIZE=250",
  "REAL_EVI_INT_SAMPLE_SEED=20260531",
  "REAL_EVI_INT_SPLIT_SEED=456",
  "REAL_EVI_INT_NITER=120",
  "REAL_EVI_INT_NBURN=20",
  "REAL_EVI_INT_L=25",
  "REAL_EVI_INT_A_LAMBDA=20",
  "REAL_EVI_INT_B_LAMBDA=1",
  "REAL_EVI_INT_GRID_N=12",
  "REAL_EVI_INT_RUN_BASELINE=false",
  "REAL_EVI_INT_SAVE_FIT=true",
  "REAL_EVI_INT_LC_MODE=dummy",
  "REAL_EVI_INT_USEFUL_ANCHOR_VAR=NIR_reflectance",
  "REAL_EVI_INT_USEFUL_VAR=GPP",
  "REAL_EVI_INT_NULL_ANCHOR_VAR=NIR_reflectance",
  "REAL_EVI_INT_NULL_VAR=view_zenith_angle"
)

results[[length(results) + 1L]] <- run_original_script(
  test_name = "MODIS interaction producer",
  source_rel = file.path(
    "scripts", "full_reproduction", "real_interaction_evi_validation.R"
  ),
  expected_files = real_int_fit,
  env = real_int_env
)

# ===========================================================================
# Summary
# ===========================================================================
summary_df <- do.call(rbind, results)

cat("\n============================================================\n")
cat("CORE RUNTIME SMOKE-TEST SUMMARY\n")
cat("============================================================\n")
print(summary_df, row.names = FALSE)

# Verify smoke tests did not modify tracked/untracked repository state.
git_after <- git_status()
repo_unchanged <- identical(git_before, git_after)

cat("\nRepository status unchanged by smoke tests:", repo_unchanged, "\n")
if (!repo_unchanged) {
  cat("\nGit status before:\n")
  cat(git_before, sep = "\n")
  cat("\n\nGit status after:\n")
  cat(git_after, sep = "\n")
  cat("\n")
}

cat("\nSmoke-test logs are under:\n", file.path(smoke_root, "logs"), "\n", sep = "")
cat("Smoke-test outputs are under:\n", smoke_root, "\n", sep = "")

n_fail <- sum(summary_df$status == "FAIL")
n_pass <- sum(summary_df$status == "PASS")
n_skip <- sum(summary_df$status == "SKIP")

cat("\nPASS:", n_pass, " FAIL:", n_fail, " SKIP:", n_skip, "\n")

if (n_fail > 0L) {
  stop(
    "CORE RUNTIME SMOKE TESTS FAILED. ",
    "Inspect the printed failing test and its log."
  )
}

if (!repo_unchanged) {
  stop(
    "Core runtime tests completed, but repository git status changed. ",
    "Inspect git status before proceeding."
  )
}

cat("\n============================================================\n")
cat("CORE RUNTIME SMOKE TESTS PASSED\n")
cat("============================================================\n")
