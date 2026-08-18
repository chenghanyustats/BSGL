# ============================================================================
# BSGL remaining full-reproduction runtime smoke tests
#
# PLACEMENT:
#   scripts/06_remaining_runtime_smoke_tests.R
#
# RUN FROM REPOSITORY ROOT, after scripts/05_runtime_smoke_tests.R:
#   source("scripts/06_remaining_runtime_smoke_tests.R")
#
# This script:
#   * reuses Phase-1 smoke fit objects for downstream/post-processing tests;
#   * runs two shortened regeneration workflows that refit by design;
#   * writes all new test outputs to tempdir();
#   * does not overwrite retained manuscript outputs;
#   * reports ggspatial as BLOCKED if it is not installed.
# ============================================================================

cat("============================================================\n")
cat("BSGL remaining runtime smoke tests\n")
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

rscript <- file.path(R.home("bin"), "Rscript")
if (!file.exists(rscript)) stop("Could not locate Rscript.")

# ---------------------------------------------------------------------------
# Locate the most recent successful Phase-1 smoke directory.
# ---------------------------------------------------------------------------
required_phase1 <- c(
  file.path("abrupt_boundary", "abrupt_boundary_mainsettings_fit.rds"),
  file.path("interaction", "main_interaction_pilot_fit.rds"),
  file.path("real_nointeraction", "real_nointeraction_scp_compare_fit.rds"),
  file.path("real_interaction", "evi_interaction_validation_fit.rds")
)

is_valid_phase1 <- function(path) {
  nzchar(path) &&
    dir.exists(path) &&
    all(file.exists(file.path(path, required_phase1)))
}

phase1_candidates <- character()

if (exists("smoke_root", inherits = TRUE)) {
  candidate <- get("smoke_root", inherits = TRUE)
  if (is.character(candidate) && length(candidate) == 1L) {
    phase1_candidates <- c(phase1_candidates, candidate)
  }
}

phase1_candidates <- c(
  phase1_candidates,
  list.dirs(tempdir(), recursive = FALSE, full.names = TRUE)
)

phase1_candidates <- unique(phase1_candidates[
  grepl("BSGL_runtime_smoke_", basename(phase1_candidates), fixed = TRUE)
])

phase1_candidates <- phase1_candidates[
  vapply(phase1_candidates, is_valid_phase1, logical(1))
]

if (!length(phase1_candidates)) {
  stop(
    "No valid Phase-1 smoke directory was found in this R session. ",
    "Run source(\"scripts/05_runtime_smoke_tests.R\") first, then rerun this script."
  )
}

mtimes <- file.info(phase1_candidates)$mtime
phase1_root <- phase1_candidates[which.max(mtimes)]

phase2_root <- file.path(
  tempdir(),
  paste0("BSGL_remaining_smoke_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(phase2_root, recursive = TRUE, showWarnings = FALSE)

cat("Repository root:", project_root, "\n")
cat("Phase-1 smoke root:", phase1_root, "\n")
cat("Phase-2 output root:", phase2_root, "\n\n")

# ---------------------------------------------------------------------------
# Git-status guard.
# ---------------------------------------------------------------------------
git_status <- function() {
  out <- tryCatch(
    system2("git", c("status", "--short"), stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  )
  sort(out)
}
git_before <- git_status()

# ---------------------------------------------------------------------------
# Text utilities.
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
  if (n != expected) {
    stop(
      "Patch mismatch for ", label, ": expected ", expected,
      " occurrence(s), found ", n, "."
    )
  }
  gsub(old, new, text, fixed = TRUE)
}

replace_one_of <- function(text, alternatives, new, label) {
  counts <- vapply(alternatives, function(z) count_fixed(text, z), integer(1))
  hit <- which(counts > 0L)
  if (length(hit) != 1L || counts[hit] != 1L) {
    stop(
      "Patch mismatch for ", label, ". Candidate counts: ",
      paste(paste0("[", counts, "] ", alternatives), collapse = " | ")
    )
  }
  gsub(alternatives[hit], new, text, fixed = TRUE)
}

r_path <- function(path) {
  gsub("\\\\", "/", normalizePath(path, mustWork = FALSE))
}

r_string <- function(path) {
  paste0('"', gsub('"', '\\"', r_path(path), fixed = TRUE), '"')
}

# ---------------------------------------------------------------------------
# Runner.
# ---------------------------------------------------------------------------
run_modified_script <- function(
    test_name,
    source_rel,
    patch_fun = identity,
    expected_files = character(),
    env = character()) {

  cat("\n------------------------------------------------------------\n")
  cat("TEST:", test_name, "\n")
  cat("SOURCE:", source_rel, "\n")
  cat("------------------------------------------------------------\n")

  source_path <- file.path(project_root, source_rel)
  if (!file.exists(source_path)) stop("Missing source script: ", source_rel)

  text <- patch_fun(read_script(source_path))

  patched_dir <- file.path(phase2_root, "patched_scripts")
  log_dir <- file.path(phase2_root, "logs")
  device_dir <- file.path(phase2_root, "graphics_devices")
  dir.create(patched_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(device_dir, recursive = TRUE, showWarnings = FALSE)

  stub <- gsub("[^A-Za-z0-9_]+", "_", test_name)
  temp_script <- file.path(patched_dir, paste0(stub, ".R"))
  log_file <- file.path(log_dir, paste0(stub, ".log"))
  device_file <- file.path(device_dir, paste0(stub, ".pdf"))

  # Some plotting scripts print grid/cowplot objects before calling ggsave().
  # Under non-interactive Rscript this would otherwise create Rplots.pdf in the
  # repository root. Open an explicit temporary device so any such incidental
  # plotting is safely redirected into the smoke-test directory.
  fallback_device_file <- file.path(
    device_dir,
    paste0(stub, "_fallback.pdf")
  )

  # Redirect BOTH the currently opened device and any later fallback default
  # device into the smoke-test directory. This prevents non-interactive
  # plotting code from creating Rplots.pdf in the repository root even if a
  # script closes/reopens graphics devices internally.
  device_preamble <- paste0(
    "options(device = function(...) grDevices::pdf(file = ",
    deparse(fallback_device_file),
    "))\n",
    "grDevices::pdf(",
    deparse(device_file),
    ")"
  )
  device_epilogue <- paste(
    "while (grDevices::dev.cur() > 1L) grDevices::dev.off()",
    sep = "\n"
  )
  text <- paste(device_preamble, text, device_epilogue, sep = "\n")

  writeLines(
    strsplit(text, "\n", fixed = TRUE)[[1]],
    temp_script,
    useBytes = TRUE
  )

  parse(file = temp_script)

  status <- system2(
    rscript,
    args = c("--vanilla", shQuote(temp_script)),
    stdout = log_file,
    stderr = log_file,
    env = env
  )
  if (is.null(status)) status <- 0L

  ok_files <- if (length(expected_files)) file.exists(expected_files) else logical()

  if (status == 0L && (length(ok_files) == 0L || all(ok_files))) {
    cat("PASS:", test_name, "\n")
    return(data.frame(
      test = test_name,
      status = "PASS",
      exit_status = status,
      note = "",
      log = log_file,
      stringsAsFactors = FALSE
    ))
  }

  cat("FAIL:", test_name, "\n")
  if (length(expected_files) && any(!ok_files)) {
    cat("Missing expected output(s):\n")
    cat(paste0("  ", expected_files[!ok_files], "\n"), sep = "")
  }

  if (file.exists(log_file)) {
    z <- readLines(log_file, warn = FALSE, encoding = "UTF-8")
    cat("\nLast log lines:\n")
    cat(tail(z, 60), sep = "\n")
    cat("\n")
  }

  data.frame(
    test = test_name,
    status = "FAIL",
    exit_status = status,
    note = "",
    log = log_file,
    stringsAsFactors = FALSE
  )
}

blocked_row <- function(test_name, note) {
  cat("\nBLOCKED:", test_name, "\n")
  cat("Reason:", note, "\n")
  data.frame(
    test = test_name,
    status = "BLOCKED",
    exit_status = NA_integer_,
    note = note,
    log = NA_character_,
    stringsAsFactors = FALSE
  )
}

results <- list()

abrupt_fit <- file.path(
  phase1_root, "abrupt_boundary", "abrupt_boundary_mainsettings_fit.rds"
)
interaction_dir <- file.path(phase1_root, "interaction")
interaction_fit <- file.path(interaction_dir, "main_interaction_pilot_fit.rds")
real_no_dir <- file.path(phase1_root, "real_nointeraction")
real_no_fit <- file.path(real_no_dir, "real_nointeraction_scp_compare_fit.rds")
real_int_dir <- file.path(phase1_root, "real_interaction")
real_int_fit <- file.path(real_int_dir, "evi_interaction_validation_fit.rds")

# ===========================================================================
# A. Pure downstream/post-processing workflows
# ===========================================================================

# 1. Abrupt-boundary plotting
abrupt_plot_out <- file.path(phase2_root, "abrupt_plots")
dir.create(abrupt_plot_out, recursive = TRUE, showWarnings = FALSE)

results[[length(results) + 1L]] <- run_modified_script(
  "abrupt-boundary plotting",
  file.path("scripts", "full_reproduction",
            "plot_abrupt_boundary_mainsettings_figures.R"),
  patch_fun = function(x) {
    x <- replace_exact(
      x,
      'fit_path <- "results/abrupt_boundary/main_simulation_settings_rerun/abrupt_boundary_mainsettings_fit.rds"',
      paste0("fit_path <- ", r_string(abrupt_fit)),
      label = "abrupt fit path"
    )
    x <- replace_one_of(
      x,
      c(
        'out_dir <- "figures/reproduction_outputs/abrupt_boundary"',
        'out_dir <- "figures/abrupt_boundary/rerun"'
      ),
      paste0("out_dir <- ", r_string(abrupt_plot_out)),
      "abrupt plot output directory"
    )
    x
  },
  expected_files = c(
    file.path(abrupt_plot_out, "beta2_BSGL_abrupt_mainsettings_n250_p5.png"),
    file.path(abrupt_plot_out, "beta5_BSGL_abrupt_mainsettings_n250_p5.png")
  )
)

# 2. Interaction SCP metrics
results[[length(results) + 1L]] <- run_modified_script(
  "interaction SCP metrics",
  file.path("scripts", "full_reproduction", "interaction_scp_metrics.R"),
  patch_fun = function(x) {
    replace_exact(
      x,
      '"results/interaction/main_effects_rerun"',
      r_string(interaction_dir),
      expected = 2L,
      label = "interaction metric input/output directory"
    )
  },
  expected_files = c(
    file.path(interaction_dir, "interaction_bsgl_gs_scp_f1_fpr.csv"),
    file.path(interaction_dir, "interaction_bsgl_scp_maps.png")
  )
)

# 3. Interaction coefficient-recovery plots
results[[length(results) + 1L]] <- run_modified_script(
  "interaction beta recovery plots",
  file.path("scripts", "full_reproduction", "plot_interaction_beta_recovery.R"),
  patch_fun = function(x) {
    replace_exact(
      x,
      '"results/interaction/main_effects_rerun"',
      r_string(interaction_dir),
      expected = 2L,
      label = "interaction beta-recovery directory"
    )
  },
  expected_files = c(
    file.path(interaction_dir, "beta6_x1x2_recovery.png"),
    file.path(interaction_dir, "beta7_x1x5_recovery.png")
  )
)

# 4. X1 contribution maps
results[[length(results) + 1L]] <- run_modified_script(
  "interaction X1 contribution plots",
  file.path("scripts", "full_reproduction", "plot_interaction_x1_contribution.R"),
  patch_fun = function(x) {
    replace_exact(
      x,
      '"results/interaction/main_effects_rerun"',
      r_string(interaction_dir),
      expected = 2L,
      label = "interaction X1 directory"
    )
  },
  expected_files = c(
    file.path(interaction_dir, "x1_contribution_maps",
              "beta12_bsgl_scp_map.png"),
    file.path(interaction_dir, "x1_contribution_maps",
              "x1_contribution_high_minus_low_x2.png")
  )
)

# 5. Organize real-data 10-variable summaries
organized_real <- file.path(phase2_root, "organized_real")
results[[length(results) + 1L]] <- run_modified_script(
  "organize MODIS 10-variable results",
  file.path("scripts", "full_reproduction", "organize_real_10variable_results.R"),
  patch_fun = function(x) {
    x <- replace_exact(
      x,
      'in_dir <- "results/real_data/nointeraction_scp_compare_n10000_rerun"',
      paste0("in_dir <- ", r_string(real_no_dir)),
      label = "organized real input"
    )
    x <- replace_exact(
      x,
      'out_dir <- "results/real_data/clean_real_data_10variable_rerun"',
      paste0("out_dir <- ", r_string(organized_real)),
      label = "organized real output"
    )
    x
  },
  expected_files = c(
    file.path(organized_real, "real_scp_10variable_wide.csv"),
    file.path(organized_real, "real_scp_10variable_bsgl_gaussian_svc.png"),
    file.path(organized_real, "real_scp_10variable_table.tex")
  )
)

# 6. MODIS 10-variable geographic SCP map
real_map_out <- file.path(phase2_root, "real_10variable_maps")
results[[length(results) + 1L]] <- run_modified_script(
  "MODIS 10-variable SCP maps",
  file.path("scripts", "full_reproduction",
            "plot_real_10variable_scp_maps_mean_lc.R"),
  patch_fun = function(x) {
    x <- replace_exact(
      x,
      'out_dir <- "results/real_data/clean_real_data_10variable_rerun"',
      paste0("out_dir <- ", r_string(real_map_out)),
      label = "real map output"
    )
    x <- replace_exact(
      x,
      'fit_obj <- readRDS("results/real_data/nointeraction_scp_compare_n10000_rerun/real_nointeraction_scp_compare_fit.rds")',
      paste0("fit_obj <- readRDS(", r_string(real_no_fit), ")"),
      label = "real map fit"
    )
    x <- replace_exact(
      x,
      "sample_size <- 10000",
      "sample_size <- 250",
      label = "real map smoke sample size"
    )
    x
  },
  expected_files = c(
    file.path(real_map_out, "real_scp_10variable_maps_bsgl_mean_lc.png"),
    file.path(real_map_out, "real_scp_10variable_maps_bsgl_mean_lc.csv")
  )
)

# 7. MODIS interaction maps
results[[length(results) + 1L]] <- run_modified_script(
  "MODIS interaction SCP maps",
  file.path("scripts", "full_reproduction",
            "plot_real_interaction_scp_maps_1x2.R"),
  patch_fun = function(x) {
    x <- replace_exact(
      x,
      'in_dir <- "results/interaction/evi_interaction_validation_n10000_lcdummy_nir_viewz_nir_gpp"',
      paste0("in_dir <- ", r_string(real_int_dir)),
      label = "real interaction map directory"
    )
    x <- replace_exact(
      x,
      "sample_size <- 10000",
      "sample_size <- 250",
      label = "real interaction map sample size"
    )
    x
  },
  expected_files = c(
    file.path(real_int_dir, "interaction_scp_maps_1x2.png"),
    file.path(real_int_dir, "interaction_scp_maps_1x2_geo.png")
  )
)

# 8. MODIS held-out predictive quantile maps
pred_out <- file.path(phase2_root, "predictive_quantiles")
results[[length(results) + 1L]] <- run_modified_script(
  "MODIS predictive quantile plots",
  file.path("scripts", "full_reproduction",
            "plot_test_evi_predictive_quantiles.R"),
  patch_fun = function(x) {
    x <- replace_one_of(
      x,
      c(
        'out_dir <- "figures/reproduction_outputs/real_data"',
        'out_dir <- "pic"'
      ),
      paste0("out_dir <- ", r_string(pred_out)),
      "predictive plot output"
    )
    x <- replace_exact(
      x,
      'idx <- if (nrow(real_data) > 10000) sample(nrow(real_data), 10000) else seq_len(nrow(real_data))',
      'idx <- if (nrow(real_data) > 250) sample(nrow(real_data), 250) else seq_len(nrow(real_data))',
      label = "predictive smoke sample size"
    )
    x <- replace_exact(
      x,
      'fit <- readRDS("results/real_data/nointeraction_scp_compare_n10000_rerun/real_nointeraction_scp_compare_fit.rds")',
      paste0("fit <- readRDS(", r_string(real_no_fit), ")"),
      label = "predictive fit path"
    )
    x
  },
  expected_files = c(
    file.path(pred_out, "test_evi_predictive_q05_mean_q95_maps.png"),
    file.path(pred_out, "test_evi_predictive_q05_mean_q95_density.png"),
    file.path(pred_out, "test_evi_predictive_q05_mean_q95_summary.csv")
  )
)

# 9. Reviewer-facing fixed four-method prediction table
real_pred_out <- file.path(phase2_root, "real_prediction_four_methods.csv")
results[[length(results) + 1L]] <- run_modified_script(
  "four-method MODIS prediction summary writer",
  file.path("scripts", "full_reproduction",
            "real_prediction_four_method_compare.R"),
  patch_fun = function(x) {
    replace_exact(
      x,
      'out_file <- file.path("results", "real_data", "real_prediction_four_methods.csv")',
      paste0("out_file <- ", r_string(real_pred_out)),
      label = "four-method summary output"
    )
  },
  expected_files = real_pred_out
)

# 10. Beta-5 retained-fit plotting path
gray_fit_candidates <- c(
  file.path(project_root, "figures", "reproduction_outputs", "simulation",
            "gray_four_method_n1000_p5_rerun", "gray_run_fit.rds"),
  file.path(project_root, "figures", "simulation",
            "gray_four_method_n1000_p5_rerun", "gray_run_fit.rds")
)
gray_fit_candidates <- gray_fit_candidates[file.exists(gray_fit_candidates)]

if (!length(gray_fit_candidates)) {
  results[[length(results) + 1L]] <- blocked_row(
    "beta5 retained-fit plotting",
    "No retained gray_run_fit.rds was found."
  )
} else {
  gray_fit <- gray_fit_candidates[1]
  beta5_out <- file.path(phase2_root, "beta5_retained")
  results[[length(results) + 1L]] <- run_modified_script(
    "beta5 retained-fit plotting",
    file.path("scripts", "full_reproduction",
              "plot_beta5_ownscale_onebar.R"),
    patch_fun = function(x) {
      x <- replace_one_of(
        x,
        c(
          'fit_file <- "figures/reproduction_outputs/simulation/gray_four_method_n1000_p5_rerun/gray_run_fit.rds"',
          'fit_file <- "figures/simulation/gray_four_method_n1000_p5_rerun/gray_run_fit.rds"'
        ),
        paste0("fit_file <- ", r_string(gray_fit)),
        "beta5 retained fit"
      )
      x <- replace_one_of(
        x,
        c(
          'out_dir <- "figures/reproduction_outputs/simulation/gray_four_method_n1000_p5_rerun"',
          'out_dir <- "figures/simulation/gray_four_method_n1000_p5_rerun"'
        ),
        paste0("out_dir <- ", r_string(beta5_out)),
        "beta5 output directory"
      )
      x
    },
    expected_files = c(
      file.path(beta5_out, "beta5_ownscale_onebar_n1000_p5.png"),
      file.path(beta5_out, "beta5_ownscale_ranges.csv")
    )
  )
}

# ===========================================================================
# B. Two remaining regeneration scripts that refit by design
# ===========================================================================

# 11. Four-method gray workflow, shortened
four_out <- file.path(phase2_root, "four_method_gray")
results[[length(results) + 1L]] <- run_modified_script(
  "four-method gray regeneration",
  file.path("scripts", "full_reproduction",
            "plot_four_method_gray_n1000_p5.R"),
  patch_fun = function(x) {
    x <- replace_exact(x, "L_range = c(25, 36),",
                       "L_range = c(25),", label = "gray CV L")
    x <- replace_exact(x, "a_lambda_range = c(15, 30, 35),",
                       "a_lambda_range = c(15),", label = "gray CV a")
    x <- replace_exact(x, "b_lambda_range = c(0.1, 1, 2),",
                       "b_lambda_range = c(1),", label = "gray CV b")
    x <- replace_exact(x, "niter = 800,",
                       "niter = 80,", label = "gray CV iterations")
    x <- replace_exact(x, "nburn = 100,",
                       "nburn = 20,", label = "gray CV burn")
    x <- replace_exact(x, "k_folds = 5",
                       "k_folds = 2", label = "gray CV folds")
    x <- replace_exact(x, "niter = 5000,",
                       "niter = 120,", expected = 2L,
                       label = "gray BSGL/GSVC iterations")
    x <- replace_exact(x, "nburn = 500,",
                       "nburn = 20,", expected = 2L,
                       label = "gray BSGL/GSVC burn")
    x <- replace_one_of(
      x,
      c(
        'out_dir <- "figures/reproduction_outputs/simulation/gray_four_method_n1000_p5_rerun"',
        'out_dir <- "figures/simulation/gray_four_method_n1000_p5_rerun"'
      ),
      paste0("out_dir <- ", r_string(four_out)),
      "gray regeneration output"
    )
    x
  },
  expected_files = c(
    file.path(four_out, "gray_run_summary.csv"),
    file.path(four_out, "gray_run_fit.rds")
  )
)

# 12. Main Figure 2 regeneration, shortened
fig2_out <- file.path(phase2_root, "fig2")
fig2_res <- file.path(phase2_root, "fig2_results")

results[[length(results) + 1L]] <- run_modified_script(
  "main Figure 2 regeneration",
  file.path("scripts", "full_reproduction",
            "regenerate_main_fig2_sim_scp_map.R"),
  patch_fun = function(x) {
    x <- replace_exact(
      x,
      'out_dir <- file.path("figures", "main")',
      paste0("out_dir <- ", r_string(fig2_out)),
      label = "fig2 output directory"
    )
    x <- replace_exact(
      x,
      'res_dir <- file.path("results", "simulation", "main_fig2")',
      paste0("res_dir <- ", r_string(fig2_res)),
      label = "fig2 result directory"
    )
    x
  },
  expected_files = c(
    file.path(fig2_out, "fig2_simulation_scp_map.png"),
    file.path(fig2_res, "main_fig2_n5000_m10_rep02_scp.csv")
  ),
  env = c(
    "FIG2_L=25",
    "FIG2_A_LAMBDA=30",
    "FIG2_B_LAMBDA=1",
    "FIG2_NITER=120",
    "FIG2_NBURN=20",
    "FIG2_SEED=123"
  )
)

# ===========================================================================
# C. plot_maps.R source dependency check
# ===========================================================================
# This is a function collection. It directly attaches ggspatial at source time.
if (!requireNamespace("ggspatial", quietly = TRUE)) {
  results[[length(results) + 1L]] <- blocked_row(
    "plot_maps.R source",
    "Optional reproduction dependency 'ggspatial' is not installed."
  )
} else {
  results[[length(results) + 1L]] <- run_modified_script(
    "plot_maps.R source",
    file.path("scripts", "full_reproduction", "plot_maps.R"),
    patch_fun = identity
  )
}

# ===========================================================================
# Summary and repository guard
# ===========================================================================
summary_df <- do.call(rbind, results)

cat("\n============================================================\n")
cat("REMAINING RUNTIME SMOKE-TEST SUMMARY\n")
cat("============================================================\n")
print(summary_df, row.names = FALSE)

git_after <- git_status()
repo_unchanged <- identical(git_before, git_after)

cat("\nRepository status unchanged by Phase-2 tests:", repo_unchanged, "\n")

n_pass <- sum(summary_df$status == "PASS")
n_fail <- sum(summary_df$status == "FAIL")
n_blocked <- sum(summary_df$status == "BLOCKED")

cat("\nPASS:", n_pass,
    " FAIL:", n_fail,
    " BLOCKED:", n_blocked, "\n")

cat("\nLogs:\n", file.path(phase2_root, "logs"), "\n", sep = "")
cat("Outputs:\n", phase2_root, "\n", sep = "")

if (n_fail > 0L) {
  stop(
    "REMAINING RUNTIME SMOKE TESTS FAILED. ",
    "Inspect the failing test log(s)."
  )
}

if (!repo_unchanged) {
  stop(
    "Smoke tests completed but git status changed. ",
    "Inspect the repository before proceeding."
  )
}

if (n_blocked > 0L) {
  cat("\n============================================================\n")
  cat("SMOKE TESTS PASSED EXCEPT FOR BLOCKED OPTIONAL DEPENDENCIES\n")
  cat("============================================================\n")
  cat("Install the reported dependency/dependencies and rerun this script.\n")
} else {
  cat("\n============================================================\n")
  cat("ALL REMAINING RUNTIME SMOKE TESTS PASSED\n")
  cat("============================================================\n")
}
