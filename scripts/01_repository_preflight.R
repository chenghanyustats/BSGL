# BSGL repository-wide reproducibility preflight.
# Run from the repository root:
#   source("scripts/01_repository_preflight.R")

source("scripts/00_project_setup.R")

cat("============================================================\n")
cat("BSGL repository preflight\n")
cat("============================================================\n\n")

required_paths <- c(
  "DESCRIPTION",
  "NAMESPACE",
  "R",
  "scripts",
  "vignettes",
  "inst/extdata/simulation_demo/n1000_p5/train_data.RData",
  "inst/extdata/simulation_demo/n1000_p5/test_data.RData",
  "inst/extdata/simulation_demo/n1000_p5/grid_data.RData",
  "inst/extdata/simulation_demo/n1000_p5/meta_info.RData",
  "data/data_cleaned_small_expanded.RData",
  "data/data_rep10/n1000_p5.RData",
  "data/data_rep10/n5000_p10.RData",
  "data/data_indep/n1000_p5/train_data.RData",
  "data/data_indep/n5000_p10/train_data.RData",
  "data/data_indep/n10000_p20/train_data.RData",
  "data/geo_naturalearth_10m/ne_10m_coastline.zip",
  "data/geo_naturalearth_10m/ne_10m_lakes.zip",
  "data/geo_naturalearth_10m/ne_10m_ocean.zip",
  "data/geo_naturalearth_10m/ne_10m_rivers_lake_centerlines.zip"
)

path_status <- data.frame(
  path = required_paths,
  exists = file.exists(project_path(required_paths)) |
           dir.exists(project_path(required_paths)),
  stringsAsFactors = FALSE
)
print(path_status, row.names = FALSE)
if (!all(path_status$exists)) {
  stop("Required repository inputs are missing.")
}
cat("\nPASS: required repository inputs are present.\n")

# Check every full-reproduction R script for parse errors.
script_files <- list.files(
  project_path("scripts", "full_reproduction"),
  pattern = "\\.R$",
  full.names = TRUE
)
parse_ok <- vapply(script_files, function(f) {
  tryCatch({ parse(file = f); TRUE }, error = function(e) {
    message("Parse failure in ", f, ": ", conditionMessage(e)); FALSE
  })
}, logical(1))
if (!all(parse_ok)) stop("At least one full-reproduction R script does not parse.")
cat("PASS: all full-reproduction R scripts parse.\n")

# Look for stale repository references that should no longer occur.
text_files <- c(
  project_path("README.md"),
  list.files(project_path("vignettes"), pattern = "\\.Rmd$", full.names = TRUE),
  script_files
)
stale_patterns <- c(
  "Qishi7/BSGL",
  "inst/extinst/extdata",
  "(^|[^[:alnum:]_])data/simulation_demo",
  "calc_beta\\(lasso_gam"
)
stale_hits <- list()
for (f in text_files) {
  z <- readLines(f, warn = FALSE, encoding = "UTF-8")
  for (pat in stale_patterns) {
    idx <- grep(pat, z, perl = TRUE)
    if (length(idx)) {
      stale_hits[[length(stale_hits) + 1L]] <- data.frame(
        file = sub(paste0("^", project_root, "/?"), "", f),
        line = idx, pattern = pat, text = z[idx],
        stringsAsFactors = FALSE
      )
    }
  }
}
if (length(stale_hits)) {
  stale_df <- do.call(rbind, stale_hits)
  print(stale_df, row.names = FALSE)
  stop("Stale path/function references remain.")
}
cat("PASS: no known stale repository references remain.\n")

# Verify the GSVC data grid used by sim_gpsvc.R exists.
gs_grid <- expand.grid(n = c(1000, 5000, 10000), p = c(5, 10, 20))
gs_dirs <- file.path(project_root, "data", "data_indep",
                     paste0("n", gs_grid$n, "_p", gs_grid$p))
if (!all(dir.exists(gs_dirs))) {
  stop("The retained GSVC paper grid is incomplete under data/data_indep/.")
}
cat("PASS: retained GSVC simulation grid is complete.\n")

# Quick non-fitting repository checks.
cat("\nRunning retained-result summary...\n")
source(project_path("scripts", "02_summarize_main_results.R"), chdir = FALSE)
cat("\nRunning figure inventory...\n")
source(project_path("scripts", "04_figure_inventory.R"), chdir = FALSE)

# Render vignettes only when requested. Their long fitting chunks are eval=FALSE.
render_vignettes <- tolower(
  Sys.getenv("BSGL_PREFLIGHT_RENDER_VIGNETTES", unset = "false")
) %in% c("1", "true", "yes")

if (render_vignettes) {
  ensure_packages(c("rmarkdown", "knitr"))
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(project_root)
  for (vf in c("vignettes/simulation_demo.Rmd",
               "vignettes/modis_real_data_demo.Rmd")) {
    cat("Rendering ", vf, "...\n", sep = "")
    rmarkdown::render(vf, quiet = FALSE, clean = TRUE)
  }
  cat("PASS: vignettes rendered.\n")
} else {
  cat("\nVignette rendering skipped. To enable it, run:\n")
  cat('Sys.setenv(BSGL_PREFLIGHT_RENDER_VIGNETTES = "true")\n')
  cat('source("scripts/01_repository_preflight.R")\n')
}

cat("\n============================================================\n")
cat("REPOSITORY PREFLIGHT PASSED\n")
cat("============================================================\n")
