get_project_root <- function() {
  find_root <- function(path) {
    path <- normalizePath(path, mustWork = TRUE)
    repeat {
      if (dir.exists(file.path(path, "R")) && dir.exists(file.path(path, "scripts"))) {
        return(path)
      }
      parent <- dirname(path)
      if (identical(parent, path)) break
      path <- parent
    }
    NULL
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
    root <- find_root(dirname(script_path))
    if (!is.null(root)) return(root)
  }

  root <- find_root(getwd())
  if (!is.null(root)) return(root)

  stop("Run scripts with Rscript from the repository root or by script file path.")
}

project_root <- get_project_root()
source_r <- function(file) source(file.path(project_root, "R", file), chdir = TRUE)
project_path <- function(...) file.path(project_root, ...)

ensure_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

ensure_packages("ggplot2")
suppressPackageStartupMessages({
  library(ggplot2)
})

load_core_bsgl <- function() {
  source_r("helpers_gplasso.R")
  source_r("helpers_gpprior.R")
  source_r("help_inclu.R")
}

load_simulation_methods <- function() {
  source_r("beta_functions.R")
  source_r("load_data.R")
  source_r("helpers_ggpgam.R")
  source_r("helpers_gplasso.R")
  source_r("helpers_gpprior.R")
  source_r("helpers_mgwr.R")
  source_r("calc_beta.R")
  source_r("beta_mse.R")
  source_r("help_inclu.R")
}
