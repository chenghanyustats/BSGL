# run_demo.R
# Standalone reviewer-facing example for the BSGL paper.
#
# This script does NOT require the BSGL GitHub repository to be cloned locally.
# It:
#   1. installs the BSGL package directly from GitHub if needed,
#   2. installs the package's required R dependencies automatically,
#   3. loads the prepared simulation demo data,
#   4. fits BSGL, GSVC, GGP-GAM, and MGWR,
#   5. writes a common test-set comparison table.
#
# Preferred long-term package layout:
#   inst/extdata/simulation_demo/n1000_p5/
# containing train_data.RData and test_data.RData.
#
# Until those files are installed with the package, this script falls back to
# downloading ONLY the two required demo-data files from the GitHub repository.

cat("===============================================\n")
cat("BSGL Package Example\n")
cat("===============================================\n\n")

# -------------------------------------------------------------------------
# 1. Install/load BSGL and required dependencies
# -------------------------------------------------------------------------
cat("1. PACKAGE: Installing/loading BSGL\n")
cat("   -------------------------------\n")

repos <- getOption("repos")
if (is.null(repos) ||
    length(repos) == 0L ||
    identical(unname(repos[["CRAN"]]), "@CRAN@") ||
    is.na(repos[["CRAN"]]) ||
    !nzchar(repos[["CRAN"]])) {
    repos["CRAN"] <- "https://cloud.r-project.org"
    options(repos = repos)
}
cat("   CRAN repository:", getOption("repos")[["CRAN"]], "\n")

# Install 'remotes' only if it is not already available.
if (!requireNamespace("remotes", quietly = TRUE)) {
    cat("   Installing installer package: remotes\n")
    install.packages("remotes")
}

# By default, install BSGL only when it is not already available.
# Set BSGL_DEMO_INSTALL=true to force a fresh GitHub installation.
force_install <- tolower(
    Sys.getenv("BSGL_DEMO_INSTALL", unset = "false")
) %in% c("1", "true", "yes")

# Treat a missing or broken local installation as unavailable rather than
# allowing namespace-load errors to abort the bootstrap before reinstalling.
bsgl_available <- function() {
    tryCatch(
        requireNamespace("BSGL", quietly = TRUE),
        error = function(e) {
            cat(
                "   Existing BSGL installation is not usable: ",
                conditionMessage(e), "\n",
                sep = ""
            )
            FALSE
        }
    )
}

if (force_install || !bsgl_available()) {
    cat("   Installing BSGL directly from GitHub...\n")
    cat("   Repository: chenghanyustats/BSGL\n")

    # We are already in the branch where BSGL is unavailable or the user
    # explicitly requested reinstallation, so force a real GitHub install.
    remotes::install_github(
        "chenghanyustats/BSGL",
        dependencies = NA,
        upgrade = "never",
        build_vignettes = FALSE,
        force = TRUE
    )
}

if (!bsgl_available()) {
    stop(
        "BSGL could not be installed or loaded from GitHub. ",
        "Check the installation messages above and .libPaths()."
    )
}

suppressPackageStartupMessages(library(BSGL))
bsgl_version <- tryCatch(
    as.character(utils::packageVersion("BSGL")),
    error = function(e) "unknown"
)
cat("   Loaded package: BSGL", bsgl_version, "\n")

metric_row <- function(method, observed, predicted) {
    if (!is.numeric(predicted)) {
        stop("Predictions are not numeric.")
    }
    if (length(predicted) != length(observed)) {
        stop(
            "Prediction length (", length(predicted),
            ") does not match observed length (", length(observed), ")."
        )
    }
    if (any(!is.finite(predicted))) {
        stop("Predictions contain NA, NaN, or infinite values.")
    }
    
    data.frame(
        Method = method,
        Status = "Success",
        Test_MSE = mean((observed - predicted)^2),
        Test_MAE = mean(abs(observed - predicted)),
        Message = "",
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

status_row <- function(method, status, message = "") {
    data.frame(
        Method = method,
        Status = status,
        Test_MSE = NA_real_,
        Test_MAE = NA_real_,
        Message = as.character(message),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

run_method <- function(method, observed, fit_and_predict) {
    cat("   Fitting ", method, "...
", sep = "")
    
    tryCatch(
        {
            predicted <- fit_and_predict()
            ans <- metric_row(method, observed, predicted)
            cat("   ", method, " completed successfully.
", sep = "")
            ans
        },
        error = function(e) {
            msg <- conditionMessage(e)
            cat("   ", method, " failed: ", msg, "
", sep = "")
            status_row(method, "Failed", msg)
        }
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

# -------------------------------------------------------------------------
# 2. Load prepared simulation data
# -------------------------------------------------------------------------
cat("\n2. DATA: Loading prepared simulation subset\n")
cat("   ---------------------------------------\n")

# Preferred route: data installed inside the package under inst/extdata/.
pkg_demo_dir <- system.file(
    "extdata", "simulation_demo", "n1000_p5",
    package = "BSGL"
)

if (nzchar(pkg_demo_dir)) {
    train_file <- file.path(pkg_demo_dir, "train_data.RData")
    test_file  <- file.path(pkg_demo_dir, "test_data.RData")
    cat("   Using demo data installed with BSGL.\n")
} else {
    # Transitional fallback:
    # download only the two files needed by this demo, not the whole repository.
    cat("   Installed package does not yet contain extdata demo files.\n")
    cat("   Downloading the two required demo-data files only...\n")
    
    demo_dir <- file.path(tempdir(), "BSGL_simulation_demo_n1000_p5")
    dir.create(demo_dir, recursive = TRUE, showWarnings = FALSE)
    
    train_file <- file.path(demo_dir, "train_data.RData")
    test_file  <- file.path(demo_dir, "test_data.RData")
    
    base_url <- paste0(
        "https://raw.githubusercontent.com/chenghanyustats/BSGL/",
        "main/inst/extdata/simulation_demo/n1000_p5/"
    )
    
    if (!file.exists(train_file)) {
        utils::download.file(
            paste0(base_url, "train_data.RData"),
            train_file,
            mode = "wb",
            quiet = FALSE
        )
    }
    
    if (!file.exists(test_file)) {
        utils::download.file(
            paste0(base_url, "test_data.RData"),
            test_file,
            mode = "wb",
            quiet = FALSE
        )
    }
}

if (!file.exists(train_file) || !file.exists(test_file)) {
    stop("Prepared simulation train/test data could not be located.")
}

load(train_file)
load(test_file)

if (!exists("train_data") || !exists("test_data")) {
    stop(
        "Demo files were found, but expected objects 'train_data' and ",
        "'test_data' were not created."
    )
}

cat("   Training data: n =", length(train_data$y),
    ", m =", ncol(train_data$X), "\n")
cat("   Test data: n =", length(test_data$y), "\n")

# -------------------------------------------------------------------------
# 3. Fit methods
# -------------------------------------------------------------------------
cat("\n3. MODELS: Fitting package methods\n")
cat("   -------------------------------\n")

L <- as.integer(Sys.getenv("DEMO_L", "25"))
niter <- as.integer(Sys.getenv("DEMO_NITER", "1000"))
nburn <- as.integer(Sys.getenv("DEMO_NBURN", "100"))
a_lam <- as.numeric(Sys.getenv("DEMO_A_LAM", "15"))
b_lam <- as.numeric(Sys.getenv("DEMO_B_LAM", "1"))
run_mgwr <- tolower(Sys.getenv("DEMO_RUN_MGWR", unset = "true")) %in%
    c("1", "true", "yes")

cat("   L =", L,
    "| niter =", niter,
    "| nburn =", nburn,
    "| a_lam =", a_lam,
    "| b_lam =", b_lam, "\n")
cat("   Run MGWR =", run_mgwr, "\n\n")

set.seed(123)

bsgl_row <- run_method(
    "BSGL",
    test_data$y,
    function() {
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
        pred_bsgl(
            bsgl_fit,
            test_data$X,
            test_data$coords
        )$mean
    }
)

gsvc_row <- run_method(
    "GSVC",
    test_data$y,
    function() {
        gsvc_fit <- fit_gs(
            y = train_data$y,
            X = train_data$X,
            coords = train_data$coords,
            L = L,
            niter = niter,
            nburn = nburn,
            verbose = FALSE
        )
        pred_gs(
            gsvc_fit,
            test_data$X,
            test_data$coords
        )$mean
    }
)

gam_row <- run_method(
    "GGP-GAM",
    test_data$y,
    function() {
        gam_fit <- fit_gam(
            train_data$X,
            train_data$y,
            train_data$coords
        )
        pred_gam(
            gam_fit,
            test_data$X,
            test_data$coords
        )
    }
)

if (run_mgwr) {
    mgwr_row <- run_method(
        "MGWR",
        test_data$y,
        function() {
            mgwr_fit <- fit_mgwr(
                train_data$X,
                train_data$y,
                train_data$coords
            )
            pred_mgwr(
                mgwr_fit,
                test_data$X,
                test_data$coords
            )
        }
    )
} else {
    mgwr_message <- paste(
        "MGWR was not run by configuration.",
        "Set DEMO_RUN_MGWR=true to run it."
    )
    cat("   MGWR not run: ", mgwr_message, "\n", sep = "")
    mgwr_row <- status_row("MGWR", "Not run", mgwr_message)
}

# -------------------------------------------------------------------------
# 4. Compare and save results
# -------------------------------------------------------------------------
cat("\n4. RESULTS: Test-set comparison\n")
cat("   ----------------------------\n")

comparison <- rbind(
    bsgl_row,
    gsvc_row,
    gam_row,
    mgwr_row
)

comparison$Test_MSE <- round(comparison$Test_MSE, 6)
comparison$Test_MAE <- round(comparison$Test_MAE, 6)

print(comparison, row.names = FALSE, na.print = "")

dir.create("results", recursive = TRUE, showWarnings = FALSE)
out_file <- file.path("results", "demo_method_comparison.csv")
write.csv(
    comparison,
    out_file,
    row.names = FALSE,
    na = ""
)

cat("\n   Saved comparison table to:", out_file, "\n")

failed_methods <- comparison$Method[comparison$Status == "Failed"]
not_run_methods <- comparison$Method[comparison$Status == "Not run"]

cat("\n===============================================\n")
if (length(failed_methods) > 0L) {
    cat("DEMO COMPLETED WITH METHOD FAILURES\n")
    cat("Failed:", paste(failed_methods, collapse = ", "), "\n")
} else {
    cat("DEMO COMPLETED SUCCESSFULLY\n")
}
if (length(not_run_methods) > 0L) {
    cat("Not run:", paste(not_run_methods, collapse = ", "), "\n")
}
cat("See the comparison table for method status and messages.\n")
cat("===============================================\n")
