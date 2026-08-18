# run_real_demo.R - Standalone MODIS EVI reviewer example.
#
# This script does NOT require the BSGL GitHub repository to be cloned locally.
# It:
#   1. installs the BSGL package directly from GitHub if needed,
#   2. installs the package's required R dependencies automatically,
#   3. loads the prepared MODIS EVI demo data,
#   4. fits the comparison methods using the existing real-data workflow,
#   5. writes one common comparison table.
#
# Preferred long-term package layout:
#   inst/extdata/real_demo/data_cleaned_small_expanded.RData
#
# Until that file is installed with the package, this script falls back to
# downloading ONLY the required demo-data file from the GitHub repository.

cat("===============================================\n")
cat("BSGL MODIS EVI Example\n")
cat("===============================================\n\n")

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

if (!requireNamespace("remotes", quietly = TRUE)) {
    cat("   Installing installer package: remotes\n")
    install.packages("remotes")
}

# By default, install BSGL only when it is not already installed.
# Set BSGL_DEMO_INSTALL=true to force a fresh GitHub installation.
force_install <- tolower(
    Sys.getenv("BSGL_DEMO_INSTALL", unset = "false")
) %in% c("1", "true", "yes")

if (force_install || !requireNamespace("BSGL", quietly = TRUE)) {
    cat("   Installing BSGL directly from GitHub...\n")
    cat("   Repository: chenghanyustats/BSGL\n")
    
    remotes::install_github(
        "chenghanyustats/BSGL",
        dependencies = NA,
        upgrade = "never",
        build_vignettes = FALSE,
        force = force_install
    )
}

if (!requireNamespace("BSGL", quietly = TRUE)) {
    stop("BSGL could not be installed or loaded.")
}

suppressPackageStartupMessages(library(BSGL))
cat("   Loaded package: BSGL",
    as.character(utils::packageVersion("BSGL")), "\n")

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

# Preferred route: data installed inside the package under inst/extdata/.
data_file <- system.file(
    "extdata",
    "real_demo",
    "data_cleaned_small_expanded.RData",
    package = "BSGL"
)

if (nzchar(data_file) && file.exists(data_file)) {
    cat("   Using MODIS demo data installed with BSGL.\n")
} else {
    # Transitional fallback:
    # download only the one file needed by this demo, not the whole repository.
    cat("   Installed package does not yet contain the real-demo extdata file.\n")
    cat("   Downloading the required MODIS demo-data file only...\n")
    
    demo_dir <- file.path(tempdir(), "BSGL_real_demo")
    dir.create(demo_dir, recursive = TRUE, showWarnings = FALSE)
    
    data_file <- file.path(
        demo_dir,
        "data_cleaned_small_expanded.RData"
    )
    
    data_url <- paste0(
        "https://raw.githubusercontent.com/chenghanyustats/BSGL/",
        "main/data/data_cleaned_small_expanded.RData"
    )
    
    if (!file.exists(data_file)) {
        utils::download.file(
            data_url,
            data_file,
            mode = "wb",
            quiet = FALSE
        )
    }
}

if (!file.exists(data_file)) {
    stop("Prepared MODIS demo data could not be located or downloaded.")
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
    "MGWR is not run by default for the real-data demo because",
    "GWmodel::gwr.multiscale can be substantially slower on real-data subsets.",
    "Set REAL_DEMO_RUN_MGWR=true to attempt the MGWR fit."
)

cat("   L =", L, "| niter =", niter, "| nburn =", nburn,
    "| a_lam =", a_lam, "| b_lam =", b_lam, "\n")
cat("   Run MGWR =", run_mgwr, "\n\n")

set.seed(20260531)

bsgl_row <- run_method(
    "BSGL",
    test_y,
    function() {
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
        pred_bsgl(bsgl_fit, test_X, test_coords)$mean
    }
)

gsvc_row <- run_method(
    "GSVC",
    test_y,
    function() {
        gsvc_fit <- fit_gs(
            y = train_y,
            X = train_X,
            coords = train_coords,
            L = L,
            niter = niter,
            nburn = nburn,
            verbose = FALSE
        )
        pred_gs(gsvc_fit, test_X, test_coords)$mean
    }
)

gam_row <- run_method(
    "GGP-GAM",
    test_y,
    function() {
        gam_fit <- fit_gam(train_X, train_y, train_coords)
        pred_gam(gam_fit, test_X, test_coords)
    }
)

if (!run_mgwr) {
    cat("   MGWR not run: ", mgwr_skip_message, "\n", sep = "")
    mgwr_row <- status_row("MGWR", "Not run", mgwr_skip_message)
} else if (length(train_y) > mgwr_max_train) {
    mgwr_message <- paste0(
        "MGWR was requested but not run because training n = ",
        length(train_y),
        " exceeds REAL_DEMO_MGWR_MAX_TRAIN = ",
        mgwr_max_train,
        ". Use a smaller REAL_DEMO_SAMPLE_SIZE or increase ",
        "REAL_DEMO_MGWR_MAX_TRAIN to attempt it."
    )
    cat("   MGWR not run: ", mgwr_message, "\n", sep = "")
    mgwr_row <- status_row("MGWR", "Not run", mgwr_message)
} else {
    mgwr_row <- run_method(
        "MGWR",
        test_y,
        function() {
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
            pred_mgwr(mgwr_fit, test_X, test_coords)
        }
    )
}

cat("\n4. RESULTS: Test-set comparison\n")
cat("   ----------------------------\n")

comparison <- rbind(
    bsgl_row,
    gsvc_row,
    gam_row,
    mgwr_row
)

comparison$Test_MSE <- round(comparison$Test_MSE, 8)
comparison$Test_MAE <- round(comparison$Test_MAE, 8)

print(comparison, row.names = FALSE, na.print = "")

dir.create("results", recursive = TRUE, showWarnings = FALSE)
out_file <- file.path("results", "real_demo_method_comparison.csv")
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
