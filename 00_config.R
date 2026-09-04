# Frozen configuration for the heterogeneous trial-specific covariate information study.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  required <- c("dplyr", "tibble", "readr", "purrr", "tidyr", "MASS", "digest")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  invisible(lapply(required, library, character.only = TRUE))
})

get_script_dir <- function() {
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(dirname(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE)))
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

if (!exists("CODE_DIR", inherits = FALSE)) CODE_DIR <- get_script_dir()
STUDY_ROOT <- normalizePath(file.path(CODE_DIR, ".."), winslash = "/", mustWork = FALSE)
DESIGN_DIR <- file.path(STUDY_ROOT, "00_design")
SIM_DATA_DIR <- file.path(STUDY_ROOT, "02_simulated_data")
SUMMARY_DIR <- file.path(STUDY_ROOT, "03_results_summary")
REPLICATION_DIR <- file.path(STUDY_ROOT, "04_results_replication")
LOG_DIR <- file.path(STUDY_ROOT, "06_logs")
VALIDATION_DIR <- file.path(STUDY_ROOT, "07_validation")
invisible(lapply(c(SIM_DATA_DIR, SUMMARY_DIR, REPLICATION_DIR, LOG_DIR, VALIDATION_DIR), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

RWD_DIR_DEFAULT <- file.path(STUDY_ROOT, "02_simulated_data", "source_data")
RWD_DIR <- Sys.getenv("EACC_RWD_DIR", RWD_DIR_DEFAULT)
RWD_FILE <- Sys.getenv("EACC_RWD_FILE", file.path(RWD_DIR, "final_analysis_dataset_CVD_and_DM.csv"))

read_int_env <- function(name, default) {
  x <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (!is.finite(x) || is.na(x)) as.integer(default) else x
}
read_num_env <- function(name, default) {
  x <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (!is.finite(x) || is.na(x)) as.numeric(default) else x
}
env_flag <- function(name, default = FALSE) {
  x <- tolower(trimws(Sys.getenv(name, ifelse(default, "true", "false"))))
  x %in% c("1", "true", "t", "yes", "y")
}

CODE_VERSION <- "heterogeneous_info_v3_shared_rematching_fullboot"
ALL_VARS <- c("age65", "htn", "male", "smoke")
VAR_LABELS <- stats::setNames(paste0("X", seq_along(ALL_VARS)), ALL_VARS)

# Frozen DGM parameters.
TRIAL_N <- 8000L
ARM_N <- 4000L
N_PAIRWISE_TRIALS <- 12L
N_NMA_AB <- 6L
N_NMA_AC <- 6L
COVARIATE_RANGE <- c(0.30, 0.70)
EM_RANGE <- c(0.10, 0.90)
CENTER_VALUE <- 0.50
BASELINE_RISK_BINARY <- 0.10
ALPHA_BINARY <- stats::qnorm(BASELINE_RISK_BINARY)
TAU <- c(A = 0, B = -0.201, C = -0.3702)
BETA_C_SCALE <- 0.50
CONTINUOUS_SD <- 1.0
TARGET_DELTAS <- c(D0 = 0.0, D0p5 = 0.5, D0p8 = 0.8, D1p2 = 1.2)
TARGET_CLIP <- c(0.02, 0.98)

# Frozen Monte Carlo and Fullboot parameters; environment overrides support smoke tests.
MC_REPS <- read_int_env("EACC_MC_REPS", 200L)
B_BOOT <- read_int_env("EACC_B_BOOT", 200L)
BOOT_CORES <- max(1L, read_int_env("EACC_BOOT_CORES", 3L))
OUTER_CHUNKS <- max(1L, read_int_env("EACC_OUTER_CHUNKS", 7L))
SAVE_IPD <- env_flag("EACC_SAVE_IPD", FALSE)
SAVE_BOOT_DRAWS <- env_flag("EACC_SAVE_BOOT_DRAWS", FALSE)
RUN_TAG <- Sys.getenv("EACC_RUN_TAG", format(Sys.time(), "RUN_%Y%m%d_%H%M%S"))
RUN_ROOT <- Sys.getenv("EACC_RUN_ROOT", file.path(REPLICATION_DIR, RUN_TAG))
dir.create(RUN_ROOT, recursive = TRUE, showWarnings = FALSE)

# Reproducible common-random-number streams.
SEED_MASTER <- read_int_env("EACC_SEED_MASTER", 202608190L)
seed_for <- function(rep_id, component, offset = 0L) {
  key <- paste(SEED_MASTER, as.integer(rep_id), as.character(component), as.integer(offset), sep = "|")
  raw <- digest::digest(key, algo = "xxhash32", serialize = FALSE)
  val <- suppressWarnings(strtoi(substr(raw, 1, 7), base = 16L))
  as.integer((val %% (.Machine$integer.max - 1L)) + 1L)
}

# Legacy shrinkage controls retained only for compatibility with archived code.
# The v3 shared-source full re-matching Fullboot path does not read either value.
ETA_SIGMA_BOOT <- read_num_env("EACC_ETA_SIGMA_BOOT", 0.8)
LAMBDA_MU_BOOT <- read_num_env("EACC_LAMBDA_MU_BOOT", 0.0)
FULLBOOT_BLUP_MODE <- "mean_scaled"
BLUP_RANDOM_SCALE <- 1.0
BLUP_DESIGN_EFFECT <- 1.0
BLUP_N_EFF_MIN <- 1.0
EACC_IMP_VAR_SCALE_MODE <- "by_imputed_vars"
FULLBOOT_EACC_IMP_VAR_SYNC <- TRUE
N_MATCH <- read_int_env("EACC_N_MATCH", TRIAL_N)

# IPF on the 16 binary covariate cells is algebraically equivalent to individual-level
# raking because all individuals in a cell receive the same calibration multiplier.
IPF_MAXIT <- 500L
IPF_TOL <- 1e-10

options(future.globals.maxSize = 16 * 1024^3, future.rng.onMisuse = "ignore")
Sys.setenv(
  OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)

clip01 <- function(x, limits = TARGET_CLIP) pmin(pmax(as.numeric(x), limits[1]), limits[2])
safe_solve <- function(A, b = NULL) {
  A <- as.matrix(A)
  tryCatch(if (is.null(b)) solve(A) else solve(A, b), error = function(e) {
    G <- MASS::ginv(A)
    if (is.null(b)) G else G %*% b
  })
}

message(sprintf(
  "[%s] MC_REPS=%d B_BOOT=%d trial_N=%d chunks=%d cores/chunk=%d run=%s",
  CODE_VERSION, MC_REPS, B_BOOT, TRIAL_N, OUTER_CHUNKS, BOOT_CORES, RUN_TAG
))
