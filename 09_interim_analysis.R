# Read-only interim analysis for a still-running simulation.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

args <- commandArgs(TRUE)
if (!length(args)) stop("Usage: Rscript 09_interim_analysis.R <run_root>")
run_root <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
files <- list.files(run_root, pattern = "performance_rows\\.csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No completed performance_rows.csv files under: ", run_root)

read_complete <- function(path) {
  x <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(x) || nrow(x) != 224L || dplyr::n_distinct(x$rep) != 1L) return(NULL)
  x
}
pieces <- lapply(files, read_complete)
pieces <- pieces[!vapply(pieces, is.null, logical(1))]
if (!length(pieces)) stop("No structurally complete replication files found")
dat <- dplyr::bind_rows(pieces)

scenario <- dat %>%
  group_by(scenario_id, outcome, architecture, target_label, target_delta, setting, contrast, method) %>%
  summarise(
    n_rep = n_distinct(rep),
    n_valid = sum(is.finite(estimate)),
    mean_truth = mean(truth, na.rm = TRUE),
    mean_estimate = mean(estimate, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    mean_abs_error = mean(abs_bias, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    coverage = mean(covered, na.rm = TRUE),
    mean_width = mean(width, na.rm = TRUE),
    empirical_sd = sd(estimate, na.rm = TRUE),
    mean_model_se = mean(se, na.rm = TRUE),
    .groups = "drop"
  )

standard_scenario <- scenario %>%
  filter(method %in% c("Meta", "NMA")) %>%
  select(scenario_id, rmse_standard = rmse, mae_standard = mean_abs_error)
scenario_cmp <- scenario %>%
  filter(method %in% c("EACC-Meta", "EACC-NMA"), n_valid > 0) %>%
  left_join(standard_scenario, by = "scenario_id") %>%
  mutate(
    rmse_ratio = rmse / rmse_standard,
    mae_reduction = mae_standard - mean_abs_error,
    better_rmse = rmse < rmse_standard,
    better_mae = mean_abs_error < mae_standard
  )

overall <- dat %>%
  group_by(outcome, setting, method) %>%
  summarise(
    rows = n(), valid_rows = sum(is.finite(estimate)),
    mean_bias = mean(bias, na.rm = TRUE),
    mean_abs_error = mean(abs_bias, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    coverage = mean(covered, na.rm = TRUE),
    mean_width = mean(width, na.rm = TRUE),
    empirical_sd = sd(estimate, na.rm = TRUE),
    mean_model_se = mean(se, na.rm = TRUE),
    .groups = "drop"
  )

standard_rows <- dat %>%
  filter(method %in% c("Meta", "NMA")) %>%
  select(rep, scenario_id, abs_bias_standard = abs_bias, sq_error_standard = sq_error,
         covered_standard = covered, width_standard = width)
paired <- dat %>%
  filter(method %in% c("EACC-Meta", "EACC-NMA"), is.finite(estimate)) %>%
  left_join(standard_rows, by = c("rep", "scenario_id"))

by_target <- paired %>%
  group_by(outcome, setting, target_label, target_delta, method) %>%
  summarise(
    n = n(),
    rmse_eacc = sqrt(mean(sq_error, na.rm = TRUE)),
    rmse_standard = sqrt(mean(sq_error_standard, na.rm = TRUE)),
    rmse_ratio = rmse_eacc / rmse_standard,
    mae_reduction = mean(abs_bias_standard - abs_bias, na.rm = TRUE),
    paired_abs_error_win_rate = mean(abs_bias < abs_bias_standard, na.rm = TRUE),
    coverage_eacc = mean(covered, na.rm = TRUE),
    coverage_standard = mean(covered_standard, na.rm = TRUE),
    width_ratio = mean(width / width_standard, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(outcome, setting, method, target_delta)

scenario_wins <- scenario_cmp %>%
  group_by(outcome, setting, method) %>%
  summarise(
    valid_scenarios = n(),
    scenarios_better_rmse = sum(better_rmse),
    proportion_better_rmse = mean(better_rmse),
    scenarios_better_mae = sum(better_mae),
    proportion_better_mae = mean(better_mae),
    median_rmse_ratio = median(rmse_ratio),
    worst_rmse_ratio = max(rmse_ratio),
    .groups = "drop"
  )

architecture_summary <- scenario_cmp %>%
  group_by(outcome, setting, architecture) %>%
  summarise(
    scenarios = n(), mean_rmse = mean(rmse),
    mean_rmse_ratio_vs_standard = mean(rmse_ratio),
    scenarios_better_rmse = sum(better_rmse),
    mean_coverage = mean(coverage),
    .groups = "drop"
  )

boot_diag <- dat %>%
  filter(method %in% c("EACC-Meta", "EACC-NMA"), is.finite(estimate)) %>%
  group_by(outcome, setting, method) %>%
  summarise(
    min_n_boot_ok = min(n_boot_ok, na.rm = TRUE),
    max_n_boot_ok = max(n_boot_ok, na.rm = TRUE),
    mean_within_W = mean(W, na.rm = TRUE),
    mean_between_Bvar = mean(Bvar, na.rm = TRUE),
    between_share_T = mean(Bvar / Tvar, na.rm = TRUE),
    .groups = "drop"
  )

target_levels <- dat %>%
  group_by(outcome, setting, target_label, target_delta, method) %>%
  summarise(
    mean_truth = mean(truth, na.rm = TRUE),
    mean_estimate = mean(estimate, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    coverage = mean(covered, na.rm = TRUE),
    mean_width = mean(width, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(outcome, setting, target_delta, method)

coverage_diag <- dat %>%
  filter(is.finite(estimate)) %>%
  group_by(scenario_id, outcome, setting, method) %>%
  summarise(
    scenario_coverage = mean(covered, na.rm = TRUE),
    empirical_error_sd = sd(bias, na.rm = TRUE),
    mean_model_se = mean(se, na.rm = TRUE),
    se_to_error_sd = mean_model_se / empirical_error_sd,
    .groups = "drop"
  ) %>%
  group_by(outcome, setting, method) %>%
  summarise(
    scenarios = n(),
    mean_coverage = mean(scenario_coverage),
    median_coverage = median(scenario_coverage),
    min_coverage = min(scenario_coverage),
    max_coverage = max(scenario_coverage),
    median_se_to_error_sd = median(se_to_error_sd, na.rm = TRUE),
    .groups = "drop"
  )

worst <- scenario_cmp %>%
  arrange(desc(rmse_ratio)) %>%
  select(scenario_id, outcome, setting, architecture, target_label, method,
         n_rep, rmse, rmse_standard, rmse_ratio, coverage, mean_bias) %>%
  slice_head(n = 20)

cat("RUN_ROOT=", run_root, "\n", sep = "")
cat("VALID_REPLICATIONS=", n_distinct(dat$rep), "\n", sep = "")
cat("REP_IDS=", paste(sort(unique(dat$rep)), collapse = ","), "\n", sep = "")
cat("TOTAL_ROWS=", nrow(dat), "\n", sep = "")
cat("STRUCTURAL_FAILURE_ROWS=", sum(dat$failed, na.rm = TRUE), "\n\n", sep = "")

options(tibble.width = Inf, dplyr.width = Inf)
cat("=== OVERALL ===\n"); print(overall, n = Inf, width = Inf)
cat("\n=== SCENARIO-LEVEL EACC SUPERIORITY COUNTS ===\n"); print(scenario_wins, n = Inf, width = Inf)
cat("\n=== EXTRAPOLATION/TARGET COMPARISON ===\n"); print(by_target, n = Inf, width = Inf)
cat("\n=== EACC RESULTS BY INFORMATION ARCHITECTURE ===\n"); print(architecture_summary, n = Inf, width = Inf)
cat("\n=== FULLBOOT DIAGNOSTICS ===\n"); print(boot_diag, n = Inf, width = Inf)
cat("\n=== TARGET-SCALE TRUTH/ESTIMATE CHECK ===\n"); print(target_levels, n = Inf, width = Inf)
cat("\n=== COVERAGE AND SE CALIBRATION ===\n"); print(coverage_diag, n = Inf, width = Inf)
cat("\n=== 20 WORST EACC/STANDARD RMSE RATIOS ===\n"); print(worst, n = Inf, width = Inf)
