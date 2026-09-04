# Aggregate replication-level performance files into scenario summaries.

CODE_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"))
source(file.path(CODE_DIR, "00_config.R"), local = FALSE)

files <- list.files(RUN_ROOT, pattern = "performance_rows\\.csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No completed replication result files under: ", RUN_ROOT)
all <- dplyr::bind_rows(lapply(files, readr::read_csv, show_col_types = FALSE))
readr::write_csv(all, file.path(SUMMARY_DIR, paste0(RUN_TAG, "_all_replication_rows.csv")))

summary <- all %>%
  dplyr::group_by(scenario_id, outcome, architecture, target_label, target_delta, setting, contrast, method) %>%
  dplyr::summarise(
    n_planned_rows = dplyr::n(), n_valid = sum(is.finite(estimate)), failure_rate = mean(failed, na.rm = TRUE),
    mean_truth = mean(truth, na.rm = TRUE), mean_estimate = mean(estimate, na.rm = TRUE),
    mean_bias = mean(.data$bias, na.rm = TRUE), mean_absolute_error = mean(abs_bias, na.rm = TRUE),
    empirical_sd = stats::sd(estimate, na.rm = TRUE), mean_model_se = mean(se, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)), coverage = mean(covered, na.rm = TRUE),
    mean_ci_width = mean(width, na.rm = TRUE),
    mcse_bias = stats::sd(.data$bias, na.rm = TRUE) / sqrt(sum(is.finite(.data$bias))),
    mcse_coverage = sqrt(coverage * (1 - coverage) / sum(is.finite(covered))),
    mcse_rmse = ifelse(rmse > 0, stats::sd(sq_error, na.rm = TRUE) / (2 * rmse * sqrt(sum(is.finite(sq_error)))), NA_real_),
    .groups = "drop"
  ) %>%
  dplyr::mutate(absolute_bias = abs(mean_bias), .after = mean_bias)

standard <- summary %>% dplyr::filter(method %in% c("Meta", "NMA")) %>%
  dplyr::select(scenario_id, rmse_standard = rmse, absolute_bias_standard = absolute_bias)
summary <- summary %>%
  dplyr::left_join(standard, by = "scenario_id") %>%
  dplyr::mutate(
    rmse_ratio_vs_standard = rmse / rmse_standard,
    absolute_bias_reduction_vs_standard = absolute_bias_standard - absolute_bias
  ) %>%
  dplyr::mutate(method_order = dplyr::if_else(method %in% c("Meta", "NMA"), 1L, 2L)) %>%
  dplyr::arrange(scenario_id, method_order) %>%
  dplyr::select(-method_order)

out <- file.path(SUMMARY_DIR, paste0(RUN_TAG, "_performance_summary.csv"))
readr::write_csv(summary, out)
message("Aggregated ", length(files), " replications: ", out)
