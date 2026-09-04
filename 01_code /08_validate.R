# Structural, estimand, shared-Fullboot, numerical, and reproducibility validation.

CODE_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"))
source(file.path(CODE_DIR, "00_config.R"), local = FALSE)
source(file.path(CODE_DIR, "01_external_data.R"), local = FALSE)
source(file.path(CODE_DIR, "02_design.R"), local = FALSE)
source(file.path(CODE_DIR, "04_eacc_methods.R"), local = FALSE)

args <- commandArgs(TRUE)
validation_run <- if (length(args)) normalizePath(args[1], winslash = "/", mustWork = TRUE) else RUN_ROOT
repro_run <- if (length(args) >= 2L) normalizePath(args[2], winslash = "/", mustWork = TRUE) else NA_character_
find_rep_dir <- function(root) {
  x <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  x <- sort(x[grepl("REP_[0-9]+$", x)])
  if (!length(x)) stop("No replication directory found under: ", root)
  x[1]
}
rep_dir <- find_rep_dir(validation_run)

perf <- readr::read_csv(file.path(rep_dir, "performance_rows.csv"), show_col_types = FALSE)
beta_tbl <- readr::read_csv(file.path(rep_dir, "effect_modifiers.csv"), show_col_types = FALSE)
manifest <- readr::read_csv(file.path(rep_dir, "replication_manifest.csv"), show_col_types = FALSE)
source_diag <- readr::read_csv(file.path(rep_dir, "fullboot_shared_external_sources.csv"), show_col_types = FALSE)
raking_diag <- readr::read_csv(file.path(rep_dir, "fullboot_trial_raking_diagnostics.csv"), show_col_types = FALSE)
expected_b <- as.integer(manifest$b_boot[1])
rep_id <- as.integer(manifest$rep[1])

checks <- list()
add_check <- function(name, pass, detail = "") {
  checks[[length(checks) + 1L]] <<- tibble::tibble(check = name, pass = isTRUE(pass), detail = as.character(detail))
}

add_check("scenario_grid_has_112_rows", nrow(make_scenario_grid()) == 112L)
add_check("replication_has_224_method_rows", nrow(perf) == 224L, nrow(perf))
add_check("replication_has_112_scenarios", dplyr::n_distinct(perf$scenario_id) == 112L)
add_check("each_scenario_has_two_methods", all(table(perf$scenario_id) == 2L), paste(range(table(perf$scenario_id)), collapse = ","))
method_counts <- table(perf$method)
expected_method_counts <- c("Meta" = 56L, "EACC-Meta" = 56L, "NMA" = 56L, "EACC-NMA" = 56L)
method_counts_ok <- setequal(names(method_counts), names(expected_method_counts)) &&
  all(as.integer(method_counts[names(expected_method_counts)]) == unname(expected_method_counts))
add_check("four_setting_specific_method_labels_each_have_56_rows", method_counts_ok, paste(names(method_counts), method_counts, collapse = ","))
add_check("pairwise_only_compares_Meta_and_EACC_Meta", setequal(unique(perf$method[perf$setting == "pairwise"]), c("Meta", "EACC-Meta")))
add_check("nma_only_compares_NMA_and_EACC_NMA", setequal(unique(perf$method[perf$setting == "nma"]), c("NMA", "EACC-NMA")))
add_check("no_structural_failures_including_A6", sum(perf$failed, na.rm = TRUE) == 0L, sum(perf$failed, na.rm = TRUE))

beta <- stats::setNames(lapply(c("A", "B", "C"), function(trt) {
  z <- beta_tbl[beta_tbl$treatment == trt, , drop = FALSE]
  stats::setNames(z$beta[match(ALL_VARS, z$variable)], ALL_VARS)
}), c("A", "B", "C"))
bB <- beta_tbl %>% dplyr::filter(treatment == "B")
bC <- beta_tbl %>% dplyr::filter(treatment == "C")
add_check("beta_B_within_0p1_0p9", all(bB$beta >= EM_RANGE[1] & bB$beta <= EM_RANGE[2]), paste(range(bB$beta), collapse = ","))
add_check("beta_B_not_fixed_at_0p5", stats::sd(bB$beta) > 0)
add_check("beta_C_is_half_beta_B", isTRUE(all.equal(bC$beta[match(bB$variable, bC$variable)], BETA_C_SCALE * bB$beta, tolerance = 1e-12)))

pool <- load_external_pool()
cell_table <- make_cell_table(pool)
binary_prob_rows <- list()
binary_truth_diffs <- numeric(0)
continuous_truth_diffs <- numeric(0)
orientation_checks <- logical(0)

empirical_binary_contrast <- function(target_row, setting, simulation_n = 10000000L) {
  margins <- stats::setNames(as.numeric(target_row[, ALL_VARS, drop = TRUE]), ALL_VARS)
  joint <- ipf_cell_prob(cell_table, margins)
  set.seed(seed_for(rep_id, paste0("validation_empirical_binary_", setting), as.integer(target_row$target_id)))
  n_cell <- as.vector(stats::rmultinom(1, simulation_n, joint$prob))
  Xc <- as.matrix(cell_table[, ALL_VARS, drop = FALSE]) - CENTER_VALUE
  eta <- vapply(c("A", "B", "C"), function(trt) {
    cell_risk <- stats::pnorm(ALPHA_BINARY + unname(TAU[trt]) + as.numeric(Xc %*% beta[[trt]]))
    events <- sum(stats::rbinom(length(n_cell), size = n_cell, prob = cell_risk))
    stats::qnorm(events / simulation_n)
  }, numeric(1))
  if (setting == "pairwise") unname(eta["B"] - eta["A"]) else unname(eta["B"] - eta["C"])
}

for (setting in c("pairwise", "nma")) {
  tm <- readr::read_csv(file.path(rep_dir, paste0("trial_target_margins_", setting, ".csv")), show_col_types = FALSE)
  om <- readr::read_csv(file.path(rep_dir, paste0("trial_observed_margins_", setting, ".csv")), show_col_types = FALSE)
  ref <- readr::read_csv(file.path(rep_dir, paste0("reference_mean_", setting, ".csv")), show_col_types = FALSE)
  tg <- readr::read_csv(file.path(rep_dir, paste0("target_populations_", setting, ".csv")), show_col_types = FALSE)
  target_values <- as.matrix(tm[, ALL_VARS])
  add_check(paste0(setting, "_trial_target_draws_within_range"), all(target_values >= COVARIATE_RANGE[1] & target_values <= COVARIATE_RANGE[2]), paste(range(target_values), collapse = ","))
  add_check(paste0(setting, "_trial_target_draws_not_fixed_0p5"), stats::sd(as.numeric(target_values)) > 0)
  add_check(paste0(setting, "_all_trial_N_8000"), all(om$N == TRIAL_N))
  calc_ref <- colMeans(as.matrix(om[, ALL_VARS]))
  stored_ref <- stats::setNames(ref$reference_mean, ref$variable)[ALL_VARS]
  add_check(paste0(setting, "_reference_is_actual_trial_mean"), isTRUE(all.equal(as.numeric(calc_ref), as.numeric(stored_ref), tolerance = 1e-12)))
  d0 <- as.numeric(tg[tg$target_delta == 0, ALL_VARS])
  add_check(paste0(setting, "_D0_equals_reference"), isTRUE(all.equal(d0, as.numeric(stored_ref), tolerance = 1e-12)))

  arm <- binary_target_arm_parameters(tg, beta, cell_table) %>% dplyr::mutate(setting = setting)
  binary_prob_rows[[setting]] <- arm
  analytic_binary <- binary_target_truth(tg, beta, setting, cell_table)
  empirical_binary <- vapply(seq_len(nrow(tg)), function(i) empirical_binary_contrast(tg[i, , drop = FALSE], setting), numeric(1))
  binary_truth_diffs <- c(binary_truth_diffs, abs(analytic_binary$truth - empirical_binary))

  old_Z <- as.matrix(tg[, ALL_VARS, drop = FALSE]) - CENTER_VALUE
  old_continuous <- if (setting == "pairwise") {
    as.numeric(TAU["B"] - TAU["A"] + old_Z %*% (beta$B - beta$A))
  } else {
    as.numeric(TAU["B"] - TAU["C"] + old_Z %*% (beta$B - beta$C))
  }
  new_continuous <- continuous_target_truth(tg, beta, setting)$truth
  continuous_truth_diffs <- c(continuous_truth_diffs, abs(old_continuous - new_continuous))

  first_arm <- arm %>% dplyr::filter(target_id == min(target_id))
  eta <- stats::setNames(first_arm$marginal_probit, first_arm$treatment)
  if (setting == "pairwise") {
    orientation_checks <- c(orientation_checks, isTRUE(all.equal(unname(fe_pairwise(eta["B"] - eta["A"], 1)["estimate"]), unname(eta["B"] - eta["A"]), tolerance = 1e-12)))
  } else {
    fit <- nma_fe(c(eta["B"] - eta["A"], eta["C"] - eta["A"]), c(1, 1), c("A", "A"), c("B", "C"))
    orientation_checks <- c(orientation_checks, isTRUE(all.equal(unname(nma_bc(fit)["estimate"]), unname(eta["B"] - eta["C"]), tolerance = 1e-12)))
  }
}

binary_probs <- dplyr::bind_rows(binary_prob_rows)
add_check("binary_target_arm_probabilities_inside_0_1_and_away_from_boundaries",
          all(is.finite(binary_probs$event_probability) & binary_probs$event_probability > 1e-6 & binary_probs$event_probability < 1 - 1e-6),
          paste(range(binary_probs$event_probability), collapse = ","))
add_check("binary_16cell_analytic_truth_matches_large_empirical_target",
          max(binary_truth_diffs) < 0.0035, max(binary_truth_diffs))
add_check("all_binary_truth_values_are_finite", all(is.finite(perf$truth[perf$outcome == "binary"])))
add_check("continuous_truth_is_exactly_unchanged", max(continuous_truth_diffs) < 1e-12, max(continuous_truth_diffs))
add_check("pairwise_and_NMA_truth_orientations_match_synthesis", all(orientation_checks), paste(orientation_checks, collapse = ","))

expected_coverage_n <- c(
  A0_Full = 4L, A1_Common75 = 3L, A2_Balanced75 = 3L,
  A3_Structured75 = 3L, A4_Common50 = 2L,
  A5_Balanced50 = 2L, A6_Structured50 = 2L
)
for (setting in c("pairwise", "nma")) {
  av <- make_availability_table(setting)
  coverage_ok <- all(vapply(seq_len(nrow(av)), function(i) {
    length(av$available_vars[[i]]) == expected_coverage_n[[av$architecture[i]]]
  }, logical(1)))
  add_check(paste0(setting, "_architecture_covariate_counts_are_4_3_3_3_2_2_2"), coverage_ok)
  a0 <- av %>% dplyr::filter(architecture == "A0_Full")
  add_check(paste0(setting, "_A0_every_trial_uses_all_four_covariates"), all(vapply(a0$available_vars, setequal, logical(1), y = ALL_VARS)))
}

eacc <- perf %>% dplyr::filter(method %in% c("EACC-Meta", "EACC-NMA"))
add_check("all_EACC_rows_have_finite_intervals", all(is.finite(eacc$estimate) & is.finite(eacc$se) & is.finite(eacc$lwr) & is.finite(eacc$upr)))
add_check("all_EACC_rows_complete_requested_fullboot_B", all(eacc$n_boot_ok == expected_b), paste(range(eacc$n_boot_ok), collapse = ","))
a6_eacc <- eacc %>% dplyr::filter(architecture == "A6_Structured50")
add_check("A6_EACC_is_fully_estimable", nrow(a6_eacc) == 16L && all(is.finite(a6_eacc$estimate) & is.finite(a6_eacc$se)), nrow(a6_eacc))

add_check("one_shared_external_source_row_per_setting_and_draw",
          nrow(source_diag) == 2L * expected_b && all(table(source_diag$setting, source_diag$draw) == 1L), nrow(source_diag))
add_check("bootstrap_external_cell_counts_change_across_draws",
          all(source_diag$changed_from_original) && all(vapply(split(source_diag$shared_source_id, source_diag$setting), dplyr::n_distinct, integer(1)) > 1L))
add_check("every_setting_draw_rerakes_all_12_trials",
          nrow(raking_diag) == 2L * expected_b * 12L && all(table(raking_diag$setting, raking_diag$draw) == 12L), nrow(raking_diag))
add_check("all_trials_within_draw_use_same_shared_source",
          all(vapply(split(raking_diag$shared_source_id, interaction(raking_diag$setting, raking_diag$draw)), dplyr::n_distinct, integer(1)) == 1L))
source_key <- source_diag %>% dplyr::select(setting, draw, expected_source_id = shared_source_id)
raking_key <- raking_diag %>% dplyr::left_join(source_key, by = c("setting", "draw"))
add_check("trial_raking_source_ids_match_shared_source_manifest", all(raking_key$shared_source_id == raking_key$expected_source_id))
add_check("trial_specific_raking_converges_to_tolerance",
          all(raking_diag$ipf_converged & raking_diag$ipf_max_abs_error <= max(1e-9, 10 * IPF_TOL)),
          max(raking_diag$ipf_max_abs_error))
add_check("every_rematched_pseudosample_has_N_MATCH", all(raking_diag$matched_n == N_MATCH), paste(range(raking_diag$matched_n), collapse = ","))

main_fullboot_files <- file.path(CODE_DIR, c("01_external_data.R", "04_eacc_methods.R", "05_run_replication.R"))
main_fullboot_text <- paste(unlist(lapply(main_fullboot_files, readLines, warn = FALSE)), collapse = "\n")
add_check("primary_fullboot_does_not_call_legacy_shrinkage_parameters",
          !grepl("ETA_SIGMA_BOOT|LAMBDA_MU_BOOT", main_fullboot_text))

if (!is.na(repro_run)) {
  repro_rep_dir <- find_rep_dir(repro_run)
  repro_perf <- readr::read_csv(file.path(repro_rep_dir, "performance_rows.csv"), show_col_types = FALSE)
  repro_source <- readr::read_csv(file.path(repro_rep_dir, "fullboot_shared_external_sources.csv"), show_col_types = FALSE)
  repro_raking <- readr::read_csv(file.path(repro_rep_dir, "fullboot_trial_raking_diagnostics.csv"), show_col_types = FALSE)
  exact_repro <- isTRUE(all.equal(perf, repro_perf, tolerance = 0, check.attributes = FALSE)) &&
    isTRUE(all.equal(source_diag, repro_source, tolerance = 0, check.attributes = FALSE)) &&
    isTRUE(all.equal(raking_diag, repro_raking, tolerance = 0, check.attributes = FALSE))
  add_check("same_seed_replication_is_exactly_reproducible", exact_repro)
} else {
  add_check("same_seed_replication_is_exactly_reproducible", FALSE, "Second run root was not supplied")
}

checks <- dplyr::bind_rows(checks)
out <- file.path(VALIDATION_DIR, paste0(basename(validation_run), "_validation_checks.csv"))
readr::write_csv(checks, out)
print(checks, n = Inf)
if (any(!checks$pass)) stop("Validation failed; see: ", out)
message("Validation passed: ", out)
