# One complete Monte Carlo replication with paired data across all 112 scenarios.

source(file.path(CODE_DIR, "00_config.R"), local = FALSE)
source(file.path(CODE_DIR, "01_external_data.R"), local = FALSE)
source(file.path(CODE_DIR, "02_design.R"), local = FALSE)
source(file.path(CODE_DIR, "03_generate_data.R"), local = FALSE)
source(file.path(CODE_DIR, "04_eacc_methods.R"), local = FALSE)

architecture_varsets <- function(availability_arch) {
  stats::setNames(availability_arch$available_vars, availability_arch$Trial)
}

varset_map_key <- function(varsets) {
  paste(vapply(sort(names(varsets)), function(tr) paste0(tr, "=", paste(varsets[[tr]], collapse = "+")), character(1)), collapse = ";")
}

label_cached_summary <- function(x, method, architecture, outcome, setting) {
  x %>% dplyr::mutate(
    method = .env$method,
    architecture = .env$architecture,
    outcome = .env$outcome,
    setting = .env$setting
  )
}

write_replication_inputs <- function(rep_dir, rep_id, beta, environments) {
  beta_tbl <- dplyr::bind_rows(lapply(names(beta), function(t) tibble::tibble(rep = rep_id, treatment = t, variable = names(beta[[t]]), beta = as.numeric(beta[[t]]))))
  readr::write_csv(beta_tbl, file.path(rep_dir, "effect_modifiers.csv"))
  for (setting in names(environments)) {
    env <- environments[[setting]]
    readr::write_csv(env$trial_targets, file.path(rep_dir, paste0("trial_target_margins_", setting, ".csv")))
    readr::write_csv(env$trial_covariate_means, file.path(rep_dir, paste0("trial_observed_margins_", setting, ".csv")))
    ref <- tibble::tibble(setting = setting, variable = names(env$targets$reference), reference_mean = as.numeric(env$targets$reference))
    readr::write_csv(ref, file.path(rep_dir, paste0("reference_mean_", setting, ".csv")))
    readr::write_csv(env$targets$targets, file.path(rep_dir, paste0("target_populations_", setting, ".csv")))
    readr::write_csv(env$diagnostics, file.path(rep_dir, paste0("generation_diagnostics_", setting, ".csv")))
    if (SAVE_IPD && !is.null(env$ipd)) {
      ipd_dir <- file.path(rep_dir, "ipd", setting); dir.create(ipd_dir, recursive = TRUE, showWarnings = FALSE)
      for (tr in names(env$ipd)) readr::write_csv(env$ipd[[tr]], file.path(ipd_dir, paste0(tr, ".csv")))
    }
  }
}

run_one_replication <- function(rep_id, pool, cell_table, b_boot = B_BOOT, write_output = TRUE) {
  started <- Sys.time()
  beta <- draw_effect_modifiers(rep_id)
  environments <- list(
    pairwise = simulate_environment(rep_id, "pairwise", beta, cell_table),
    nma = simulate_environment(rep_id, "nma", beta, cell_table)
  )
  availability <- list(pairwise = make_availability_table("pairwise"), nma = make_availability_table("nma"))
  scenario_grid <- make_scenario_grid()
  result_rows <- list(); boot_draw_rows <- list()
  external_source_diagnostics <- list(); external_raking_diagnostics <- list()

  for (setting in c("pairwise", "nma")) {
    env <- environments[[setting]]
    targets <- env$targets$targets
    # One shared external bootstrap source per draw and setting.  Trial-specific
    # re-raking is completed once here and reused across outcomes/architectures.
    shared_external <- prepare_shared_fullboot_external(
      cell_table = cell_table, trial_covariate_means = env$trial_covariate_means,
      rep_id = rep_id, setting = setting, b_boot = b_boot
    )
    external_source_diagnostics[[setting]] <- shared_external$source_diagnostics
    external_raking_diagnostics[[setting]] <- shared_external$raking_diagnostics
    for (outcome in c("binary", "continuous")) {
      agd <- env$agd[[outcome]]
      cache <- new.env(parent = emptyenv())
      for (arch in arch_levels) {
        avail_arch <- availability[[setting]] %>% dplyr::filter(architecture == arch)

        # Standard synthesis is architecture-invariant but is labelled for every scenario.
        result_rows[[length(result_rows) + 1L]] <- run_standard_synthesis(agd, targets, arch, outcome, setting)

        # EACC uses exactly the covariates available in each trial under this
        # architecture.  Architecture is the information-structure dimension;
        # it is not crossed with an additional M1/M2/M3 strategy dimension.
        method <- if (setting == "pairwise") "EACC-Meta" else "EACC-NMA"
        varsets <- architecture_varsets(avail_arch)
        key <- varset_map_key(varsets)
        if (!exists(key, envir = cache, inherits = FALSE)) {
          fb <- run_fullboot_eacc(
            agd = agd, shared_external_draws = shared_external$draws, targets = targets,
            varsets_by_trial = varsets, method = method, architecture = arch,
            outcome = outcome, setting = setting, rep_id = rep_id, b_boot = b_boot
          )
          assign(key, fb, envir = cache)
        }
        fb <- get(key, envir = cache, inherits = FALSE)
        result_rows[[length(result_rows) + 1L]] <- label_cached_summary(fb$summary, method, arch, outcome, setting)
        if (SAVE_BOOT_DRAWS && !is.null(fb$draws)) {
          boot_draw_rows[[length(boot_draw_rows) + 1L]] <- fb$draws %>%
            dplyr::mutate(rep = rep_id, method = method, architecture = arch, outcome = outcome, setting = setting, .before = 1)
        }
      }
    }
  }

  results <- dplyr::bind_rows(result_rows) %>%
    dplyr::mutate(target_label_join = names(TARGET_DELTAS)[target_id]) %>%
    dplyr::left_join(
      scenario_grid %>% dplyr::select(scenario_id, outcome, architecture, target_label, target_delta, setting),
      by = c("outcome", "architecture", "setting", "target_label_join" = "target_label")
    ) %>%
    dplyr::select(-target_label_join, -target_delta) %>%
    dplyr::left_join(
      dplyr::bind_rows(lapply(names(environments), function(setting) {
        dplyr::bind_rows(lapply(c("binary", "continuous"), function(outcome) {
          truth_original_scale(
            environments[[setting]]$targets$targets, beta, setting,
            outcome = outcome, cell_table = cell_table
          ) %>% dplyr::mutate(setting = setting, outcome = outcome)
        }))
      })) %>% dplyr::select(setting, outcome, target_id, target_label, target_delta, contrast, truth),
      by = c("setting", "outcome", "target_id")
    ) %>%
    dplyr::mutate(
      rep = rep_id, bias = estimate - truth, abs_bias = abs(bias), sq_error = bias^2,
      covered = dplyr::if_else(is.finite(lwr) & is.finite(upr), as.integer(lwr <= truth & truth <= upr), NA_integer_),
      width = upr - lwr, failed = as.integer(!is.finite(estimate) | !is.finite(se)),
      .before = 1
    ) %>%
    dplyr::mutate(method_order = dplyr::if_else(method %in% c("Meta", "NMA"), 1L, 2L)) %>%
    dplyr::arrange(scenario_id, method_order) %>%
    dplyr::select(-method_order)

  rep_dir <- file.path(RUN_ROOT, sprintf("REP_%04d", rep_id))
  if (isTRUE(write_output)) {
    dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
    write_replication_inputs(rep_dir, rep_id, beta, environments)
    readr::write_csv(results, file.path(rep_dir, "performance_rows.csv"))
    readr::write_csv(dplyr::bind_rows(external_source_diagnostics), file.path(rep_dir, "fullboot_shared_external_sources.csv"))
    readr::write_csv(dplyr::bind_rows(external_raking_diagnostics), file.path(rep_dir, "fullboot_trial_raking_diagnostics.csv"))
    if (SAVE_BOOT_DRAWS && length(boot_draw_rows)) readr::write_csv(dplyr::bind_rows(boot_draw_rows), file.path(rep_dir, "fullboot_draws.csv"))
    manifest <- tibble::tibble(
      rep = rep_id, code_version = CODE_VERSION, b_boot = b_boot, trial_n = TRIAL_N,
      started_at = format(started, "%Y-%m-%d %H:%M:%S"),
      finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
      result_rows = nrow(results), failed_rows = sum(results$failed, na.rm = TRUE)
    )
    readr::write_csv(manifest, file.path(rep_dir, "replication_manifest.csv"))
  }
  list(results = results, beta = beta, environments = environments,
       external_source_diagnostics = dplyr::bind_rows(external_source_diagnostics),
       external_raking_diagnostics = dplyr::bind_rows(external_raking_diagnostics),
       boot_draws = if (length(boot_draw_rows)) dplyr::bind_rows(boot_draw_rows) else NULL)
}
