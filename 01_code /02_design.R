# Trial structures, information architectures, random DGM parameters, and targets.

source(file.path(CODE_DIR, "00_config.R"), local = FALSE)

make_trial_table <- function(setting = c("pairwise", "nma")) {
  setting <- match.arg(setting)
  if (setting == "pairwise") {
    tibble::tibble(
      setting = setting,
      trial_index = seq_len(N_PAIRWISE_TRIALS),
      Trial = sprintf("PW-%02d", seq_len(N_PAIRWISE_TRIALS)),
      Edge = "AB", Trt1 = "A", Trt2 = "B", N = TRIAL_N
    )
  } else {
    dplyr::bind_rows(
      tibble::tibble(setting = setting, trial_index = seq_len(N_NMA_AB), Trial = sprintf("AB-%02d", seq_len(N_NMA_AB)), Edge = "AB", Trt1 = "A", Trt2 = "B", N = TRIAL_N),
      tibble::tibble(setting = setting, trial_index = N_NMA_AB + seq_len(N_NMA_AC), Trial = sprintf("AC-%02d", seq_len(N_NMA_AC)), Edge = "AC", Trt1 = "A", Trt2 = "C", N = TRIAL_N)
    )
  }
}

arch_levels <- c("A0_Full", "A1_Common75", "A2_Balanced75", "A3_Structured75",
                 "A4_Common50", "A5_Balanced50", "A6_Structured50")

setv <- function(...) intersect(ALL_VARS, c(...))

availability_for_pairwise <- function(arch, trial_index) {
  switch(arch,
    A0_Full = ALL_VARS,
    A1_Common75 = setv("age65", "htn", "male"),
    A2_Balanced75 = if (trial_index <= 4L) setv("age65", "htn", "male") else if (trial_index <= 8L) setv("age65", "htn", "smoke") else setv("age65", "male", "smoke"),
    A3_Structured75 = if (trial_index <= 6L) setv("age65", "htn", "male") else setv("age65", "male", "smoke"),
    A4_Common50 = setv("age65", "htn"),
    A5_Balanced50 = if (trial_index <= 4L) setv("age65", "htn") else if (trial_index <= 8L) setv("age65", "male") else setv("age65", "smoke"),
    A6_Structured50 = if (trial_index <= 6L) setv("age65", "htn") else setv("male", "smoke"),
    stop("Unknown architecture: ", arch)
  )
}

availability_for_nma <- function(arch, edge, within_edge_index) {
  switch(arch,
    A0_Full = ALL_VARS,
    A1_Common75 = setv("age65", "htn", "male"),
    A2_Balanced75 = {
      group <- ceiling(within_edge_index / 2)
      list(setv("age65", "htn", "male"), setv("age65", "htn", "smoke"), setv("age65", "male", "smoke"))[[group]]
    },
    A3_Structured75 = if (edge == "AB") setv("age65", "htn", "male") else setv("age65", "male", "smoke"),
    A4_Common50 = setv("age65", "htn"),
    A5_Balanced50 = {
      group <- ceiling(within_edge_index / 2)
      list(setv("age65", "htn"), setv("age65", "male"), setv("age65", "smoke"))[[group]]
    },
    A6_Structured50 = if (edge == "AB") setv("age65", "htn") else setv("male", "smoke"),
    stop("Unknown architecture: ", arch)
  )
}

make_availability_table <- function(setting = c("pairwise", "nma")) {
  setting <- match.arg(setting)
  trials <- make_trial_table(setting)
  purrr::map_dfr(arch_levels, function(arch) {
    out <- trials
    if (setting == "pairwise") {
      out$available_vars <- lapply(out$trial_index, function(i) availability_for_pairwise(arch, i))
    } else {
      within <- ave(out$trial_index, out$Edge, FUN = seq_along)
      out$available_vars <- Map(function(e, i) availability_for_nma(arch, e, i), out$Edge, within)
    }
    out$architecture <- arch
    out$coverage_n <- lengths(out$available_vars)
    out$available_label <- vapply(out$available_vars, paste, collapse = "|", FUN.VALUE = character(1))
    out
  })
}

common_vars_for_architecture <- function(availability_arch) {
  if (!nrow(availability_arch)) return(character(0))
  Reduce(intersect, availability_arch$available_vars)
}

make_scenario_grid <- function() {
  g <- tidyr::crossing(
    outcome = c("binary", "continuous"),
    architecture = arch_levels,
    target_label = names(TARGET_DELTAS),
    setting = c("pairwise", "nma")
  ) %>%
    dplyr::mutate(
      target_delta = unname(TARGET_DELTAS[target_label]),
      outcome_order = match(outcome, c("binary", "continuous")),
      arch_order = match(architecture, arch_levels),
      target_order = match(target_label, names(TARGET_DELTAS)),
      setting_order = match(setting, c("pairwise", "nma"))
    ) %>%
    dplyr::arrange(outcome_order, arch_order, target_order, setting_order) %>%
    dplyr::mutate(scenario_id = sprintf("SC%03d", dplyr::row_number()), .before = 1) %>%
    dplyr::select(-dplyr::ends_with("_order"))
  stopifnot(nrow(g) == 112L)
  g
}

draw_effect_modifiers <- function(rep_id) {
  set.seed(seed_for(rep_id, "effect_modifiers"))
  beta_B <- stats::setNames(stats::runif(length(ALL_VARS), EM_RANGE[1], EM_RANGE[2]), ALL_VARS)
  list(
    A = stats::setNames(rep(0, length(ALL_VARS)), ALL_VARS),
    B = beta_B,
    C = BETA_C_SCALE * beta_B
  )
}

draw_trial_margins <- function(rep_id, setting) {
  trials <- make_trial_table(setting)
  set.seed(seed_for(rep_id, paste0("trial_margins_", setting)))
  P <- matrix(stats::runif(nrow(trials) * length(ALL_VARS), COVARIATE_RANGE[1], COVARIATE_RANGE[2]),
              nrow = nrow(trials), ncol = length(ALL_VARS), byrow = TRUE,
              dimnames = list(trials$Trial, ALL_VARS))
  dplyr::bind_cols(trials, tibble::as_tibble(P))
}

make_targets_from_trials <- function(trial_covariate_means) {
  stopifnot(all(c("Trial", "N", ALL_VARS) %in% names(trial_covariate_means)))
  ref <- vapply(ALL_VARS, function(v) stats::weighted.mean(trial_covariate_means[[v]], trial_covariate_means$N), numeric(1))
  z0 <- stats::qnorm(clip01(ref, c(1e-6, 1 - 1e-6)))
  direction <- rep(1, length(ALL_VARS)); direction <- direction / sqrt(sum(direction^2))
  rows <- lapply(seq_along(TARGET_DELTAS), function(i) {
    p <- stats::pnorm(z0 + unname(TARGET_DELTAS[i]) * direction)
    p <- clip01(p, TARGET_CLIP)
    tibble::tibble(target_id = i, target_label = names(TARGET_DELTAS)[i], target_delta = unname(TARGET_DELTAS[i]),
                   !!!stats::setNames(as.list(p), ALL_VARS))
  })
  list(reference = stats::setNames(ref, ALL_VARS), targets = dplyr::bind_rows(rows))
}

continuous_target_truth <- function(targets, beta, setting) {
  Z <- as.matrix(targets[, ALL_VARS, drop = FALSE]) - CENTER_VALUE
  if (setting == "pairwise") {
    truth <- as.numeric(TAU["B"] - TAU["A"] + Z %*% (beta$B - beta$A))
    contrast <- "B-A"
  } else {
    truth <- as.numeric(TAU["B"] - TAU["C"] + Z %*% (beta$B - beta$C))
    contrast <- "B-C"
  }
  targets %>% dplyr::transmute(target_id, target_label, target_delta, contrast = contrast, truth = truth)
}

# Marginal arm risks and marginal probits for the target joint distribution.
# Target joint cell probabilities are obtained by IPF from the same external
# 16-cell framework used elsewhere in the simulation.
binary_target_arm_parameters <- function(targets, beta, cell_table) {
  if (missing(cell_table) || is.null(cell_table)) stop("cell_table is required for binary marginal truth")
  Xc <- as.matrix(cell_table[, ALL_VARS, drop = FALSE]) - CENTER_VALUE
  treatments <- c("A", "B", "C")
  purrr::map_dfr(seq_len(nrow(targets)), function(i) {
    margins <- stats::setNames(as.numeric(targets[i, ALL_VARS, drop = TRUE]), ALL_VARS)
    target_joint <- ipf_cell_prob(cell_table, margins)
    purrr::map_dfr(treatments, function(trt) {
      cell_risk <- stats::pnorm(ALPHA_BINARY + unname(TAU[trt]) + as.numeric(Xc %*% beta[[trt]]))
      p_target <- sum(target_joint$prob * cell_risk)
      tibble::tibble(
        target_id = targets$target_id[i], target_label = targets$target_label[i],
        target_delta = targets$target_delta[i], treatment = trt,
        event_probability = p_target, marginal_probit = stats::qnorm(p_target),
        target_ipf_converged = target_joint$converged,
        target_ipf_max_abs_error = target_joint$max_abs_error
      )
    })
  })
}

binary_target_truth <- function(targets, beta, setting, cell_table) {
  arms <- binary_target_arm_parameters(targets, beta, cell_table)
  wide <- arms %>%
    dplyr::select(target_id, target_label, target_delta, treatment, marginal_probit) %>%
    tidyr::pivot_wider(names_from = treatment, values_from = marginal_probit)
  if (setting == "pairwise") {
    wide %>% dplyr::transmute(target_id, target_label, target_delta, contrast = "B-A", truth = B - A)
  } else {
    # nma_bc() returns d_BA - d_CA, which is the B-C contrast.
    wide %>% dplyr::transmute(target_id, target_label, target_delta, contrast = "B-C", truth = B - C)
  }
}

truth_original_scale <- function(targets, beta, setting, outcome, cell_table = NULL) {
  outcome <- match.arg(outcome, c("binary", "continuous"))
  if (outcome == "binary") binary_target_truth(targets, beta, setting, cell_table) else continuous_target_truth(targets, beta, setting)
}
