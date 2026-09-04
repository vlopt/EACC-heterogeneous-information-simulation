# Generate trial IPD, outcomes, matched external samples, and complete AgD.

source(file.path(CODE_DIR, "00_config.R"), local = FALSE)
source(file.path(CODE_DIR, "01_external_data.R"), local = FALSE)
source(file.path(CODE_DIR, "02_design.R"), local = FALSE)

te_binary_from_counts <- function(r_t, n_t, r_a, n_a, cc = 0.5, eps = 1e-6) {
  if (any(!is.finite(c(r_t, n_t, r_a, n_a))) || n_t <= 0 || n_a <= 0) return(c(TE = NA_real_, se = NA_real_))
  p_t <- (r_t + cc) / (n_t + 2 * cc)
  p_a <- (r_a + cc) / (n_a + 2 * cc)
  p_t <- pmin(pmax(p_t, eps), 1 - eps); p_a <- pmin(pmax(p_a, eps), 1 - eps)
  z_t <- stats::qnorm(p_t); z_a <- stats::qnorm(p_a)
  v_t <- p_t * (1 - p_t) / (n_t * stats::dnorm(z_t)^2)
  v_a <- p_a * (1 - p_a) / (n_a * stats::dnorm(z_a)^2)
  c(TE = z_t - z_a, se = sqrt(v_t + v_a))
}

te_continuous_from_values <- function(y_t, y_a) {
  n_t <- length(y_t); n_a <- length(y_a)
  if (n_t < 2L || n_a < 2L) return(c(TE = NA_real_, se = NA_real_))
  c(TE = mean(y_t) - mean(y_a), se = sqrt(stats::var(y_t) / n_t + stats::var(y_a) / n_a))
}

generate_outcomes <- function(dat, beta) {
  z <- as.matrix(dat[, ALL_VARS, drop = FALSE]) - CENTER_VALUE
  beta_by_row <- t(vapply(dat$Trt, function(t) beta[[t]], numeric(length(ALL_VARS))))
  interaction_term <- rowSums(z * beta_by_row)
  treatment_term <- unname(TAU[dat$Trt]) + interaction_term
  eta_binary <- ALPHA_BINARY + treatment_term
  eps <- stats::rnorm(nrow(dat), 0, 1)
  dat$latent_error <- eps
  dat$eta_binary <- eta_binary
  dat$y_binary <- as.integer(eta_binary + eps > 0)
  dat$y_continuous <- treatment_term + CONTINUOUS_SD * eps
  dat
}

effect_row <- function(dat, outcome, type, known_var = NULL, known_value = NULL) {
  trt1 <- unique(dat$Trt1); trt2 <- unique(dat$Trt2)
  if (length(trt1) != 1L || length(trt2) != 1L) stop("Invalid treatment labels in trial data")
  a <- dat[dat$Trt == trt1, , drop = FALSE]
  t <- dat[dat$Trt == trt2, , drop = FALSE]
  n_a <- nrow(a); n_t <- nrow(t)
  if (outcome == "binary") {
    r_a <- sum(a$y_binary); r_t <- sum(t$y_binary)
    est <- te_binary_from_counts(r_t, n_t, r_a, n_a)
    mean_a <- if (n_a) r_a / n_a else NA_real_; mean_t <- if (n_t) r_t / n_t else NA_real_
    sd_a <- sd_t <- NA_real_
  } else {
    est <- te_continuous_from_values(t$y_continuous, a$y_continuous)
    mean_a <- mean(a$y_continuous); mean_t <- mean(t$y_continuous)
    sd_a <- stats::sd(a$y_continuous); sd_t <- stats::sd(t$y_continuous)
    r_a <- r_t <- NA_real_
  }
  x <- stats::setNames(rep(NA_real_, length(ALL_VARS)), ALL_VARS)
  if (identical(type, "overall")) x <- vapply(ALL_VARS, function(v) mean(dat[[v]]), numeric(1))
  if (!is.null(known_var)) x[[known_var]] <- as.numeric(known_value)
  tibble::tibble(
    outcome = outcome, Study = unique(dat$Trial), Edge = unique(dat$Edge),
    Trt1 = trt1, Trt2 = trt2, type = type,
    TE = unname(est["TE"]), se = unname(est["se"]),
    n1 = n_a, n2 = n_t, n_eff = n_a + n_t,
    event1 = r_a, event2 = r_t, mean1 = mean_a, mean2 = mean_t, sd1 = sd_a, sd2 = sd_t,
    !!!as.list(x)
  )
}

build_complete_agd <- function(dat, outcome = c("binary", "continuous")) {
  outcome <- match.arg(outcome)
  rows <- list(effect_row(dat, outcome, "overall"))
  for (v in ALL_VARS) {
    for (value in c(1L, 0L)) {
      sub <- dat[dat[[v]] == value, , drop = FALSE]
      rows[[length(rows) + 1L]] <- effect_row(sub, outcome, paste0(v, "=", value), v, value)
    }
  }
  out <- dplyr::bind_rows(rows)
  out$.agd_row_id <- paste(out$Study, out$outcome, out$type, sep = "|")
  out
}

simulate_one_trial <- function(trial_row, target_margins, beta, cell_table, rep_id, setting) {
  trial_id <- as.character(trial_row$Trial)
  trial_seed <- seed_for(rep_id, paste0("trial_ipd_", setting, "_", trial_id))
  rr <- raked_binary_sample(cell_table, target_margins, n = TRIAL_N, seed = trial_seed)
  dat <- rr$sample
  set.seed(seed_for(rep_id, paste0("randomize_outcome_", setting, "_", trial_id)))
  dat <- dat[sample.int(nrow(dat), nrow(dat), replace = FALSE), , drop = FALSE]
  trt1 <- as.character(trial_row$Trt1); trt2 <- as.character(trial_row$Trt2)
  dat$Trt <- c(rep(trt1, ARM_N), rep(trt2, ARM_N))
  dat$Trial <- trial_id; dat$Edge <- as.character(trial_row$Edge); dat$Trt1 <- trt1; dat$Trt2 <- trt2
  dat <- generate_outcomes(dat, beta)

  observed <- vapply(ALL_VARS, function(v) mean(dat[[v]]), numeric(1))
  ext <- raked_binary_sample(
    cell_table, observed, n = N_MATCH,
    seed = seed_for(rep_id, paste0("external_match_", setting, "_", trial_id))
  )
  ext_moments <- external_moments(ext$sample)

  diag <- tibble::tibble(
    setting = setting, Trial = trial_id, Edge = as.character(trial_row$Edge), N = TRIAL_N,
    ipf_converged = rr$ipf_converged, ipf_iterations = rr$ipf_iterations,
    ipf_max_abs_error = rr$ipf_max_abs_error,
    max_sample_margin_deviation = max(abs(rr$observed - rr$target)),
    external_max_margin_deviation = max(abs(ext$observed - observed)),
    binary_events = sum(dat$y_binary), binary_risk = mean(dat$y_binary)
  )
  list(
    ipd = dat,
    agd_binary = build_complete_agd(dat, "binary"),
    agd_continuous = build_complete_agd(dat, "continuous"),
    external = ext_moments,
    target_margins = stats::setNames(as.numeric(target_margins[ALL_VARS]), ALL_VARS),
    observed_margins = observed,
    diagnostics = diag
  )
}

simulate_environment <- function(rep_id, setting, beta, cell_table) {
  trial_targets <- draw_trial_margins(rep_id, setting)
  trials <- make_trial_table(setting)
  out <- vector("list", nrow(trials)); names(out) <- trials$Trial
  for (i in seq_len(nrow(trials))) {
    tr <- trials[i, , drop = FALSE]
    margins <- stats::setNames(as.numeric(trial_targets[i, ALL_VARS, drop = TRUE]), ALL_VARS)
    out[[i]] <- simulate_one_trial(tr, margins, beta, cell_table, rep_id, setting)
  }
  cov_means <- dplyr::bind_rows(lapply(seq_along(out), function(i) {
    tibble::tibble(Trial = trials$Trial[i], Edge = trials$Edge[i], N = TRIAL_N,
                   !!!as.list(out[[i]]$observed_margins))
  }))
  list(
    setting = setting,
    trial_table = trials,
    trial_targets = trial_targets,
    trial_covariate_means = cov_means,
    targets = make_targets_from_trials(cov_means),
    agd = list(
      binary = dplyr::bind_rows(lapply(out, `[[`, "agd_binary")),
      continuous = dplyr::bind_rows(lapply(out, `[[`, "agd_continuous"))
    ),
    external = lapply(out, `[[`, "external"),
    ipd = if (SAVE_IPD) lapply(out, `[[`, "ipd") else NULL,
    diagnostics = dplyr::bind_rows(lapply(out, `[[`, "diagnostics"))
  )
}
