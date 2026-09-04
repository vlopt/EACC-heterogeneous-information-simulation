# Standard synthesis and RSM-style Fullboot-EACC with trial-specific covariate sets.

source(file.path(CODE_DIR, "00_config.R"), local = FALSE)
source(file.path(CODE_DIR, "02_design.R"), local = FALSE)

subset_agd_for_vars <- function(agd_trial, vars) {
  wanted <- c("overall", unlist(lapply(vars, function(v) paste0(v, "=", c(1, 0)))))
  agd_trial %>% dplyr::filter(type %in% wanted) %>% dplyr::arrange(match(type, wanted))
}

conditional_covariance <- function(mu, Sigma, x, vars) {
  mu <- mu[vars]; S <- Sigma[vars, vars, drop = FALSE]
  obs <- which(is.finite(x)); mis <- which(!is.finite(x))
  if (!length(mis)) return(list(mean = numeric(0), C = matrix(numeric(0), 0, 0), missing = character(0)))
  missing <- vars[mis]
  if (!length(obs)) {
    C <- regularize_psd(S[mis, mis, drop = FALSE])
    rownames(C) <- colnames(C) <- missing
    return(list(mean = mu[mis], C = C, missing = missing))
  }
  Soo <- S[obs, obs, drop = FALSE]; Smo <- S[mis, obs, drop = FALSE]
  Smm <- S[mis, mis, drop = FALSE]
  diag(Soo) <- pmax(diag(Soo), 1e-10)
  cm <- as.numeric(mu[mis] + Smo %*% safe_solve(Soo, x[obs] - mu[obs])); names(cm) <- missing
  C <- Smm - Smo %*% safe_solve(Soo, t(Smo))
  rownames(C) <- colnames(C) <- missing
  list(mean = cm, C = regularize_psd(C, 1e-12), missing = missing)
}

blup_impute_selected <- function(agd, mu, Sigma, vars, random = TRUE, seed = NULL) {
  if (!length(vars)) return(list(data = agd, info = tibble::tibble()))
  if (!is.null(seed)) set.seed(as.integer(seed))
  out <- agd
  info <- vector("list", nrow(out))
  for (i in seq_len(nrow(out))) {
    x <- vapply(vars, function(v) suppressWarnings(as.numeric(out[[v]][i])), numeric(1)); names(x) <- vars
    cond <- conditional_covariance(mu, Sigma, x, vars)
    if (length(cond$missing)) {
      values <- cond$mean
      if (isTRUE(random)) {
        ne <- suppressWarnings(as.numeric(out$n_eff[i]))
        if (!is.finite(ne) || ne <= 0) ne <- BLUP_N_EFF_MIN
        Cdraw <- cond$C * BLUP_RANDOM_SCALE^2 * BLUP_DESIGN_EFFECT / max(ne, BLUP_N_EFF_MIN)
        values <- tryCatch(as.numeric(MASS::mvrnorm(1, cond$mean, regularize_psd(Cdraw, 1e-12))), error = function(e) cond$mean)
        names(values) <- cond$missing
      }
      for (v in cond$missing) out[[v]][i] <- pmin(pmax(values[[v]], 0), 1)
    }
    info[[i]] <- tibble::tibble(
      .agd_row_id = as.character(out$.agd_row_id[i]),
      missing_vars = list(cond$missing), C = list(cond$C)
    )
  }
  list(data = out, info = dplyr::bind_rows(info))
}

eacc_imp_lambda <- function(k_imp) {
  k_imp <- max(1L, as.integer(k_imp))
  if (identical(EACC_IMP_VAR_SCALE_MODE, "by_imputed_vars")) 1 / k_imp else 1
}

eacc_version_a_variance <- function(dat, X, y, w, coefficients, imp_info) {
  p <- length(coefficients); coef_names <- names(coefficients)
  zero <- matrix(0, p, p, dimnames = list(coef_names, coef_names))
  if (is.null(imp_info) || !nrow(imp_info)) return(list(raw = zero, scaled = zero, k_imp = 0L, lambda = 1))
  A_inv <- safe_solve(t(X) %*% (X * as.numeric(w)))
  rownames(A_inv) <- colnames(A_inv) <- coef_names
  resid <- as.numeric(y - X %*% coefficients)
  Vraw <- zero; used <- character(0)
  row_match <- match(as.character(dat$.agd_row_id), as.character(imp_info$.agd_row_id))
  for (i in seq_len(nrow(dat))) {
    ii <- row_match[i]
    if (is.na(ii)) next
    missing <- imp_info$missing_vars[[ii]]
    model_missing <- intersect(missing, setdiff(coef_names, "(Intercept)"))
    if (!length(model_missing)) next
    C <- imp_info$C[[ii]]
    if (is.null(C) || !length(C)) next
    C <- C[model_missing, model_missing, drop = FALSE]
    k <- length(model_missing)
    Sr <- matrix(0, p, k, dimnames = list(coef_names, model_missing))
    for (v in model_missing) Sr[v, v] <- 1
    beta_m <- as.numeric(coefficients[model_missing])
    xr <- as.numeric(X[i, coef_names, drop = TRUE])
    Br <- resid[i] * Sr - tcrossprod(xr, beta_m)
    Gr <- A_inv %*% (as.numeric(w[i]) * Br)
    add <- Gr %*% regularize_psd(C, 1e-12) %*% t(Gr)
    add[!is.finite(add)] <- 0
    Vraw <- Vraw + (add + t(add)) / 2
    used <- unique(c(used, model_missing))
  }
  k_imp <- length(used); lambda <- eacc_imp_lambda(k_imp)
  list(raw = Vraw, scaled = lambda * Vraw, k_imp = k_imp, lambda = lambda)
}

wls_trial_predict <- function(agd_imp, vars, targets, imp_info = NULL) {
  dat <- agd_imp %>%
    dplyr::mutate(dplyr::across(c(TE, se, dplyr::all_of(vars)), ~ suppressWarnings(as.numeric(.)))) %>%
    dplyr::filter(dplyr::if_all(c(TE, se, dplyr::all_of(vars)), is.finite), se > 0)
  p <- 1L + length(vars)
  if (nrow(dat) < p) return(tibble::tibble())
  X <- cbind(`(Intercept)` = 1, as.matrix(dat[, vars, drop = FALSE]))
  y <- dat$TE
  w0 <- 1 / dat$se^2
  w <- pmin(w0, as.numeric(stats::quantile(w0, 0.99, na.rm = TRUE)))
  XtWX <- t(X) %*% (X * w)
  coef <- as.numeric(safe_solve(XtWX, t(X) %*% (w * y)))
  names(coef) <- colnames(X)
  resid <- as.numeric(y - X %*% coef)
  rank <- qr(X * sqrt(w))$rank
  sigma2 <- sum(w * resid^2) / max(1L, nrow(X) - rank)
  Vreg <- sigma2 * safe_solve(XtWX)
  rownames(Vreg) <- colnames(Vreg) <- names(coef)
  imp <- if (isTRUE(FULLBOOT_EACC_IMP_VAR_SYNC)) eacc_version_a_variance(dat, X, y, w, coef, imp_info) else list(raw = 0 * Vreg, scaled = 0 * Vreg, k_imp = 0L, lambda = 1)
  Vtotal <- regularize_psd(Vreg + imp$scaled, 1e-14)

  purrr::map_dfr(seq_len(nrow(targets)), function(i) {
    xv <- c(`(Intercept)` = 1, stats::setNames(as.numeric(targets[i, vars, drop = TRUE]), vars))
    est <- sum(xv[names(coef)] * coef)
    vreg <- as.numeric(t(xv) %*% Vreg %*% xv)
    vimp <- as.numeric(t(xv) %*% imp$scaled %*% xv)
    vtot <- as.numeric(t(xv) %*% Vtotal %*% xv)
    tibble::tibble(
      target_id = targets$target_id[i], estimate = est,
      within_var = pmax(vtot, 0), within_se = sqrt(pmax(vtot, 0)),
      var_reg = pmax(vreg, 0), var_imp = pmax(vimp, 0),
      imp_var_k = imp$k_imp, imp_var_lambda = imp$lambda
    )
  })
}

fe_pairwise <- function(te, se) {
  ok <- is.finite(te) & is.finite(se) & se > 0
  if (!any(ok)) return(c(estimate = NA_real_, se = NA_real_))
  w <- 1 / se[ok]^2
  c(estimate = sum(w * te[ok]) / sum(w), se = sqrt(1 / sum(w)))
}

nma_fe <- function(y, se, t1, t2, ref = "A") {
  ok <- is.finite(y) & is.finite(se) & se > 0
  y <- y[ok]; se <- se[ok]; t1 <- as.character(t1[ok]); t2 <- as.character(t2[ok])
  trts <- sort(unique(c(t1, t2))); nonref <- setdiff(trts, ref)
  if (!length(y) || !length(nonref)) return(NULL)
  X <- matrix(0, length(y), length(nonref), dimnames = list(NULL, nonref))
  for (i in seq_along(y)) {
    if (t2[i] != ref) X[i, t2[i]] <- 1
    if (t1[i] != ref) X[i, t1[i]] <- X[i, t1[i]] - 1
  }
  XtWX <- t(X) %*% (X / se^2); XtWy <- t(X) %*% (y / se^2)
  d <- as.numeric(safe_solve(XtWX, XtWy)); names(d) <- nonref
  V <- safe_solve(XtWX); rownames(V) <- colnames(V) <- nonref
  list(d = d, V = V, ref = ref)
}

nma_bc <- function(fit) {
  if (is.null(fit) || !all(c("B", "C") %in% names(fit$d))) return(c(estimate = NA_real_, se = NA_real_))
  est <- fit$d["B"] - fit$d["C"]
  v <- fit$V["B", "B"] + fit$V["C", "C"] - 2 * fit$V["B", "C"]
  c(estimate = as.numeric(est), se = sqrt(pmax(0, v)))
}

synthesize_trial_predictions <- function(preds, setting) {
  purrr::map_dfr(sort(unique(preds$target_id)), function(tid) {
    x <- preds[preds$target_id == tid, , drop = FALSE]
    pooled <- if (setting == "pairwise") fe_pairwise(x$estimate, sqrt(x$within_var)) else nma_bc(nma_fe(x$estimate, sqrt(x$within_var), x$Trt1, x$Trt2))
    tibble::tibble(target_id = tid, estimate = unname(pooled["estimate"]), within_var = unname(pooled["se"])^2)
  })
}

fullboot_one_draw <- function(b, agd, shared_external_draws, targets, varsets_by_trial, setting, rep_id) {
  trial_rows <- list()
  trials <- names(varsets_by_trial)
  external_draw <- shared_external_draws[[b]]
  if (is.null(external_draw) || is.null(external_draw$trial_moments)) stop("Missing shared Fullboot external draw: ", b)
  for (trial in trials) {
    vars <- varsets_by_trial[[trial]]
    if (!length(vars)) next
    raw <- agd[agd$Study == trial, , drop = FALSE]
    raw <- subset_agd_for_vars(raw, vars)
    # The same shared external source was bootstrapped once for this draw, then
    # independently re-raked to every trial's observed margins.  These moments
    # contain no mean/covariance shrinkage toward the original match.
    mom <- external_draw$trial_moments[[trial]]
    if (is.null(mom)) stop("Missing re-matched external moments for trial: ", trial)
    imp <- blup_impute_selected(
      raw, mom$mu, mom$Sigma, vars, random = TRUE,
      seed = seed_for(rep_id, paste0("fullboot_blup_", setting, "_", trial, "_", paste(vars, collapse = "-")), b)
    )
    pr <- wls_trial_predict(imp$data, vars, targets, imp$info)
    if (!nrow(pr)) next
    map <- raw[1, c("Study", "Edge", "Trt1", "Trt2"), drop = FALSE]
    trial_rows[[length(trial_rows) + 1L]] <- pr %>%
      dplyr::mutate(Study = trial, Edge = map$Edge, Trt1 = map$Trt1, Trt2 = map$Trt2)
  }
  preds <- dplyr::bind_rows(trial_rows)
  if (!nrow(preds)) return(tibble::tibble())
  synthesize_trial_predictions(preds, setting) %>% dplyr::mutate(draw = b, .before = 1)
}

combine_fullboot <- function(draws, method, architecture, outcome, setting) {
  if (!nrow(draws)) return(tibble::tibble())
  draws %>%
    dplyr::filter(is.finite(estimate), is.finite(within_var)) %>%
    dplyr::group_by(target_id) %>%
    dplyr::summarise(
      Bvar = if (dplyr::n() >= 2L) stats::var(.data$estimate) else 0,
      W = mean(.data$within_var), n_boot_ok = dplyr::n(),
      estimate = mean(.data$estimate), .groups = "drop"
    ) %>%
    dplyr::mutate(
      Tvar = W + (1 + 1 / pmax(1, n_boot_ok)) * Bvar,
      se = sqrt(pmax(Tvar, 0)), lwr = estimate - 1.96 * se, upr = estimate + 1.96 * se,
      method = method, architecture = architecture, outcome = outcome, setting = setting,
      .before = 1
    )
}

run_fullboot_eacc <- function(agd, shared_external_draws, targets, varsets_by_trial, method, architecture, outcome, setting, rep_id, b_boot = B_BOOT) {
  if (length(shared_external_draws) < b_boot) stop("Insufficient shared external Fullboot draws")
  draws <- if (BOOT_CORES > 1L && requireNamespace("future.apply", quietly = TRUE)) {
    future.apply::future_lapply(seq_len(b_boot), function(b) fullboot_one_draw(b, agd, shared_external_draws, targets, varsets_by_trial, setting, rep_id), future.seed = TRUE)
  } else {
    lapply(seq_len(b_boot), function(b) fullboot_one_draw(b, agd, shared_external_draws, targets, varsets_by_trial, setting, rep_id))
  }
  draw_tbl <- dplyr::bind_rows(draws)
  ans <- combine_fullboot(draw_tbl, method, architecture, outcome, setting)
  list(summary = ans, draws = if (SAVE_BOOT_DRAWS) draw_tbl else NULL)
}

run_standard_synthesis <- function(agd, targets, architecture, outcome, setting) {
  overall <- agd %>% dplyr::filter(type == "overall")
  pooled <- if (setting == "pairwise") fe_pairwise(overall$TE, overall$se) else nma_bc(nma_fe(overall$TE, overall$se, overall$Trt1, overall$Trt2))
  method_label <- if (setting == "pairwise") "Meta" else "NMA"
  tibble::tibble(
    method = method_label, architecture = architecture, outcome = outcome, setting = setting,
    target_id = targets$target_id, estimate = unname(pooled["estimate"]),
    se = unname(pooled["se"]), lwr = estimate - 1.96 * se, upr = estimate + 1.96 * se,
    W = se^2, Bvar = 0, Tvar = se^2, n_boot_ok = NA_integer_
  )
}
