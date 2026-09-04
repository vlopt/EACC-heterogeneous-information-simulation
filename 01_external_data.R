# External binary-covariate pool and efficient IPF/raking helpers.

source(file.path(CODE_DIR, "00_config.R"), local = FALSE)

first_matching_name <- function(nms, patterns) {
  low <- tolower(nms)
  for (p in patterns) {
    hit <- grep(p, low, perl = TRUE)
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

as_binary01 <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) return(as.integer(x))
  if (is.character(x)) return(as.integer(tolower(trimws(x)) %in% c("1", "yes", "y", "true", "male", "m", "current", "ever", "男")))
  x <- suppressWarnings(as.numeric(x))
  if (all(stats::na.omit(x) %in% c(0, 1))) return(as.integer(x))
  if (all(stats::na.omit(x) %in% c(1, 2))) return(as.integer(x == 1))
  as.integer(x >= 0.5)
}

standardize_external_covariates <- function(df) {
  nms <- names(df)
  age65_name <- first_matching_name(nms, c("age65plus", "age.?65", "65plus"))
  age_name <- first_matching_name(nms, c("^age_final$", "^age$", "age_year"))
  htn_name <- first_matching_name(nms, c("hypertension", "^htn$"))
  male_name <- first_matching_name(nms, c("^male$", "sex", "gender"))
  smoke_name <- first_matching_name(nms, c("smoking", "smok"))
  if (is.na(htn_name) || is.na(male_name) || is.na(smoke_name) || (is.na(age65_name) && is.na(age_name))) {
    stop("External database does not contain the required four covariates")
  }
  age65 <- if (!is.na(age65_name)) as_binary01(df[[age65_name]]) else as.integer(as.numeric(df[[age_name]]) >= 65)
  out <- tibble::tibble(
    age65 = age65,
    htn = as_binary01(df[[htn_name]]),
    male = as_binary01(df[[male_name]]),
    smoke = as_binary01(df[[smoke_name]])
  )
  out <- stats::na.omit(out)
  for (v in ALL_VARS) out[[v]] <- as.integer(out[[v]] == 1L)
  out
}

load_external_pool <- function(path = RWD_FILE) {
  if (!file.exists(path)) stop("External database not found: ", path)
  raw <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  x <- standardize_external_covariates(raw)
  if (nrow(x) < 1000L) stop("External pool is unexpectedly small after standardization")
  attr(x, "source_path") <- normalizePath(path, winslash = "/", mustWork = TRUE)
  x
}

make_cell_table <- function(pool, pseudo_count = 0.5) {
  grid <- expand.grid(rep(list(0:1), length(ALL_VARS)))
  names(grid) <- ALL_VARS
  key <- do.call(interaction, c(as.list(pool[, ALL_VARS, drop = FALSE]), list(drop = FALSE, lex.order = TRUE)))
  observed <- as.data.frame(table(key), stringsAsFactors = FALSE)
  names(observed) <- c("key", "Freq")
  grid$key <- do.call(interaction, c(as.list(grid[, ALL_VARS, drop = FALSE]), list(drop = FALSE, lex.order = TRUE)))
  counts <- observed$Freq[match(as.character(grid$key), as.character(observed$key))]
  counts[is.na(counts)] <- 0
  grid$source_count <- as.numeric(counts)
  grid$count <- grid$source_count + pseudo_count
  grid$base_prob <- grid$count / sum(grid$count)
  grid$key <- NULL
  out <- tibble::as_tibble(grid)
  attr(out, "source_n") <- sum(out$source_count)
  attr(out, "pseudo_count") <- pseudo_count
  out
}

ipf_cell_prob <- function(cell_table, targets, maxit = IPF_MAXIT, tol = IPF_TOL) {
  targets <- stats::setNames(as.numeric(targets[ALL_VARS]), ALL_VARS)
  if (any(!is.finite(targets)) || any(targets <= 0) || any(targets >= 1)) stop("Invalid IPF target margins")
  p <- as.numeric(cell_table$base_prob)
  converged <- FALSE
  iter <- 0L
  for (it in seq_len(maxit)) {
    iter <- it
    for (v in ALL_VARS) {
      one <- cell_table[[v]] == 1L
      cur <- sum(p[one])
      if (!is.finite(cur) || cur <= 0 || cur >= 1) stop("Degenerate IPF cell margin for ", v)
      p[one] <- p[one] * targets[[v]] / cur
      p[!one] <- p[!one] * (1 - targets[[v]]) / (1 - cur)
      p <- p / sum(p)
    }
    achieved <- vapply(ALL_VARS, function(v) sum(p[cell_table[[v]] == 1L]), numeric(1))
    if (max(abs(achieved - targets)) < tol) {
      converged <- TRUE
      break
    }
  }
  achieved <- vapply(ALL_VARS, function(v) sum(p[cell_table[[v]] == 1L]), numeric(1))
  list(prob = p / sum(p), achieved = achieved, target = targets, converged = converged, iterations = iter,
       max_abs_error = max(abs(achieved - targets)))
}

sample_from_cells <- function(cell_table, prob, n, seed) {
  set.seed(as.integer(seed))
  counts <- as.vector(stats::rmultinom(1, size = as.integer(n), prob = prob))
  idx <- rep.int(seq_len(nrow(cell_table)), counts)
  out <- cell_table[idx, ALL_VARS, drop = FALSE]
  out <- out[sample.int(nrow(out), nrow(out), replace = FALSE), , drop = FALSE]
  rownames(out) <- NULL
  tibble::as_tibble(out)
}

raked_binary_sample <- function(cell_table, targets, n, seed) {
  fit <- ipf_cell_prob(cell_table, targets)
  sample <- sample_from_cells(cell_table, fit$prob, n = n, seed = seed)
  observed <- vapply(ALL_VARS, function(v) mean(sample[[v]]), numeric(1))
  list(
    sample = sample,
    target = stats::setNames(as.numeric(targets[ALL_VARS]), ALL_VARS),
    expected = fit$achieved,
    observed = observed,
    deviations = observed - as.numeric(targets[ALL_VARS]),
    ipf_converged = fit$converged,
    ipf_iterations = fit$iterations,
    ipf_max_abs_error = fit$max_abs_error,
    cell_prob = fit$prob
  )
}

regularize_psd <- function(S, eps = 1e-10) {
  S <- as.matrix(S)
  rn <- rownames(S); cn <- colnames(S)
  S <- (S + t(S)) / 2
  S[!is.finite(S)] <- 0
  ee <- eigen(S, symmetric = TRUE)
  out <- ee$vectors %*% diag(pmax(ee$values, eps), nrow = length(ee$values)) %*% t(ee$vectors)
  out <- (out + t(out)) / 2
  rownames(out) <- rn; colnames(out) <- cn
  out
}

external_moments <- function(x, keep_sample = TRUE) {
  X <- as.matrix(x[, ALL_VARS, drop = FALSE])
  mu <- colMeans(X)
  S <- stats::cov(X)
  rownames(S) <- colnames(S) <- ALL_VARS
  list(mu = mu, Sigma = regularize_psd(S), sample = if (isTRUE(keep_sample)) tibble::as_tibble(X) else NULL)
}

# Exact cell-level equivalent of drawing N_E individuals with replacement from
# the original external source.  No covariance or mean shrinkage is applied.
bootstrap_external_cell_table <- function(cell_table, rep_id, setting, b) {
  source_counts <- as.numeric(cell_table$source_count)
  source_n <- as.integer(sum(source_counts))
  if (!is.finite(source_n) || source_n <= 0L) stop("Invalid external source cell counts")
  source_prob <- source_counts / source_n
  set.seed(seed_for(rep_id, paste0("fullboot_shared_external_", setting), b))
  boot_counts <- as.vector(stats::rmultinom(1, size = source_n, prob = source_prob))
  out <- cell_table
  out$source_count <- boot_counts
  out$count <- boot_counts
  out$base_prob <- boot_counts / sum(boot_counts)
  source_id <- digest::digest(paste(boot_counts, collapse = "|"), algo = "xxhash64", serialize = FALSE)
  list(
    cell_table = out,
    source_id = source_id,
    counts = boot_counts,
    source_n = source_n,
    changed_from_original = !identical(as.numeric(boot_counts), as.numeric(source_counts))
  )
}

# Prepare the external part of Fullboot exactly once per setting and bootstrap
# draw.  Every trial is re-raked to its own observed margins from the same
# shared bootstrap source, and the resulting replicate-specific moments are
# reused by all outcomes and information architectures.
prepare_shared_fullboot_external <- function(cell_table, trial_covariate_means, rep_id, setting, b_boot) {
  stopifnot(all(c("Trial", ALL_VARS) %in% names(trial_covariate_means)))
  draws <- vector("list", b_boot)
  source_rows <- vector("list", b_boot)
  raking_rows <- vector("list", b_boot)
  cell_labels <- apply(as.matrix(cell_table[, ALL_VARS, drop = FALSE]), 1, paste0, collapse = "")
  count_names <- paste0("cell_", cell_labels, "_count")

  for (b in seq_len(b_boot)) {
    shared <- bootstrap_external_cell_table(cell_table, rep_id, setting, b)
    moments_by_trial <- vector("list", nrow(trial_covariate_means))
    names(moments_by_trial) <- trial_covariate_means$Trial
    trial_diag <- vector("list", nrow(trial_covariate_means))

    for (i in seq_len(nrow(trial_covariate_means))) {
      trial <- as.character(trial_covariate_means$Trial[i])
      margins <- stats::setNames(as.numeric(trial_covariate_means[i, ALL_VARS, drop = TRUE]), ALL_VARS)
      rematched <- raked_binary_sample(
        shared$cell_table, margins, n = N_MATCH,
        seed = seed_for(rep_id, paste0("fullboot_rematch_", setting, "_", trial), b)
      )
      moments_by_trial[[trial]] <- external_moments(rematched$sample, keep_sample = FALSE)
      trial_diag[[i]] <- tibble::tibble(
        setting = setting, draw = b, Trial = trial,
        shared_source_id = shared$source_id,
        matched_n = nrow(rematched$sample),
        ipf_converged = rematched$ipf_converged,
        ipf_iterations = rematched$ipf_iterations,
        ipf_max_abs_error = rematched$ipf_max_abs_error,
        matched_max_margin_deviation = max(abs(rematched$deviations))
      )
    }

    draws[[b]] <- list(
      draw = b, shared_source_id = shared$source_id,
      source_counts = shared$counts, trial_moments = moments_by_trial
    )
    source_rows[[b]] <- tibble::tibble(
      setting = setting, draw = b, shared_source_id = shared$source_id,
      source_n = shared$source_n, changed_from_original = shared$changed_from_original,
      !!!stats::setNames(as.list(shared$counts), count_names)
    )
    raking_rows[[b]] <- dplyr::bind_rows(trial_diag)
  }

  list(
    draws = draws,
    source_diagnostics = dplyr::bind_rows(source_rows),
    raking_diagnostics = dplyr::bind_rows(raking_rows)
  )
}
