# Run a resumable consecutive replication chunk.

CODE_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/"))
source(file.path(CODE_DIR, "00_config.R"), local = FALSE)
source(file.path(CODE_DIR, "05_run_replication.R"), local = FALSE)

parse_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), commandArgs(TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

rep_start <- as.integer(parse_arg("rep_start", 1L))
rep_end <- as.integer(parse_arg("rep_end", MC_REPS))
chunk_label <- parse_arg("chunk_label", sprintf("rep_%04d_%04d", rep_start, rep_end))
if (!is.finite(rep_start) || !is.finite(rep_end) || rep_start < 1L || rep_end < rep_start) stop("Invalid replication range")

if (BOOT_CORES > 1L && requireNamespace("future", quietly = TRUE)) {
  future::plan(future::multisession, workers = BOOT_CORES)
  on.exit(future::plan(future::sequential), add = TRUE)
}

pool <- load_external_pool()
cell_table <- make_cell_table(pool)
chunk_log <- file.path(LOG_DIR, paste0(RUN_TAG, "_", chunk_label, "_status.csv"))
status <- list()

for (rep_id in seq.int(rep_start, rep_end)) {
  rep_dir <- file.path(RUN_ROOT, sprintf("REP_%04d", rep_id))
  result_file <- file.path(rep_dir, "performance_rows.csv")
  if (file.exists(result_file)) {
    message(sprintf("[%s] REP %d already complete; skipped", chunk_label, rep_id))
    status[[length(status) + 1L]] <- tibble::tibble(rep = rep_id, status = "skipped_complete", message = "", timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    next
  }
  message(sprintf("[%s] Starting REP %d", chunk_label, rep_id))
  started <- Sys.time()
  ans <- tryCatch({
    run_one_replication(rep_id, pool, cell_table, b_boot = B_BOOT, write_output = TRUE)
    list(ok = TRUE, message = "")
  }, error = function(e) list(ok = FALSE, message = conditionMessage(e)))
  status[[length(status) + 1L]] <- tibble::tibble(
    rep = rep_id, status = if (ans$ok) "complete" else "failed", message = ans$message,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  readr::write_csv(dplyr::bind_rows(status), chunk_log)
  if (!ans$ok) message(sprintf("[%s] REP %d FAILED: %s", chunk_label, rep_id, ans$message))
}

message(sprintf("[%s] Chunk finished", chunk_label))
