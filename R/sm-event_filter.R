#' Internal: split events into analysis-group row indices
#'
#' @param events Data frame of events, already chr-normalized and
#'   count-filtered by [prepare_event_map()].
#' @param schema One of the `event_schema_*` lists.
#' @param opts Result of [splicing_options()].
#' @return Named list of integer index vectors, one per requested group.
#' @keywords internal
#' @noRd
filter_events <- function(events, schema, opts) {

  if (identical(schema$group_set, "Single")) {
    return(list(Single = seq_len(nrow(events))))
  }

  event_fdr       <- opts$event_fdr
  control_pval    <- opts$control_pval
  psi_cutoff      <- opts$psi_cutoff
  psi_control_max <- opts$psi_control_max

  fd <- events$FDR
  pv <- events$PValue
  dp <- events$IncLevelDifference

  sig_mask     <- fd < event_fdr    & pv < event_fdr
  control_mask <- pv > control_pval & fd > control_pval

  groups_all <- list(
    Negative = which(sig_mask     & dp < psi_cutoff[1]),
    Positive = which(sig_mask     & dp > psi_cutoff[2]),
    Control  = which(control_mask & abs(dp) < psi_control_max)
  )

  requested <- intersect(c("Negative", "Positive", "Control"), opts$groups)
  out       <- groups_all[requested]

  for (g in requested) {
    if (length(out[[g]]) == 0L) {
      cli::cli_warn(c(
        "{.val {g}} group has no events after filtering.",
        "i" = "Check {.arg event_fdr} / {.arg control_pval} / {.arg psi_cutoff} / {.arg psi_control_max}."
      ))
    }
  }

  out
}


.apply_count_filter <- function(events, opts) {
  if (!isTRUE(opts$min_count > 0L)) return(events)

  count_cols <- c("IJC_SAMPLE_1", "SJC_SAMPLE_1",
                  "IJC_SAMPLE_2", "SJC_SAMPLE_2")
  missing_cc <- setdiff(count_cols, colnames(events))
  if (length(missing_cc) > 0L) {
    abort_invalid_arg(c(
      "{.arg events} is missing junction-count column{?s}: {.field {missing_cc}}.",
      "i" = "Set {.arg min_count = 0} in {.fn splicing_options} to skip count filtering."
    ))
  }

  # Sum comma-separated per-replicate counts (e.g. "12,7,9") to one total per
  # row, vectorised across the whole column in a single rowsum() pass.
  sum_replicate_counts <- function(column) {
    replicates_per_row <- strsplit(as.character(column), ",", fixed = TRUE)
    counts_per_row     <- lengths(replicates_per_row)
    all_replicate_counts <- suppressWarnings(
      as.numeric(unlist(replicates_per_row, use.names = FALSE))
    )
    all_replicate_counts[is.na(all_replicate_counts)] <- 0
    row_of_each_replicate <- rep.int(seq_along(replicates_per_row),
                                     counts_per_row)
    row_totals <- rowsum(all_replicate_counts, row_of_each_replicate,
                          reorder = FALSE)
    totals_per_row <- numeric(length(replicates_per_row))
    totals_per_row[as.integer(rownames(row_totals))] <- row_totals[, 1L]
    totals_per_row
  }

  inclusion_sample_1 <- sum_replicate_counts(events[["IJC_SAMPLE_1"]])
  skipping_sample_1  <- sum_replicate_counts(events[["SJC_SAMPLE_1"]])
  inclusion_sample_2 <- sum_replicate_counts(events[["IJC_SAMPLE_2"]])
  skipping_sample_2  <- sum_replicate_counts(events[["SJC_SAMPLE_2"]])

  keep <- (inclusion_sample_1 + skipping_sample_1) > opts$min_count &
          (inclusion_sample_2 + skipping_sample_2) > opts$min_count

  if (isTRUE(opts$verbose)) {
    cli::cli_inform(
      "min_count = {opts$min_count}: kept {sum(keep)}/{length(keep)} event{?s}."
    )
  }
  events[keep, , drop = FALSE]
}
