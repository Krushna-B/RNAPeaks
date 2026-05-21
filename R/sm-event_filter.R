#' Internal: drop low-coverage events and split into analysis groups
#'
#' First applies the junction-count coverage filter (skipped when
#' `opts$min_count == 0`), then partitions the surviving rows into
#' Negative / Positive / Control by FDR / PValue / IncLevelDifference.
#'
#' @param events Data frame of events, already chr-normalized.
#' @param schema One of the `event_schema_*` lists.
#' @param opts Result of [splicing_options()].
#' @return List with `events` (the count-filtered data frame) and
#'   `groups_idx` (named list of row-index vectors, one per requested group,
#'   indexing into `events`).
#' @keywords internal
#' @noRd
filter_events <- function(events, schema, opts) {
  #Apply count filtering (if count filtering is on)
  events <- .apply_count_filter(events, opts)

  #If only 1 group then return
  if (identical(schema$group_set, "Single")) {
    return(list(events     = events,
                groups_idx = list(Single = seq_len(nrow(events)))))
  }

  #Filter out and create Positive, Negative, and Controls Groups
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

  requested  <- intersect(c("Negative", "Positive", "Control"), opts$groups)
  groups_idx <- groups_all[requested]

  for (g in requested) {
    if (length(groups_idx[[g]]) == 0L) {
      cli::cli_warn(c(
        "{.val {g}} group has no events after filtering.",
        "i" = "Check {.arg event_fdr} / {.arg control_pval} / {.arg psi_cutoff} / {.arg psi_control_max}."
      ))
    }
  }

  list(events = events, groups_idx = groups_idx)
}


#Filtering based on min count
.apply_count_filter <- function(events, opts) {
  #If min count > 0 then it is on
  if (!isTRUE(opts$min_count > 0L)) return(events)

  #Designated Count cols, show error if not present
  count_cols <- c("IJC_SAMPLE_1", "SJC_SAMPLE_1",
                  "IJC_SAMPLE_2", "SJC_SAMPLE_2")
  missing_cc <- setdiff(count_cols, colnames(events))
  if (length(missing_cc) > 0L) {
    abort_invalid_arg(c(
      "{.arg events} is missing junction-count column{?s}: {.field {missing_cc}}.",
      "i" = "Set {.arg min_count = 0} in {.fn splicing_options} to skip count filtering."
    ))
  }

  #Helper function that does a summation of of all comma separated values in each column
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

  #Min Count filtering
  inclusion_sample_1 <- sum_replicate_counts(events[["IJC_SAMPLE_1"]])
  skipping_sample_1  <- sum_replicate_counts(events[["SJC_SAMPLE_1"]])
  inclusion_sample_2 <- sum_replicate_counts(events[["IJC_SAMPLE_2"]])
  skipping_sample_2  <- sum_replicate_counts(events[["SJC_SAMPLE_2"]])

  keep <- (inclusion_sample_1 + skipping_sample_1) > opts$min_count &
          (inclusion_sample_2 + skipping_sample_2) > opts$min_count

  #Information in Verbose mode
  if (isTRUE(opts$verbose)) {
    cli::cli_inform(
      "min_count = {opts$min_count}: kept {sum(keep)}/{length(keep)} event{?s}."
    )
  }
  events[keep, , drop = FALSE]
}
