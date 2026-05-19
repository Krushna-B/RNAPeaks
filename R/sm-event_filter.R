#' Internal: split events table into analysis groups
#'
#' @param events Data frame of events.
#' @param schema One of the `event_schema_*` lists from `R/event_schema.R`.
#' @param opts Result of [splicing_options()].
#' @return Named list of data frames, one per requested group.
#' @keywords internal
#' @noRd
filter_events <- function(events, schema, opts) {

  check_data_frame(events, "events", required = schema$required_cols)
  events$chr <- normalize_chr(events$chr)
  events     <- .apply_count_filter(events, opts)

  # Single-distribution event types (e.g. UTR).
  if (identical(schema$group_set, "Single")) {
    events$group <- "Single"
    return(list(Single = events))
  }

  # Neg/Pos require both PValue and FDR filter
  # Control requires both PValue and FDR filter
  event_fdr       <- opts$event_fdr
  control_pval    <- opts$control_pval
  psi_cutoff      <- opts$psi_cutoff
  psi_control_max <- opts$psi_control_max

  fd <- events$FDR
  pv <- events$PValue
  dp <- events$IncLevelDifference

  sig_mask     <- fd < event_fdr    & pv < event_fdr
  control_mask <- pv > control_pval & fd > control_pval

  #filter
  groups_all <- list(
    Negative = events[sig_mask     & dp < psi_cutoff[1],        , drop = FALSE],
    Positive = events[sig_mask     & dp > psi_cutoff[2],        , drop = FALSE],
    Control  = events[control_mask & abs(dp) < psi_control_max, , drop = FALSE]
  )

  requested <- intersect(c("Negative", "Positive", "Control"), opts$groups)
  out       <- groups_all[requested]

  #Check groups aren't empty
  for (g in requested) {
    if (nrow(out[[g]]) == 0L) {
      cli::cli_warn(c(
        "{.val {g}} group has no events after filtering.",
        "i" = "Check {.arg event_fdr} / {.arg control_pval} / {.arg psi_cutoff} / {.arg psi_control_max}."
      ))
    } else {
      out[[g]]$group <- g
    }
  }

  out
}


# Drop events with low junction-read coverage.
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

  parse_count <- function(x) {
    vapply(strsplit(as.character(x), ",", fixed = TRUE), function(v) {
      sum(suppressWarnings(as.numeric(v)), na.rm = TRUE)
    }, numeric(1L))
  }
  in_1 <- parse_count(events[["IJC_SAMPLE_1"]])
  sk_1 <- parse_count(events[["SJC_SAMPLE_1"]])
  in_2 <- parse_count(events[["IJC_SAMPLE_2"]])
  sk_2 <- parse_count(events[["SJC_SAMPLE_2"]])
  keep <- (in_1 + sk_1) > opts$min_count &
          (in_2 + sk_2) > opts$min_count

  if (isTRUE(opts$verbose)) {
    cli::cli_inform(
      "min_count = {opts$min_count}: kept {sum(keep)}/{length(keep)} event{?s}."
    )
  }
  events[keep, , drop = FALSE]
}
