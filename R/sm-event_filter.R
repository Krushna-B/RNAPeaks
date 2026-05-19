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

  # Single-distribution event types (e.g. UTR).
  if (identical(schema$group_set, "Single")) {
    events$group <- "Single"
    return(list(Single = events))
  }

  # Neg/Pos require BOTH PValue and FDR below `event_fdr`.
  # Control requires BOTH PValue and FDR above `control_pval`.
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
    Negative = events[sig_mask     & dp <= psi_cutoff[1],        , drop = FALSE],
    Positive = events[sig_mask     & dp >= psi_cutoff[2],        , drop = FALSE],
    Control  = events[control_mask & abs(dp) <  psi_control_max, , drop = FALSE]
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
