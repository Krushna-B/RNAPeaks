#' Internal: per-position significance for splicing / sequence maps
#'
#' Compares Positive / Negative per-position hit counts against Control at
#' each plot position. One-sided test for enrichment.
#'
#' @param group_counts Integer vector of per-position hit counts for the
#'   tested group.
#' @param n_group Integer number of events in the tested group.
#' @param control_counts Integer vector of per-position hit counts for
#'   Control.
#' @param n_control Integer number of Control events.
#' @param opts Result of [splicing_options()].
#'
#' @keywords internal
#' @name event_significance


#' @rdname event_significance
#' @noRd
test_per_position <- function(group_counts, n_group,
                              control_counts, n_control, opts) {
  #Choose which significance test
  pvalues <- switch(
    opts$stat_test,
    fisher   = .fisher_per_position(group_counts, n_group,
                                    control_counts, n_control),
    binomial = .binomial_per_position(group_counts, n_group,
                                      control_counts, n_control),
    abort_invalid_arg(
      "Unknown {.arg stat_test}: {.val {opts$stat_test}}."
    )
  )

  #Adjust P value (FDR Correction)
  pvalue_adj <- if (isTRUE(opts$use_fdr)) {
    stats::p.adjust(pvalues, method = "BH")
  } else {
    pvalues
  }

  data.frame(
    position    = seq_along(pvalues),
    pvalue      = pvalues,
    pvalue_adj  = pvalue_adj,
    significant = pvalue_adj < opts$fdr_threshold,
    stringsAsFactors = FALSE
  )
}

#Fisher Test
.fisher_per_position <- function(group_hits, n_group, control_hits, n_control) {
  total_hits   <- group_hits + control_hits
  total_misses <- (n_group + n_control) - total_hits
  stats::phyper(
    q           = group_hits - 1L,
    m           = total_hits,
    n           = total_misses,
    k           = n_group,
    lower.tail  = FALSE
  )
}

#Binomial Test
.binomial_per_position <- function(group_hits, n_group, control_hits, n_control) {
  control_rate <- control_hits / n_control
  control_rate[!is.finite(control_rate)] <- 0
  stats::pbinom(
    q          = group_hits - 1L,
    size       = n_group,
    prob       = control_rate,
    lower.tail = FALSE
  )
}
