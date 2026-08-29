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
#' @param control_stats Result of `bootstrap_control()` (per-position Control
#'   mean/SD plus the bootstrap `sample_size`). Only used by
#'   `"fisher-bootstrap"`. May be `NULL`.
#' @param opts Result of [splicing_options()].
#'
#' @keywords internal
#' @name event_significance


#' @rdname event_significance
#' @noRd
test_per_position <- function(group_counts, n_group,
                              control_counts, n_control,
                              control_stats = NULL, opts) {
  #Both tests are the same hypergeometric; they differ only in the Control
  #numbers fed to it. fisher-bootstrap substitutes the resampled Control mean
  #(rounded to a whole count) and the bootstrap sample size.
  if (identical(opts$stat_test, "fisher-bootstrap")) {
    boot <- .bootstrap_control_cell(control_stats, control_counts, n_control)
    control_counts <- boot$control_counts
    n_control      <- boot$n_control
  }

  pvalues <- .fisher_per_position(group_counts, n_group,
                                  control_counts, n_control)

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

#Reconstruct the Control cell for fisher-bootstrap. The bootstrap stores a
#per-position hit *frequency*; Fisher needs whole counts, so multiply by the
#bootstrap sample size and round to nearest. Falls back to the raw Control
#pool when the bootstrap is degenerate (sample_size < 1).
.bootstrap_control_cell <- function(control_stats, control_counts, n_control) {
  sample_size <- if (!is.null(control_stats)) control_stats$sample_size else 0L
  if (is.null(control_stats) || length(sample_size) != 1L || sample_size < 1L) {
    cli::cli_warn(c(
      "{.val fisher-bootstrap}: no usable bootstrap sample; using raw Control counts.",
      "i" = "Check {.arg control_multiplier} / Control pool size in {.fn splicing_options}."
    ))
    return(list(control_counts = control_counts, n_control = n_control))
  }
  list(
    control_counts = round(control_stats$mean_per_position * sample_size),
    n_control      = as.integer(sample_size)
  )
}
