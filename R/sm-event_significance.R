#' Internal: per-position significance for splicing / sequence maps
#'
#' Compares a Positive / Negative group score matrix against the Control
#' score matrix at each plot position. One-sided test for enrichment
#'
#' @keywords internal
#' @name event_significance


#' @rdname event_significance
#' @noRd
test_per_position <- function(group_score_matrix,
                              control_score_matrix,
                              opts) {

  group_hits   <- colSums(group_score_matrix)
  control_hits <- colSums(control_score_matrix)
  n_group      <- nrow(group_score_matrix)
  n_control    <- nrow(control_score_matrix)

  #Pick Significance Test
  pvalues <- switch(
    opts$stat_test,
    fisher   = .fisher_per_position(group_hits, n_group, control_hits, n_control),
    binomial = .binomial_per_position(group_hits, n_group, control_hits, n_control),
    abort_invalid_arg(
      "Unknown {.arg stat_test}: {.val {opts$stat_test}}."
    )
  )

  #FDR Correction
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


# Fisher's exact test, one-sided.
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

# One-sided binomial test.
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
