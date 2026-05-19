#' Analysis options for splicing and sequence maps
#'
#' Returns a validated list of options controlling event filtering, region
#' geometry, control bootstrap, moving average, and significance testing for
#' splicing-map / sequence-map entry points.
#'
#' @param width_exon Integer. Bp to extend into exons from each junction.
#' @param width_intron Integer. Bp to extend into introns from each junction.
#' @param event_fdr Maximum value of the rMATS `FDR` column for an event to
#'   be eligible as Positive or Negative.
#' @param control_pval Minimum value of the rMATS `PValue` column for an
#'   event to be eligible as a Control.
#' @param psi_cutoff Length-2 numeric `c(neg, pos)` giving the
#'   `IncLevelDifference` thresholds. Events with `ΔΨ <= neg` are Negative,
#'   `ΔΨ >= pos` are Positive.
#' @param psi_control_max Maximum `abs(IncLevelDifference)` allowed for an
#'   event to enter the Control pool.
#' @param groups Character subset of `c("Negative", "Positive", "Control")`
#'   specifying which event groups to compute.
#' @param control_multiplier Numeric multiplier for control sample size.
#'   Controls drawn per iteration is
#'   `(n_positive + n_negative) * control_multiplier`.
#' @param control_iterations Integer number of bootstrap iterations for the
#'   control distribution.
#' @param moving_average Window size for moving-average smoothing.
#'   `NULL` or `0` disables smoothing.
#' @param stat_test Per-position significance test for Positive / Negative
#'   vs. Control. One of `"fisher"` (Fisher's exact test) or `"binomial"`.
#' @param use_fdr If `TRUE`, apply Benjamini-Hochberg FDR correction to the
#'   per-position p-values before thresholding. If `FALSE`, threshold the raw
#'   p-values.
#' @param fdr_threshold Significance threshold applied to the (FDR-corrected
#'   or raw) p-values.
#' @param verbose If `TRUE`, print progress messages.
#'
#' @return A named list of validated analysis options.
#' @family params
#' @export
splicing_options <- function(
  #Region geometry
  width_exon          = 50,
  width_intron        = 250,

  # Event filtering
  event_fdr           = 0.05,
  control_pval        = 0.95,
  psi_cutoff          = c(-0.1, 0.1),
  psi_control_max     = 0.005,

  # Groups
  groups              = c("Negative", "Positive", "Control"),

  #Control bootstrap
  control_multiplier  = 2.0,
  control_iterations  = 20,

  # Smoothing
  moving_average      = 10,

  # Significance
  stat_test           = "fisher",
  use_fdr             = TRUE,
  fdr_threshold       = 0.05,

  # Diagnostics
  verbose             = TRUE
) {

  # Region geometry
  check_scalar_int(width_exon,   "width_exon",   min = 1)
  check_scalar_int(width_intron, "width_intron", min = 1)

  # Event filtering
  check_unit_interval(event_fdr,    "event_fdr")
  check_unit_interval(control_pval, "control_pval")

  if (!is.numeric(psi_cutoff) || length(psi_cutoff) != 2L || anyNA(psi_cutoff)) {
    abort_invalid_arg(c(
      "{.arg psi_cutoff} must be a length-2 numeric with no NAs.",
      "x" = "Got {.cls {class(psi_cutoff)[1]}} of length {length(psi_cutoff)}."
    ))
  }
  if (any(psi_cutoff < -1) || any(psi_cutoff > 1)) {
    abort_invalid_arg(c(
      "{.arg psi_cutoff} values must lie in {.code [-1, 1]}.",
      "x" = "Got {.code c({psi_cutoff[1]}, {psi_cutoff[2]})}."
    ))
  }
  if (psi_cutoff[1] >= psi_cutoff[2]) {
    abort_invalid_arg(c(
      "{.arg psi_cutoff} must satisfy {.code psi_cutoff[1] < psi_cutoff[2]} (neg, pos).",
      "x" = "Got {.code c({psi_cutoff[1]}, {psi_cutoff[2]})}."
    ))
  }

  check_scalar_number(psi_control_max, "psi_control_max", min = 0, max = 1)
  if (psi_control_max >= min(abs(psi_cutoff))) {
    abort_invalid_arg(c(
      "{.arg psi_control_max} must be strictly less than {.code min(abs(psi_cutoff))}.",
      "x" = "Got {.val {psi_control_max}}; {.code min(abs(psi_cutoff))} is {.val {min(abs(psi_cutoff))}}.",
      "i" = "Otherwise the Control pool would overlap with Positive / Negative."
    ))
  }

  # Groups
  valid_groups <- c("Negative", "Positive", "Control")
  if (!is.character(groups) || length(groups) == 0L || anyNA(groups) ||
      anyDuplicated(groups) || !all(groups %in% valid_groups)) {
    abort_invalid_arg(c(
      "{.arg groups} must be a non-empty, non-duplicated character subset of {.or {.val {valid_groups}}}.",
      "x" = "Got {.val {groups}}."
    ))
  }

  # Control bootstrap
  check_scalar_number(control_multiplier, "control_multiplier", min = 0)
  if (control_multiplier == 0) {
    abort_invalid_arg("{.arg control_multiplier} must be greater than {.val 0}.")
  }
  check_scalar_int(control_iterations, "control_iterations", min = 1)

  # Smoothing: NULL or non-negative int (0 disables)
  if (is.null(moving_average)) {
    moving_average <- 0L
  }
  check_scalar_int(moving_average, "moving_average", min = 0)

  # Significance
  check_string(stat_test, "stat_test", choices = c("fisher", "binomial"))
  check_flag(use_fdr, "use_fdr")
  check_unit_interval(fdr_threshold, "fdr_threshold")

  # Diagnostics
  check_flag(verbose, "verbose")

  as.list(environment())
}
