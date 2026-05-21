#' Internal: bootstrap-resample the Control pool
#'
#' Builds a sparse `[n_control x n_positions]` hit matrix from `ctrl_hits`,
#' resamples controls with replacement, and reports per-position mean and SD.
#'
#' @param ctrl_hits `list(event_id, col_idx)` from the scorer for Control,
#'   or `NULL` when Control is absent / empty.
#' @param n_control Integer number of Control events (rows in the implied
#'   matrix).
#' @param n_positions Integer total positions per event
#'   (`n_regions * region_width`).
#' @param n_pos,n_neg Integer counts of Positive / Negative events.
#' @param opts Result of [splicing_options()].
#'
#' @keywords internal
#' @noRd
bootstrap_control <- function(ctrl_hits, n_control, n_positions,
                              n_pos, n_neg, opts) {
  if (!is.numeric(n_pos) || length(n_pos) != 1L ||
      is.na(n_pos) || n_pos < 0) {
    abort_invalid_arg("{.arg n_pos} must be a single non-negative number.")
  }
  if (!is.numeric(n_neg) || length(n_neg) != 1L ||
      is.na(n_neg) || n_neg < 0) {
    abort_invalid_arg("{.arg n_neg} must be a single non-negative number.")
  }

  zeros <- list(
    mean_per_position = rep(0, n_positions),
    sd_per_position   = rep(0, n_positions)
  )

  if (n_control == 0L) {
    cli::cli_warn(c(
      "No Control events available; SD ribbon will be flat at zero.",
      "i" = "Check {.arg control_pval} / {.arg psi_control_max} in {.fn splicing_options}."
    ))
    return(zeros)
  }
  if (n_positions == 0L) return(zeros)

  sample_size  <- as.integer(round((n_pos + n_neg) * opts$control_multiplier))
  n_iterations <- as.integer(opts$control_iterations)

  ctrl_counts <- if (!is.null(ctrl_hits) && length(ctrl_hits$col_idx) > 0L) {
    tabulate(ctrl_hits$col_idx, nbins = n_positions)
  } else {
    integer(n_positions)
  }

  if (sample_size < 1L) {
    cli::cli_warn(
      "Bootstrap sample size is zero ({.code (n_pos + n_neg) * control_multiplier}); using raw Control mean with zero SD."
    )
    return(list(
      mean_per_position = ctrl_counts / n_control,
      sd_per_position   = rep(0, n_positions)
    ))
  }

  if (is.null(ctrl_hits) || length(ctrl_hits$col_idx) == 0L) {
    return(zeros)
  }

  M <- Matrix::sparseMatrix(
    i = ctrl_hits$event_id,
    j = ctrl_hits$col_idx,
    x = 1,
    dims = c(n_control, n_positions)
  )

  sampled <- sample.int(n_control, sample_size * n_iterations, replace = TRUE)
  iter_id <- rep(seq_len(n_iterations), each = sample_size)
  S <- Matrix::sparseMatrix(
    i = sampled,
    j = iter_id,
    x = 1,
    dims = c(n_control, n_iterations)
  )

  iter_means <- as.matrix(Matrix::crossprod(S, M)) / sample_size

  mean_per_position <- colMeans(iter_means)
  if (n_iterations > 1L) {
    centered_sq     <- sweep(iter_means, 2, mean_per_position, FUN = "-")^2
    sd_per_position <- sqrt(colSums(centered_sq) / (n_iterations - 1L))
  } else {
    sd_per_position <- rep(0, n_positions)
  }

  list(
    mean_per_position = mean_per_position,
    sd_per_position   = sd_per_position
  )
}
