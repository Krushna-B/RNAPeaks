#' Internal: bootstrap-resample the Control pool
#'
#' @keywords internal
#' @noRd
bootstrap_control <- function(control_score_matrix, n_pos, n_neg, opts) {
  #Validate Params
  if (!is.matrix(control_score_matrix)) {
    abort_invalid_arg(
      "{.arg control_score_matrix} must be a matrix."
    )
  }
  if (!is.numeric(n_pos) || length(n_pos) != 1L ||
      is.na(n_pos) || n_pos < 0) {
    abort_invalid_arg(
      "{.arg n_pos} must be a single non-negative number."
    )
  }
  if (!is.numeric(n_neg) || length(n_neg) != 1L ||
      is.na(n_neg) || n_neg < 0) {
    abort_invalid_arg(
      "{.arg n_neg} must be a single non-negative number."
    )
  }

  n_control       <- nrow(control_score_matrix)
  total_positions <- ncol(control_score_matrix)

  zeros <- list(
    mean_per_position = rep(0, total_positions),
    sd_per_position   = rep(0, total_positions)
  )

  if (n_control == 0L) {
    cli::cli_warn(c(
      "No Control events available; SD ribbon will be flat at zero.",
      "i" = "Check {.arg control_pval} / {.arg psi_control_max} in {.fn splicing_options}."
    ))
    return(zeros)
  }
  if (total_positions == 0L) return(zeros)

  sample_size  <- as.integer(round((n_pos + n_neg) * opts$control_multiplier))
  n_iterations <- as.integer(opts$control_iterations)

  if (sample_size < 1L) {
    cli::cli_warn(
      "Bootstrap sample size is zero ({.code (n_pos + n_neg) * control_multiplier}); using raw Control mean with zero SD."
    )
    return(list(
      mean_per_position = colMeans(control_score_matrix),
      sd_per_position   = rep(0, total_positions)
    ))
  }

  # Sampling-count matrix S [n_control x n_iterations]
  sampled  <- sample.int(n_control, sample_size * n_iterations, replace = TRUE)
  iter_id  <- rep(seq_len(n_iterations), each = sample_size)
  flat_idx <- (iter_id - 1L) * n_control + sampled
  S <- matrix(
    tabulate(flat_idx, nbins = n_control * n_iterations),
    nrow = n_control, ncol = n_iterations
  )

  #Computer Means per positions
  iter_means <- crossprod(S, control_score_matrix) / sample_size

  mean_per_position <- colMeans(iter_means)
  if (n_iterations > 1L) {
    centered_sq     <- sweep(iter_means, 2, mean_per_position, FUN = "-")^2
    sd_per_position <- sqrt(colSums(centered_sq) / (n_iterations - 1L))
  } else {
    sd_per_position <- rep(0, total_positions)
  }

  list(
    mean_per_position = mean_per_position,
    sd_per_position   = sd_per_position
  )
}
