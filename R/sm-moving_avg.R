#' Internal: moving-average smoothing for plotting
#'
#' @param x Numeric vector of length `n_regions * region_width`.
#' @param window Integer window size. `NULL` or `<= 1` returns `x` unchanged.
#' @param n_regions Number of regions concatenated end-to-end in `x`.
#' @param region_width Width of each region in positions.
#' @return Numeric vector the same length as `x`.
#' @keywords internal
#' @noRd
apply_moving_average <- function(x, window, n_regions, region_width) {
  #Validate Params
  if (!is.numeric(x)) {
    abort_invalid_arg("{.arg x} must be a numeric vector.")
  }
  expected_len <- n_regions * region_width
  if (length(x) != expected_len) {
    abort_invalid_arg(c(
      "{.arg x} length must equal {.code n_regions * region_width} = {expected_len}.",
      "x" = "Got length {length(x)}."
    ))
  }

  if (is.null(window)) return(x)
  if (!is.numeric(window) || length(window) != 1L || is.na(window)) {
    abort_invalid_arg("{.arg window} must be {.code NULL} or a single non-NA number.")
  }
  if (window <= 1L) return(x)

  #Moving Average
  N    <- length(x)
  w    <- as.integer(window)
  half <- w %/% 2L

  region_idx   <- ((seq_len(N) - 1L) %/% region_width) + 1L
  region_start <- (region_idx - 1L) * region_width + 1L
  region_end   <- region_idx * region_width

  gs     <- c(0, cumsum(x))
  starts <- pmax(seq_len(N) - half - 1L,        region_start - 1L)
  ends   <- pmin(seq_len(N) + (w - half - 1L),  region_end)

  (gs[ends + 1L] - gs[starts + 1L]) / (ends - starts)
}
