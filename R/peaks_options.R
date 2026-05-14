#' BED peak processing options for [plot_gene()] and [plot_region()]
#'
#' Returns a validated list of options controlling how the input BED is
#' filtered, grouped into tracks, ordered, and capped before plotting.
#'
#' @param split_by Optional column name in the BED used to split peaks into
#'   per-target tracks. Must not be one of the canonical columns
#'   (`chr`, `start`, `end`, `strand`).
#' @param omit Optional character vector of target names to drop.
#' @param order_by `"Count"` (most peaks first) or `"Alphabetical"`. Ignored
#'   when `order_in` is supplied.
#' @param order_in Optional character vector giving an explicit target order.
#' @param collapse Min gap (bp) for merging adjacent peaks within a target.
#'   `0` disables merging.
#' @param max_proteins Cap on the number of peak tracks displayed.
#'
#' @return A named list of validated peak-processing options.
#' @export
peaks_options <- function(split_by     = NULL,
                          omit         = NULL,
                          order_by     = "Count",
                          order_in     = NULL,
                          collapse     = 0,
                          max_proteins = 100) {

  #Validate Params

  # split_by: NULL or single non-empty string
  if (!is.null(split_by)) {
    check_string(split_by, "split_by")
    if (!nzchar(split_by)) {
      abort_invalid_arg("{.arg split_by} must be a non-empty string or {.code NULL}.")
    }
  }

  # omit: NULL or character vector with no NAs
  if (!is.null(omit) && (!is.character(omit) || anyNA(omit))) {
    abort_invalid_arg("{.arg omit} must be a character vector with no NAs, or {.code NULL}.")
  }

  # order_by: one of two values
  check_string(order_by, "order_by", choices = c("Count", "Alphabetical"))

  # order_in: NULL or character vector with no NAs
  if (!is.null(order_in) && (!is.character(order_in) || anyNA(order_in))) {
    abort_invalid_arg("{.arg order_in} must be a character vector with no NAs, or {.code NULL}.")
  }

  # collapse: non-negative scalar
  check_scalar_number(collapse, "collapse", min = 0)

  # max_proteins: positive integer
  check_scalar_int(max_proteins, "max_proteins", min = 1)

  as.list(environment())
}
