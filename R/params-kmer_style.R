#' Visual styling for the k-mer enrichment plots
#'
#' Returns a validated list of every visual setting used by the scatter and
#' rank plots produced by [kmer_enrichment()].
#'
#' @param point_color Fill color for the k-mer points.
#' @param point_size Size of the k-mer points.
#' @param point_alpha Opacity of the k-mer points (0-1).
#' @param label_size Font size of the top-`n` k-mer text labels on the
#'   scatter plot.
#' @param label_check_overlap If `TRUE`, ggplot2 hides labels that would
#'   overlap already-drawn ones rather than plotting them on top of each
#'   other.
#' @param ref_line_color Color of the reference lines (the `y = x` diagonal
#'   on the scatter and the zero line on the rank plot).
#' @param ref_line_type Line type of the reference lines. One of `"solid"`,
#'   `"dashed"`, `"dotted"`, `"dotdash"`, `"longdash"`, or `"twodash"`.
#' @param title_size,title_color Plot title appearance.
#' @param axis_text_size Font size for axis tick labels.
#'
#' @return A named list of validated styling parameters.
#' @family params
#' @export
kmer_style <- function(
  # Points
  point_color         = "#3B6FB6",
  point_size          = 1.4,
  point_alpha         = 0.6,

  # k-mer labels
  label_size          = 3,
  label_check_overlap = TRUE,

  # Reference lines
  ref_line_color      = "grey60",
  ref_line_type       = "dashed",

  # Title
  title_size          = 14,
  title_color         = "black",

  # Axis
  axis_text_size      = 11
) {
  # Points
  check_color(point_color, "point_color", allow_na = FALSE)
  check_scalar_number(point_size, "point_size", min = 0)
  check_unit_interval(point_alpha, "point_alpha")

  # Labels
  check_scalar_number(label_size, "label_size", min = 0)
  check_flag(label_check_overlap, "label_check_overlap")

  # Reference lines
  check_color(ref_line_color, "ref_line_color", allow_na = FALSE)
  check_string(ref_line_type, "ref_line_type",
               choices = c("solid", "dashed", "dotted", "dotdash",
                           "longdash", "twodash"))

  # Title
  check_scalar_number(title_size, "title_size", min = 0)
  check_color(title_color, "title_color", allow_na = FALSE)

  # Axis
  check_scalar_number(axis_text_size, "axis_text_size", min = 0)

  # Return only the declared parameters, not the helper locals used for validation.
  mget(names(formals(kmer_style)))
}
