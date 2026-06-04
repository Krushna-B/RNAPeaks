#' Visual styling for splicing and sequence maps
#'
#' Returns a validated list of every visual setting used when rendering
#' splicing-map / sequence-map outputs
#'
#' @param group_colors Named list mapping group name to line color. Must
#'   contain entries `Positive`, `Negative`, and `Control`.
#' @param line_width Width of the frequency lines.
#' @param line_alpha Opacity of the frequency lines (0-1).
#' @param ribbon_alpha Opacity of the SD ribbon around the Control line (0-1).
#' @param show_significance If `TRUE`, draw colored bars above the plot
#'   indicating positions where Positive / Negative differ significantly from
#'   Control.
#' @param title_size,title_color Plot title appearance.
#' @param axis_text_size Font size for axis tick labels.
#' @param ylab Y-axis label.
#' @param boundary_col Color of dashed vertical boundary lines between
#'   regions.
#' @param exon_col Fill color for exon boxes in the bottom schematic.
#' @param isoform_label_size Font size of the "Long isoform" / "Short
#'   isoform" labels in the a5ss / a3ss schematic.
#' @param isoform_label_nudge_x,isoform_label_nudge_y Horizontal / vertical
#'   shift applied to those labels (x in plot bp units, y in data units), so
#'   they can be repositioned. Both default to `0`.
#' @param legend_position One of `"bottom"`, `"top"`, `"left"`, `"right"`,
#'   or `"none"`.
#'
#' @return A named list of validated styling parameters.
#' @family params
#' @export
splicing_style <- function(
  # Group colors
  group_colors      = list(Positive = "blue",
                           Negative = "red",
                           Control  = "black"),

  # Line / ribbon
  line_width        = 0.8,
  line_alpha        = 1,
  ribbon_alpha      = 0.3,

  # Significance bars
  show_significance = TRUE,

  # Title
  title_size        = 20,
  title_color       = "black",

  # Axis
  axis_text_size    = 11,
  ylab              = "Frequency",

  # Schematic
  boundary_col      = "gray70",
  exon_col          = "navy",

  # Isoform labels (a5ss / a3ss schematic)
  isoform_label_size    = 2.6,
  isoform_label_nudge_x = 10,
  isoform_label_nudge_y = 0,

  # Legend
  legend_position   = "bottom"
) {

  # Group colors
  required_groups <- c("Positive", "Negative", "Control")
  if (!is.list(group_colors) || is.null(names(group_colors)) ||
      !all(required_groups %in% names(group_colors))) {
    abort_invalid_arg(c(
      "{.arg group_colors} must be a named list with entries {.or {.val {required_groups}}}.",
      "x" = "Got names {.val {names(group_colors)}}."
    ))
  }
  for (nm in required_groups) {
    check_color(group_colors[[nm]], paste0("group_colors$", nm), allow_na = FALSE)
  }

  # Line / ribbon
  check_scalar_number(line_width, "line_width", min = 0)
  check_unit_interval(line_alpha,   "line_alpha")
  check_unit_interval(ribbon_alpha, "ribbon_alpha")

  # Significance bars
  check_flag(show_significance, "show_significance")

  # Title
  check_scalar_number(title_size, "title_size", min = 0)
  check_color(title_color, "title_color", allow_na = FALSE)

  # Axis
  check_scalar_number(axis_text_size, "axis_text_size", min = 0)
  check_string(ylab, "ylab")

  # Schematic
  check_color(boundary_col, "boundary_col", allow_na = FALSE)
  check_color(exon_col,     "exon_col",     allow_na = FALSE)

  # Isoform labels
  check_scalar_number(isoform_label_size,    "isoform_label_size", min = 0)
  check_scalar_number(isoform_label_nudge_x, "isoform_label_nudge_x")
  check_scalar_number(isoform_label_nudge_y, "isoform_label_nudge_y")

  # Legend
  check_string(legend_position, "legend_position",
               choices = c("bottom", "top", "left", "right", "none"))

  as.list(environment())
}
