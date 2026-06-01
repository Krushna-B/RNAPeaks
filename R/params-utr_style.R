#' Visual styling for the UTR plot
#'
#' Returns a validated list of every visual setting used by
#' [plot_utr_binding()].
#'
#' @param line_width Width of the per-track frequency lines.
#' @param line_alpha Opacity of the per-track frequency lines (0-1).
#' @param utr_fill Fill color for the thin UTR box in the schematic.
#' @param cds_fill Fill color for the adjacent (taller) CDS box in the
#'   schematic.
#' @param schematic_height Height of the schematic boxes as a fraction of
#'   the maximum plotted frequency. Larger values make the gene diagram
#'   taller relative to the curve.
#' @param label_size Font size of the `5' UTR` / `3' UTR` / `CDS` part
#'   labels in the schematic.
#' @param label_position Where to place the part labels: `"ends"` (at the
#'   outer end of each box, clear of the percentage markers) or `"center"`
#'   (centered under each box).
#' @param pct_label_size Font size of the `0%`-`100%` position markers.
#' @param pct_label_color Color of the position markers (ticks and text).
#' @param palette Character vector of colors for the curves when more than
#'   one BED track is plotted. If unnamed it is applied positionally
#'   (recycled when there are more tracks than colors). If *named* by BED
#'   track name, each track gets its mapped color regardless of order.
#' @param single_track_color Color used when `track_colors` is `NULL` and
#'   exactly one BED track is plotted.
#' @param title_size,title_color Plot title appearance.
#' @param axis_text_size Font size for axis tick labels.
#' @param ylab Y-axis label.
#' @param legend_position One of `"bottom"`, `"top"`, `"left"`, `"right"`,
#'   or `"none"`.
#'
#' @return A named list of validated styling parameters.
#' @family params
#' @export
utr_style <- function(
  # Line
  line_width         = 0.8,
  line_alpha         = 1,

  # Schematic
  utr_fill           = "lightgray",
  cds_fill           = "navy",
  schematic_height   = 0.06,
  label_size         = 3.2,
  label_position     = "ends",
  pct_label_size     = 3,
  pct_label_color    = "grey30",

  # Track colors (fallback when track_colors arg is NULL)
  palette            = c("black", "red", "blue", "darkgreen",
                         "purple", "orange", "brown", "magenta"),
  single_track_color = "blue",

  # Title
  title_size         = 20,
  title_color        = "black",

  # Axis
  axis_text_size     = 11,
  ylab               = "Frequency",

  # Legend
  legend_position    = "bottom"
) {
  # Line
  check_scalar_number(line_width, "line_width", min = 0)
  check_unit_interval(line_alpha, "line_alpha")

  # Schematic
  check_color(utr_fill, "utr_fill", allow_na = FALSE)
  check_color(cds_fill, "cds_fill", allow_na = FALSE)
  check_scalar_number(schematic_height, "schematic_height", min = 0)
  check_scalar_number(label_size,       "label_size",       min = 0)
  check_string(label_position, "label_position",
               choices = c("ends", "center"))
  check_scalar_number(pct_label_size,   "pct_label_size",   min = 0)
  check_color(pct_label_color, "pct_label_color", allow_na = FALSE)

  # Palette
  if (!is.character(palette) || length(palette) == 0L || anyNA(palette)) {
    abort_invalid_arg(c(
      "{.arg palette} must be a non-empty character vector of colors.",
      "x" = "Got {.cls {class(palette)[1]}} of length {length(palette)}."
    ))
  }
  for (i in seq_along(palette)) {
    check_color(palette[i], paste0("palette[", i, "]"), allow_na = FALSE)
  }
  check_color(single_track_color, "single_track_color", allow_na = FALSE)

  # Title
  check_scalar_number(title_size, "title_size", min = 0)
  check_color(title_color, "title_color", allow_na = FALSE)

  # Axis
  check_scalar_number(axis_text_size, "axis_text_size", min = 0)
  check_string(ylab, "ylab")

  # Legend
  check_string(legend_position, "legend_position",
               choices = c("bottom", "top", "left", "right", "none"))

  as.list(environment())
}
