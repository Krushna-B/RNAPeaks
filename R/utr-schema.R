#' Internal: UTR plot schema
#'
#' Two-region single-distribution layout: 5' UTR (left panel) and 3' UTR
#' (right panel), each 100 bins wide. The
#' plotter draws one density curve per BED track.
#'
#' @keywords internal
#' @noRd
event_schema_utr <- list(
  n_bins = 100L,

  # x-axis layout: two 100-wide panels separated by a gap.
  plot_layout = function(n_bins = 100L, gap = 80L) {
    bin_width <- as.integer(n_bins)
    region_starts <- c(0L, bin_width + gap)
    list(
      region_starts    = region_starts,
      bin_width        = bin_width,
      gap              = gap,
      boundary_offsets = integer(0),
      x_max            = region_starts[2L] + bin_width
    )
  },

  # Two UTR rectangles joined by a thin black baseline (the // sits on
  # the baseline). 5' / 3' labels go beneath each rectangle.
  build_schematic_layers = function(layout, style, y_min, exon_height) {
    utr_fill     <- style$utr_fill     %||% "lightgray"
    baseline_col <- style$baseline_col %||% "black"

    left_xmin  <- layout$region_starts[1L]
    left_xmax  <- layout$region_starts[1L] + layout$bin_width
    right_xmin <- layout$region_starts[2L]
    right_xmax <- layout$region_starts[2L] + layout$bin_width

    rect_df <- data.frame(
      xmin = c(left_xmin,  right_xmin),
      xmax = c(left_xmax,  right_xmax),
      ymin = rep(y_min - exon_height, 2L),
      ymax = rep(y_min, 2L),
      fill = c(utr_fill, utr_fill),
      stringsAsFactors = FALSE
    )

    mid_y <- y_min - exon_height / 2
    seg_df <- data.frame(
      x    = left_xmax,
      xend = right_xmin,
      y    = mid_y,
      yend = mid_y
    )

    break_x <- left_xmax + layout$gap / 2

    label_df <- data.frame(
      x     = c((left_xmin + left_xmax) / 2,
                (right_xmin + right_xmax) / 2),
      y     = y_min - exon_height * 1.6,
      label = c("5'", "3'"),
      stringsAsFactors = FALSE
    )

    list(
      ggplot2::geom_segment(
        data    = seg_df,
        mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        color = baseline_col, linewidth = 1.5, inherit.aes = FALSE
      ),
      ggplot2::geom_rect(
        data    = rect_df,
        mapping = ggplot2::aes(xmin = xmin, xmax = xmax,
                                ymin = ymin, ymax = ymax, fill = fill),
        color = "black", linewidth = 0.5, inherit.aes = FALSE
      ),
      ggplot2::annotate("text", x = break_x, y = mid_y,
                        label = "//", size = 7,
                        fontface = "bold", vjust = 0.5),
      ggplot2::geom_text(
        data    = label_df,
        mapping = ggplot2::aes(x = x, y = y, label = label),
        size = 5, fontface = "bold", inherit.aes = FALSE
      )
    )
  }
)
