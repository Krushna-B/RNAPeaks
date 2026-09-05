#' Internal: UTR plot schema
#'
#' Single-region layout used for one UTR side at a time. The density curve
#' covers a 100-bin window (`[0, bin_width]`); a CDS block is drawn
#' adjacent to convey orientation: to the right of a 5' UTR (UTR -> start
#' codon) and to the left of a 3' UTR (stop codon -> UTR). Spliced
#' coordinates from the scorer already run 5' -> 3', so `0%` is the
#' transcription-/translation-start-proximal end and `100%` the distal end.
#'
#' @keywords internal
#' @noRd
event_schema_utr <- list(
  n_bins = 100L,

  # Single UTR region [0, bin_width] plus an adjacent CDS block of
  # `cds_width` bp. Side decides which end the CDS sits on (set later).
  plot_layout = function(n_bins = 100L, cds_width = 25L) {
    bin_width <- as.integer(n_bins)
    list(
      region_start = 0L,
      bin_width    = bin_width,
      cds_width    = as.integer(cds_width)
    )
  },

  # Gene-anatomy schematic for one side. A thin UTR box spans the window
  # and a taller CDS box abuts it. Percentage markers (0-100%) run
  # along the UTR box. `side` is "utr5" or "utr3".
  build_schematic_layers = function(layout, style, y_min, exon_height, side) {
    utr_fill <- style$utr_fill %||% "lightgray"
    cds_fill <- style$cds_fill %||% "white"

    bw <- layout$bin_width
    cw <- layout$cds_width
    is5 <- identical(side, "utr5")

    if (is5) {
      cds_xmin <- bw;  cds_xmax <- bw + cw
    } else {
      cds_xmin <- -cw; cds_xmax <- 0L
    }

    # CDS box is full height; UTR box is thinner and vertically centered on
    # the CDS box so the two read as one gene segment.
    base_top    <- y_min
    base_bottom <- y_min - exon_height
    mid_y       <- y_min - exon_height / 2
    utr_half    <- exon_height * 0.32

    cds_df <- data.frame(
      xmin = cds_xmin, xmax = cds_xmax,
      ymin = base_bottom, ymax = base_top,
      fill = cds_fill, stringsAsFactors = FALSE
    )
    utr_df <- data.frame(
      xmin = 0L, xmax = bw,
      ymin = mid_y - utr_half, ymax = mid_y + utr_half,
      fill = utr_fill, stringsAsFactors = FALSE
    )

    # Part-label placement. "ends" pushes each label to the outer end of
    # its box (clear of the centred percentage markers); "center" sits it
    # under the middle of the box.
    if (identical(style$label_position %||% "ends", "center")) {
      utr_lab_x <- bw / 2;        utr_hjust <- 0.5
      cds_lab_x <- (cds_xmin + cds_xmax) / 2; cds_hjust <- 0.5
    } else if (is5) {
      utr_lab_x <- 0L;      utr_hjust <- 0
      cds_lab_x <- bw + cw; cds_hjust <- 1
    } else {
      utr_lab_x <- bw;      utr_hjust <- 1
      cds_lab_x <- -cw;     cds_hjust <- 0
    }
    label_df <- data.frame(
      x     = c(utr_lab_x, cds_lab_x),
      y     = base_bottom - exon_height * 0.55,
      label = c(if (is5) "5' UTR" else "3' UTR", "CDS"),
      hjust = c(utr_hjust, cds_hjust),
      stringsAsFactors = FALSE
    )

    label_size <- style$label_size      %||% 3.2
    pct_size   <- style$pct_label_size  %||% 3
    pct_col    <- style$pct_label_color %||% "grey30"

    # Percentage markers along the UTR box.
    pct      <- c(0, 0.25, 0.5, 0.75, 1)
    tick_x   <- bw * pct
    tick_top <- base_bottom - exon_height * 1.0
    tick_bot <- base_bottom - exon_height * 1.3
    tick_df  <- data.frame(x = tick_x, y = tick_bot, yend = tick_top)
    pct_df   <- data.frame(
      x     = tick_x,
      y     = base_bottom - exon_height * 1.9,
      label = paste0(pct * 100, "%"),
      stringsAsFactors = FALSE
    )

    list(
      ggplot2::geom_rect(
        data    = cds_df,
        mapping = ggplot2::aes(xmin = xmin, xmax = xmax,
                                ymin = ymin, ymax = ymax, fill = fill),
        color = "black", linewidth = 0.5, inherit.aes = FALSE
      ),
      ggplot2::geom_rect(
        data    = utr_df,
        mapping = ggplot2::aes(xmin = xmin, xmax = xmax,
                                ymin = ymin, ymax = ymax, fill = fill),
        color = "black", linewidth = 0.5, inherit.aes = FALSE
      ),
      ggplot2::geom_text(
        data    = label_df,
        mapping = ggplot2::aes(x = x, y = y, label = label),
        hjust = label_df$hjust, size = label_size,
        fontface = "bold", inherit.aes = FALSE
      )
      #ggplot2::geom_segment(
      #  data    = tick_df,
      #  mapping = ggplot2::aes(x = x, xend = x, y = y, yend = yend),
      #  color = pct_col, linewidth = 0.4, inherit.aes = FALSE
      #),
      #ggplot2::geom_text(
      #  data    = pct_df,
      #  mapping = ggplot2::aes(x = x, y = y, label = label),
      #  size = pct_size, color = pct_col, inherit.aes = FALSE
      #)
    )
  }
)
