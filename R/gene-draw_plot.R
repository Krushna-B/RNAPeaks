#' Render the peaks-on-region ggplot
#'
#' Internal plotter used by both the single-transcript and multi-gene
#' region workflows. All visual settings live in `style`
#' The structural input is a single `region` data frame.
#'
#' @param region Data frame of structural rows (one transcript or many)
#'   with `start`, `end`, `strand`, `seqnames`, `type`, `y_start`, `y_end`,
#'   and (optionally) `gene_name`, `gene_id`, `transcript_id`. `type` is
#'   used to split rows into exons / UTRs / introns.
#' @param peaks Prepared peaks with `start`, `end`, `y_start`, `y_end`,
#'   `group_name`, `col`, `xpos`.
#' @param is_region `TRUE` for multi-gene region plots; switches title to
#'   coordinates and enables per-gene labels.
#' @param style List of style settings; see [peaks_plot_style()].
#'
#' @return A ggplot object.
#' @noRd
#' @family gene
draw_plot <- function(region,
                      peaks,
                      is_region = FALSE,
                      style     = peaks_plot_style()) {
  # Validate Params
  check_data_frame(region, "region",
                   required = c("start", "end", "strand", "seqnames",
                                "type", "y_start", "y_end"))
  check_data_frame(peaks, "peaks",
                   required = c("start", "end", "y_start", "y_end",
                                "group_name"))
  if (is_region) {
    check_data_frame(region, "region", required = "gene_id")
  }

  # Identify Strand
  strand <- region$strand[1]
  if (is.na(strand) || !strand %in% c("+", "-")) {
    abort_invalid_arg(c(
      "{.code region$strand[1]} must be {.val +} or {.val -}.",
      "x" = "Got {.val {strand}}."
    ))
  }

  # Flip View; Only Applies to (-) Strand
  flip <- style$five_to_three
  if (flip && strand == "+") {
    cli::cli_inform("{.code five_to_three} only applies on {.val -} strand; ignored.")
    flip <- FALSE
  }

  # Calculate Axis Range
  xlim  <- c(min(region$start), max(region$end))
  x_min <- xlim[1]; x_max <- xlim[2]
  pad   <- style$axis_pad_bp
  x_lo  <- x_min - pad
  x_hi  <- x_max + pad

  # Split Region into different features
  exons   <- region[region$type %in% c("CDS", "exon"), , drop = FALSE]
  utrs    <- region[region$type %in% c("five_prime_utr",
                                       "three_prime_utr"), , drop = FALSE]
  introns <- region[region$type == "intron", , drop = FALSE]

  if (nrow(introns)) {
    introns$mid_y     <- (introns$y_start + introns$y_end) / 2
    introns$dir_start <- ifelse(introns$strand == "+", introns$start, introns$end)
    introns$dir_end   <- ifelse(introns$strand == "+", introns$end,   introns$start)
  }

  # arrows <- build_intron_arrows(introns, style)

  # left margin: fit longest protein label,
  # and in region mode also the longest gene/transcript label drawn on the left.
  left_margin <- style$plot_left_margin
  if (is.null(left_margin)) {
    margins <- find_left_margin(peaks$group_name, style$protein_label_size)
    if (is_region) {
      margins <- c(margins,
                   find_left_margin(unique(gene_label_text(region)),
                                    style$gene_label_size))
    }
    left_margin <- max(margins)
  }

  # Protein Label Position Calculation
  # Flipped view: labels go to the right of the gene (3' end visually on the left).
  # Default + region: labels go to the left of the gene/region.
  peaks$xpos <- if (flip) {
    x_max + style$protein_label_x_offset_bp
  } else {
    x_min - style$protein_label_x_offset_bp
  }

  # Background Bands for the peaks
  bg <- background_table(peaks, xlim)
  bg$fill <- ifelse(seq_len(nrow(bg)) %% 2 == 0,
                    style$band_even_fill,
                    style$band_odd_fill)

  # Build 5' and 3' Tags
  tags <- make_strand_labels(region)

  # Build Title
  coord_str <- paste0("Chr ", as.character(region$seqnames[1]),
                      style$subtitle_sep,
                      scales::comma(x_min), "-", scales::comma(x_max), " bp")
  plot_title <- if (is_region) coord_str else gene_title(region)
  plot_sub   <- if (is_region) NULL      else coord_str

  # Build Scale: Account for flip
  breaks <- seq(x_min, x_max, length.out = style$axis_breaks_n)
  scale_x <- if (flip) {
    ggplot2::scale_x_reverse(
      name   = "Genomic position (bp)",
      limits = c(x_hi, x_lo),
      breaks = breaks,
      labels = scales::label_comma(accuracy = 1),
      expand = c(0, 0)
    )
  } else {
    ggplot2::scale_x_continuous(
      name   = "Genomic position (bp)",
      limits = c(x_lo, x_hi),
      breaks = breaks,
      labels = scales::label_comma(accuracy = 1),
      expand = c(0, 0)
    )
  }

  # Build Plot
  g <- ggplot2::ggplot() +
    ggplot2::theme_classic() +

    # Intron baselin
    ggplot2::geom_segment(
      data    = introns,
      mapping = ggplot2::aes(x = dir_start, xend = dir_end,
                             y = mid_y,     yend = mid_y),
      color     = style$intron_color,
      linewidth = style$intron_linewidth,
      lineend   = "butt"
    ) +

    # UTRs
    ggplot2::geom_rect(
      data    = utrs,
      mapping = ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5,
                             ymin = y_start,    ymax = y_end),
      fill = style$utr_color, color = NA
    ) +

    # Direction arrows along introns
    # ggplot2::geom_segment(
    #   data    = arrows,
    #   mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    #   arrow = grid::arrow(
    #     type   = "open",
    #     length = grid::unit(style$intron_arrow_len_in, "inches")
    #   ),
    #   color     = style$intron_color,
    #   linewidth = style$intron_linewidth
    # ) +

    # Exons
    ggplot2::geom_rect(
      data    = exons,
      mapping = ggplot2::aes(xmin = start - 0.5, xmax = end + 0.5,
                             ymin = y_start,    ymax = y_end),
      fill = style$exon_color, color = NA
    ) +

    scale_x +
    ggplot2::coord_cartesian(clip = "off") +

    # Alternating background bands + thin separators
    ggplot2::geom_rect(
      data    = bg,
      mapping = ggplot2::aes(xmin = x_start, xmax = x_end,
                             ymin = y_start, ymax = y_end),
      fill = bg$fill, color = NA
    ) +
    ggplot2::geom_segment(
      data    = bg,
      mapping = ggplot2::aes(x = x_start, xend = x_end,
                             y = y_end,   yend = y_end),
      color     = style$band_sep_color,
      linewidth = style$band_sep_linewidth
    ) +

    # Peaks
    ggplot2::geom_rect(
      data    = peaks,
      mapping = ggplot2::aes(xmin = start,   xmax = end,
                             ymin = y_start, ymax = y_end),
      inherit.aes = FALSE,
      fill        = style$peak_color,
      alpha       = style$peak_alpha,
      color       = style$peak_border_color,
      linewidth   = style$peak_border_linewidth
    ) +

    # Protein row labels
    ggplot2::geom_text(
      data        = transform(peaks[!duplicated(peaks$group_name), ],
                              label_y = (y_start + y_end) / 2),
      mapping     = ggplot2::aes(label = group_name, x = xpos, y = label_y),
      inherit.aes = FALSE,
      hjust       = 1,
      size        = style$protein_label_size,
      color       = style$protein_label_color
    ) +

    # 5'/3' tags
    ggplot2::geom_text(
      data    = tags$left,
      mapping = ggplot2::aes(label = Label, x = X, y = Y),
      hjust = 1, size = style$strand_label_size, color = style$strand_label_color
    ) +
    ggplot2::geom_text(
      data    = tags$right,
      mapping = ggplot2::aes(label = Label, x = X, y = Y),
      hjust = 0, size = style$strand_label_size, color = style$strand_label_color
    ) +

    ggplot2::ggtitle(plot_title, subtitle = plot_sub) +
    plot_theme(style, left_margin)

  # Region Plot: Transcript labels on side
  if (is_region) g <- g + region_gene_labels(region, xlim, style)

  # Highlighting Region
  g <- add_highlight_band(g, c(x_lo, x_hi), style)

  # Junction Lines
  if (isTRUE(style$show_junctions)) {
    jx <- unique(c(exons$start, exons$end, utrs$start, utrs$end))
    if (length(jx)) {
      g <- g + ggplot2::geom_vline(
        xintercept = jx,
        linetype   = style$junction_linetype,
        color      = style$junction_color,
        linewidth  = style$junction_linewidth,
        alpha      = style$junction_alpha
      )
    }
  }

  g
}



#Helper Functions for gene-draw_plot
#ggplot theme block
plot_theme <- function(style, left_margin) {
  ggplot2::theme(
    plot.margin   = ggplot2::margin(style$plot_top_margin,
                                    style$plot_right_margin,
                                    style$plot_bottom_margin,
                                    left_margin),
    plot.title    = ggplot2::element_text(
      hjust = 0.5, size = style$title_size,
      color = style$title_color, face = "bold.italic"
    ),
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5, size = style$subtitle_size,
      color = style$subtitle_color, face = "bold",
      margin = ggplot2::margin(t = 2, b = 8)
    ),
    axis.title.x = ggplot2::element_text(
      size = style$axis_title_size,
      color = style$title_color,
      face = "bold"
    ),
    axis.text.x  = ggplot2::element_text(size = style$axis_text_size),
    axis.ticks.x = ggplot2::element_line(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.y  = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    panel.grid       = ggplot2::element_blank(),
    panel.border     = ggplot2::element_blank(),
    axis.line        = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
    strip.background = ggplot2::element_blank(),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.text      = ggplot2::element_text(size = 5)
  )
}

# Title built from gene name's + transcript
gene_title <- function(region) {
  names_ <- unique(stats::na.omit(region$gene_name))
  if (!length(names_)) names_ <- unique(region$gene_id)
  out <- paste(names_, collapse = ", ")

  tx <- unique(stats::na.omit(region$transcript_id))
  if (length(tx)) out <- paste0(out, " (", tx[1], ")")
  out
}

# Resolve the label text transcripts when plot region
gene_label_text <- function(region) {
  name_col <- as.character(region$gene_name)
  id_col   <- as.character(region$gene_id)
  use_id   <- is.na(name_col) | !nzchar(name_col)
  ifelse(use_id, id_col, name_col)
}

# Per-transcript labels for a region plot
region_gene_labels <- function(region, xlim, style) {
  agg <- stats::aggregate(y_start ~ gene_id, region, mean)

  lab_map <- unique(data.frame(
    gene_id = as.character(region$gene_id),
    label   = gene_label_text(region),
    stringsAsFactors = FALSE
  ))
  df <- merge(agg, lab_map, by = "gene_id", all.x = TRUE)
  df$x <- xlim[1] - style$axis_pad_bp * style$gene_label_x_offset

  ggplot2::geom_text(
    data    = df,
    mapping = ggplot2::aes(x = x, y = y_start, label = label),
    hjust = 1,
    size  = style$gene_label_size,
    color = style$gene_label_color
  )
}

# One band per unique (y_start, y_end) and spans x_min to x_max.
background_table <- function(peaks, xlim) {
  bg <- unique(peaks[c("y_start", "y_end")])
  bg <- bg[order(bg$y_start), , drop = FALSE]
  bg$x_start <- xlim[1]
  bg$x_end   <- xlim[2]
  bg
}

# Width for the left margin based on the longest protein and gene name.
find_left_margin <- function(labels, label_size, pad_pt = 8) {
  if (!length(labels)) return(pad_pt)

  pts_per_mm <- 72.27 / 25.4
  longest    <- labels[which.max(nchar(labels))]
  w_pt <- grid::convertWidth(
    grid::grobWidth(grid::textGrob(
      longest,
      gp = grid::gpar(fontsize = label_size * pts_per_mm)
    )),
    "pt", valueOnly = TRUE
  )
  w_pt + pad_pt
}

# Overlay a highlight rectangle, clipped to the plotted x range.
add_highlight_band <- function(g, xlim_padded, style) {
  if (is.null(style$highlight)) return(g)

  xmin <- max(xlim_padded[1], style$highlight[1])
  xmax <- min(xlim_padded[2], style$highlight[2])
  if (xmin >= xmax) return(g)

  g + ggplot2::annotate(
    "rect",
    xmin = xmin, xmax = xmax,
    ymin = -Inf, ymax = Inf,
    fill  = style$highlight_color,
    alpha = style$highlight_opacity
  )
}
