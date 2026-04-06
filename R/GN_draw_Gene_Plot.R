
Draw_Gene_Plot <- function(Gene_s,
                          Exons,
                           UTRs,
                           Intron_s,
                           Arrow_df,
                           bed,
                           gene,
                           peak_col,
                           peaks_width,

            # styling knobs (defaults = are set)
                          exon_fill = "navy",
                          utr_fill  = "lightgray",

                          band_even_fill = "#F7F8FA",
                          band_odd_fill  = "#FFFFFF",
                          band_sep_color = "#E5E7EB",
                          band_sep_linewidth = 0.4,

                          intron_color = "gray60",
                          intron_linewidth = 0.9,
                          intron_arrow_len_in = 0.15,       # inches

                          peak_alpha = 0.95,
                          peak_border_color = NA,
                          peak_border_linewidth = 0.4,

                          label_size = 5,                    # protein labels
                          label_color = "black",
                          strand_label_size = 5,             # 5′/3′ labels
                          strand_label_color = "black",

                          title_size = 25,
                          title_color = "black",
                          subtitle_size = 12,
                          subtitle_color = "black",
                          axis_title_size = 11,
                          axis_text_size  = 9,

                          x_lims = NULL,
                          axis_pad_bp   = 500,               # ± window (bp)
                          axis_breaks_n = 5,
                          subtitle_sep = ": ",            # ASCII-safe separator

                          plot_right_margin = 50,
                          plot_top_margin = 30,
                          plot_bottom_margin = 30,
                          plot_left_margin = NULL,

                          is_region_plot = FALSE,  #Flag for is this a region plot
                          gene_label_x_offset = 0.25,     # fraction of axis_pad_bp
                          gene_label_size     = 5,
                          gene_label_color    = "black",
                          max_proteins = 40,
                          protein_label_x_offset = 100,

                          highlighted_region_start = NULL,
                          highlighted_region_stop = NULL,
                          highlighted_region_color = "pink",
                          highlighted_region_opacity = 0.30,

                          five_to_three = FALSE,
                          gene_strand = "+",

                          show_junctions = FALSE,
                          junction_color = "gray40",
                          junction_linetype = "dashed",
                          junction_linewidth = 0.4,
                          junction_alpha = 0.7
                          ){




  # Compute a left plot margin large enough to fit the longest protein label
  if(is.null(plot_left_margin)){
    left_margin_pt <- Finding_Left_Margin(bed)
  } else{
    left_margin_pt <- plot_left_margin
  }


  #Base Pair Locations (5 evenly spaced out)
  if(is.null(x_lims)){
    x_min <- min(gene$start)
    x_max <- max(gene$end)
  } else{
    x_min <- x_lims[1]
    x_max <- x_lims[2]
  }

  # Determine if we need to flip x-axis (5' to 3' orientation for negative strand)
  should_flip <- isTRUE(five_to_three) && gene_strand == "-"

  #Protein limiting
  all_proteins <- unique(bed$group_name)
  n_proteins   <- length(all_proteins)
  if (n_proteins > max_proteins) {

    warning(sprintf(
      "Plot has %d protein tracks; showing only the first %d.",
      n_proteins, max_proteins
    ))

    keep_names <- all_proteins[1:max_proteins]
    bed <- bed[bed$group_name %in% keep_names, , drop = FALSE]
  }

  #Setting Protein label offset for plotting region
  if (is_region_plot) {
    bed$xpos <- x_min - protein_label_x_offset
  }

  # When flipping, move labels to the right side (which appears on left after flip)
  # Position labels just past x_max, similar offset as normal case uses past x_min
  if (should_flip) {
    bed$xpos <- x_max + protein_label_x_offset
    label_hjust <- 1  # Right-align so text extends left (into margin)
  } else {
    label_hjust <- 1
  }
  # Keep margins the same - left margin always holds label space
  final_left_margin <- left_margin_pt
  final_right_margin <- plot_right_margin
  # Arrow head size
  arrow_head <- grid::unit(intron_arrow_len_in, "inches")

  #Sub Title
  chr <- as.character(Gene_s$seqnames[1])
  sub_txt <- paste0(
    "Chr ", chr, subtitle_sep,
    scales::comma(x_min), "-", scales::comma(x_max), " bp"
  )



  # Build alternating background stripes so peak rows are easy to read
  Bg_tab<-Background_Table(df=bed,Start=min(gene$start),End=max(gene$end))
  Bg_tab$band_id <- seq_len(nrow(Bg_tab))
  Bg_tab$fill <- ifelse(Bg_tab$band_id %% 2 == 0, band_even_fill, band_odd_fill)

  # 5′ / 3′ tag positions
  labs <- Make_Strand_Labels(Gene_s, offset = 100)

  # Direction for intron baseline arrows:
  Intron_s$dir_start <- ifelse(Intron_s$strand == "+", Intron_s$start, Intron_s$end)
  Intron_s$dir_end   <- ifelse(Intron_s$strand == "+", Intron_s$end, Intron_s$start)


  # ---- title text ----
  gene_names <- unique(stats::na.omit(gene$gene_name))
  if (length(gene_names) == 0) {
    # fall back to gene_id(s)
    gene_names <- unique(gene$gene_id)
  }
  title_txt <- paste(gene_names, collapse = ", ")

  # transcript (only if present)
  trans_ids <- unique(stats::na.omit(gene$transcript_id))
  if (length(trans_ids) > 0) {
    title_txt <- paste0(title_txt, " (", trans_ids[1], ")")
  }

  #If region plot then only show positions
  if (isTRUE(is_region_plot)) {
    plot_title <- sub_txt
    plot_sub   <- NULL
  } else {
    plot_title <- title_txt
    plot_sub   <- sub_txt
  }



#--------PLOT----------
  g <- ggplot2::ggplot() +
    # Clean base theme and allow drawing outside panel
    ggplot2::theme_classic() +

    ggplot2::theme(
        plot.margin = ggplot2::margin(plot_top_margin, final_right_margin, plot_bottom_margin, final_left_margin)  # top, right, bottom, LEFT
    ) +

    # Intron baselines with small directional arrows along each intron row
    ggplot2::geom_segment(data=Intron_s,
                 ggplot2::aes(x=dir_start, xend=dir_end, y=mid_y, yend=mid_y),
                 color=intron_color, linewidth = intron_linewidth,
                 lineend = "butt",
                 # arrow=arrow(type="open", length=unit(intron_arrow_len_in,"inches"))
                 ) +
    ggplot2::geom_rect(data=UTRs,
              ggplot2::aes(xmin=start-0.5, xmax=end+0.5, ymin=y_start, ymax=y_end),
              fill=utr_fill, color=NA)+

    ggplot2::geom_segment(data=Arrow_df,
                 ggplot2::aes(x=x, xend=xend, y=y, yend=yend),
                 arrow=grid::arrow(type="open", length=arrow_head),
                 color=intron_color, linewidth=intron_linewidth)+

    # Exons (CDS) and UTRs as horizontal blocks on the gene row
    ggplot2::geom_rect(data=Exons,
              ggplot2::aes(xmin=start-0.5, xmax=end+0.5, ymin=y_start, ymax=y_end),
              fill=exon_fill, color=NA) +


    #Bottom row with spaced out bp markers
    # xlim(x_min,x_max)+
    # Use scale_x_reverse for 5' to 3' orientation on negative strand
    {
      if (should_flip) {
        ggplot2::scale_x_reverse(
          name   = "Genomic position (bp)",
          limits = c(x_max + axis_pad_bp, x_min - axis_pad_bp),
          breaks = seq(x_min, x_max, length.out = axis_breaks_n),
          labels = scales::label_comma(accuracy = 1),
          expand = c(0, 0)
        )
      } else {
        ggplot2::scale_x_continuous(
          name   = "Genomic position (bp)",
          limits = c(x_min - axis_pad_bp, x_max + axis_pad_bp),
          breaks = seq(x_min, x_max, length.out = axis_breaks_n),
          labels = scales::label_comma(accuracy = 1),
          expand = c(0, 0)
        )
      }
    } +
    ggplot2::coord_cartesian( clip = "off") +
    # Title + general theme cleanup
    ggplot2::theme(strip.background = ggplot2::element_blank()) +

    ggplot2::ggtitle(
      plot_title,
      subtitle = plot_sub
    ) +

    # Alternating background bands behind each protein row
    ggplot2::geom_rect(
      data = Bg_tab,
      ggplot2::aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end),
      fill = Bg_tab$fill, color = NA
    ) +

    # Thin separators at the top of each band
    ggplot2::geom_segment(
      data = Bg_tab,
      ggplot2::aes(x = x_start, xend = x_end, y = y_end, yend = y_end),
      color = band_sep_color, linewidth = band_sep_linewidth
    )+


    ggplot2::theme(
      # Global style for axes and panel (hide y-axis; transparent panel/plot)

      legend.position="bottom",
      legend.direction = "horizontal",
      legend.text      = ggplot2::element_text(size = 5),

      axis.title.x        = ggplot2::element_text(size = axis_title_size,color = title_color,   face = "bold"),  # bottom
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "bold",    # center + bold subtitle
                                   size = subtitle_size,
                                   color = subtitle_color,
                                   margin = ggplot2::margin(t = 2, b = 8)), # tuck closer under title

      axis.text.x         = ggplot2::element_text(size = axis_text_size),
      axis.ticks.x        = ggplot2::element_line(),
      axis.title.x.top    = ggplot2::element_text(size = axis_title_size, face = "bold"),  # top (relative)
      axis.text.x.top     = ggplot2::element_text(size = axis_text_size),
      axis.ticks.x.top    = ggplot2::element_line(),
      axis.title.y=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank(),

      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),

      panel.background = ggplot2::element_rect(fill = "transparent",colour = NA), # or theme_blank()
      plot.background = ggplot2::element_rect(fill = "transparent",colour = NA),
      plot.title = ggplot2::element_text(hjust=0.5,size=title_size,,color = title_color,   face = "bold.italic")
      )+

    # Peak rectangles for each protein row.
    # Uses per-row color found in bed$col; outline hidden ("transparent").
    ggplot2::geom_rect(inherit.aes = FALSE,
              data=bed,
              mapping=ggplot2::aes(xmin=start,xmax=end,ymin=y_start,ymax=y_end),
              fill=bed$col,
              alpha= peak_alpha,
              color = peak_border_color,
              linewidth = peak_border_linewidth)+


    # Protein labels
    ggplot2::geom_text(
      inherit.aes = FALSE,
      data = transform(bed[!duplicated(bed$group_name), ],
                       label_y = (y_start + y_end) / 2),
      mapping = ggplot2::aes(label = group_name, x = xpos, y = label_y),
      hjust = label_hjust,
      size = label_size,
      color = label_color,
    ) +



    # 5′ / 3′ strand tags
    ggplot2::geom_text(data = labs$left,ggplot2::aes(label=Label,x=X,y=Y,hjust=1),size=strand_label_size,color = strand_label_color)+
    ggplot2::geom_text(data = labs$right,ggplot2::aes(label=Label,x=X,y=Y,hjust=0),size=strand_label_size,color = strand_label_color)

  #Gene Labels Only Plot if Region Plot flagged is enabled
  if (is_region_plot) {
    # 1) pick one row per gene
    gl <- Gene_s[!duplicated(Gene_s$gene_id), ]

    # 2) force both to character (THIS is the fix)
    gl$gene_name <- as.character(gl$gene_name)
    gl$gene_id   <- as.character(gl$gene_id)

    # bad if NA, empty
    bad <- is.na(gl$gene_name) | gl$gene_name == ""

    # 3) choose label
    gl$label <- ifelse(bad, gl$gene_id, gl$gene_name)

    # 4) position
    gl$label_y <- (gl$y_start + gl$y_end) / 2
    gl$label_x <- x_min - axis_pad_bp * gene_label_x_offset

    # Add gene labels to plot
    g <- g + ggplot2::geom_text(
      data = gl,
      ggplot2::aes(x = label_x, y = label_y, label = label),
      hjust = 1,
      size  = gene_label_size,
      color = gene_label_color
    )
  }
  g <- add_highlight_band(
    g,
    hstart  = highlighted_region_start,
    hstop   = highlighted_region_stop,
    color   = highlighted_region_color,
    opacity = highlighted_region_opacity,
    x_min   = x_min,
    x_max   = x_max,
    axis_pad_bp = axis_pad_bp
  )

  if (isTRUE(show_junctions)) {
    junction_x <- unique(c(Exons$start, Exons$end, UTRs$start, UTRs$end))
    if (length(junction_x) > 0) {
      g <- g + ggplot2::geom_vline(
        xintercept = junction_x,
        linetype   = junction_linetype,
        color      = junction_color,
        linewidth  = junction_linewidth,
        alpha      = junction_alpha
      )
    }
  }


  return(g)
}

# Draw a single BAM coverage track as a ggplot panel.
#
# Returns a minimal ggplot (geom_area) that is meant to be stacked above the
# gene plot via patchwork.  The BAM label is shown as the y-axis title so it
# sits on the left-hand side of the panel, matching the gene plot layout.
#
# @param cov_df    data.frame with columns `pos` and `coverage`.
# @param label     Character label displayed on the left (y-axis title).
# @param x_min     Numeric. Left bound of the gene (matches gene plot).
# @param x_max     Numeric. Right bound of the gene (matches gene plot).
# @param axis_pad_bp Numeric. Same padding used in the gene plot.
# @param fill_col  Fill colour for the area. Default "steelblue".
# @param fill_alpha Opacity of the fill. Default 0.75.
# @param line_col  Colour of the top line. Default same as fill_col.
# @param five_to_three Logical. Reverse x-axis for negative-strand genes.
# @param gene_strand Character. "+" or "-".
#
# @return A ggplot2 object.
Draw_BAM_Track <- function(cov_df,
                           label,
                           x_min,
                           x_max,
                           axis_pad_bp    = 500,
                           fill_col       = "navy",
                           fill_alpha     = 0.75,
                           line_col       = NULL,
                           label_size     = 9,
                           axis_text_size = 8,
                           ylim           = NULL,
                           five_to_three  = FALSE,
                           gene_strand    = "+") {

  if (is.null(line_col)) line_col <- fill_col
  should_flip <- isTRUE(five_to_three) && gene_strand == "-"

  # y_display_max = the labelled top value; axis limit gets 10% headroom above it
  if (!is.null(ylim)) {
    y_min         <- ylim[1]
    y_display_max <- ylim[2]
  } else {
    y_min         <- 0
    y_display_max <- max(cov_df$coverage, na.rm = TRUE)
    if (y_display_max == 0) y_display_max <- 1
  }
  y_axis_max <- y_display_max * 1.1   # headroom so top label stays inside panel

  g <- ggplot2::ggplot(cov_df, ggplot2::aes(x = pos, y = coverage)) +
    ggplot2::geom_area(fill = fill_col, colour = line_col,
                       alpha = fill_alpha, linewidth = 0.3) +
    ggplot2::scale_y_continuous(
      limits = c(y_min, y_axis_max),
      breaks = c(y_min, y_display_max),
      labels = c(as.character(y_min), scales::comma(y_display_max)),
      expand = c(0, 0)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.x      = ggplot2::element_blank(),
      axis.text.x       = ggplot2::element_blank(),
      axis.ticks.x      = ggplot2::element_blank(),
      axis.line.x       = ggplot2::element_blank(),
      axis.title.y      = ggplot2::element_text(size = label_size, angle = 0,
                                                 vjust = 0.5, hjust = 1,
                                                 face = "bold"),
      axis.text.y       = ggplot2::element_text(size = axis_text_size),
      axis.ticks.y      = ggplot2::element_blank(),
      axis.line.y       = ggplot2::element_blank(),
      panel.grid        = ggplot2::element_blank(),
      panel.border      = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(0, 5, 0, 5)
    ) +
    ggplot2::labs(y = label)

  if (should_flip) {
    g <- g + ggplot2::scale_x_reverse(
      limits = c(x_max + axis_pad_bp, x_min - axis_pad_bp),
      expand = c(0, 0)
    )
  } else {
    g <- g + ggplot2::scale_x_continuous(
      limits = c(x_min - axis_pad_bp, x_max + axis_pad_bp),
      expand = c(0, 0)
    )
  }

  g
}


#---------Helper Functions---------
Finding_Left_Margin <- function(bed,
                                label_size_mm = 5,              # font size (mm) used to measure labels
                                pad_pt       = 8,              # extra left padding (points)
                                pts_per_mm   = 72.27 / 25.4){

  #Converts Label Size from mm -> pt
  label_size_pt <- label_size_mm * pts_per_mm

  # pick the longest label in bed file by character count
  longest_label <- bed$group_name[ which.max(nchar(bed$group_name)) ]

  # measure width of the longest label (in points)
  w_pt <- grid::convertWidth(
    grid::grobWidth(grid::textGrob(longest_label, gp = grid::gpar(fontsize = label_size_pt))),
    "pt", valueOnly = TRUE
  )

  # add a small pad so text doesn't hug the plot area
  return (w_pt + pad_pt)
}



add_highlight_band <- function(g,
                               hstart = NULL,
                               hstop  = NULL,
                               color  = "yellow",
                               opacity = 0.15,
                               x_min = NULL,
                               x_max = NULL,
                               axis_pad_bp = 0) {
  if (is.null(hstart) || is.null(hstop) || is.null(x_min) || is.null(x_max)){
    return(g)
  }

  hstart <- as.numeric(hstart)
  hstop <- as.numeric(hstop)


  if (is.na(hstart) || is.na(hstop)){
    return(g)
  }
  if (hstart > hstop) { tmp <- hstart; hstart <- hstop; hstop <- tmp }


  alpha <- max(0, min(1, as.numeric(opacity)))

  clip_min <- x_min - axis_pad_bp
  clip_max <- x_max + axis_pad_bp
  xmin <- max(clip_min, hstart)
  xmax <- min(clip_max, hstop)

  if (xmin < xmax) {
    g <- g + ggplot2::annotate(
      "rect",
      xmin = xmin, xmax = xmax,
      ymin = -Inf, ymax = Inf,
      fill = color, alpha = alpha
    )
  }
  g
}
