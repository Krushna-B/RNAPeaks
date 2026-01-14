
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
                          highlighted_region_opacity = 0.30
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
  gene_names <- unique(na.omit(gene$gene_name))
  if (length(gene_names) == 0) {
    # fall back to gene_id(s)
    gene_names <- unique(gene$gene_id)
  }
  title_txt <- paste(gene_names, collapse = ", ")

  # transcript (only if present)
  trans_ids <- unique(na.omit(gene$transcript_id))
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
  g <- ggplot() +
    # Clean base theme and allow drawing outside panel
    theme_classic() +

    theme(
        plot.margin = margin(plot_top_margin, plot_right_margin, plot_bottom_margin, left_margin_pt)  # top, right, bottom, LEFT
    ) +

    # Intron baselines with small directional arrows along each intron row
    geom_segment(data=Intron_s,
                 aes(x=dir_start, xend=dir_end, y=mid_y, yend=mid_y),
                 color=intron_color, linewidth = intron_linewidth,
                 lineend = "butt",
                 # arrow=arrow(type="open", length=unit(intron_arrow_len_in,"inches"))
                 ) +
    geom_rect(data=UTRs,
              aes(xmin=start-0.5, xmax=end+0.5, ymin=y_start, ymax=y_end),
              fill=utr_fill, color=NA)+

    geom_segment(data=Arrow_df,
                 aes(x=x, xend=xend, y=y, yend=yend),
                 arrow=arrow(type="open", length=arrow_head),
                 color=intron_color, linewidth=intron_linewidth)+

    # Exons (CDS) and UTRs as horizontal blocks on the gene row
    geom_rect(data=Exons,
              aes(xmin=start-0.5, xmax=end+0.5, ymin=y_start, ymax=y_end),
              fill=exon_fill, color=NA) +


    #Bottom row with spaced out bp markers
    # xlim(x_min,x_max)+
    scale_x_continuous(
      name   = "Genomic position (bp)",
      limits = c(x_min - axis_pad_bp, x_max + axis_pad_bp),
      breaks = seq(x_min,x_max,length.out = axis_breaks_n),
      labels = scales::label_comma(accuracy = 1),
      expand = c(0, 0),
    )+
    coord_cartesian( clip = "off") +
    # Title + general theme cleanup
    theme(strip.background = element_blank()) +

    ggtitle(
      plot_title,
      subtitle = plot_sub
    ) +

    # Alternating background bands behind each protein row
    geom_rect(
      data = Bg_tab,
      aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end),
      fill = Bg_tab$fill, color = NA
    ) +

    # Thin separators at the top of each band
    geom_segment(
      data = Bg_tab,
      aes(x = x_start, xend = x_end, y = y_end, yend = y_end),
      color = band_sep_color, linewidth = band_sep_linewidth
    )+


    theme(
      # Global style for axes and panel (hide y-axis; transparent panel/plot)

      legend.position="bottom",
      legend.direction = "horizontal",
      legend.text      = element_text(size = 5),

      axis.title.x        = element_text(size = axis_title_size,color = title_color,   face = "bold"),  # bottom
      plot.subtitle = element_text(hjust = 0.5, face = "bold",    # center + bold subtitle
                                   size = subtitle_size,
                                   color = subtitle_color,
                                   margin = margin(t = 2, b = 8)), # tuck closer under title

      axis.text.x         = element_text(size = axis_text_size),
      axis.ticks.x        = element_line(),
      axis.title.x.top    = element_text(size = axis_title_size, face = "bold"),  # top (relative)
      axis.text.x.top     = element_text(size = axis_text_size),
      axis.ticks.x.top    = element_line(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),

      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),

      panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
      plot.background = element_rect(fill = "transparent",colour = NA),
      plot.title = element_text(hjust=0.5,size=title_size,,color = title_color,   face = "bold.italic")
      )+

    # Peak rectangles for each protein row.
    # Uses per-row color found in bed$col; outline hidden ("transparent").
    geom_rect(inherit.aes = FALSE,
              data=bed,
              mapping=aes(xmin=start,xmax=end,ymin=y_start,ymax=y_end),
              fill=bed$col,
              alpha= peak_alpha,
              color = peak_border_color,
              linewidth = peak_border_linewidth)+


    # Protein labels
    geom_text(
      inherit.aes = FALSE,
      data = transform(bed[!duplicated(bed$group_name), ],
                       label_y = (y_start + y_end) / 2),
      mapping = aes(label = group_name, x = xpos, y = label_y),
      hjust = 1,
      size = label_size,
      color = label_color,
    ) +



    # 5′ / 3′ strand tags
    geom_text(data = labs$left,aes(label=Label,x=X,y=Y,hjust=1),size=strand_label_size,color = strand_label_color)+
    geom_text(data = labs$right,aes(label=Label,x=X,y=Y,hjust=0),size=strand_label_size,color = strand_label_color)

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
    g <- g + geom_text(
      data = gl,
      aes(x = label_x, y = label_y, label = label),
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


  return(g)
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
  w_pt <- convertWidth(
    grobWidth(textGrob(longest_label, gp = gpar(fontsize = label_size_pt))),
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
