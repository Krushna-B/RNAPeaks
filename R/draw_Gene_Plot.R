
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

                          title_size = 30,
                          title_color = "black",
                          subtitle_size = 12,
                          subtitle_color = "black",
                          axis_title_size = 11,
                          axis_text_size  = 9,

                          axis_pad_bp   = 500,               # ± window (bp)
                          axis_breaks_n = 5,
                          subtitle_sep = " ~ "               # ASCII-safe separator
                          ){

  # Compute a left plot margin large enough to fit the longest protein label
  left_margin_pt <- Finding_Left_Margin(bed)




  #Base Pair Locations (5 evenly spaced out)
  x_min <- min(gene$start)
  x_max <- max(gene$end)
  span_bp <- x_max - x_min
  # axis_pad_bp <- min(max(round(span_bp * 0.06), 10), round(span_bp * 0.45))
  lims   <- c(x_min , x_max )

  #Arrow head size
  arrow_head <- grid::unit(intron_arrow_len_in, "inches")



  #Sub Title
  chr <- as.character(Gene_s$seqnames[1])
  sub_txt <- paste0(
    "Chr ", chr, subtitle_sep,
    scales::comma(x_min), "–", scales::comma(x_max), " bp"
  )




  # Build alternating background stripes so peak rows are easy to read
  Bg_tab<-Background_Table(df=bed,Start=min(gene$start),End=max(gene$end))
  Bg_tab$band_id <- seq_len(nrow(Bg_tab))
  Bg_tab$fill <- ifelse(Bg_tab$band_id %% 2 == 0, band_even_fill, band_odd_fill)

  # 5′ / 3′ tag positions
  labs <- Make_Strand_Labels(Gene_s, offset = 100)

  # Direction for intron baseline arrows:
  Intron_s$dir_start <- ifelse(Gene_s$strand[1] == "+", Intron_s$start, Intron_s$end)
  Intron_s$dir_end   <- ifelse(Gene_s$strand[1] == "+", Intron_s$end, Intron_s$start)


#--------PLOT----------
  g <- ggplot() +
    # Clean base theme and allow drawing outside panel
    theme_classic() +

    theme(
        plot.margin = margin(10, 20, 10, left_margin_pt)  # top, right, bottom, LEFT
    ) +

    # Intron baselines with small directional arrows along each intron row
    geom_segment(data=Intron_s,
                 aes(x=dir_start, xend=dir_end, y=mid_y, yend=mid_y),
                 color=intron_color, linewidth = intron_linewidth,
                 lineend = "butt",
                 # arrow=arrow(type="open", length=unit(intron_arrow_len_in,"inches"))
                 ) +

    geom_segment(data=Arrow_df,
                 aes(x=x, xend=xend, y=y, yend=yend),
                 arrow=arrow(type="open", length=arrow_head),
                 color=intron_color, linewidth=intron_linewidth)+

    # Exons (CDS) and UTRs as horizontal blocks on the gene row
    geom_rect(data=Exons,
              aes(xmin=start-0.5, xmax=end+0.5, ymin=y_start, ymax=y_end),
              fill=exon_fill, color=NA) +
    geom_rect(data=UTRs,
              aes(xmin=start-0.5, xmax=end+0.5, ymin=y_start, ymax=y_end),
              fill=utr_fill, color=NA)+

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
      paste(gene$gene_name[1],"(",gene$transcript_id[1],")",sep=""),
      subtitle = sub_txt
      )+

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
      plot.title = element_text(hjust=0.5,size=title_size,,color = title_color,   face = "bold")
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
    )+


    # 5′ / 3′ strand tags
    geom_text(data = labs$left,aes(label=Label,x=X,y=Y,hjust=1),size=strand_label_size,color = strand_label_color)+
    geom_text(data = labs$right,aes(label=Label,x=X,y=Y,hjust=0),size=strand_label_size,color = strand_label_color)

  return(g)

}
