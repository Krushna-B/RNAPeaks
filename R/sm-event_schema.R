#' Internal: event-type schemas
#'
#' Schema slots:
#' \describe{
#'   \item{name}{Human-readable label.}
#'   \item{required_cols}{Character vector of columns the input table must
#'     contain. Universal columns plus event-type anchors.}
#'   \item{n_regions}{Number of regions per event.}
#'   \item{region_width}{`function(width_exon, width_intron)` returning the
#'     nominal bp width of each region. Same for every region within a
#'     schema.}
#'   \item{build_regions}{`function(events, width_exon, width_intron)`
#'     returning a `GRanges` of `nrow(events) * n_regions` ranges in
#'     event-then-region order, carrying `event_id` and `region_idx` in
#'     `mcols`.}
#'   \item{region_labels}{Length `n_regions` character vector of facet labels
#'     (transcript order, 5' -> 3').}
#'   \item{schematic}{String key for the plotter's bottom-strip
#'     drawing}
#'   \item{group_set}{`c("Negative","Positive","Control")` for trio-mode
#'     event types, or `"Single"` for single-distribution modes (e.g. UTR).}
#' }
#'
#' @keywords internal
#' @noRd
#' @name event_schema

# Columns every event-table schema must provide.
.event_required_base <- c(
  "chr", "strand", "GeneID",
  "PValue", "FDR", "IncLevelDifference"
)

# Long-only extension width for the a5ss / a3ss cartoon: the median real
# extension (to scale, capped to the intron side) when event stats exist,
# else a small fixed cue.
.schematic_ext <- function(schematic_data, layout) {
  intron_w <- layout$boundary_offsets[2L]
  if (!is.null(schematic_data) &&
      is.finite(schematic_data$median_ext) &&
      schematic_data$median_ext > 0) {
    min(as.integer(round(schematic_data$median_ext)), intron_w)
  } else {
    max(8L, round(layout$bin_width * 0.12))
  }
}

#' @rdname event_schema
#' @noRd
event_schema_se <- list(
  name          = "Skipped Exon",

  required_cols = c(
    .event_required_base,
    "exonStart_0base", "exonEnd",
    "upstreamES", "upstreamEE",
    "downstreamES", "downstreamEE"
  ),

  n_regions     = 4L,

  region_width  = function(width_exon, width_intron) {
    as.integer(width_exon + width_intron + 1L)
  },

  #How to make regions for specific schema
  build_regions = function(events, width_exon, width_intron) {
    n <- nrow(events)
    if (n == 0L) {
      return(GenomicRanges::GRanges())
    }

    upS <- events$upstreamES
    upE <- events$upstreamEE
    eS  <- events$exonStart_0base
    eE  <- events$exonEnd
    dnS <- events$downstreamES
    dnE <- events$downstreamEE

    # R1: upstream exon end into 1st intron
    s1 <- pmax(upE - width_exon,   upS); e1 <- pmin(upE + width_intron, eS)
    # R2: exon start into intron
    s2 <- pmax(eS  - width_intron, upE); e2 <- pmin(eS  + width_exon,   eE)
    # R3: exon End into exon
    s3 <- pmax(eE  - width_exon,   eS);  e3 <- pmin(eE  + width_intron, dnS)
    # R4: downstream exon start into intron
    s4 <- pmax(dnS - width_intron, eE);  e4 <- pmin(dnS + width_exon,   dnE)

    starts <- cbind(s1, s2, s3, s4)
    ends   <- cbind(e1, e2, e3, e4)

    GenomicRanges::GRanges(
      seqnames   = rep(events$chr, each = 4L),
      ranges     = IRanges::IRanges(start = as.vector(t(starts)),
                                    end   = as.vector(t(ends))),
      strand     = rep(events$strand, each = 4L),
      event_id   = rep(seq_len(n), each = 4L),
      region_idx = rep(seq_len(4L), times = n)
    )
  },

  #Labels for each region
  region_labels = c(
    "Upstream exon | intron",
    "Intron | Skipped exon",
    "Skipped exon | intron",
    "Intron | downstream exon"
  ),

  schematic     = "se",

  group_set     = c("Negative", "Positive", "Control"),

  # x-axis layout for the plotter.
  plot_layout = function(width_exon, width_intron) {
    bin_width <- as.integer(width_exon + width_intron)
    gap       <- 80L
    region_starts <- c(
      0L,
      bin_width + gap,
      2L * bin_width + gap,
      3L * bin_width + 2L * gap
    )
    list(
      region_starts    = region_starts,
      bin_width        = bin_width,
      gap              = gap,
      boundary_offsets = c(width_exon, width_intron, width_exon, width_intron),
      x_max            = region_starts[4L] + bin_width
    )
  },

  # Schematic layers for the splicing event plot
  build_schematic_layers = function(layout, style, y_min, exon_height,
                                    schematic_data = NULL) {
    bs <- layout$region_starts + layout$boundary_offsets
    exon_df <- data.frame(
      xmin = c(layout$region_starts[1L], bs[2L], bs[4L]),
      xmax = c(bs[1L], bs[3L], layout$region_starts[4L] + layout$bin_width),
      ymin = rep(y_min - exon_height, 3L),
      ymax = rep(y_min, 3L),
      fill = c("white", style$exon_col, "white"),
      stringsAsFactors = FALSE
    )
    intron_y <- y_min - exon_height / 2
    intron_df <- data.frame(
      x    = c(bs[1L], bs[3L]),
      xend = c(bs[2L], bs[4L]),
      y    = rep(intron_y, 2L),
      yend = rep(intron_y, 2L)
    )
    break_x <- c(
      layout$bin_width + layout$gap / 2,
      3L * layout$bin_width + layout$gap + layout$gap / 2
    )
    list(
      ggplot2::geom_rect(
        data = exon_df,
        ggplot2::aes(xmin = xmin, xmax = xmax,
                     ymin = ymin, ymax = ymax, fill = fill),
        color = "black", linewidth = 0.5, inherit.aes = FALSE
      ),
      ggplot2::geom_segment(
        data = intron_df,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        color = "black", linewidth = 1.5, inherit.aes = FALSE
      ),
      ggplot2::annotate("text", x = break_x, y = y_min,
                        label = "//", size = 8, fontface = "bold", vjust = 1)
    )
  }
)


#' @rdname event_schema
#' @noRd
event_schema_ri <- list(
  name          = "Retained Intron",

  required_cols = c(
    .event_required_base,
    "upstreamES", "upstreamEE",
    "downstreamES", "downstreamEE"
  ),

  n_regions     = 2L,

  region_width  = function(width_exon, width_intron) {
    as.integer(width_exon + width_intron + 1L)
  },

  #How to make regions for specific schema
  build_regions = function(events, width_exon, width_intron) {
    n <- nrow(events)
    if (n == 0L) {
      return(GenomicRanges::GRanges())
    }

    upS <- events$upstreamES
    upE <- events$upstreamEE
    dnS <- events$downstreamES
    dnE <- events$downstreamEE

    # R1: upstream exon end into retained intron
    s1 <- pmax(upE - width_exon,   upS); e1 <- pmin(upE + width_intron, dnS)
    # R2: retained intron end into downstream exon
    s2 <- pmax(dnS - width_intron, upE); e2 <- pmin(dnS + width_exon,   dnE)

    starts <- cbind(s1, s2)
    ends   <- cbind(e1, e2)

    GenomicRanges::GRanges(
      seqnames   = rep(events$chr, each = 2L),
      ranges     = IRanges::IRanges(start = as.vector(t(starts)),
                                    end   = as.vector(t(ends))),
      strand     = rep(events$strand, each = 2L),
      event_id   = rep(seq_len(n), each = 2L),
      region_idx = rep(seq_len(2L), times = n)
    )
  },

  #Regions Labels
  region_labels = c(
    "Upstream exon | intron",
    "Intron | downstream exon"
  ),

  schematic     = "ri",

  group_set     = c("Negative", "Positive", "Control"),

  # x-axis layout for the plotter.
  plot_layout = function(width_exon, width_intron) {
    bin_width <- as.integer(width_exon + width_intron)
    gap       <- 80L
    region_starts <- c(0L, bin_width + gap)
    list(
      region_starts    = region_starts,
      bin_width        = bin_width,
      gap              = gap,
      boundary_offsets = c(width_exon, width_intron),
      x_max            = region_starts[2L] + bin_width
    )
  },

  # Schematic layers for the splicing event plot
  build_schematic_layers = function(layout, style, y_min, exon_height,
                                    schematic_data = NULL) {
    bs <- layout$region_starts + layout$boundary_offsets
    exon_df <- data.frame(
      xmin = c(layout$region_starts[1L], bs[2L]),
      xmax = c(bs[1L], layout$region_starts[2L] + layout$bin_width),
      ymin = rep(y_min - exon_height, 2L),
      ymax = rep(y_min, 2L),
      fill = c("white", "white"),
      stringsAsFactors = FALSE
    )
    intron_y <- y_min - exon_height / 2
    intron_df <- data.frame(
      x    = bs[1L],
      xend = bs[2L],
      y    = intron_y,
      yend = intron_y
    )
    list(
      ggplot2::geom_rect(
        data = exon_df,
        ggplot2::aes(xmin = xmin, xmax = xmax,
                     ymin = ymin, ymax = ymax, fill = fill),
        color = "black", linewidth = 0.5, inherit.aes = FALSE
      ),
      ggplot2::geom_segment(
        data = intron_df,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        color = "black", linewidth = 1.5, inherit.aes = FALSE
      )
    )
  }
)

#' @rdname event_schema
#' @noRd
event_schema_a5ss <- list(
  name = "Alternative 5' Splice Site",

  required_cols = c(
    .event_required_base,
    "longExonStart_0base", "longExonEnd",
    "shortES", "shortEE",
    "flankingES", "flankingEE"
  ),

  n_regions     = 2L,

  region_width  = function(width_exon, width_intron) {
    as.integer(width_exon + width_intron + 1L)
  },

  #How to make regions for specific schema
  build_regions = function(events, width_exon, width_intron) {
    n <- nrow(events)
    if (n == 0L) {
      return(GenomicRanges::GRanges())
    }

    plus <- events$strand == "+"

    # Anchor on the short-isoform donor; the long donor falls inside the
    # window on the intron side. The flanking edge faces the intron.
    alt_prox   <- ifelse(plus, events$shortEE, events$shortES)
    flank_edge <- ifelse(plus, events$flankingES, events$flankingEE)

    # Genomic left/right order (the scorer reorients - strand). On + the alt
    # exon is left of the flanking exon; on - it is reversed.
    left_anchor  <- ifelse(plus, alt_prox,   flank_edge)
    right_anchor <- ifelse(plus, flank_edge, alt_prox)

    # Far exon edges to clamp against; the alt exon's far edge is its shared
    # splice site (+ shares the start, - shares the end).
    left_far  <- ifelse(plus, events$shortES,    events$flankingES)
    right_far <- ifelse(plus, events$flankingEE, events$shortEE)

    s1 <- pmax(left_anchor  - width_exon,   left_far)
    e1 <- pmin(left_anchor  + width_intron, right_anchor)
    s2 <- pmax(right_anchor - width_intron, left_anchor)
    e2 <- pmin(right_anchor + width_exon,   right_far)

    starts <- cbind(s1, s2)
    ends   <- cbind(e1, e2)

    GenomicRanges::GRanges(
      seqnames   = rep(events$chr, each = 2L),
      ranges     = IRanges::IRanges(start = as.vector(t(starts)),
                                    end   = as.vector(t(ends))),
      strand     = rep(events$strand, each = 2L),
      event_id   = rep(seq_len(n), each = 2L),
      region_idx = rep(seq_len(2L), times = n)
    )
  },

  region_labels = c(
    "Alternative 5' splice site",
    "Flanking exon"
  ),

  schematic = "a5ss",

  schematic_stat = function(events) {
    ext <- (events$longExonEnd - events$longExonStart_0base) -
           (events$shortEE - events$shortES)
    list(median_ext = stats::median(ext, na.rm = TRUE))
  },

  group_set     = c("Negative", "Positive", "Control"),

  # x-axis layout for the plotter.
  plot_layout = function(width_exon, width_intron) {
    bin_width <- as.integer(width_exon + width_intron)
    gap       <- 80L
    region_starts <- c(0L, bin_width + gap)
    list(
      region_starts    = region_starts,
      bin_width        = bin_width,
      gap              = gap,
      boundary_offsets = c(width_exon, width_intron),
      x_max            = region_starts[2L] + bin_width
    )
  },

  # Schematic layers for the splicing event plot
  build_schematic_layers = function(layout, style, y_min, exon_height,
                                    schematic_data = NULL) {
    bs <- layout$region_starts + layout$boundary_offsets

    ext <- .schematic_ext(schematic_data, layout)

    # R1 = alternative exon (solid = short isoform), R2 = flanking exon.
    exon_df <- data.frame(
      xmin = c(layout$region_starts[1L], bs[2L]),
      xmax = c(bs[1L], layout$region_starts[2L] + layout$bin_width),
      ymin = rep(y_min - exon_height, 2L),
      ymax = rep(y_min, 2L),
      fill = c(style$exon_col, "white"),
      stringsAsFactors = FALSE
    )

    extension_df <- data.frame(
      xmin = bs[1L], xmax = bs[1L] + ext,
      ymin = y_min - exon_height, ymax = y_min
    )

    intron_y  <- y_min - exon_height / 2
    intron_df <- data.frame(x = bs[1L] + ext, xend = bs[2L],
                            y = intron_y, yend = intron_y)

    label_df <- data.frame(
      x     = c(bs[1L], bs[1L] + ext) + style$isoform_label_nudge_x,
      label = c("Short isoform", "Long isoform"),
      stringsAsFactors = FALSE
    )
    label_y <- y_min - exon_height * 1.45 + style$isoform_label_nudge_y

    list(
      ggplot2::geom_rect(
        data = exon_df,
        ggplot2::aes(xmin = xmin, xmax = xmax,
                     ymin = ymin, ymax = ymax, fill = fill),
        color = "black", linewidth = 0.5, inherit.aes = FALSE
      ),
      ggplot2::geom_rect(
        data = extension_df,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "white", color = "black", linetype = "dotted",
        linewidth = 0.6, inherit.aes = FALSE
      ),
      ggplot2::geom_segment(
        data = intron_df,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        color = "black", linewidth = 1.5, inherit.aes = FALSE
      ),
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = x, y = label_y, label = label),
        hjust = c(1, 0), vjust = 1, size = style$isoform_label_size,
        fontface = "italic", inherit.aes = FALSE
      )
    )
  }
)


#' @rdname event_schema
#' @noRd
event_schema_a3ss <- list(
  name = "Alternative 3' Splice Site",

  required_cols = c(
    .event_required_base,
    "longExonStart_0base", "longExonEnd",
    "shortES", "shortEE",
    "flankingES", "flankingEE"
  ),

  n_regions     = 2L,

  region_width  = function(width_exon, width_intron) {
    as.integer(width_exon + width_intron + 1L)
  },

  #How to make regions for specific schema
  build_regions = function(events, width_exon, width_intron) {
    n <- nrow(events)
    if (n == 0L) {
      return(GenomicRanges::GRanges())
    }

    plus <- events$strand == "+"

    # Anchor on the short-isoform acceptor; the long acceptor falls inside the
    # window on the intron side. The flanking edge faces the intron.
    alt_prox   <- ifelse(plus, events$shortES, events$shortEE)
    flank_edge <- ifelse(plus, events$flankingEE, events$flankingES)

    # Genomic left/right order (the scorer reorients - strand). On + the
    # flanking exon is left of the alt exon; on - it is reversed.
    left_anchor  <- ifelse(plus, flank_edge, alt_prox)
    right_anchor <- ifelse(plus, alt_prox,   flank_edge)

    # Far exon edges to clamp against; the alt exon's far edge is its shared
    # splice site (+ shares the end, - shares the start).
    left_far  <- ifelse(plus, events$flankingES, events$shortES)
    right_far <- ifelse(plus, events$shortEE,    events$flankingEE)

    s1 <- pmax(left_anchor  - width_exon,   left_far)
    e1 <- pmin(left_anchor  + width_intron, right_anchor)
    s2 <- pmax(right_anchor - width_intron, left_anchor)
    e2 <- pmin(right_anchor + width_exon,   right_far)

    starts <- cbind(s1, s2)
    ends   <- cbind(e1, e2)

    GenomicRanges::GRanges(
      seqnames   = rep(events$chr, each = 2L),
      ranges     = IRanges::IRanges(start = as.vector(t(starts)),
                                    end   = as.vector(t(ends))),
      strand     = rep(events$strand, each = 2L),
      event_id   = rep(seq_len(n), each = 2L),
      region_idx = rep(seq_len(2L), times = n)
    )
  },

  region_labels = c(
    "Flanking exon",
    "Alternative 3' splice site"
  ),

  schematic = "a3ss",

  schematic_stat = function(events) {
    ext <- (events$longExonEnd - events$longExonStart_0base) -
           (events$shortEE - events$shortES)
    list(median_ext = stats::median(ext, na.rm = TRUE))
  },

  group_set     = c("Negative", "Positive", "Control"),

  # x-axis layout for the plotter.
  plot_layout = function(width_exon, width_intron) {
    bin_width <- as.integer(width_exon + width_intron)
    gap       <- 80L
    region_starts <- c(0L, bin_width + gap)
    list(
      region_starts    = region_starts,
      bin_width        = bin_width,
      gap              = gap,
      boundary_offsets = c(width_exon, width_intron),
      x_max            = region_starts[2L] + bin_width
    )
  },

  # Schematic layers for the splicing event plot
  build_schematic_layers = function(layout, style, y_min, exon_height,
                                    schematic_data = NULL) {
    bs <- layout$region_starts + layout$boundary_offsets

    ext <- .schematic_ext(schematic_data, layout)

    # R1 = flanking exon, R2 = alternative exon (solid = short isoform).
    exon_df <- data.frame(
      xmin = c(layout$region_starts[1L], bs[2L]),
      xmax = c(bs[1L], layout$region_starts[2L] + layout$bin_width),
      ymin = rep(y_min - exon_height, 2L),
      ymax = rep(y_min, 2L),
      fill = c("white", style$exon_col),
      stringsAsFactors = FALSE
    )

    extension_df <- data.frame(
      xmin = bs[2L] - ext, xmax = bs[2L],
      ymin = y_min - exon_height, ymax = y_min
    )

    intron_y  <- y_min - exon_height / 2
    intron_df <- data.frame(x = bs[1L], xend = bs[2L] - ext,
                            y = intron_y, yend = intron_y)

    label_df <- data.frame(
      x     = c(bs[2L] - ext, bs[2L]) + style$isoform_label_nudge_x,
      label = c("Long isoform", "Short isoform"),
      stringsAsFactors = FALSE
    )
    label_y <- y_min - exon_height * 1.45 + style$isoform_label_nudge_y

    list(
      ggplot2::geom_rect(
        data = exon_df,
        ggplot2::aes(xmin = xmin, xmax = xmax,
                     ymin = ymin, ymax = ymax, fill = fill),
        color = "black", linewidth = 0.5, inherit.aes = FALSE
      ),
      ggplot2::geom_rect(
        data = extension_df,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "white", color = "black", linetype = "dotted",
        linewidth = 0.6, inherit.aes = FALSE
      ),
      ggplot2::geom_segment(
        data = intron_df,
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        color = "black", linewidth = 1.5, inherit.aes = FALSE
      ),
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = x, y = label_y, label = label),
        hjust = c(1, 0), vjust = 1, size = style$isoform_label_size,
        fontface = "italic", inherit.aes = FALSE
      )
    )
  }
)
















