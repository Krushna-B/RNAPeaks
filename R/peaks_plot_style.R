#' Default styling for Peak Plots
#' Returns a named list of every visual setting used by [draw_plot()].
#'
#' @param exon_color,exon_height Exon/CDS rectangles.
#' @param utr_color,utr_height UTR rectangles.
#' @param intron_color,intron_linewidth,intron_arrow_len_in Intron line and arrow geometry.
#' @param total_arrows,max_per_intron Density controls for intron arrows.
#' @param peak_color,peak_height,peak_alpha,peak_border_color,peak_border_linewidth Peak rectangles.
#' @param band_even_fill,band_odd_fill,band_sep_color,band_sep_linewidth
#'   Alternating background bands behind protein rows.
#' @param protein_label_size,protein_label_color,protein_label_x_offset_bp Protein label appearance.
#' @param strand_label_size,strand_label_color 5'/3' tag appearance.
#' @param gene_label_size,gene_label_color,gene_label_x_offset
#'   Per-gene labels (region plots).
#' @param title_size,title_color,subtitle_size,subtitle_color,subtitle_sep
#'   Plot title / subtitle.
#' @param axis_title_size,axis_text_size,axis_pad_bp,axis_breaks_n Axis.
#' @param five_to_three If `TRUE` and the region's strand is `"-"`, the
#'   x-axis is flipped so 5' is on the left. Ignored on `"+"` strand.
#' @param plot_top_margin,plot_right_margin,plot_bottom_margin,plot_left_margin
#'   Plot margins in points; `plot_left_margin = NULL` auto-fits to the
#'   longest protein label.
#' @param bam_fill_color,bam_fill_alpha BAM coverage track fill.
#' @param bam_label_size,bam_axis_text_size BAM track label and y-axis text.
#' @param bam_ylim Optional `c(min, max)` shared by all BAM tracks; `NULL`
#'   derives a common scale from the maximum coverage across tracks.
#' @param bam_track_height Relative height of each BAM panel (gene panel = 4).
#' @param highlight Optional `c(start, stop)` for an overlay band.
#' @param highlight_color,highlight_opacity Highlight band appearance.
#' @param show_junctions,junction_color,junction_linetype,junction_linewidth,junction_alpha
#'   Optional dashed vertical lines at exon/UTR boundaries.
#'
#' @return A named list of styling parameters.
#' @export
peaks_plot_style <- function(
  # region features
  exon_color             = "navy",
  exon_height            = 0.5,
  utr_color              = "lightgray",
  utr_height             = 0.3,
  intron_color           = "gray60",
  intron_linewidth       = 0.9,
  intron_arrow_len_in    = 0.15,
  total_arrows           = 6,
  max_per_intron         = 2,

  # peaks
  peak_color             = "purple",
  peak_height            = 0.3,
  peak_alpha             = 0.95,
  peak_border_color      = NA,
  peak_border_linewidth  = 0.4,

  # background bands
  band_even_fill         = "#F7F8FA",
  band_odd_fill          = "#FFFFFF",
  band_sep_color         = "#E5E7EB",
  band_sep_linewidth     = 0.4,

  # protein labels
  protein_label_size        = 5,
  protein_label_color       = "black",
  protein_label_x_offset_bp = 100,
  strand_label_size      = 5,
  strand_label_color     = "black",
  gene_label_size        = 5,
  gene_label_color       = "black",
  gene_label_x_offset    = 0.25,

  # title and subtitle
  title_size             = 25,
  title_color            = "black",
  subtitle_size          = 12,
  subtitle_color         = "black",
  subtitle_sep           = ": ",

  # Axis
  axis_title_size        = 11,
  axis_text_size         = 9,
  axis_pad_bp            = 500,
  axis_breaks_n          = 5,
  five_to_three          = FALSE,

  # plot margins
  plot_top_margin        = 30,
  plot_right_margin      = 50,
  plot_bottom_margin     = 30,
  plot_left_margin       = NULL,

  # BAM Coverage Tracks
  bam_fill_color         = "navy",
  bam_fill_alpha         = 0.75,
  bam_label_size         = 9,
  bam_axis_text_size     = 8,
  bam_ylim               = NULL,
  bam_track_height       = 1,

  # Highlighting Region
  highlight              = NULL,
  highlight_color        = "pink",
  highlight_opacity      = 0.30,

  # Junction Lines
  show_junctions         = FALSE,
  junction_color         = "gray40",
  junction_linetype      = "dashed",
  junction_linewidth     = 0.4,
  junction_alpha         = 0.7
) {
  # region features
  check_color(exon_color, "exon_color")
  check_scalar_number(exon_height, "exon_height", min = 0)
  check_color(utr_color, "utr_color")
  check_scalar_number(utr_height, "utr_height", min = 0)
  check_color(intron_color, "intron_color")
  check_scalar_number(intron_linewidth, "intron_linewidth", min = 0)
  check_scalar_number(intron_arrow_len_in, "intron_arrow_len_in", min = 0)
  check_scalar_int(total_arrows, "total_arrows", min = 0)
  check_scalar_int(max_per_intron, "max_per_intron", min = 0)

  # peaks
  check_color(peak_color, "peak_color")
  check_scalar_number(peak_height, "peak_height", min = 0)
  check_unit_interval(peak_alpha, "peak_alpha")
  check_color(peak_border_color, "peak_border_color")
  check_scalar_number(peak_border_linewidth, "peak_border_linewidth", min = 0)

  # background bands
  check_color(band_even_fill, "band_even_fill")
  check_color(band_odd_fill,  "band_odd_fill")
  check_color(band_sep_color, "band_sep_color")
  check_scalar_number(band_sep_linewidth, "band_sep_linewidth", min = 0)

  # protein labels
  check_scalar_number(protein_label_size, "protein_label_size", min = 0)
  check_color(protein_label_color, "protein_label_color")
  check_scalar_number(protein_label_x_offset_bp, "protein_label_x_offset_bp", min = 0)
  check_scalar_number(strand_label_size, "strand_label_size", min = 0)
  check_color(strand_label_color, "strand_label_color")
  check_scalar_number(gene_label_size, "gene_label_size", min = 0)
  check_color(gene_label_color, "gene_label_color")
  check_scalar_number(gene_label_x_offset, "gene_label_x_offset", min = 0)

  # title and subtitle
  check_scalar_number(title_size, "title_size", min = 0)
  check_color(title_color, "title_color")
  check_scalar_number(subtitle_size, "subtitle_size", min = 0)
  check_color(subtitle_color, "subtitle_color")
  check_string(subtitle_sep, "subtitle_sep")

  # Axis
  check_scalar_number(axis_title_size, "axis_title_size", min = 0)
  check_scalar_number(axis_text_size, "axis_text_size", min = 0)
  check_scalar_number(axis_pad_bp, "axis_pad_bp", min = 0)
  check_scalar_int(axis_breaks_n, "axis_breaks_n", min = 2)
  check_flag(five_to_three, "five_to_three")

  # plot margins
  check_scalar_number(plot_top_margin,    "plot_top_margin",    min = 0)
  check_scalar_number(plot_right_margin,  "plot_right_margin",  min = 0)
  check_scalar_number(plot_bottom_margin, "plot_bottom_margin", min = 0)
  check_scalar_number_or_null(plot_left_margin, "plot_left_margin", min = 0)

  # BAM Coverage Tracks
  check_color(bam_fill_color, "bam_fill_color")
  check_unit_interval(bam_fill_alpha, "bam_fill_alpha")
  check_scalar_number(bam_label_size, "bam_label_size", min = 0)
  check_scalar_number(bam_axis_text_size, "bam_axis_text_size", min = 0)
  check_range_or_null(bam_ylim, "bam_ylim")
  check_scalar_number(bam_track_height, "bam_track_height", min = 0)

  # Highlighting Region
  check_range_or_null(highlight, "highlight")
  check_color(highlight_color, "highlight_color")
  check_unit_interval(highlight_opacity, "highlight_opacity")

  # Junction Lines
  check_flag(show_junctions, "show_junctions")
  check_color(junction_color, "junction_color")
  check_string(junction_linetype, "junction_linetype",
               choices = c("solid", "dashed", "dotted", "dotdash",
                           "longdash", "twodash", "blank"))
  check_scalar_number(junction_linewidth, "junction_linewidth", min = 0)
  check_unit_interval(junction_alpha, "junction_alpha")

  as.list(environment())
}
