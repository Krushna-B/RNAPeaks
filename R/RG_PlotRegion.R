#' Plot RNA-Binding Protein Peaks on a Genomic Region
#'
#' Creates a publication-quality visualization of RNA-binding protein peaks
#' overlaid on multiple gene structures within a specified genomic region.
#' Unlike `PlotGene()`, this function can display multiple genes that fall
#' within the specified coordinates.
#'
#' @param Chr Chromosome name (e.g., "1", "X", "chr1").
#' @param Start Start position of the genomic region (bp).
#' @param End End position of the genomic region (bp).
#' @param Strand Strand to display ("+" or "-").
#' @param geneID Optional gene identifier to focus on a specific gene.
#' @param gtf Optional pre-loaded GTF annotation data frame. If NULL, annotations
#'   are loaded from AnnotationHub based on species.
#' @param bed A data frame containing BED-format peak data.
#' @param species Species for annotation lookup. Either "Human" or "Mouse".
#' @param TxID Optional transcript ID to plot a specific transcript isoform.
#' @param Target_col Column name in bed containing the protein/target identifiers.
#' @param omit Character vector of target names to exclude from the plot.
#' @param order_by Method for ordering protein tracks: "Count" (default),
#'   "Target", or "Region".
#' @param order_in Optional character vector specifying exact order of targets.
#' @param merge Minimum gap width (bp) for merging nearby peaks.
#' @param peaks_width Vertical height of each peak track row.
#' @param utr_col Color for UTR regions.
#' @param peak_col Color for peak rectangles.
#' @param exon_width Vertical height of exon rectangles.
#' @param utr_width Vertical height of UTR rectangles.
#' @param exon_col Color for exon/CDS regions.
#' @param RNA_Peaks_File_Path File path to save the output PDF plot.
#' @param Bed_File_Path File path to save the filtered BED data as CSV.
#' @param ... Additional styling arguments passed to internal plotting functions.
#'   See Styling Parameters section below.
#'
#' @section Styling Parameters:
#' The following parameters can be passed via `...` to customize the plot appearance:
#'
#' \strong{Gene Structure Colors:}
#' \describe{
#'   \item{exon_fill}{Fill color for exon/CDS regions. Default: "navy"}
#'   \item{utr_fill}{Fill color for UTR regions. Default: "lightgray"}
#'   \item{intron_color}{Color for intron lines. Default: "gray60"}
#'   \item{intron_linewidth}{Line width for introns. Default: 0.9}
#'   \item{intron_arrow_len_in}{Length of intron direction arrows in inches. Default: 0.15}
#' }
#'
#' \strong{Peak Styling:}
#' \describe{
#'   \item{peak_alpha}{Opacity of peak rectangles. Default: 0.95}
#'   \item{peak_border_color}{Border color for peaks. Default: NA (no border)}
#'   \item{peak_border_linewidth}{Border line width for peaks. Default: 0.4}
#' }
#'
#' \strong{Background Bands:}
#' \describe{
#'   \item{band_even_fill}{Fill color for even-numbered protein track bands. Default: "#F7F8FA"}
#'   \item{band_odd_fill}{Fill color for odd-numbered protein track bands. Default: "#FFFFFF"}
#'   \item{band_sep_color}{Color for band separator lines. Default: "#E5E7EB"}
#'   \item{band_sep_linewidth}{Line width for band separators. Default: 0.4}
#' }
#'
#' \strong{Labels:}
#' \describe{
#'   \item{label_size}{Font size for protein labels. Default: 5}
#'   \item{label_color}{Color for protein labels. Default: "black"}
#'   \item{strand_label_size}{Font size for 5'/3' strand labels. Default: 5}
#'   \item{strand_label_color}{Color for strand labels. Default: "black"}
#'   \item{protein_label_x_offset}{Horizontal offset for protein labels in bp. Default: 100}
#' }
#'
#' \strong{Title and Axes:}
#' \describe{
#'   \item{title_size}{Font size for plot title. Default: 25}
#'   \item{title_color}{Color for plot title. Default: "black"}
#'   \item{subtitle_size}{Font size for plot subtitle. Default: 12}
#'   \item{subtitle_color}{Color for plot subtitle. Default: "black"}
#'   \item{subtitle_sep}{Separator between gene name and coordinates in subtitle. Default: ": "}
#'   \item{axis_title_size}{Font size for axis titles. Default: 11}
#'   \item{axis_text_size}{Font size for axis tick labels. Default: 9}
#' }
#'
#' \strong{Axis and Layout:}
#' \describe{
#'   \item{x_lims}{Custom x-axis limits as c(min, max). Default: NULL (auto)}
#'   \item{axis_pad_bp}{Padding in base pairs added to each side of the plot. Default: 500}
#'   \item{axis_breaks_n}{Number of axis tick breaks. Default: 5}
#'   \item{max_proteins}{Maximum number of protein tracks to display. Default: 40}
#' }
#'
#' \strong{Plot Margins:}
#' \describe{
#'   \item{plot_right_margin}{Right margin in points. Default: 50}
#'   \item{plot_top_margin}{Top margin in points. Default: 30}
#'   \item{plot_bottom_margin}{Bottom margin in points. Default: 30}
#'   \item{plot_left_margin}{Left margin in points. Default: NULL (auto)}
#' }
#'
#' \strong{Region Plot Specific:}
#' \describe{
#'   \item{gene_label_x_offset}{Horizontal offset for gene labels as fraction of axis_pad_bp. Default: 0.25}
#'   \item{gene_label_size}{Font size for gene name labels. Default: 5}
#'   \item{gene_label_color}{Color for gene name labels. Default: "black"}
#' }
#'
#' \strong{Highlighted Region:}
#' \describe{
#'   \item{highlighted_region_start}{Start position of region to highlight. Default: NULL}
#'   \item{highlighted_region_stop}{End position of region to highlight. Default: NULL}
#'   \item{highlighted_region_color}{Color for highlighted region. Default: "pink"}
#'   \item{highlighted_region_opacity}{Opacity for highlighted region. Default: 0.30}
#' }
#'
#' @return A named list containing:
#' \describe{
#'   \item{plot}{A ggplot2 object of the peak visualization}
#'   \item{csv}{The filtered BED data frame used for plotting}
#' }
#' Access with \code{result$plot} and \code{result$csv}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Load GTF annotation (do this once, takes time on first call)
#'   gtf <- LoadGTF(species = "Human")
#'
#'   # ----- Using included sample data -----
#'   # sample_bed is included with the package and ready to use
#'   result <- PlotRegion(
#'     bed = sample_bed,
#'     gtf = gtf,
#'     Chr = "12",
#'     Start = 56000000,
#'     End = 56050000,
#'     Strand = "+"
#'   )
#'
#'   # Access results
#'   result$plot
#'   result$csv
#'
#'   # ----- Using your own BED file -----
#'   # 1. Read your BED file
#'   my_bed <- read.table("my_peaks.bed", header = FALSE, sep = "\t")
#'
#'   # 2. Check Bed file
#'   my_bed <- checkBed(my_bed)
#'
#'   # 3. Plot peaks on region
#'   result <- PlotRegion(
#'     bed = my_bed,
#'     gtf = gtf,
#'     Chr = "12",
#'     Start = 56000000,
#'     End = 56050000,
#'     Strand = "+"
#'   )
#'
#'   # Access results
#'   result$plot
#'   result$csv
#' }
PlotRegion <- function(Chr = NULL,
                       Start = NULL,
                       End = NULL,
                       Strand = NULL,
                       geneID = NULL,
                       gtf = NULL,
                       bed = NULL,
                       species = "Human",
                       TxID = NA,
                       Target_col = NULL,
                       omit = c(),
                       order_by = "Count",
                       order_in = NULL,
                       merge = 0,
                       peaks_width = 0.3,
                       utr_col = "dark gray",
                       peak_col = "Blue",
                       exon_width = 0.5,
                       utr_width = 0.3,
                       exon_col = "black",
                       RNA_Peaks_File_Path = "~/Desktop/RNAPeaks.pdf",
                       Bed_File_Path = "~/Desktop/BEDFILE_PEAKS.csv",
                       ...) {

  if (is.null(Chr) | is.null(Start) | is.null(End) | is.null(Strand)) {
    stop("Need to provide Chr, Start, Stop and Strand parameters to visualize region.")
  }

  clipped_region <- Build_Region_Structure(
    gtf = gtf,
    Chr = Chr,
    Start = Start,
    End = End,
    Strand = Strand,
    geneID = geneID,
    TxID = TxID
  )
  clipped_region <- clipped_region[order(clipped_region$start, clipped_region$end), ]

  # Check the bed file
  colnames(bed)[which(colnames(bed) == Target_col)] <- "target"
  bed <- checkBed(bed)

  # Filter bed
  bed <- FilterBed(
    bed = bed,
    chr = clipped_region$seqnames[1],
    start = min(clipped_region$start),
    end = max(clipped_region$end),
    strand = clipped_region$strand[1],
    omit = omit,
    collapse = merge
  )

  # Option to give an order to rank proteins in
  if (!is.null(order_in)) {
    Target_rank <- order_in
  } else {
    Target_rank <- OrderPeak(bed = bed, order_by = order_by)
  }

  Plot <- Get_Multi_Plot_by_Region(
    gtf = gtf,
    bed = bed,
    Chr = Chr,
    Start = Start,
    End = End,
    region_with_multiple_genes = clipped_region,
    Strand = Strand,
    geneID = geneID,
    rank_ = Target_rank,
    gene = clipped_region,
    TxID = TxID,
    peaks_width = peaks_width,
    exon_width = exon_width,
    utr_width = utr_width,
    exon_col = exon_col,
    utr_col = utr_col,
    ...
  )

  # Save to Desktop
  ggplot2::ggsave(RNA_Peaks_File_Path, Plot, height = 12, width = 16)
  utils::write.csv(bed, Bed_File_Path, row.names = FALSE)

  Plot_and_Peaks <- list(plot = Plot, csv = bed)
  return(Plot_and_Peaks)
}
