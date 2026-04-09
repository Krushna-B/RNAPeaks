#' Plot RNA-Binding Protein Peaks on a Genomic Region
#'
#' Creates a quality visualization of RNA-binding protein peaks
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
#' @param total_arrows Total number of directional arrows drawn across all
#'   introns to indicate transcription direction. Default is 12.
#' @param max_per_intron Maximum number of directional arrows drawn per intron.
#'   Default is 5.
#' @param five_to_three Logical. If TRUE and Strand is "-", flips the x-axis so
#'   5' is on the left. Default FALSE.
#' @param bam_files Optional. A named character vector of BAM file paths to
#'   display as coverage tracks above the gene structure. Names are used as
#'   track labels on the left-hand side of each panel. If unnamed, the filename
#'   (without extension) is used as the label. BAM files must be sorted and
#'   indexed (a `.bai` file must exist alongside each BAM).
#'   Example: \code{c("Sample A" = "/path/to/a.bam", "Sample B" = "/path/to/b.bam")}
#' @param bam_fill_col Fill color for BAM coverage tracks. A single color
#'   applied to all tracks, or a character vector the same length as
#'   \code{bam_files} for per-track colours. Default \code{"navy"}.
#' @param bam_fill_alpha Opacity of BAM track fill. Default \code{0.75}.
#' @param bam_label_size Font size of the BAM track name label on the left. Default \code{9}.
#' @param bam_axis_text_size Font size of the 0 and max coverage values. Default \code{8}.
#' @param bam_ylim Optional global y-axis limits \code{c(min, max)} applied to
#'   all BAM tracks. If \code{NULL}, all tracks share a common scale derived
#'   from the maximum coverage across all BAMs.
#' @param bam_track_height Relative height of each BAM panel compared to the
#'   gene plot panel (which is always 4 units). Default \code{1}.
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
#' \strong{Junction Lines:}
#' \describe{
#'   \item{show_junctions}{Logical. If TRUE, draws vertical dashed lines at exon/intron boundaries. Default: FALSE}
#'   \item{junction_color}{Color for junction lines. Default: "gray40"}
#'   \item{junction_linetype}{Line type for junction lines. Default: "dashed"}
#'   \item{junction_linewidth}{Line width for junction lines. Default: 0.4}
#'   \item{junction_alpha}{Opacity for junction lines. Default: 0.7}
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
                       utr_col = "darkgray",
                       peak_col = "Blue",
                       exon_width = 0.5,
                       utr_width = 0.3,
                       exon_col = "navy",
                       total_arrows = 12,
                       max_per_intron = 5,
                       five_to_three = FALSE,
                       bam_files = NULL,
                       bam_fill_col = "navy",
                       bam_fill_alpha = 0.75,
                       bam_label_size = 9,
                       bam_axis_text_size = 8,
                       bam_ylim = NULL,
                       bam_track_height = 1,
                       RNA_Peaks_File_Path = "~/Desktop/RNAPeaks.pdf",
                       Bed_File_Path = "~/Desktop/BEDFILE_PEAKS.csv",
                       ...) {

  #Ensure params are passed
  if (is.null(Chr) | is.null(Start) | is.null(End) | is.null(Strand)) {
    stop("Need to provide Chr, Start, Stop and Strand parameters to visualize region.")
  }

  #Build region area
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
    peak_col = peak_col,
    total_arrows = total_arrows,
    max_per_intron = max_per_intron,
    five_to_three = five_to_three,
    ...
  )

  #BAM coverage tracks
  if (!is.null(bam_files)) {

    # Ensure names exist; fall back to filename without extension
    if (is.null(names(bam_files))) {
      names(bam_files) <- vapply(bam_files, BAM_Label_From_Path, character(1))
    } else {
      missing_names <- names(bam_files) == "" | is.na(names(bam_files))
      names(bam_files)[missing_names] <- vapply(
        bam_files[missing_names], BAM_Label_From_Path, character(1)
      )
    }

    n_bam     <- length(bam_files)
    fill_cols <- rep_len(bam_fill_col, n_bam)

    # axis_pad_bp may be overridden via ...; use same default as Draw_Gene_Plot
    dots        <- list(...)
    axis_pad_bp <- if ("axis_pad_bp" %in% names(dots)) dots$axis_pad_bp else 500

    # Pre-compute coverage for all BAMs so global limits can be derived
    cov_list <- vector("list", n_bam)
    for (i in seq_len(n_bam)) {
      cov_list[[i]] <- Compute_BAM_Coverage(
        bam_path = bam_files[[i]],
        chr      = as.character(Chr),
        start    = Start,
        end      = End
      )
    }

    # Global y limits: user-supplied or derived from all tracks together
    if (!is.null(bam_ylim)) {
      resolved_ylim <- bam_ylim
    } else {
      global_max <- max(vapply(cov_list, function(d) max(d$coverage, na.rm = TRUE), numeric(1)))
      resolved_ylim <- c(0, if (global_max == 0) 1 else global_max)
    }

    bam_track_plots <- vector("list", n_bam)
    for (i in seq_len(n_bam)) {
      bam_track_plots[[i]] <- Draw_BAM_Track(
        cov_df         = cov_list[[i]],
        label          = names(bam_files)[i],
        x_min          = Start,
        x_max          = End,
        axis_pad_bp    = axis_pad_bp,
        fill_col       = fill_cols[i],
        fill_alpha     = bam_fill_alpha,
        label_size     = bam_label_size,
        axis_text_size = bam_axis_text_size,
        ylim           = resolved_ylim,
        five_to_three  = five_to_three,
        gene_strand    = Strand
      )
    }

    # Lift title from gene plot and move to top of combined figure
    plot_title <- Plot$labels$title
    plot_sub   <- Plot$labels$subtitle

    existing_m <- Plot$theme$plot.margin
    Plot <- Plot + ggplot2::theme(
      plot.title    = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.margin   = ggplot2::margin(
        t = 0,
        r = grid::convertUnit(existing_m[2], "pt", valueOnly = TRUE),
        b = grid::convertUnit(existing_m[3], "pt", valueOnly = TRUE),
        l = grid::convertUnit(existing_m[4], "pt", valueOnly = TRUE)
      )
    )

    all_panels <- c(bam_track_plots, list(Plot))
    heights    <- c(rep(bam_track_height, n_bam), 4)
    Plot       <- patchwork::wrap_plots(all_panels, ncol = 1, heights = heights) +
      patchwork::plot_annotation(
        title    = plot_title,
        subtitle = plot_sub,
        theme    = ggplot2::theme(
          plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold.italic",
                                                size = 25, color = "black"),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                size = 12, color = "black",
                                                margin = ggplot2::margin(t = 2, b = 8))
        )
      )
  }

  # Save files if paths are provided
  if (!is.null(RNA_Peaks_File_Path)) {
    ggplot2::ggsave(RNA_Peaks_File_Path, Plot, height = 12, width = 16)
  }
  if (!is.null(Bed_File_Path)) {
    utils::write.csv(bed, Bed_File_Path, row.names = FALSE)
  }

  Plot_and_Peaks <- list(plot = Plot, csv = bed)
  return(Plot_and_Peaks)
}
