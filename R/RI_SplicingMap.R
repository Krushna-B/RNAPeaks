#' Create Retained Intron Splicing Map
#'
#' Analyzes protein binding frequency across retained intron junction regions.
#' Uses a 2-region structure to show where protein binding sites appear
#' relative to the upstream exon/intron and intron/downstream exon boundaries.
#' Filters events into Retained, Excluded, and Control groups.
#'
#' @param bed_file Either a file path to a BED file or a data frame containing
#'   BED data with columns: chr, start, end, tag, score, strand
#' @param RIMATS  A data frame containing retained intron outputs.
#' @param moving_average Integer specifying the window size for moving average
#'   smoothing. Set to NULL or 0 to disable smoothing. Default is 50.
#' @param WidthIntoExon Integer specifying how many bp to extend into exons.
#'   Default is 50.
#' @param WidthIntoIntron Integer specifying how many bp to extend into introns.
#'   Default is 300.
#' @param p_valueRetainedAndExclusion P-value threshold for retained/excluded events.
#'   Default is 0.05.
#' @param p_valueControls P-value threshold for control events. Default is 0.95.
#' @param retained_IncLevelDifference Inclusion level difference threshold for
#'   retained events. Default is 0.1.
#' @param exclusion_IncLevelDifference Inclusion level difference threshold for
#'   excluded events. Default is -0.1.
#' @param Min_Count Minimum read count threshold. Default is 50. Set to 0 or NULL to skip filtering.
#' @param read_count_cols Column names for junction read counts, in order
#'   c(IJC_s1, SJC_s1, IJC_s2, SJC_s2). Defaults to the standard rMATS names.
#'   Only used when \code{Min_Count > 0}.
#' @param groups Character vector specifying which event groups to process.
#'   Options are "Retained", "Excluded", and/or "Control". Default is
#'   c("Retained", "Excluded", "Control") to process all groups.
#' @param control_multiplier Numeric multiplier for control sample size. The
#'   number of control events sampled per iteration is
#'   (n_retained + n_excluded) * control_multiplier. Default is 2.0.
#' @param control_iterations Integer number for sampling iterations for control
#'   sampling. The final control frequency is the mean across iterations, with
#'   standard deviation shown as a shaded band. Default is 20.
#' @param z_threshold Z-score threshold for significance testing. Default is 1.96.
#'   Only used when use_fdr = FALSE.
#' @param min_consecutive Minimum number of consecutive significant positions
#'   required to form a significant region. Default is 10.
#' @param one_sided Logical. If TRUE (default), only test for enrichment.
#' @param use_fdr Logical. If TRUE, use FDR-corrected p-values. Default is TRUE.
#' @param fdr_threshold FDR threshold for significance when use_fdr = TRUE.
#'   Default is 0.05.
#' @param show_significance Logical. If TRUE (default), displays colored bars above
#'   the plot indicating regions where Retained/Excluded differ significantly
#'   from Control based on z-test.
#' @param return_data Logical. If TRUE, returns the frequency data frame instead
#'   of a plot. Default is FALSE.
#' @param return_diagnostics Logical. If TRUE, returns a list containing the
#'   frequency data and bootstrap diagnostics. Default is FALSE.
#' @param verbose Logical. If TRUE, prints progress messages. Default is TRUE.
#' @param progress_callback Optional function to report progress. Default is NULL.
#' @param title Character string for the plot title. Default is "".
#' @param retained_col Color for the Retained group line. Default is "blue".
#' @param excluded_col Color for the Excluded group line. Default is "red".
#' @param control_col Color for the Control group line. Default is "black".
#' @param line_width Numeric line width for the frequency lines. Default is 0.8.
#' @param line_alpha Numeric alpha for the frequency lines. Default is 0.7.
#' @param ribbon_alpha Numeric alpha for the SD ribbon around Control. Default is 0.3.
#' @param title_size Numeric font size for the plot title. Default is 20.
#' @param title_color Color for the plot title text. Default is "black".
#' @param axis_text_size Numeric font size for y-axis tick labels. Default is 11.
#' @param boundary_col Color for the dashed vertical boundary lines. Default is "gray70".
#' @param exon_col Unused parameter kept for API consistency. Default is "navy".
#' @param legend_position Position of the legend. Default is "bottom".
#' @param ylab Label for the y-axis. Default is "Frequency".
#'
#' @return A ggplot object showing protein binding frequency across the 2 regions
#'   for Retained, Excluded, and Control groups. The bottom schematic shows two
#'   exon boxes connected by a single intron line. Returns a data frame if
#'   return_data = TRUE.
#'
#' @details
#' The function divides each retained intron event into 2 regions of
#' (WidthIntoExon + WidthIntoIntron) bp each:
#' \itemize{
#'   \item Region 1: Upstream exon end to retained intron start
#'   \item Region 2: Retained intron end to downstream exon start
#' }
#'
#' @examples
#' \dontrun{
#' # Load BED file and RI.MATS data
#' bed <- checkBed("peaks.bed")
#' rimats <- read.table("RI.MATS.JC.txt", header = TRUE)
#'
#' # Basic usage
#' createRetainedIntronSplicingMap(bed_file = bed, RIMATS = rimats)
#'
#' # Return data instead of plot
#' freq_data <- createRetainedIntronSplicingMap(bed_file = bed, RIMATS = rimats,
#'                                       return_data = TRUE)
#' }
#'
#' @export
createRetainedIntronSplicingMap <- function(bed_file,
                                     RIMATS,
                                     moving_average = 50,
                                     WidthIntoExon = 50,
                                     WidthIntoIntron = 300,
                                     p_valueRetainedAndExclusion = 0.05,
                                     p_valueControls = 0.95,
                                     retained_IncLevelDifference = 0.1,
                                     exclusion_IncLevelDifference = -0.1,
                                     Min_Count = 50,
                                     read_count_cols = c("IJC_SAMPLE_1", "SJC_SAMPLE_1",
                                                         "IJC_SAMPLE_2", "SJC_SAMPLE_2"),
                                     groups = c("Retained", "Excluded", "Control"),
                                     control_multiplier = 2.0,
                                     control_iterations = 20,
                                     z_threshold = 1.96,
                                     min_consecutive = 10,
                                     one_sided = TRUE,
                                     use_fdr = TRUE,
                                     fdr_threshold = 0.05,
                                     show_significance = TRUE,
                                     return_data = FALSE,
                                     return_diagnostics = FALSE,
                                     verbose = TRUE,
                                     progress_callback = NULL,
                                     title = "",
                                     retained_col = "blue",
                                     excluded_col = "red",
                                     control_col = "black",
                                     line_width = 0.8,
                                     line_alpha = 0.7,
                                     ribbon_alpha = 0.3,
                                     title_size = 20,
                                     title_color = "black",
                                     axis_text_size = 11,
                                     boundary_col = "gray70",
                                     exon_col = "navy",
                                     legend_position = "bottom",
                                     ylab = "Frequency") {

  required_cols <- c("chr", "strand",
                     "upstreamES", "upstreamEE",
                     "downstreamES", "downstreamEE",
                     "GeneID", "PValue", "FDR", "IncLevelDifference",
                     "IncLevel1", "IncLevel2")
  missing_cols <- setdiff(required_cols, colnames(RIMATS))
  if (length(missing_cols) > 0) {
    stop("RIMATS is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  .splicing_map_worker(
    mats_data                    = RIMATS,
    bed_file                     = bed_file,
    bins_fn                      = make_ri_bins_matrix,
    n_bins                       = 2L,
    plot_fn                      = plot_retained_intron_map,
    moving_average               = moving_average,
    WidthIntoExon                = WidthIntoExon,
    WidthIntoIntron              = WidthIntoIntron,
    p_valueRetainedAndExclusion  = p_valueRetainedAndExclusion,
    p_valueControls              = p_valueControls,
    retained_IncLevelDifference  = retained_IncLevelDifference,
    exclusion_IncLevelDifference = exclusion_IncLevelDifference,
    Min_Count                    = Min_Count,
    read_count_cols              = read_count_cols,
    groups                       = groups,
    control_multiplier           = control_multiplier,
    control_iterations           = control_iterations,
    z_threshold                  = z_threshold,
    min_consecutive              = min_consecutive,
    one_sided                    = one_sided,
    use_fdr                      = use_fdr,
    fdr_threshold                = fdr_threshold,
    show_significance            = show_significance,
    return_data                  = return_data,
    return_diagnostics           = return_diagnostics,
    verbose                      = verbose,
    progress_callback            = progress_callback,
    title                        = title,
    retained_col                 = retained_col,
    excluded_col                 = excluded_col,
    control_col                  = control_col,
    line_width                   = line_width,
    line_alpha                   = line_alpha,
    ribbon_alpha                 = ribbon_alpha,
    title_size                   = title_size,
    title_color                  = title_color,
    axis_text_size               = axis_text_size,
    boundary_col                 = boundary_col,
    exon_col                     = exon_col,
    legend_position              = legend_position,
    ylab                         = ylab
  )
}
