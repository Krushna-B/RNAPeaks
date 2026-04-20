#' Create Retained Intron Sequence Map
#'
#' Analyzes the frequency of a target sequence motif across retained intron
#' junction regions. Compares motif frequency between Retained, Excluded, and
#' Control events to identify position-specific enrichment patterns around the
#' upstream exon/intron and intron/downstream exon boundaries.
#'
#' @param RIMATS A data frame containing retained intron outputs.
#' @param sequence Character string or character vector of sequence motifs to search
#'   for (e.g., \code{"YCAY"} or \code{c("YCAY", "CCCC")}). Supports IUPAC ambiguity
#'   codes. When multiple motifs are provided, behavior depends on \code{motif_mode}.
#' @param motif_mode How to handle multiple motifs. \code{"combined"} (default) treats
#'   all motifs as a single hit set, a position counts if any motif matches there
#'   and returns one plot. \code{"individual"} runs the full analysis independently for
#'   each motif and returns a named list of plots (one per motif). Ignored when
#'   \code{sequence} is a single motif.
#' @param genome A BSgenome object. Default uses BSgenome.Hsapiens.UCSC.hg38.
#' @param moving_average Integer specifying the window size for moving average
#'   smoothing. Set to NULL or 0 to disable smoothing. Default is 40.
#' @param WidthIntoExon Integer specifying how many bp to extend into exons.
#'   Default is 50.
#' @param WidthIntoIntron Integer specifying how many bp to extend into introns.
#'   Default is 250.
#' @param p_valueRetainedAndExclusion P-value threshold for retained/excluded events.
#'   Default is 0.05.
#' @param p_valueControls P-value threshold for control events. Default is 0.95.
#' @param retained_IncLevelDifference Inclusion level difference threshold for
#'   retained events. Default is 0.1.
#' @param exclusion_IncLevelDifference Inclusion level difference threshold for
#'   excluded events. Default is -0.1.
#' @param Min_Count Minimum read count threshold. Default is 50.
#' @param groups Character vector specifying which event groups to process.
#'   Options are "Retained", "Excluded", and/or "Control". Default is
#'   c("Retained", "Excluded", "Control") to process all groups.
#' @param control_multiplier Numeric multiplier for control sample size.
#'   Default is 2.0.
#' @param control_iterations Integer number for sampling iterations for control
#'   sampling. Default is 20.
#' @param z_threshold Z-score threshold for significance testing. Default is 1.96.
#'   Only used when use_fdr = FALSE.
#' @param min_consecutive Minimum number of consecutive significant positions
#'   required to form a significant region. Default is 10.
#' @param one_sided Logical. If TRUE (default), only test for enrichment.
#' @param use_fdr Logical. If TRUE, use FDR-corrected p-values. Default is TRUE.
#' @param fdr_threshold FDR threshold for significance when use_fdr = TRUE.
#'   Default is 0.05.
#' @param show_significance Logical. If TRUE (default), displays colored bars above
#'   the plot indicating significant regions.
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
#' @return A ggplot object showing sequence motif frequency across the 2 regions
#'   for Retained, Excluded, and Control groups. The bottom schematic shows two
#'   exon boxes connected by a single intron line. Returns a data frame if
#'   \code{return_data = TRUE}. When \code{motif_mode = "individual"} and multiple
#'   motifs are supplied, returns a named list of ggplot objects (or data frames),
#'   one entry per motif.
#'
#' @details
#' The function divides each retained intron event into 2 regions of
#' (WidthIntoExon + WidthIntoIntron) bp each:
#' \itemize{
#'   \item Region 1: Upstream exon end to retained intron start
#'   \item Region 2: Retained intron end to downstream exon start
#' }
#'
#' At each position, the function checks if the target sequence starts there.
#' The frequency is calculated as: (events with motif at position) / (total events)
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' rimats <- read.table("RI.MATS.JC.txt", header = TRUE)
#'
#' # Basic usage
#' createRetainedIntronSequenceMap(RIMATS = rimats, sequence = "CCCC")
#'
#' # Search for YCAY motif (Y = C or T)
#' createRetainedIntronSequenceMap(RIMATS = rimats, sequence = "YCAY")
#'
#' # Return data instead of plot
#' freq_data <- createRetainedIntronSequenceMap(RIMATS = rimats,
#'                                               sequence = "GGGG",
#'                                               return_data = TRUE)
#' }
#'
#' @export
createRetainedIntronSequenceMap <- function(RIMATS,
                                             sequence,
                                             motif_mode = c("combined", "individual"),
                                             genome = NULL,
                                             moving_average = 40,
                                             WidthIntoExon = 50,
                                             WidthIntoIntron = 250,
                                             p_valueRetainedAndExclusion = 0.05,
                                             p_valueControls = 0.95,
                                             retained_IncLevelDifference = 0.1,
                                             exclusion_IncLevelDifference = -0.1,
                                             Min_Count = 50,
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
                     "IJC_SAMPLE_1", "SJC_SAMPLE_1",
                     "IJC_SAMPLE_2", "SJC_SAMPLE_2")
  missing_cols <- setdiff(required_cols, colnames(RIMATS))
  if (length(missing_cols) > 0) {
    stop("RIMATS is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (missing(sequence) || !is.character(sequence) ||
      length(sequence) == 0 || any(nchar(trimws(sequence)) == 0)) {
    stop("A valid sequence motif (or character vector of motifs) must be provided")
  }
  sequence   <- toupper(trimws(sequence))
  motif_mode <- match.arg(motif_mode)

  if (motif_mode == "individual" && length(sequence) > 1) {
    plot_list <- lapply(sequence, function(motif) {
      motif_title <- if (nchar(title) == 0) motif else paste0(title, " \u2014 ", motif)
      createRetainedIntronSequenceMap(
        RIMATS                       = RIMATS,
        sequence                     = motif,
        motif_mode                   = "combined",
        genome                       = genome,
        moving_average               = moving_average,
        WidthIntoExon                = WidthIntoExon,
        WidthIntoIntron              = WidthIntoIntron,
        p_valueRetainedAndExclusion  = p_valueRetainedAndExclusion,
        p_valueControls              = p_valueControls,
        retained_IncLevelDifference  = retained_IncLevelDifference,
        exclusion_IncLevelDifference = exclusion_IncLevelDifference,
        Min_Count                    = Min_Count,
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
        progress_callback            = NULL,
        title                        = motif_title,
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
    })
    names(plot_list) <- sequence
    return(plot_list)
  }

  .sequence_map_worker(
    mats_data                    = RIMATS,
    bins_fn                      = make_ri_bins_matrix,
    n_bins                       = 2L,
    plot_fn                      = plot_retained_intron_map,
    genome                       = genome,
    sequence                     = sequence,
    moving_average               = moving_average,
    WidthIntoExon                = WidthIntoExon,
    WidthIntoIntron              = WidthIntoIntron,
    p_valueRetainedAndExclusion  = p_valueRetainedAndExclusion,
    p_valueControls              = p_valueControls,
    retained_IncLevelDifference  = retained_IncLevelDifference,
    exclusion_IncLevelDifference = exclusion_IncLevelDifference,
    Min_Count                    = Min_Count,
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
