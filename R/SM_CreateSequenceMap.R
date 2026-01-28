#' Create Sequence Map
#'
#' Analyzes the frequency of a target sequence motif across splicing junction
#' regions.Filters events into Retained, Excluded, and Control groups.
#'
#' @param SEMATS A data frame containing SE.MATS output with columns:
#'   chr, strand, upstreamES, upstreamEE, exonStart_0base, exonEnd,
#'   downstreamES, downstreamEE, GeneID, PValue, FDR, IncLevelDifference,
#'   IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2, IncLevel1, IncLevel2
#' @param sequence Character string of the target sequence motif to search for
#'   (e.g., "CCCC", "YGCY"). Supports IUPAC ambiguity codes.
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
#'   c("Retained", "Excluded", "Control") to process all groups. Use
#'   c("Retained", "Excluded") to skip the Control group (which can be large).
#' @param return_data Logical. If TRUE, returns the frequency data frame instead
#'   of a plot. Default is FALSE.
#'
#' @return A ggplot object showing sequence frequency across the 4 regions
#'   for Retained (blue), Excluded (red), and Control (black) groups,
#'   or a data frame if return_data = TRUE.
#'
#' @details
#' The function divides each splicing event into 4 regions of (WidthIntoExon +
#' WidthIntoIntron) bp each:
#' \itemize{
#'   \item Region 1: Upstream exon end to first intron
#'   \item Region 2: First intron end to middle (skipped) exon start
#'   \item Region 3: Middle exon end to second intron
#'   \item Region 4: Second intron end to downstream exon start
#' }
#'
#' Events are filtered into three groups:
#' \itemize{
#'   \item Retained: Significant events (PValue < threshold) with positive inclusion
#'   \item Excluded: Significant events (PValue < threshold) with negative inclusion
#'   \item Control: Non-significant events
#' }
#'
#' At each position, the function checks if the target sequence starts there.
#' The frequency is calculated as: (events with motif at position) / (total events)
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Basic usage
#' createSequenceMap(SEMATS = sample_se.mats, sequence = "CCCC")
#'
#' # Search for YCAY motif (Y = C or T)
#' createSequenceMap(SEMATS = sample_se.mats, sequence = "YCAY")
#'
#' # Return data instead of plot
#' freq_data <- createSequenceMap(SEMATS = sample_se.mats,
#'                                 sequence = "GGGG",
#'                                 return_data = TRUE)
#' }
#'
#' @export
createSequenceMap <- function(SEMATS,
                               sequence,
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
                               cores = 1,
                               return_data = FALSE) {

  # Load default genome if not provided
  if (is.null(genome)) {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      stop("BSgenome.Hsapiens.UCSC.hg38 is required. Install with:\n",
           "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')")
    }
    genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
  }

  # Validate sequence input
  if (missing(sequence) || !is.character(sequence) || nchar(sequence) == 0) {
    stop("A valid sequence motif must be provided")
  }
  sequence <- toupper(sequence)

  # Validate groups parameter
  valid_groups <- c("Retained", "Excluded", "Control")
  groups <- match.arg(groups, valid_groups, several.ok = TRUE)

  if (length(groups) == 0) {
    stop("At least one group must be specified")
  }

  message(paste0("Processing groups: ", paste(groups, collapse = ", ")))

  # Cap cores at max available - 1
  max_cores <- parallel::detectCores() - 1
  if (is.na(max_cores) || max_cores < 1) max_cores <- 1
  cores <- min(cores, max_cores)
  cores <- max(cores, 1)
  options(future.globals.maxSize = 8 * 1024^3)

  # If using parallel, warm up workers early while we do other setup
  warmup_future <- NULL
  if (cores > 1) {
    message(sprintf("Starting %d parallel workers in background...", cores))
    future::plan(future::multisession, workers = cores)

    # workers will spawn and load packages
    warmup_future <- future::future({
      # library(Biostrings)
      # library(IRanges)
      # library(GenomicRanges)
      # library(GenomeInfoDb)
      TRUE
    }, seed = TRUE)
  }


  # Filter SEMATS into Controls, Retained, and Excluded
  filtered_events <- filter_SEMATS_events(
    SEMATS,
    p_valueRetainedAndExclusion = p_valueRetainedAndExclusion,
    p_valueControls = p_valueControls,
    retained_IncLevelDifference = retained_IncLevelDifference,
    exclusion_IncLevelDifference = exclusion_IncLevelDifference,
    Min_Count = Min_Count
  )

  bin_width <- WidthIntoExon + WidthIntoIntron

  # Process each group
  process_group <- function(data, group_name) {
    if (nrow(data) == 0) {
      message(paste0("No events found for group: ", group_name))
      return(data.frame(
        global_position = 1:(4 * bin_width),
        frequency = 0,
        moving_avg = 0,
        group = group_name
      ))
    }

    # Build bins matrix
    bins_gr <- make_bins_matrix(data,
                                 WidthIntoExon = WidthIntoExon,
                                 WidthIntoIntron = WidthIntoIntron)

    # Calculate sequence frequency
    freq_data <- calculate_sequence_frequency(bins_gr,
                                               sequence,
                                               bsgenome_obj = genome,
                                               bin_width,
                                              cores = cores)

    total_events <- length(unique(bins_gr$event_id))
    freq_data$frequency <- freq_data$match_count / total_events

    # Apply moving average if specified
    if (!is.null(moving_average) && moving_average > 0) {
      freq_data <- freq_data %>%
        dplyr::mutate(
          bin = dplyr::case_when(
            global_position <= bin_width ~ 1,
            global_position <= 2 * bin_width ~ 2,
            global_position <= 3 * bin_width ~ 3,
            TRUE ~ 4
          )
        ) %>%
        dplyr::arrange(global_position) %>%
        dplyr::group_by(bin) %>%
        dplyr::mutate(moving_avg = slider::slide_dbl(frequency,
                                                       mean,
                                                       .before = floor((moving_average - 1) / 2),
                                                       .after = floor((moving_average - 1) / 2),
                                                       .complete = FALSE)) %>%
        dplyr::ungroup()
    } else {
      freq_data <- freq_data %>%
        dplyr::mutate(
          bin = dplyr::case_when(
            global_position <= bin_width ~ 1,
            global_position <= 2 * bin_width ~ 2,
            global_position <= 3 * bin_width ~ 3,
            TRUE ~ 4
          ),
          moving_avg = frequency
        )
    }

    freq_data$group <- group_name
    return(freq_data)
  }

  # Wait for parallel workers to be ready
  if (!is.null(warmup_future)) {
    invisible(future::value(warmup_future))
  }

  # Process only selected groups
  results_list <- list()

  if ("Retained" %in% groups) {
    message("Processing Retained events...")
    results_list$Retained <- process_group(filtered_events$Retained, "Retained")
  }

  if ("Excluded" %in% groups) {
    message("Processing Excluded events...")
    results_list$Excluded <- process_group(filtered_events$Excluded, "Excluded")
  }

  if ("Control" %in% groups) {
    message("Processing Control events...")
    results_list$Control <- process_group(filtered_events$Control, "Control")
  }

  # Combine selected groups
  combined_data <- do.call(rbind, results_list)

  # Clean up parallel workers
  if (cores > 1) {
    future::plan(future::sequential)
  }


  # Plot using the shared plotting function
  plot_splicing_sequence_map(combined_data,
                              WidthIntoExon = WidthIntoExon,
                              WidthIntoIntron = WidthIntoIntron,
                              title = paste0("Sequence Frequency: ", sequence))
}
