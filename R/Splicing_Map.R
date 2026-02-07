#' Create Splicing Map
#'
#' Analyzes protein binding frequency across splicing junction regions.
#' Uses a 4-region structure to show where protein binding sites appear
#' relative to exon/intron boundaries. Filters events into Retained,
#' Excluded, and Control groups.
#'
#' @param bed_file Either a file path to a BED file or a data frame containing
#'   BED data with columns: chr, start, end, tag, score, strand
#' @param SEMATS A data frame containing SE.MATS output with columns:
#'   chr, strand, upstreamES, upstreamEE, exonStart_0base, exonEnd,
#'   downstreamES, downstreamEE, GeneID, PValue, FDR, IncLevelDifference,
#'   IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2, IncLevel1, IncLevel2
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
#' @param Min_Count Minimum read count threshold. Default is 50.
#' @param groups Character vector specifying which event groups to process.
#'   Options are "Retained", "Excluded", and/or "Control". Default is
#'   c("Retained", "Excluded", "Control") to process all groups. Use
#'   c("Retained", "Excluded") to skip the Control group (which can be large).
#' @param cores Number of cores for parallel processing. Default is 1 (sequential).
#' @param return_data Logical. If TRUE, returns the frequency data frame instead
#'   of a plot. Default is FALSE.
#'
#' @return A ggplot object showing protein binding frequency across the 4 regions
#'   for Retained (blue), Excluded (red), and Control (black) groups,
#'   or a data frame if return_data = TRUE.
#'
#' @details
#' The function divides each splicing event into 4 regions of (WidthIntoExon +
#' WidthIntoIntron) bp each:
#' \itemize{
#'   \item Region 1 (UE-UI5): Upstream exon end to first intron
#'   \item Region 2 (UI3-EX3): First intron end to middle (skipped) exon start
#'   \item Region 3 (EX5-DI5): Middle exon end to second intron
#'   \item Region 4 (DI3-DE): Second intron end to downstream exon start
#' }
#'
#' Events are filtered into three groups:
#' \itemize{
#'   \item Retained: Significant events (PValue < threshold) with negative IncLevelDifference
#'   \item Excluded: Significant events (PValue < threshold) with positive IncLevelDifference
#'   \item Control: Non-significant events with stable inclusion levels
#' }
#'
#' @examples
#' \dontrun{
#' # Load BED file and SE.MATS data
#' bed <- read.table("peaks.bed")
#' semats <- read.table("SE.MATS.JC.txt", header = TRUE)
#'
#' # Basic usage
#' createSplicingMap(bed_file = bed, SEMATS = semats)
#'
#' # Use parallel processing
#' createSplicingMap(bed_file = bed, SEMATS = semats, cores = 4)
#'
#' # Return data instead of plot
#' freq_data <- createSplicingMap(bed_file = bed,
#'                                 SEMATS = semats,
#'                                 return_data = TRUE)
#' }
#'
#' @export
createSplicingMap <- function(bed_file,
                               SEMATS,
                               moving_average = 50,
                               WidthIntoExon = 50,
                               WidthIntoIntron = 300,
                               p_valueRetainedAndExclusion = 0.05,
                               p_valueControls = 0.95,
                               retained_IncLevelDifference = 0.1,
                               exclusion_IncLevelDifference = -0.1,
                               Min_Count = 50,
                               groups = c("Retained", "Excluded", "Control"),
                               cores = 1,
                               return_data = FALSE) {

  # Load BED file if path provided
  if (is.character(bed_file)) {
    bed_data <- utils::read.table(bed_file)
  } else {
    bed_data <- bed_file
  }

  # Normalize chromosome names (handle lowercase x, y, m and chr prefix mismatches)
  bed_data$chr <- toupper(bed_data$chr)
  SEMATS$chr <- sub("^chr", "", SEMATS$chr)

  # Convert BED to GRanges and reduce overlapping peaks

  buckets <- GenomicRanges::makeGRangesFromDataFrame(
    bed_data,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )
  buckets <- GenomicRanges::reduce(buckets)

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

  # Filter SEMATS into Controls, Retained, and Excluded using helper function
  filtered_events <- filter_SEMATS_events(
    SEMATS,
    p_valueRetainedAndExclusion = p_valueRetainedAndExclusion,
    p_valueControls = p_valueControls,
    retained_IncLevelDifference = retained_IncLevelDifference,
    exclusion_IncLevelDifference = exclusion_IncLevelDifference,
    Min_Count = Min_Count
  )

  bin_width <- WidthIntoExon + WidthIntoIntron + 1

  # Combine all selected groups into one data frame with group labels
  combined_events <- do.call(rbind, lapply(groups, function(g) {
    data <- filtered_events[[g]]
    if (nrow(data) > 0) {
      data$group <- g
    }
    data
  }))

  # Count events per group for frequency calculation (convert to named vector for safe subsetting)
  event_counts <- table(combined_events$group)
  events_per_group <- setNames(as.integer(event_counts), names(event_counts))

  if (nrow(combined_events) == 0) {
    message("No events found in any selected group")
    combined_data <- do.call(rbind, lapply(groups, function(g) {
      data.frame(
        global_position = 1:(4 * bin_width),
        frequency = 0,
        moving_avg = 0,
        group = g
      )
    }))
  } else {
    message(sprintf("Processing %d total events across %d groups in single pass...",
                    nrow(combined_events), length(unique(combined_events$group))))

    # Build bins matrix for ALL events at once (group column is preserved)
    bins_gr <- make_bins_matrix(combined_events,
                                 WidthIntoExon = WidthIntoExon,
                                 WidthIntoIntron = WidthIntoIntron)

    # Calculate binding frequency for all groups in one vectorized pass
    freq_data <- calculate_binding_frequency(bins_gr,
                                              buckets,
                                              bin_width,
                                              cores = cores)

    # Calculate per-group frequency using per-group event counts
    freq_data$frequency <- freq_data$overlap_count / as.numeric(events_per_group[freq_data$group])

    # Apply moving average per group
    combined_data <- freq_data %>%
      dplyr::group_by(group) %>%
      dplyr::group_modify(~ calculate_moving_average(.x, moving_average, bins = bin_width)) %>%
      dplyr::ungroup()
  }

  # Handle any groups with no events (add zero rows)
  missing_groups <- setdiff(groups, unique(combined_data$group))
  if (length(missing_groups) > 0) {
    message(paste0("No events found for groups: ", paste(missing_groups, collapse = ", ")))
    zero_data <- do.call(rbind, lapply(missing_groups, function(g) {
      data.frame(
        global_position = 1:(4 * bin_width),
        overlap_count = 0L,
        frequency = 0,
        moving_avg = 0,
        group = g
      )
    }))
    combined_data <- rbind(combined_data, zero_data)
  }

  if (return_data) {
    return(combined_data)
  }

  # Plot using the shared plotting function
  plot_splicing_sequence_map(combined_data,
                              WidthIntoExon = WidthIntoExon,
                              WidthIntoIntron = WidthIntoIntron,
                              title = paste0("Splicing Map Peaks"))
}

