#' Create Retained Intron Splicing Map
#'
#' Analyzes protein binding frequency across retained intron junction regions.
#' Uses a 2-region structure to show where protein binding sites appear
#' relative to the upstream exon/intron and intron/downstream exon boundaries.
#' Filters events into Retained, Excluded, and Control groups.
#'
#' @param bed_file Either a file path to a BED file or a data frame containing
#'   BED data with columns: chr, start, end, tag, score, strand
#' @param RIMATS A data frame containing rMATS output with columns:
#'   chr, strand, upstreamES, upstreamEE, downstreamES, downstreamEE,
#'   GeneID, PValue, FDR, IncLevelDifference, IJC_SAMPLE_1, SJC_SAMPLE_1,
#'   IJC_SAMPLE_2, SJC_SAMPLE_2, IncLevel1, IncLevel2
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
#'   c("Retained", "Excluded", "Control") to process all groups.
#' @param control_multiplier Numeric multiplier for control sample size. The
#'   number of control events sampled per iteration is
#'   (n_retained + n_excluded) * control_multiplier. Default is 2.0.
#' @param control_iterations Integer number for sampling iterations for control
#'   sampling. The final control frequency is the mean across iterations, with
#'   standard deviation shown as a shaded band. Default is 20.
#' @param cores Number of cores for parallel processing. Default is 1 (sequential).
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
#'   \item Region 1 (UE-RI5): Upstream exon end to retained intron
#'   \item Region 2 (RI3-DE): Retained intron end to downstream exon start
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
#' # Load BED file and RI.MATS data
#' bed <- checkBed("peaks.bed")
#' rimats <- read.table("RI.MATS.JC.txt", header = TRUE)
#'
#' # Basic usage
#' createRetainedIntronMap(bed_file = bed, RIMATS = rimats)
#'
#' # Return data instead of plot
#' freq_data <- createRetainedIntronMap(bed_file = bed, RIMATS = rimats,
#'                                       return_data = TRUE)
#' }
#'
#' @export
createRetainedIntronMap <- function(bed_file,
                                     RIMATS,
                                     moving_average = 50,
                                     WidthIntoExon = 50,
                                     WidthIntoIntron = 300,
                                     p_valueRetainedAndExclusion = 0.05,
                                     p_valueControls = 0.95,
                                     retained_IncLevelDifference = 0.1,
                                     exclusion_IncLevelDifference = -0.1,
                                     Min_Count = 50,
                                     groups = c("Retained", "Excluded", "Control"),
                                     control_multiplier = 2.0,
                                     control_iterations = 20,
                                     cores = 1,
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

  # Validate required RI.MATS columns
  required_cols <- c("chr", "strand",
                      "upstreamES", "upstreamEE",
                      "exonStart_0base", "exonEnd",
                      "downstreamES", "downstreamEE",
                      "GeneID", "PValue", "FDR", "IncLevelDifference")
  missing_cols <- setdiff(required_cols, colnames(RIMATS))
  if (length(missing_cols) > 0) {
    stop("RIMATS is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Load BED file if path provided
  if (is.character(bed_file)) {
    bed_data <- utils::read.table(bed_file)
  } else {
    bed_data <- bed_file
  }

  # Normalize chromosome names
  bed_data$chr <- toupper(bed_data$chr)
  RIMATS$chr <- sub("^chr", "", RIMATS$chr)

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

  if (verbose) message("Processing groups: ", paste(groups, collapse = ", "))

  .report_progress <- function(current, total, detail = NULL) {
    if (is.function(progress_callback)) {
      try(progress_callback(current, total, detail), silent = TRUE)
    }
  }

  # Cap cores at max available
  max_cores <- parallel::detectCores() - 1
  if (is.na(max_cores) || max_cores < 1) max_cores <- 1
  cores <- min(cores, max_cores)
  cores <- max(cores, 1)

  # Filter events using the shared SE.MATS filter (same statistical columns)
  filtered_events <- filter_SEMATS_events(
    RIMATS,
    p_valueRetainedAndExclusion = p_valueRetainedAndExclusion,
    p_valueControls = p_valueControls,
    retained_IncLevelDifference = retained_IncLevelDifference,
    exclusion_IncLevelDifference = exclusion_IncLevelDifference,
    Min_Count = Min_Count
  )

  bin_width <- WidthIntoExon + WidthIntoIntron + 1
  n_bins <- 2

  # Helper function to process a group
  process_group <- function(data, group_name, all_data_n = NULL) {
    if (nrow(data) == 0) {
      if (verbose) message("No events found for group: ", group_name)
      return(data.frame(
        global_position = 1:(n_bins * bin_width),
        overlap_count = 0L,
        frequency = 0,
        bin = rep(1:n_bins, each = bin_width),
        moving_avg = 0,
        group = group_name
      ))
    }

    data$group <- group_name
    bins_gr <- make_ri_bins_matrix(data,
                                    WidthIntoExon = WidthIntoExon,
                                    WidthIntoIntron = WidthIntoIntron)

    freq_data <- calculate_binding_frequency(bins_gr,
                                              buckets,
                                              bin_width,
                                              cores = cores,
                                              n_bins = n_bins)

    total_events <- if (!is.null(all_data_n)) all_data_n else nrow(data)
    freq_data$frequency <- freq_data$overlap_count / total_events

    freq_data <- calculate_moving_average(freq_data, moving_average, bins = bin_width)

    freq_data$group <- group_name
    return(freq_data)
  }

  results_list <- list()

  if ("Retained" %in% groups) {
    if (verbose) message("Processing Retained events...")
    .report_progress(1, 100, "Processing Retained events...")
    results_list$Retained <- process_group(filtered_events$Retained, "Retained")
    results_list$Retained$moving_avg_sd <- 0
  }

  if ("Excluded" %in% groups) {
    if (verbose) message("Processing Excluded events...")
    .report_progress(5, 100, "Processing Excluded events...")
    results_list$Excluded <- process_group(filtered_events$Excluded, "Excluded")
    results_list$Excluded$moving_avg_sd <- 0
  }

  if ("Control" %in% groups) {
    if (verbose) message("Processing Control events with sampling...")
    .report_progress(10, 100, "Preparing control sampling...")
    control_data <- filtered_events$Control
    n_controls <- nrow(control_data)

    n_retained <- nrow(filtered_events$Retained)
    n_excluded <- nrow(filtered_events$Excluded)
    sample_size <- round((n_retained + n_excluded) * control_multiplier)

    if (verbose) {
      message(sprintf("  Retained: %d, Excluded: %d, Sample size: %d",
                      n_retained, n_excluded, sample_size))
      message(sprintf("  Controls: %d, Iterations: %d",
                      n_controls, control_iterations))
    }

    if (n_controls == 0) {
      if (verbose) message("No control events found")
      results_list$Control <- data.frame(
        global_position = 1:(n_bins * bin_width),
        overlap_count = 0L,
        frequency = 0,
        bin = rep(1:n_bins, each = bin_width),
        moving_avg = 0,
        group = "Control",
        moving_avg_sd = 0
      )
    } else if (sample_size >= n_controls || sample_size == 0) {
      if (verbose) message("Using all controls without bootstrap")
      results_list$Control <- process_group(control_data, "Control")
      results_list$Control$moving_avg_sd <- 0
    } else {
      iteration_results <- vector("list", control_iterations)

      apply_moving_avg <- function(freq_vec) {
        if (is.null(moving_average) || moving_average <= 0) return(freq_vec)
        result <- numeric(length(freq_vec))
        half_window <- floor((moving_average - 1) / 2)
        for (b in 1:n_bins) {
          bin_start <- (b - 1) * bin_width + 1
          bin_end <- b * bin_width
          bin_vals <- freq_vec[bin_start:bin_end]
          smoothed <- slider::slide_dbl(bin_vals, mean,
                                         .before = half_window,
                                         .after = half_window,
                                         .complete = FALSE)
          result[bin_start:bin_end] <- smoothed
        }
        return(result)
      }

      if (verbose) message("  Pre-computing binding cache for all control events...")
      .report_progress(10, 100, "Building control binding cache...")

      cache_matrix <- precompute_binding_cache(
        events_data = control_data,
        protein = buckets,
        WidthIntoExon = WidthIntoExon,
        WidthIntoIntron = WidthIntoIntron,
        verbose = verbose,
        n_bins = n_bins,
        make_bins_fn = make_ri_bins_matrix
      )

      if (verbose) message("  Cache built. Starting bootstrap iterations...")

      all_sampled_indices <- lapply(seq_len(control_iterations), function(i) {
        sample(n_controls, sample_size, replace = FALSE)
      })

      pb <- progress::progress_bar$new(
        format = "  Sampling iterations [:bar] :current/:total (:percent) eta::eta",
        total = control_iterations, clear = FALSE, width = 80
      )

      loop_start <- 20
      loop_end <- 90

      for (iter in seq_len(control_iterations)) {
        pb$tick()
        progress_pct <- loop_start + (iter / control_iterations) * (loop_end - loop_start)
        .report_progress(
          progress_pct,
          100,
          sprintf("Control sampling iteration %d/%d", iter, control_iterations)
        )

        sampled_ids <- all_sampled_indices[[iter]]
        overlap_counts <- rowSums(cache_matrix[, sampled_ids, drop = FALSE])
        freq_vec <- overlap_counts / sample_size
        iteration_results[[iter]] <- apply_moving_avg(freq_vec)
      }

      freq_matrix <- do.call(cbind, iteration_results)
      mean_freq <- rowMeans(freq_matrix, na.rm = TRUE)
      sd_freq <- apply(freq_matrix, 1, stats::sd, na.rm = TRUE)

      results_list$Control <- data.frame(
        global_position = 1:(n_bins * bin_width),
        overlap_count = NA_integer_,
        frequency = mean_freq,
        bin = rep(1:n_bins, each = bin_width),
        moving_avg = mean_freq,
        group = "Control",
        moving_avg_sd = sd_freq
      )

      bootstrap_matrix <- freq_matrix
    }
  }

  combined_data <- dplyr::bind_rows(results_list)

  missing_groups <- setdiff(groups, unique(combined_data$group))
  if (length(missing_groups) > 0) {
    if (verbose) message("No events for: ", paste(missing_groups, collapse = ", "))
    zero_data <- do.call(rbind, lapply(missing_groups, function(g) {
      data.frame(
        global_position = 1:(n_bins * bin_width),
        overlap_count = 0L,
        frequency = 0,
        bin = rep(1:n_bins, each = bin_width),
        moving_avg = 0,
        group = g,
        moving_avg_sd = 0
      )
    }))
    combined_data <- dplyr::bind_rows(combined_data, zero_data)
  }

  if (return_data) return(combined_data)

  if (return_diagnostics) {
    diagnostics <- list(
      data = combined_data,
      bootstrap_matrix = if (exists("bootstrap_matrix")) bootstrap_matrix else NULL,
      n_iterations = if (exists("bootstrap_matrix")) control_iterations else 0,
      sample_size = if (exists("sample_size")) sample_size else NA,
      n_controls = if (exists("n_controls")) n_controls else NA
    )
    return(diagnostics)
  }

  .report_progress(92, 100, "Combining results...")

  sig_regions <- NULL
  if (show_significance && "Control" %in% groups) {
    .report_progress(96, 100, "Calculating significance...")
    control_has_sd <- any(combined_data$moving_avg_sd[combined_data$group == "Control"] > 0,
                          na.rm = TRUE)
    if (control_has_sd) {
      if (verbose) message("Calculating significance...")
      sig_result <- calculate_significance(
        combined_data,
        z_threshold = z_threshold,
        min_consecutive = min_consecutive,
        compare_to = "Control",
        one_sided = one_sided,
        use_fdr = use_fdr,
        fdr_threshold = fdr_threshold
      )
      sig_regions <- sig_result$significant_regions

      if (verbose) {
        if (!is.null(sig_regions) && nrow(sig_regions) > 0) {
          message(sprintf("Found %d significant regions", nrow(sig_regions)))
        } else {
          message("No significant regions found")
        }
      }
    } else if (verbose) {
      message("Skipping significance: Control SD is zero")
    }
  }

  .report_progress(100, 100, "Complete")
  plot_retained_intron_map(combined_data,
                            WidthIntoExon = WidthIntoExon,
                            WidthIntoIntron = WidthIntoIntron,
                            title = title,
                            sig_regions = sig_regions,
                            retained_cutoff = retained_IncLevelDifference,
                            excluded_cutoff = exclusion_IncLevelDifference,
                            retained_col = retained_col,
                            excluded_col = excluded_col,
                            control_col = control_col,
                            line_width = line_width,
                            line_alpha = line_alpha,
                            ribbon_alpha = ribbon_alpha,
                            title_size = title_size,
                            title_color = title_color,
                            axis_text_size = axis_text_size,
                            boundary_col = boundary_col,
                            exon_col = exon_col,
                            legend_position = legend_position,
                            ylab = ylab)
}
