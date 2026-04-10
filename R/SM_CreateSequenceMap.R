
#'  Create Sequence Map

#' Analyzes the frequency of a target sequence motif across splicing junction
#' regions. Compares motif frequency between Retained, Excluded, and Control
#' splicing events to identify position-specific enrichment patterns.
#'
#' @param SEMATS A data frame containing SE.MATS output with columns:
#'   chr, strand, upstreamES, upstreamEE, exonStart_0base, exonEnd,
#'   downstreamES, downstreamEE, GeneID, PValue, FDR, IncLevelDifference,
#'   IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2, IncLevel1, IncLevel2
#' @param sequence Character string or character vector of sequence motifs to search
#'   for (e.g., \code{"YCAY"} or \code{c("YCAY", "CCCC")}). Supports IUPAC ambiguity
#'   codes. When multiple motifs are provided, behaviour depends on \code{motif_mode}.
#' @param motif_mode How to handle multiple motifs. \code{"combined"} (default) treats
#'   all motifs as a single hit set — a position counts if any motif matches there —
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
#' @param control_multiplier Numeric multiplier for control sample size. The
#'   number of control events sampled per iteration is
#'   (n_retained + n_excluded) * control_multiplier. Default is 2.0.
#' @param control_iterations Integer number for sampling iterations for control
#'   sampling. The final control frequency is the mean across iterations, with
#'   standard deviation shown as a shaded band. Default is 20.
#' @param z_threshold Z-score threshold for significance testing. Default is 1.96
#' Only used when use_fdr = FALSE.
#' @param min_consecutive Minimum number of consecutive significant positions
#'   required to form a significant region. Default is 10. Helps reduce false
#'   positives from noise.
#' @param one_sided Logical. If TRUE (default), only test for enrichment
#'   (frequency > control). If FALSE, test for both enrichment and depletion.
#' @param use_fdr Logical. If TRUE (default), use FDR-corrected p-values (Benjamini-Hochberg)
#'   for significance testing. If FALSE, use z_threshold directly.
#' @param fdr_threshold FDR threshold for significance when use_fdr = TRUE.
#'   Default is 0.05.
#' @param show_significance Logical. If TRUE (default), displays colored bars above
#'   the plot indicating regions where Retained/Excluded differ significantly
#'   from Control based on z-test.
#' @param return_data Logical. If TRUE, returns the frequency data frame instead
#'   of a plot. Default is FALSE.
#' @param return_diagnostics Logical. If TRUE, returns a list containing the
#'   frequency data, raw bootstrap iteration results (for normality testing),
#'   and significance results. Useful for validating bootstrap assumptions.
#'   Default is FALSE.
#' @param verbose Logical. If TRUE, prints progress messages. Default is TRUE.
#' @param progress_callback Optional function to report progress. Called with two
#'   arguments: current iteration number and total iterations. Default is NULL.
#' @param title Character string for the plot title. Default is "" (no title).
#' @param retained_col Color for the Retained group line. Default is "blue".
#' @param excluded_col Color for the Excluded group line. Default is "red".
#' @param control_col Color for the Control group line. Default is "black".
#' @param line_width Numeric line width for the frequency lines. Default is 0.8.
#' @param line_alpha Numeric alpha (opacity) for the frequency lines. Default is 0.7.
#' @param ribbon_alpha Numeric alpha for the SD ribbon around Control. Default is 0.3.
#' @param title_size Numeric font size for the plot title. Default is 20.
#' @param title_color Color for the plot title text. Default is "black".
#' @param axis_text_size Numeric font size for y-axis tick labels. Default is 11.
#' @param boundary_col Color for the dashed vertical boundary lines. Default is "gray70".
#' @param exon_col Fill color for the skipped (middle) exon in the schematic. Default is "navy".
#' @param legend_position Position of the legend. Default is "bottom".
#' @param ylab Label for the y-axis. Default is "Frequency".
#'
#' @return A ggplot object showing sequence frequency across the 4 regions
#'   for Retained, Excluded, and Control groups. Significant regions (z-test vs
#'   Control) are shown as colored bars above the plot. Returns a data frame if
#'   \code{return_data = TRUE}. When \code{motif_mode = "individual"} and multiple
#'   motifs are supplied, returns a named list of ggplot objects (or data frames),
#'   one entry per motif.
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
#' # Single motif — basic usage (unchanged)
#' createSequenceMap(SEMATS = sample_se.mats, sequence = "YCAY")
#'
#' # Multiple motifs — combined: one plot, hit if any motif matches
#' createSequenceMap(SEMATS = sample_se.mats,
#'                   sequence = c("YCAY", "CCCC"),
#'                   motif_mode = "combined")
#'
#' # Multiple motifs — individual: named list of plots, one per motif
#' plots <- createSequenceMap(SEMATS = sample_se.mats,
#'                             sequence = c("YCAY", "CCCC", "GGGG"),
#'                             motif_mode = "individual")
#' plots[["YCAY"]]
#' plots[["CCCC"]]
#' }
#'
#' @export
createSequenceMap <- function(SEMATS,
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

  # Load default genome if not provided (done once here so individual mode
  # doesn't re-load it on each recursive call)
  if (is.null(genome)) {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      stop("BSgenome.Hsapiens.UCSC.hg38 is required. Install with:\n",
           "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')")
    }
    genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
  }

  # Validate sequence input
  if (missing(sequence) || !is.character(sequence) ||
      length(sequence) == 0 || any(nchar(trimws(sequence)) == 0)) {
    stop("A valid sequence motif (or character vector of motifs) must be provided")
  }
  sequence <- toupper(trimws(sequence))
  motif_mode <- match.arg(motif_mode)

  # Individual mode: run full pipeline once per motif, return named list
  if (motif_mode == "individual" && length(sequence) > 1) {
    plot_list <- lapply(sequence, function(motif) {
      motif_title <- if (nchar(title) == 0) motif else paste0(title, " \u2014 ", motif)
      createSequenceMap(
        SEMATS = SEMATS,
        sequence = motif,
        motif_mode = "combined",
        genome = genome,
        moving_average = moving_average,
        WidthIntoExon = WidthIntoExon,
        WidthIntoIntron = WidthIntoIntron,
        p_valueRetainedAndExclusion = p_valueRetainedAndExclusion,
        p_valueControls = p_valueControls,
        retained_IncLevelDifference = retained_IncLevelDifference,
        exclusion_IncLevelDifference = exclusion_IncLevelDifference,
        Min_Count = Min_Count,
        groups = groups,
        control_multiplier = control_multiplier,
        control_iterations = control_iterations,
        z_threshold = z_threshold,
        min_consecutive = min_consecutive,
        one_sided = one_sided,
        use_fdr = use_fdr,
        fdr_threshold = fdr_threshold,
        show_significance = show_significance,
        return_data = return_data,
        return_diagnostics = return_diagnostics,
        verbose = verbose,
        progress_callback = NULL,
        title = motif_title,
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
        ylab = ylab
      )
    })
    names(plot_list) <- sequence
    return(plot_list)
  }

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

  # Filter SEMATS into Controls, Retained, and Excluded
  filtered_events <- filter_SEMATS_events(
    SEMATS,
    p_valueRetainedAndExclusion = p_valueRetainedAndExclusion,
    p_valueControls = p_valueControls,
    retained_IncLevelDifference = retained_IncLevelDifference,
    exclusion_IncLevelDifference = exclusion_IncLevelDifference,
    Min_Count = Min_Count
  )

  bin_width <- WidthIntoExon + WidthIntoIntron + 1

  # Process each group
  process_group <- function(data, group_name) {
    if (nrow(data) == 0) {
      if (verbose) message("No events found for group: ", group_name)
      return(data.frame(
        global_position = 1:(4 * bin_width),
        match_count = 0,
        frequency = 0,
        bin = rep(1:4, each = bin_width),
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
                                               bin_width
                                              )

    total_events <- length(unique(bins_gr$event_id))
    freq_data$frequency <- freq_data$match_count / total_events

    # Apply moving average using helper function
    freq_data <- calculate_moving_average(freq_data, moving_average, bins = bin_width)

    freq_data$group <- group_name
    return(freq_data)
  }

  # Process only selected groups
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

    # Calculate sample size based on retained + excluded counts
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
        global_position = 1:(4 * bin_width),
        match_count = 0,
        frequency = 0,
        bin = rep(1:4, each = bin_width),
        moving_avg = 0,
        group = "Control",
        moving_avg_sd = 0
      )
    } else if (sample_size >= n_controls || sample_size == 0) {
      # If sample size >= available controls, just use all controls once
      if (verbose) message("Using all controls without bootstrap")
      results_list$Control <- process_group(control_data, "Control")
      results_list$Control$moving_avg_sd <- 0
    } else {
      #Pre-compute motif matches for all control events once
      if (verbose) message("  Pre-computing sequence cache for all control events...")
      .report_progress(12, 100, "Building sequence cache...")

      cache_matrix <- precompute_sequence_cache(
        events_data = control_data,
        sequence = sequence,
        bsgenome_obj = genome,
        WidthIntoExon = WidthIntoExon,
        WidthIntoIntron = WidthIntoIntron,
        verbose = verbose
      )

      if (verbose) message("  Cache built. Running bootstrap iterations...")
      .report_progress(20, 100, "Running bootstrap iterations...")

      # Pre-generate all random samples
      all_sampled_indices <- lapply(seq_len(control_iterations), function(i) {
        sample(n_controls, sample_size, replace = FALSE)
      })

      # Helper to apply moving average to a frequency vector
      apply_moving_avg <- function(freq_vec) {
        if (is.null(moving_average) || moving_average <= 0) {
          return(freq_vec)
        }
        # Apply moving average per bin
        result <- numeric(length(freq_vec))
        for (b in 1:4) {
          start_idx <- (b - 1) * bin_width + 1
          end_idx <- b * bin_width
          bin_data <- freq_vec[start_idx:end_idx]
          result[start_idx:end_idx] <- slider::slide_dbl(
            bin_data,
            mean,
            .before = floor(moving_average / 2),
            .after = ceiling(moving_average / 2) - 1,
            .complete = FALSE
          )
        }
        result
      }

      # Bootstrap sampling using cached matrix
      iteration_results <- vector("list", control_iterations)

      pb <- progress::progress_bar$new(
        format = "  Sampling iterations [:bar] :current/:total (:percent) eta::eta",
        total = control_iterations, clear = FALSE, width = 80
      )

      loop_start <- 20
      loop_end <- 90

      for (iter in seq_len(control_iterations)) {
        pb$tick()

        sampled_ids <- all_sampled_indices[[iter]]
        match_counts <- rowSums(cache_matrix[, sampled_ids, drop = FALSE])
        freq_vec <- match_counts / sample_size
        iteration_results[[iter]] <- apply_moving_avg(freq_vec)

        progress_pct <- loop_start + (iter / control_iterations) * (loop_end - loop_start)
        .report_progress(
          progress_pct,
          100,
          sprintf("Control sampling iteration %d/%d", iter, control_iterations)
        )
      }

      # Combine results and calculate mean/sd
      freq_matrix <- do.call(cbind, iteration_results)
      mean_freq <- rowMeans(freq_matrix, na.rm = TRUE)
      sd_freq <- apply(freq_matrix, 1, stats::sd, na.rm = TRUE)

      results_list$Control <- data.frame(
        global_position = 1:(4 * bin_width),
        match_count = NA,
        frequency = mean_freq,
        bin = rep(1:4, each = bin_width),
        moving_avg = mean_freq,
        group = "Control",
        moving_avg_sd = sd_freq
      )

      # Store raw bootstrap matrix for diagnostics
      bootstrap_matrix <- freq_matrix
    }
  }

  # Combine selected groups
  combined_data <- dplyr::bind_rows(results_list)

  # Return data if requested
  if (return_data) {
    return(combined_data)
  }

  # Return diagnostics if requested (for normality testing)
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

  # Calculate significance if Control group is present and has SD
  sig_regions <- NULL
  if (show_significance && "Control" %in% groups) {
    control_has_sd <- any(combined_data$moving_avg_sd[combined_data$group == "Control"] > 0,
                          na.rm = TRUE)
    if (control_has_sd) {
      if (verbose) message("Calculating significance...")
      .report_progress(96, 100, "Calculating significance...")
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

  # Plot using the shared plotting function
  .report_progress(100, 100, "Complete")
  plot_splicing_sequence_map(combined_data,
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
