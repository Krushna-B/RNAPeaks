
#' Analyzes the frequency of a target sequence motif across splicing junction
#' regions. Compares motif frequency between Retained, Excluded, and Control
#' splicing events to identify position-specific enrichment patterns.
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
#' @param control_multiplier Numeric multiplier for control sample size. The
#'   number of control events sampled per iteration is
#'   (n_retained + n_excluded) * control_multiplier. Default is 1.0.
#' @param control_iterations Integer number for sampling iterations for control
#'   sampling. The final control frequency is the mean across iterations, with
#'   standard deviation shown as a shaded band. Default is 20.
#' @param cores Integer number of cores for parallel processing. Default is 1
#'   (sequential). Set higher for faster processing on multi-core systems.
#' @param z_threshold Z-score threshold for significance testing. Default is 1.96
#'   (corresponds to p < 0.05 two-tailed).
#' @param min_consecutive Minimum number of consecutive significant positions
#'   required to form a significant region. Default is 10. Helps reduce false
#'   positives from noise.
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
#'
#' @return A ggplot object showing sequence frequency across the 4 regions
#'   for Retained (blue), Excluded (red), and Control (black) groups.
#'   Significant regions (z-test vs Control) are shown as colored bars above
#'   the plot. Returns a data frame if return_data = TRUE.
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
                               control_multiplier = 2.0,
                               control_iterations = 20,
                               cores = 1,
                               z_threshold = 1.96,
                               min_consecutive = 10,
                               show_significance = TRUE,
                               return_data = FALSE,
                               return_diagnostics = FALSE,
                               verbose = TRUE) {

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

  if (verbose) message("Processing groups: ", paste(groups, collapse = ", "))

  # Cap cores at max available - 1
  max_cores <- parallel::detectCores() - 1
  if (is.na(max_cores) || max_cores < 1) max_cores <- 1
  cores <- min(cores, max_cores)
  cores <- max(cores, 1)
  options(future.globals.maxSize = 8 * 1024^3)
  # If using parallel, warm up workers early while we do other setup
  warmup_future <- NULL
  if (cores > 1) {
    if (verbose) message(sprintf("Starting %d parallel workers...", cores))
    future::plan(future::multisession, workers = cores)

    # workers will spawn and load packages
    warmup_future <- future::future({
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
                                               bin_width,
                                              cores = cores)

    total_events <- length(unique(bins_gr$event_id))
    freq_data$frequency <- freq_data$match_count / total_events

    # Apply moving average using helper function
    freq_data <- calculate_moving_average(freq_data, moving_average, bins = bin_width)

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
    if (verbose) message("Processing Retained events...")
    results_list$Retained <- process_group(filtered_events$Retained, "Retained")
    results_list$Retained$moving_avg_sd <- 0
  }

  if ("Excluded" %in% groups) {
    if (verbose) message("Processing Excluded events...")
    results_list$Excluded <- process_group(filtered_events$Excluded, "Excluded")
    results_list$Excluded$moving_avg_sd <- 0
  }

  if ("Control" %in% groups) {
    if (verbose) message("Processing Control events with sampling...")

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
      # sampling
      iteration_results <- vector("list", control_iterations)

      pb <- progress::progress_bar$new(
        format = "  Sampling iterations [:bar] :current/:total (:percent) eta::eta",
        total = control_iterations, clear = FALSE, width = 80
      )

      for (iter in seq_len(control_iterations)) {
        pb$tick()

        # Random sample of controls
        sampled_indices <- sample(n_controls, sample_size, replace = FALSE)
        sampled_controls <- control_data[sampled_indices, ]

        # Process this sample
        iter_result <- process_group(sampled_controls, "Control")
        iteration_results[[iter]] <- iter_result$moving_avg
      }

      # Combine results and calculate mean/sd
      freq_matrix <- do.call(cbind, iteration_results)
      mean_freq <- rowMeans(freq_matrix, na.rm = TRUE)
      sd_freq <- apply(freq_matrix, 1, sd, na.rm = TRUE)

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

  # Clean up parallel workers
  if (cores > 1) {
    future::plan(future::sequential)
  }

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

  # Calculate significance if Control group is present and has SD
  sig_regions <- NULL
  if (show_significance && "Control" %in% groups) {
    control_has_sd <- any(combined_data$moving_avg_sd[combined_data$group == "Control"] > 0,
                          na.rm = TRUE)
    if (control_has_sd) {
      if (verbose) message("Calculating significance...")
      sig_result <- calculate_significance(
        combined_data,
        z_threshold = z_threshold,
        min_consecutive = min_consecutive,
        compare_to = "Control"
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
  plot_splicing_sequence_map(combined_data,
                             WidthIntoExon = WidthIntoExon,
                             WidthIntoIntron = WidthIntoIntron,
                             title = paste0("Sequence Frequency: ", sequence),
                             sig_regions = sig_regions)
}
