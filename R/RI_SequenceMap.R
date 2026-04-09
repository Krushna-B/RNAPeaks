#' Create Retained Intron Sequence Map
#'
#' Analyzes the frequency of a target sequence motif across retained intron
#' junction regions. Compares motif frequency between Retained, Excluded, and
#' Control events to identify position-specific enrichment patterns around the
#' upstream exon/intron and intron/downstream exon boundaries.
#'
#' @param RIMATS A data frame containing rMATS output with columns:
#'   chr, strand, upstreamES, upstreamEE, downstreamES, downstreamEE,
#'   GeneID, PValue, FDR, IncLevelDifference, IJC_SAMPLE_1, SJC_SAMPLE_1,
#'   IJC_SAMPLE_2, SJC_SAMPLE_2, IncLevel1, IncLevel2
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
#'   c("Retained", "Excluded", "Control") to process all groups.
#' @param control_multiplier Numeric multiplier for control sample size.
#'   Default is 2.0.
#' @param control_iterations Integer number for sampling iterations for control
#'   sampling. Default is 20.
#' @param cores Integer number of cores for parallel processing. Default is 1.
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

  .report_progress <- function(current, total, detail = NULL) {
    if (is.function(progress_callback)) {
      try(progress_callback(current, total, detail), silent = TRUE)
    }
  }

  # Cap cores at max available - 1
  max_cores <- parallel::detectCores() - 1
  if (is.na(max_cores) || max_cores < 1) max_cores <- 1
  cores <- min(cores, max_cores)
  cores <- max(cores, 1)
  options(future.globals.maxSize = 8 * 1024^3)

  warmup_future <- NULL
  if (cores > 1) {
    if (verbose) message(sprintf("Starting %d parallel workers...", cores))
    future::plan(future::multisession, workers = cores)
    warmup_future <- future::future({ TRUE }, seed = TRUE)
  }

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

  process_group <- function(data, group_name) {
    if (nrow(data) == 0) {
      if (verbose) message("No events found for group: ", group_name)
      return(data.frame(
        global_position = 1:(n_bins * bin_width),
        match_count = 0,
        frequency = 0,
        bin = rep(1:n_bins, each = bin_width),
        moving_avg = 0,
        group = group_name
      ))
    }

    bins_gr <- make_ri_bins_matrix(data,
                                    WidthIntoExon = WidthIntoExon,
                                    WidthIntoIntron = WidthIntoIntron)

    freq_data <- calculate_sequence_frequency(bins_gr,
                                               sequence,
                                               bsgenome_obj = genome,
                                               bin_width,
                                               n_bins = n_bins)

    total_events <- length(unique(bins_gr$event_id))
    freq_data$frequency <- freq_data$match_count / total_events

    freq_data <- calculate_moving_average(freq_data, moving_average, bins = bin_width)

    freq_data$group <- group_name
    return(freq_data)
  }

  if (!is.null(warmup_future)) {
    invisible(future::value(warmup_future))
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
        match_count = 0,
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
      if (verbose) message("  Pre-computing sequence cache for all control events...")
      .report_progress(12, 100, "Building sequence cache...")

      cache_matrix <- precompute_sequence_cache(
        events_data = control_data,
        sequence = sequence,
        bsgenome_obj = genome,
        WidthIntoExon = WidthIntoExon,
        WidthIntoIntron = WidthIntoIntron,
        verbose = verbose,
        n_bins = n_bins,
        make_bins_fn = make_ri_bins_matrix
      )

      if (verbose) message("  Cache built. Running bootstrap iterations...")
      .report_progress(20, 100, "Running bootstrap iterations...")

      all_sampled_indices <- lapply(seq_len(control_iterations), function(i) {
        sample(n_controls, sample_size, replace = FALSE)
      })

      apply_moving_avg <- function(freq_vec) {
        if (is.null(moving_average) || moving_average <= 0) return(freq_vec)
        result <- numeric(length(freq_vec))
        for (b in 1:n_bins) {
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

      iteration_results <- vector("list", control_iterations)
      loop_start <- 20
      loop_end <- 90

      for (iter in seq_len(control_iterations)) {
        sampled_ids <- all_sampled_indices[[iter]]
        match_counts <- rowSums(cache_matrix[, sampled_ids, drop = FALSE])
        freq_vec <- match_counts / sample_size
        iteration_results[[iter]] <- apply_moving_avg(freq_vec)

        if (iter %% 10 == 0 || iter == control_iterations) {
          progress_pct <- loop_start + (iter / control_iterations) * (loop_end - loop_start)
          .report_progress(
            progress_pct,
            100,
            sprintf("Bootstrap iteration %d/%d", iter, control_iterations)
          )
        }
      }

      freq_matrix <- do.call(cbind, iteration_results)
      mean_freq <- rowMeans(freq_matrix, na.rm = TRUE)
      sd_freq <- apply(freq_matrix, 1, stats::sd, na.rm = TRUE)

      results_list$Control <- data.frame(
        global_position = 1:(n_bins * bin_width),
        match_count = NA,
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

  # Clean up parallel workers
  if (cores > 1) {
    future::plan(future::sequential)
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
