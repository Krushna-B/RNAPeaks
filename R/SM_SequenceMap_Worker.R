.sequence_map_worker <- function(mats_data,
                                  bins_fn,
                                  n_bins,
                                  plot_fn,
                                  genome,
                                  sequence,
                                  moving_average,
                                  WidthIntoExon,
                                  WidthIntoIntron,
                                  p_valueRetainedAndExclusion,
                                  p_valueControls,
                                  retained_IncLevelDifference,
                                  exclusion_IncLevelDifference,
                                  Min_Count,
                                  groups,
                                  control_multiplier,
                                  control_iterations,
                                  z_threshold,
                                  min_consecutive,
                                  one_sided,
                                  use_fdr,
                                  fdr_threshold,
                                  show_significance,
                                  return_data,
                                  return_diagnostics,
                                  verbose,
                                  progress_callback,
                                  title,
                                  retained_col,
                                  excluded_col,
                                  control_col,
                                  line_width,
                                  line_alpha,
                                  ribbon_alpha,
                                  title_size,
                                  title_color,
                                  axis_text_size,
                                  boundary_col,
                                  exon_col,
                                  legend_position,
                                  ylab) {

  if (is.null(genome)) {
    if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
      stop("BSgenome.Hsapiens.UCSC.hg38 is required. Install with:\n",
           "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')")
    }
    genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
  }

  valid_groups <- c("Retained", "Excluded", "Control")
  groups <- match.arg(groups, valid_groups, several.ok = TRUE)
  if (length(groups) == 0) stop("At least one group must be specified")

  if (verbose) message("Processing groups: ", paste(groups, collapse = ", "))

  .report_progress <- function(current, total, detail = NULL) {
    if (is.function(progress_callback)) {
      try(progress_callback(current, total, detail), silent = TRUE)
    }
  }

  filtered_events <- filter_SEMATS_events(
    mats_data,
    p_valueRetainedAndExclusion  = p_valueRetainedAndExclusion,
    p_valueControls              = p_valueControls,
    retained_IncLevelDifference  = retained_IncLevelDifference,
    exclusion_IncLevelDifference = exclusion_IncLevelDifference,
    Min_Count                    = Min_Count
  )

  bin_width <- WidthIntoExon + WidthIntoIntron + 1

  process_group <- function(data, group_name) {
    if (nrow(data) == 0) {
      if (verbose) message("No events found for group: ", group_name)
      return(data.frame(
        global_position = seq_len(n_bins * bin_width),
        match_count     = 0,
        frequency       = 0,
        bin             = rep(seq_len(n_bins), each = bin_width),
        moving_avg      = 0,
        group           = group_name,
        n_events        = 0L
      ))
    }
    n_events_val <- nrow(data)
    bins_gr   <- bins_fn(data, WidthIntoExon = WidthIntoExon, WidthIntoIntron = WidthIntoIntron)
    freq_data <- calculate_sequence_frequency(bins_gr, sequence,
                                               bsgenome_obj = genome,
                                               bin_width,
                                               n_bins = n_bins)
    freq_data$frequency <- freq_data$match_count / length(unique(bins_gr$event_id))
    freq_data <- calculate_moving_average(freq_data, moving_average, bins = bin_width)
    freq_data$group    <- group_name
    freq_data$n_events <- n_events_val
    freq_data
  }

  apply_moving_avg <- function(freq_vec) {
    if (is.null(moving_average) || moving_average <= 0) return(freq_vec)
    result <- numeric(length(freq_vec))
    for (b in seq_len(n_bins)) {
      start_idx <- (b - 1) * bin_width + 1
      end_idx   <- b * bin_width
      result[start_idx:end_idx] <- slider::slide_dbl(
        freq_vec[start_idx:end_idx], mean,
        .before = floor(moving_average / 2),
        .after  = ceiling(moving_average / 2) - 1,
        .complete = FALSE
      )
    }
    result
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
    n_controls   <- nrow(control_data)
    n_retained   <- nrow(filtered_events$Retained)
    n_excluded   <- nrow(filtered_events$Excluded)
    sample_size  <- round((n_retained + n_excluded) * control_multiplier)

    if (verbose) {
      message(sprintf("  Retained: %d, Excluded: %d, Sample size: %d",
                      n_retained, n_excluded, sample_size))
      message(sprintf("  Controls: %d, Iterations: %d", n_controls, control_iterations))
    }

    if (n_controls == 0) {
      if (verbose) message("No control events found")
      results_list$Control <- data.frame(
        global_position = seq_len(n_bins * bin_width),
        match_count     = 0,
        frequency       = 0,
        bin             = rep(seq_len(n_bins), each = bin_width),
        moving_avg      = 0,
        group           = "Control",
        moving_avg_sd   = 0,
        n_events        = 0L
      )
    } else if (sample_size >= n_controls || sample_size == 0) {
      if (verbose) message("Using all controls without bootstrap")
      results_list$Control <- process_group(control_data, "Control")
      results_list$Control$moving_avg_sd <- 0
    } else {
      if (verbose) message("  Pre-computing sequence cache for all control events...")
      .report_progress(12, 100, "Building sequence cache...")

      cache_matrix <- precompute_sequence_cache(
        events_data     = control_data,
        sequence        = sequence,
        bsgenome_obj    = genome,
        WidthIntoExon   = WidthIntoExon,
        WidthIntoIntron = WidthIntoIntron,
        verbose         = verbose,
        n_bins          = n_bins,
        make_bins_fn    = bins_fn
      )

      if (verbose) message("  Cache built. Running bootstrap iterations...")
      .report_progress(20, 100, "Running bootstrap iterations...")

      all_sampled_indices <- lapply(seq_len(control_iterations), function(i) {
        sample(n_controls, sample_size, replace = FALSE)
      })

      iteration_results <- vector("list", control_iterations)

      pb <- progress::progress_bar$new(
        format = "  Sampling iterations [:bar] :current/:total (:percent) eta::eta",
        total = control_iterations, clear = FALSE, width = 80
      )

      loop_start <- 20
      loop_end   <- 90

      for (iter in seq_len(control_iterations)) {
        pb$tick()
        sampled_ids            <- all_sampled_indices[[iter]]
        match_counts           <- rowSums(cache_matrix[, sampled_ids, drop = FALSE])
        iteration_results[[iter]] <- apply_moving_avg(match_counts / sample_size)
        .report_progress(
          loop_start + (iter / control_iterations) * (loop_end - loop_start),
          100,
          sprintf("Control sampling iteration %d/%d", iter, control_iterations)
        )
      }

      freq_matrix <- do.call(cbind, iteration_results)
      mean_freq   <- rowMeans(freq_matrix, na.rm = TRUE)
      sd_freq     <- apply(freq_matrix, 1, stats::sd, na.rm = TRUE)

      results_list$Control <- data.frame(
        global_position = seq_len(n_bins * bin_width),
        match_count     = NA,
        frequency       = mean_freq,
        bin             = rep(seq_len(n_bins), each = bin_width),
        moving_avg      = mean_freq,
        group           = "Control",
        moving_avg_sd   = sd_freq,
        n_events        = n_controls
      )

      bootstrap_matrix <- freq_matrix
    }
  }

  combined_data <- dplyr::bind_rows(results_list)

  if (return_data) return(combined_data)

  if (return_diagnostics) {
    return(list(
      data             = combined_data,
      bootstrap_matrix = if (exists("bootstrap_matrix")) bootstrap_matrix else NULL,
      n_iterations     = if (exists("bootstrap_matrix")) control_iterations else 0,
      sample_size      = if (exists("sample_size")) sample_size else NA,
      n_controls       = if (exists("n_controls")) n_controls else NA
    ))
  }

  .report_progress(92, 100, "Combining results...")

  sig_regions <- NULL
  if (show_significance && "Control" %in% groups) {
    control_has_sd <- any(
      combined_data$moving_avg_sd[combined_data$group == "Control"] > 0,
      na.rm = TRUE
    )
    if (control_has_sd) {
      if (verbose) message("Calculating significance...")
      .report_progress(96, 100, "Calculating significance...")
      sig_result  <- calculate_significance(
        combined_data,
        z_threshold     = z_threshold,
        min_consecutive = min_consecutive,
        compare_to      = "Control",
        one_sided       = one_sided,
        use_fdr         = use_fdr,
        fdr_threshold   = fdr_threshold
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
  plot_fn(combined_data,
          WidthIntoExon                = WidthIntoExon,
          WidthIntoIntron              = WidthIntoIntron,
          title                        = title,
          sig_regions                  = sig_regions,
          retained_cutoff              = retained_IncLevelDifference,
          excluded_cutoff              = exclusion_IncLevelDifference,
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
          ylab                         = ylab)
}
