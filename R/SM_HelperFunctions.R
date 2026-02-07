# SM_HelperFunctions.R
# Helper functions shared between Splicing Map and Sequence Map analyses

#' Filter SEMATS data into Controls, Retained, and Excluded events
#'
#' @param SEMATS A data frame containing SE.MATS output
#' @param p_valueRetainedAndExclusion P-value threshold for retained/excluded (default 0.05)
#' @param p_valueControls P-value threshold for controls (default 0.95)
#' @param retained_IncLevelDifference Threshold for retained events (default 0.1)
#' @param exclusion_IncLevelDifference Threshold for excluded events (default -0.1)
#' @param Min_Count Minimum read count threshold (default 50)
#'
#' @return A list with three data frames: Retained, Excluded, Control
#' @keywords internal
filter_SEMATS_events <- function(SEMATS,
                                  p_valueRetainedAndExclusion = 0.05,
                                  p_valueControls = 0.95,
                                  retained_IncLevelDifference = -0.1,
                                  exclusion_IncLevelDifference = 0.1,
                                  Min_Count = 50) {

  # Calculate counts from junction counts
  SEMATS$IN_Count_1 <- sapply(SEMATS$IJC_SAMPLE_1, function(x) {
    sum(as.numeric(strsplit(as.character(x), ",")[[1]]))
  })
  SEMATS$SK_Count_1 <- sapply(SEMATS$SJC_SAMPLE_1, function(x) {
    sum(as.numeric(strsplit(as.character(x), ",")[[1]]))
  })
  SEMATS$IN_Count_2 <- sapply(SEMATS$IJC_SAMPLE_2, function(x) {
    sum(as.numeric(strsplit(as.character(x), ",")[[1]]))
  })
  SEMATS$SK_Count_2 <- sapply(SEMATS$SJC_SAMPLE_2, function(x) {
    sum(as.numeric(strsplit(as.character(x), ",")[[1]]))
  })

  SEMATS$Total_Count_1 <- SEMATS$IN_Count_1 + SEMATS$SK_Count_1
  SEMATS$Total_Count_2 <- SEMATS$IN_Count_2 + SEMATS$SK_Count_2

  # Calculate inclusion levels
  SEMATS$Inc_1 <- sapply(SEMATS$IncLevel1, function(x) {
    mean(as.numeric(strsplit(as.character(x), ",")[[1]]), na.rm = TRUE)
  })
  SEMATS$Inc_2 <- sapply(SEMATS$IncLevel2, function(x) {
    mean(as.numeric(strsplit(as.character(x), ",")[[1]]), na.rm = TRUE)
  })

  # Filter by minimum count
  SEMATS <- SEMATS[which(SEMATS$Total_Count_1 > Min_Count &
                           SEMATS$Total_Count_2 > Min_Count), ]

  # Generate event categories
  Excluded <- SEMATS %>%
    dplyr::filter(PValue < p_valueRetainedAndExclusion,
                  FDR < p_valueRetainedAndExclusion,
                  IncLevelDifference > exclusion_IncLevelDifference)

  Retained <- SEMATS %>%
    dplyr::filter(PValue < p_valueRetainedAndExclusion,
                  FDR < p_valueRetainedAndExclusion,
                  IncLevelDifference < retained_IncLevelDifference)

  Control <- SEMATS %>%
    dplyr::filter(PValue > p_valueControls,
                  FDR > p_valueControls,
                  Inc_1 > 0.9,
                  abs(IncLevelDifference) < 0.005)

  return(list(
    Retained = Retained,
    Excluded = Excluded,
    Control = Control
  ))
}


#' Create bins matrix for splicing events
#' Divides each splicing event into 4 bins around exon/intron boundaries.
#' @param MAT A data frame with SE.MATS columns. If a 'group' column exists,
#'   it will be preserved in the output GRanges.
#' @param WidthIntoExon Width to extend into exons (default 50)
#' @param WidthIntoIntron Width to extend into introns (default 250)
#'
#' @return A GRanges object with bins for all events
#' @keywords internal
make_bins_matrix <- function(MAT, WidthIntoExon, WidthIntoIntron) {
  n <- nrow(MAT)

  # Renaming variables for clarity
  exonStart <- MAT$exonStart_0base
  exonEnd   <- MAT$exonEnd
  upS <- MAT$upstreamES
  upE <- MAT$upstreamEE
  downS <- MAT$downstreamES
  downE <- MAT$downstreamEE

  # Bin1: Upstream exon end to first intron
  bin1_start <- pmax(upE - WidthIntoExon, upS)
  bin1_end <- pmin(upE + WidthIntoIntron, exonStart)

  # Bin2: First intron end to middle (skipped) exon start
  bin2_start <- pmax(exonStart - WidthIntoIntron, upE)
  bin2_end <- pmin(exonStart + WidthIntoExon, exonEnd)

  # Bin3: Middle exon end to second intron
  bin3_start <- pmax(exonEnd - WidthIntoExon, exonStart)
  bin3_end <- pmin(exonEnd + WidthIntoIntron, downS)

  # Bin4: Second intron end to downstream exon start
  bin4_start <- pmax(downS - WidthIntoIntron, exonEnd)
  bin4_end <- pmin(downS + WidthIntoExon, downE)

  starts <- cbind(bin1_start, bin2_start, bin3_start, bin4_start)
  ends   <- cbind(bin1_end,   bin2_end,   bin3_end,   bin4_end)

  # Build GRanges
  gr <- GenomicRanges::GRanges(
    seqnames = rep(MAT$chr, each = 4),
    ranges   = IRanges::IRanges(start = as.vector(t(starts)),
                                end   = as.vector(t(ends))),
    strand   = rep(MAT$strand, each = 4),
    geneID = rep(MAT$GeneID, each = 4),
    event_id = rep(seq_len(nrow(MAT)), each = 4)
  )

  # Preserve group column if present
  if ("group" %in% colnames(MAT)) {
    GenomicRanges::mcols(gr)$group <- rep(MAT$group, each = 4)
  }
  gr
}


#' Calculate binding frequency using vectorized batch
#'
#' Calculates protein binding frequency across all bins using
#' vectorized operations and a single findOverlaps call.
#'
#' @param bins_gr GRanges object from make_bins_matrix. If it has a 'group'
#'   metadata column, results are returned per group.
#' @param protein GRanges object of protein binding sites
#' @param bin_width Width of each bin region
#' @param cores Number of cores for parallel processing (default 1, Sync)
#'
#' @return Data frame with global_position and overlap_count columns.
#'   If bins_gr has a 'group' column, includes group column.
#' @keywords internal
calculate_binding_frequency <- function(bins_gr, protein, bin_width, cores = 1) {

  n_bins <- length(bins_gr)
  total_positions <- 4 * bin_width

  # Check for group metadata
  has_groups <- "group" %in% names(GenomicRanges::mcols(bins_gr))
  if (has_groups) {
    bins_groups <- as.character(GenomicRanges::mcols(bins_gr)$group)
    unique_groups <- unique(bins_groups)
    group_to_idx <- setNames(seq_along(unique_groups), unique_groups)
  }

  # Pre-filter bins that overlap protein ranges
  overlapping_idx <- which(IRanges::overlapsAny(bins_gr, protein))

  if (length(overlapping_idx) == 0) {
    message("No bins overlap protein binding sites")
    if (has_groups) {
      return(do.call(rbind, lapply(unique_groups, function(g) {
        data.frame(global_position = seq_len(total_positions), overlap_count = 0L, group = g)
      })))
    }
    return(data.frame(global_position = seq_len(total_positions), overlap_count = 0L))
  }

  # Subset to overlapping bins only
  overlapping_bins <- bins_gr[overlapping_idx]
  n_overlapping <- length(overlapping_bins)

  message(sprintf("Processing %d bins that overlap protein sites...", n_overlapping))

  # Pre-compute bin metadata

  bins_chr <- as.character(GenomicRanges::seqnames(overlapping_bins))
  bins_start <- GenomicRanges::start(overlapping_bins)
  bins_end <- GenomicRanges::end(overlapping_bins)
  bins_strand <- as.character(GenomicRanges::strand(overlapping_bins))
  bins_event_id <- GenomicRanges::mcols(overlapping_bins)$event_id

  # Calculate bin indices (1-4) based on original position in bins_gr
  original_indices <- overlapping_idx
  bin_indices <- ((original_indices - 1) %% 4) + 1

  # Filter out invalid bins (start > end)
  valid_mask <- bins_start <= bins_end
  if (!all(valid_mask)) {
    overlapping_bins <- overlapping_bins[valid_mask]
    bins_chr <- bins_chr[valid_mask]
    bins_start <- bins_start[valid_mask]
    bins_end <- bins_end[valid_mask]
    bins_strand <- bins_strand[valid_mask]
    bins_event_id <- bins_event_id[valid_mask]
    bin_indices <- bin_indices[valid_mask]
    if (has_groups) {
      bins_groups_subset <- bins_groups[overlapping_idx][valid_mask]
    }
    n_overlapping <- length(overlapping_bins)
  } else if (has_groups) {
    bins_groups_subset <- bins_groups[overlapping_idx]
  }

  if (n_overlapping == 0) {
    if (has_groups) {
      return(do.call(rbind, lapply(unique_groups, function(g) {
        data.frame(global_position = seq_len(total_positions), overlap_count = 0L, group = g)
      })))
    }
    return(data.frame(global_position = seq_len(total_positions), overlap_count = 0L))
  }

  # Calculate bin widths and total positions to expand

  bin_widths <- bins_end - bins_start + 1L
  total_1bp_positions <- sum(bin_widths)

  message(sprintf("Expanding %d bins to %d 1bp positions...", n_overlapping, total_1bp_positions))

  # Vectorized expansion: create all 1bp positions at once
  # For each bin, we need positions from start to end
  bin_rep_idx <- rep(seq_len(n_overlapping), times = bin_widths)

  # Calculate genomic positions using vectorized cumsum trick
  # Create sequence 1,2,3,...,w1, 1,2,3,...,w2, etc.
  position_in_bin <- sequence(bin_widths)
  genomic_positions <- bins_start[bin_rep_idx] + position_in_bin - 1L

  # Create 1bp GRanges for all positions at once
  all_positions_gr <- GenomicRanges::GRanges(
    seqnames = bins_chr[bin_rep_idx],
    ranges = IRanges::IRanges(start = genomic_positions, width = 1L),
    strand = bins_strand[bin_rep_idx]
  )

  message("Finding overlaps with protein binding sites...")

  # Single findOverlaps call for ALL positions
  hits <- GenomicRanges::findOverlaps(all_positions_gr, protein, ignore.strand = FALSE)

  # Get indices of positions that have overlaps
  hit_query_idx <- S4Vectors::queryHits(hits)

  if (length(hit_query_idx) == 0) {
    message("No position overlaps found")
    if (has_groups) {
      return(do.call(rbind, lapply(unique_groups, function(g) {
        data.frame(global_position = seq_len(total_positions), overlap_count = 0L, group = g)
      })))
    }
    return(data.frame(global_position = seq_len(total_positions), overlap_count = 0L))
  }

  # Get unique hit positions (a position might overlap multiple protein regions)
  hit_query_idx <- unique(hit_query_idx)

  message(sprintf("Found %d positions with overlaps, calculating global positions...", length(hit_query_idx)))

  # For hit positions, calculate global positions
  hit_bin_idx <- bin_rep_idx[hit_query_idx]
  hit_position_in_bin <- position_in_bin[hit_query_idx]
  hit_bin_indices <- bin_indices[hit_bin_idx]
  hit_strands <- bins_strand[hit_bin_idx]

  # Calculate global positions (strand-aware)
  # Plus strand: (bin_index - 1) * bin_width + position_in_bin
  # Minus strand: (4 - bin_index) * bin_width + (bin_width - position_in_bin + 1)
  is_minus <- hit_strands == "-"
  global_pos <- integer(length(hit_query_idx))
  global_pos[!is_minus] <- (hit_bin_indices[!is_minus] - 1L) * bin_width + hit_position_in_bin[!is_minus]
  global_pos[is_minus] <- (4L - hit_bin_indices[is_minus]) * bin_width + (bin_width - hit_position_in_bin[is_minus] + 1L)

  # Clamp to valid range

  global_pos <- pmax(1L, pmin(global_pos, total_positions))

  # Aggregate results
  if (has_groups) {
    hit_groups <- bins_groups_subset[hit_bin_idx]
    hit_group_indices <- group_to_idx[hit_groups]
    n_groups <- length(unique_groups)

    # Combined index for tabulation: (group_idx - 1) * total_positions + global_pos
    combined_idx <- (hit_group_indices - 1L) * total_positions + global_pos
    counts <- tabulate(combined_idx, nbins = n_groups * total_positions)

    result <- data.frame(
      global_position = rep(seq_len(total_positions), n_groups),
      overlap_count = counts,
      group = rep(unique_groups, each = total_positions)
    )
  } else {
    counts <- tabulate(global_pos, nbins = total_positions)
    result <- data.frame(
      global_position = seq_len(total_positions),
      overlap_count = counts
    )
  }

  return(result)
}


#' Find overlaps between bins and protein binding sites (legacy)
#' @param bins GRanges object of bins
#' @param protein GRanges object of protein binding sites
#'
#' @return Data frame with overlap information for each position
#' @keywords internal
find_overlaps <- function(bins, protein) {
  # Pre-filter - only check bins that actually overlap protein ranges
  overlapping_bins <- subsetByOverlaps(bins, protein)

  if (length(overlapping_bins) == 0) {
    return(data.frame(
      event_id = integer(),
      bin_index = integer(),
      position = integer(),
      has_overlap = logical()
    ))
  }

  # Progress Bar
  n <- length(overlapping_bins)
  pb <- progress::progress_bar$new(
    format = "  finding overlaps [:bar] :current/:total (:percent) eta::eta",
    total = n, clear = FALSE, width = 120
  )

  results <- lapply(seq_along(overlapping_bins), function(i) {
    pb$tick()
    bin <- overlapping_bins[i]
    if (GenomicRanges::start(bin) > GenomicRanges::end(bin)) {
      return(NULL)
    }
    chr <- as.character(GenomicRanges::seqnames(bin))
    chr_protein <- protein[GenomicRanges::seqnames(protein) == chr]
    positions <- GenomicRanges::start(bin):GenomicRanges::end(bin)

    # Every 1 bp position
    pos_gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(positions, width = 1),
                                      strand = GenomicRanges::strand(bin))
    overlaps <- GenomicRanges::countOverlaps(pos_gr, chr_protein, ignore.strand = FALSE) > 0
    data.frame(
      event_id = GenomicRanges::mcols(bin)$event_id,
      geneID = GenomicRanges::mcols(bin)$geneID,
      seqnames = chr,
      strand = as.character(GenomicRanges::strand(bin)),
      bin_index = ((i - 1) %% 4) + 1,  # 1-4 for each event
      position = seq_along(positions),
      genomic_position = positions,
      has_overlap = overlaps
    )
  })
  do.call(rbind, results)
}


#' Calculate overlap frequency at each position (legacy)
#'
#' @param overlap_df Data frame from find_overlaps
#' @param total_events Total number of events
#' @param bin_width Width of each bin
#'
#' @return Data frame with global_position and frequency columns
#' @keywords internal
calculate_overlap_frequency <- function(overlap_df, total_events, bin_width) {

  # Flips Global Position for minus strand
  overlap_df <- overlap_df %>%
    dplyr::group_by(event_id) %>%
    dplyr::mutate(
      global_position = (bin_index - 1) * bin_width + position,
      global_position = ifelse(
        strand == "-",
        (4 * bin_width) - global_position + 1,
        global_position
      )
    ) %>%
    dplyr::ungroup()

  # Count overlaps at each global position
  overlap_counts <- stats::aggregate(
    has_overlap ~ global_position,
    data = overlap_df,
    FUN = sum
  )
  # Calculate frequency (overlaps / TOTAL events in dataset)

  overlap_counts$frequency <- overlap_counts$has_overlap / total_events

  # Ensure all positions are represented (fill missing with 0)
  all_positions <- data.frame(global_position = 1:(4 * bin_width))
  result <- merge(all_positions, overlap_counts[, c("global_position", "frequency")],
                  by = "global_position", all.x = TRUE)
  result$frequency[is.na(result$frequency)] <- 0

  return(result)
}


#' Calculate moving average of frequency data
#' @param freq_data Data frame with global_position and frequency columns
#' @param window_size Window size for moving average (NULL or 0 to disable)
#' @param bins Width of each bin region
#' @return Data frame with additional moving_avg column
#' @keywords internal
calculate_moving_average <- function(freq_data, window_size = NULL, bins = NULL) {

  # Add bin column
  freq_data <- freq_data %>%
    dplyr::mutate(
      bin = dplyr::case_when(
        global_position <= bins ~ 1,
        global_position <= bins * 2 ~ 2,
        global_position <= bins * 3 ~ 3,
        TRUE ~ 4
      )
    ) %>%
    dplyr::arrange(global_position)

  # Apply moving average only if window_size is provided
  if (!is.null(window_size) && window_size > 0) {

    freq_data <- freq_data %>%
      dplyr::arrange(bin, global_position) %>%
      dplyr::group_by(bin) %>%
      dplyr::mutate(moving_avg = slider::slide_dbl(frequency,
                                                    mean,
                                                    .before = floor((window_size - 1) / 2),
                                                    .after = floor((window_size - 1) / 2),
                                                    .complete = FALSE)) %>%
      dplyr::ungroup()
  } else {
    # No moving average - just copy frequency to moving_avg
    freq_data <- freq_data %>%
      dplyr::mutate(moving_avg = frequency)
  }
  return(freq_data)
}


#' Calculate sequence motif frequency across bins
#'
#' @param bins_gr GRanges object from make_bins_matrix
#' @param sequence Target sequence motif
#' @param bsgenome_obj BSgenome object
#' @param bin_width Width of each bin
#' @param cores Number of cores for parallel processing (default 1 = sequential).
#'   Caller is responsible for validating/capping this value.
#'
#' @return Data frame with global_position and match_count columns
#' @keywords internal
calculate_sequence_frequency <- function(bins_gr, sequence, bsgenome_obj, bin_width, cores = 1) {

  seq_length <- nchar(sequence)

  # Convert sequence to DNAString for matching (handles IUPAC codes)
  pattern <- Biostrings::DNAString(sequence)

  # Force evaluation to avoid scoping issues
  force(bsgenome_obj)
  n <- length(bins_gr)

  # Pre-compute all metadata vectors
  bins_chr <- as.character(GenomicRanges::seqnames(bins_gr))
  bins_start <- GenomicRanges::start(bins_gr)
  bins_end <- GenomicRanges::end(bins_gr)
  bins_strand <- as.character(GenomicRanges::strand(bins_gr))
  bin_indices <- ((seq_len(n) - 1) %% 4) + 1
  bin_lengths <- bins_end - bins_start + 1

  # Calculate extended ends for sequence extraction
  # (extended to allow matching at the end of each bin)
  chr_lengths <- GenomeInfoDb::seqlengths(bsgenome_obj)
  extended_ends <- pmin(bins_end + seq_length - 1, chr_lengths[bins_chr])

  # Filter out bins with invalid chromosomes (NA in chr_lengths)
  valid_idx <- which(!is.na(extended_ends) & extended_ends >= bins_start)

  if (length(valid_idx) == 0) {
    warning("No sequences could be extracted. Check chromosome naming convention.")
    return(data.frame(global_position = 1:(4 * bin_width), match_count = 0))
  }

  # Create GRanges for BATCH sequence extraction
  extended_gr <- GenomicRanges::GRanges(
    seqnames = bins_chr[valid_idx],
    ranges = IRanges::IRanges(start = bins_start[valid_idx], end = extended_ends[valid_idx]),
    strand = bins_strand[valid_idx]
  )

  # Batch getSeq
  message("Extracting ", length(extended_gr), " sequences in batch...")
  all_seqs <- tryCatch({
    Biostrings::getSeq(bsgenome_obj, extended_gr)
  }, error = function(e) {
    warning("Batch sequence extraction failed: ", e$message)
    return(NULL)
  })

  if (is.null(all_seqs)) {
    warning("No sequences could be extracted. Check chromosome naming convention.")
    return(data.frame(global_position = 1:(4 * bin_width), match_count = 0))
  }


  # Subset metadata to valid indices
  valid_bin_indices <- bin_indices[valid_idx]
  valid_bin_lengths <- bin_lengths[valid_idx]
  valid_strands <- bins_strand[valid_idx]

  # Get all pattern matches
  message("Finding pattern matches...")
  n_valid <- length(all_seqs)

  hits_all <- Biostrings::vmatchPattern(pattern, all_seqs, fixed = FALSE)

  starts_all <- BiocGenerics::start(hits_all)
  # Convert S4 IntegerList to regular R list for parallel serialization
  starts_all <- as.list(starts_all)

  total_positions <- 4 * bin_width

  # Process and collect global positions of matches
  if (cores > 1) {
    # Parallel processing - each chunk returns global positions of matches
    # Dynamic chunk size
    chunk_size <- max(500, ceiling(n_valid / (cores * 2)))
    idx_chunks <- split(seq_len(n_valid), ceiling(seq_len(n_valid) / chunk_size))
    n_chunks <- length(idx_chunks)

    message(sprintf("Processing %d bins in %d chunks using %d cores...", n_valid, n_chunks, cores))

    # Use progressr for parallel progress reporting
    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")

    all_match_positions <- progressr::with_progress({
      p <- progressr::progressor(steps = 100)

      chunk_results <- future.apply::future_lapply(idx_chunks, function(chunk_indices) {
        match_positions <- vector("list", length(chunk_indices))

        for (pos in seq_along(chunk_indices)) {
          idx <- chunk_indices[pos]
          bin_length <- valid_bin_lengths[idx]
          starts <- starts_all[[idx]]
          starts <- starts[starts <= bin_length]

          if (length(starts) > 0) {
            bin_idx <- valid_bin_indices[idx]
            strand <- valid_strands[idx]

            # Calculate global positions (strand-aware)
            if (strand == "+") {
              global_pos <- (bin_idx - 1) * bin_width + starts
            } else {
              global_pos <- (4 - bin_idx) * bin_width + (bin_width - starts + 1)
            }
            match_positions[[pos]] <- global_pos
          }
        }

        p()  # Signal chunk completed
        unlist(match_positions)
      }, future.seed = TRUE)

      unlist(chunk_results)
    })

    message("Done.")

  } else {
    # Sequential processing with progress bar
    pb <- progress::progress_bar$new(
      format = "  processing [:bar] :current/:total (:percent) eta::eta",
      total = n_valid, clear = FALSE, width = 120
    )

    match_positions <- vector("list", n_valid)

    for (i in seq_len(n_valid)) {
      pb$tick()

      bin_length <- valid_bin_lengths[i]
      starts <- starts_all[[i]]
      starts <- starts[starts <= bin_length]

      if (length(starts) > 0) {
        bin_idx <- valid_bin_indices[i]
        strand <- valid_strands[i]

        # Calculate global positions (strand-aware)
        if (strand == "+") {
          global_pos <- (bin_idx - 1) * bin_width + starts
        } else {
          global_pos <- (4 - bin_idx) * bin_width + (bin_width - starts + 1)
        }
        match_positions[[i]] <- global_pos
      }
    }

    all_match_positions <- unlist(match_positions)
  }

  # Use tabulate to count matches at each position
  if (length(all_match_positions) > 0) {
    match_count <- tabulate(all_match_positions, nbins = total_positions)
  } else {
    match_count <- integer(total_positions)
  }

  # Create final data frame only once at the end
  data.frame(
    global_position = seq_len(total_positions),
    match_count = match_count
  )
}


#' Calculate Z-scores and Find Significant Regions
#'
#' Compares Retained and Excluded frequencies against Control using z-scores.
#' Identifies positions where the difference is statistically significant.
#'
#' @param freq_data Data frame with columns: global_position, moving_avg, group,
#'   and moving_avg_sd (for Control group from bootstrap)
#' @param z_threshold Z-score threshold for significance (default 1.96 for p < 0.05)
#' @param min_consecutive Minimum consecutive significant positions to form a region
#'   (default 10). Helps reduce false positives from noise.
#' @param compare_to Which group to compare against. Default "Control".
#'
#' @return A list with:
#'   \item{z_scores}{Data frame with z-scores for each position and comparison}
#'   \item{significant_regions}{Data frame with start/end of significant regions}
#'
#' @export
calculate_significance <- function(freq_data,
                                    z_threshold = 1.96,
                                    min_consecutive = 10,
                                    compare_to = "Control") {

  # Get Control data (must have SD from bootstrap)
  control_data <- freq_data %>%
    dplyr::filter(group == compare_to) %>%
    dplyr::select(global_position, moving_avg, moving_avg_sd) %>%
    dplyr::rename(control_mean = moving_avg, control_sd = moving_avg_sd)

  if (nrow(control_data) == 0) {
    stop("No Control group found in freq_data")
  }

  if (all(control_data$control_sd == 0, na.rm = TRUE)) {
    warning("Control SD is all zeros. Cannot calculate z-scores. ",
            "Make sure bootstrap iterations > 1 in createSequenceMap.")
    return(list(z_scores = NULL, significant_regions = NULL))
  }

  # Get other groups to compare
  other_groups <- unique(freq_data$group[freq_data$group != compare_to])

  z_scores_list <- list()
  sig_regions_list <- list()

  for (grp in other_groups) {
    grp_data <- freq_data %>%
      dplyr::filter(group == grp) %>%
      dplyr::select(global_position, moving_avg) %>%
      dplyr::rename(grp_freq = moving_avg)

    # Merge with control
    merged <- dplyr::left_join(grp_data, control_data, by = "global_position")

    # Calculate z-scores (handle zero SD)
    merged <- merged %>%
      dplyr::mutate(
        z_score = dplyr::if_else(
          control_sd > 0,
          (grp_freq - control_mean) / control_sd,
          0
        ),
        significant = abs(z_score) >= z_threshold,
        direction = dplyr::case_when(
          z_score >= z_threshold ~ "enriched",
          z_score <= -z_threshold ~ "depleted",
          TRUE ~ "ns"
        ),
        comparison = paste0(grp, "_vs_", compare_to)
      )

    z_scores_list[[grp]] <- merged

    # Find consecutive significant regions
    if (any(merged$significant, na.rm = TRUE)) {
      # Use run-length encoding to find consecutive runs
      rle_sig <- rle(merged$significant)
      ends <- cumsum(rle_sig$lengths)
      starts <- c(1, ends[-length(ends)] + 1)

      # Filter for significant runs that meet minimum length
      sig_runs <- which(rle_sig$values == TRUE & rle_sig$lengths >= min_consecutive)

      if (length(sig_runs) > 0) {
        regions <- data.frame(
          start_pos = merged$global_position[starts[sig_runs]],
          end_pos = merged$global_position[ends[sig_runs]],
          length = rle_sig$lengths[sig_runs],
          comparison = paste0(grp, "_vs_", compare_to),
          group = grp
        )

        # Add mean z-score for each region
        regions$mean_z <- sapply(seq_len(nrow(regions)), function(i) {
          idx <- merged$global_position >= regions$start_pos[i] &
                 merged$global_position <= regions$end_pos[i]
          mean(merged$z_score[idx], na.rm = TRUE)
        })

        regions$direction <- ifelse(regions$mean_z > 0, "enriched", "depleted")

        sig_regions_list[[grp]] <- regions
      }
    }
  }

  # Combine results
  z_scores <- dplyr::bind_rows(z_scores_list)
  significant_regions <- dplyr::bind_rows(sig_regions_list)

  return(list(
    z_scores = z_scores,
    significant_regions = significant_regions
  ))
}


#' Add region labels to combined data
#'
#' @param Combined Data frame with Pos column
#' @param WidthIntoExon Width into exon
#' @param WidthIntoIntron Width into intron
#'
#' @return Data frame with Region column added
#' @keywords internal
add_regions <- function(Combined, WidthIntoExon = 50, WidthIntoIntron = 300) {
  bin_width <- WidthIntoExon + WidthIntoIntron
  Combined$Region <- NA
  Combined$Region[which(Combined$Pos <= WidthIntoExon)] <- "UE"
  Combined$Region[which(Combined$Pos >= (WidthIntoExon + 1) &
                          Combined$Pos <= bin_width)] <- "UI5"
  Combined$Region[which(Combined$Pos >= (bin_width + 1) &
                          Combined$Pos <= (bin_width + WidthIntoIntron))] <- "UI3"
  Combined$Region[which(Combined$Pos >= (bin_width + WidthIntoIntron + 1) &
                          Combined$Pos <= (2 * bin_width))] <- "EX3"
  Combined$Region[which(Combined$Pos >= (2 * bin_width + 1) &
                          Combined$Pos <= (2 * bin_width + WidthIntoExon))] <- "EX5"
  Combined$Region[which(Combined$Pos >= (2 * bin_width + WidthIntoExon + 1) &
                          Combined$Pos <= (3 * bin_width))] <- "DI5"
  Combined$Region[which(Combined$Pos >= (3 * bin_width + 1) &
                          Combined$Pos <= (3 * bin_width + WidthIntoIntron))] <- "DI3"
  Combined$Region[which(Combined$Pos >= (3 * bin_width + WidthIntoIntron + 1) &
                          Combined$Pos <= (4 * bin_width))] <- "DE"
  return(Combined)
}


#' Plot Splicing/Sequence Map
#'
#' Creates a visualization of frequency data across splicing regions
#' for Retained, Excluded, and Control groups. Used by both createSplicingMap
#' and createSequenceMap.
#'
#' @param freq_data Data frame with global_position, frequency, moving_avg, and group
#' @param WidthIntoExon Width into exon (default 50)
#' @param WidthIntoIntron Width into intron (default 250)
#' @param title Plot title
#' @param sig_regions Data frame with significant regions from calculate_significance.
#'   Should have columns: start_pos, end_pos, group. If NULL, no significance bars shown.
#'
#' @return A ggplot object
#' @keywords internal
plot_splicing_sequence_map <- function(freq_data,
                                        WidthIntoExon = 50,
                                        WidthIntoIntron = 250,
                                        title = NULL,
                                        sig_regions = NULL) {

  bin_width <- WidthIntoExon + WidthIntoIntron
  gap <- 80

  # Ensure bin column exists
  if (!"bin" %in% names(freq_data)) {
    freq_data <- freq_data %>%
      dplyr::mutate(
        bin = dplyr::case_when(
          global_position <= bin_width ~ 1,
          global_position <= 2 * bin_width ~ 2,
          global_position <= 3 * bin_width ~ 3,
          TRUE ~ 4
        )
      )
  }

  # Handle empty or all-NA data
  valid_avg <- freq_data$moving_avg[!is.na(freq_data$moving_avg)]
  if (length(valid_avg) == 0) {
    warning("No valid moving_avg values. Using frequency instead.")
    freq_data$moving_avg <- freq_data$frequency
    valid_avg <- freq_data$moving_avg[!is.na(freq_data$moving_avg)]
  }

  # Height for exons - position at bottom of plot
  y_max <- if (length(valid_avg) > 0) max(valid_avg) else 0.01
  data_min <- if (length(valid_avg) > 0) min(valid_avg) else 0
  y_range <- y_max - data_min
  if (y_range == 0) y_range <- y_max * 0.1
  exon_height <- y_range * 0.08

  # Place exon schematic below the data
  y_min <- min(0, data_min) - y_range * 0.05

  # Calculate region starts
  region_starts <- c(0, bin_width + gap, 2 * bin_width + gap, 3 * bin_width + 2 * gap)

  boundary1 <- region_starts[1] + WidthIntoExon
  boundary2 <- region_starts[2] + WidthIntoIntron
  boundary3 <- region_starts[3] + WidthIntoExon
  boundary4 <- region_starts[4] + WidthIntoIntron

  # Dashed boundary lines locations
  boundary_lines <- data.frame(
    xintercept = c(boundary1, boundary2, boundary3, boundary4)
  )

  exon_regions <- data.frame(
    xmin = c(region_starts[1], boundary2, boundary4),
    xmax = c(boundary1, boundary3, region_starts[4] + bin_width),
    ymin = rep(y_min - exon_height, 3),
    ymax = rep(y_min, 3),
    fill = c("white", "navy", "white")
  )

  # Intron lines
  intron_y <- y_min - exon_height / 2
  intron_segments <- data.frame(
    x = c(boundary1, boundary3),
    xend = c(boundary2, boundary4),
    y = rep(intron_y, 2),
    yend = rep(intron_y, 2),
    linetype = c("solid", "solid")
  )

  # Breaks
  break_x <- c(bin_width + gap / 2, 3 * bin_width + gap + gap / 2)

  # Align data to proper region
  freq_data <- freq_data %>%
    dplyr::mutate(
      position_in_bin = global_position - (bin - 1) * (bin_width + 1),
      schematic_position = dplyr::case_when(
        bin == 1 ~ position_in_bin,
        bin == 2 ~ bin_width + gap + position_in_bin,
        bin == 3 ~ 2 * bin_width + gap + position_in_bin,
        bin == 4 ~ 3 * bin_width + 2 * gap + position_in_bin
      )
    )

  # Set factor levels for group (only include groups present in data)
  present_groups <- intersect(c("Retained", "Excluded", "Control"), unique(freq_data$group))
  freq_data$group <- factor(freq_data$group, levels = present_groups)

  # Add ribbon data for groups with standard deviation
  if ("moving_avg_sd" %in% names(freq_data)) {
    freq_data <- freq_data %>%
      dplyr::mutate(
        ymin = moving_avg - moving_avg_sd,
        ymax = moving_avg + moving_avg_sd
      )
  }

  plot <- ggplot2::ggplot(freq_data,
                          ggplot2::aes(x = schematic_position,
                                       y = moving_avg,
                                       color = group,
                                       group = interaction(bin, group)))

  # Add ribbon for standard deviation if available
  if ("moving_avg_sd" %in% names(freq_data)) {
    ribbon_data <- freq_data %>%
      dplyr::filter(moving_avg_sd > 0)
    if (nrow(ribbon_data) > 0) {
      # Map group names to colors for ribbon fill
      group_colors <- c("Retained" = "blue", "Excluded" = "red", "Control" = "gray50")
      ribbon_data$ribbon_fill <- group_colors[as.character(ribbon_data$group)]

      plot <- plot +
        ggplot2::geom_ribbon(data = ribbon_data,
                             ggplot2::aes(x = schematic_position,
                                          ymin = ymin, ymax = ymax,
                                          group = interaction(bin, group),
                                          fill = ribbon_fill),
                             alpha = 0.3, color = NA, inherit.aes = FALSE)
    }
  }

  plot <- plot +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = c("Retained" = "blue", "Excluded" = "red", "Control" = "black"),
      name = "Event Type"
    ) +

    # Vertical dashed lines at boundaries
    ggplot2::geom_vline(data = boundary_lines,
                        ggplot2::aes(xintercept = xintercept),
                        linetype = "dashed", color = "gray70", linewidth = 0.5) +

    # Break symbols
    ggplot2::annotate("text", x = break_x, y = y_min,
                      label = "//", size = 8, fontface = "bold",vjust = 1) +

    # Zero line
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +

    # Exon boxes
    ggplot2::geom_rect(data = exon_regions,
                       ggplot2::aes(xmin = xmin, xmax = xmax,
                                    ymin = ymin, ymax = ymax, fill = fill),
                       color = "black", linewidth = 0.5, inherit.aes = FALSE) +
    ggplot2::scale_fill_identity() +

    # Intron lines
    ggplot2::geom_segment(data = intron_segments,
                          ggplot2::aes(x = x, xend = xend,
                                       y = y, yend = yend,
                                       linetype = linetype),
                          color = "black", linewidth = 1.5, inherit.aes = FALSE) +
    ggplot2::scale_linetype_identity() +

    ggplot2::scale_y_continuous(limits = c(y_min - exon_height * 1.5,
                                           y_max * 1.05 + y_range * 0.15)) +

    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 11, color = "black"),
      axis.ticks.y = ggplot2::element_line(color = "black"),
      ggplot2::labs(x = NULL, y = "Frequency", title = title),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 20,
                                          color = "black", face = "bold.italic"),
      legend.position = "bottom"
    )

  # Add significance bars if provided
  if (!is.null(sig_regions) && nrow(sig_regions) > 0) {
    # Helper to get bin number from global position
    get_bin <- function(pos) {
      dplyr::case_when(
        pos <= bin_width ~ 1L,
        pos <= 2 * bin_width ~ 2L,
        pos <= 3 * bin_width ~ 3L,
        TRUE ~ 4L
      )
    }

    # Helper function to convert global position to schematic position
    global_to_schematic <- function(pos) {
      bin <- get_bin(pos)
      position_in_bin <- pos - (bin - 1) * bin_width
      schematic_pos <- dplyr::case_when(
        bin == 1 ~ position_in_bin,
        bin == 2 ~ bin_width + gap + position_in_bin,
        bin == 3 ~ 2 * bin_width + gap + position_in_bin,
        bin == 4 ~ 3 * bin_width + 2 * gap + position_in_bin
      )
      return(schematic_pos)
    }

    # Calculate max y value per bin for positioning bars just above data
    bin_max_y <- freq_data %>%
      dplyr::group_by(bin) %>%
      dplyr::summarise(max_y = max(moving_avg, na.rm = TRUE), .groups = "drop")

    # Split significance regions at bin boundaries
    split_regions <- list()
    for (i in seq_len(nrow(sig_regions))) {
      reg <- sig_regions[i, ]
      start_bin <- get_bin(reg$start_pos)
      end_bin <- get_bin(reg$end_pos)

      if (start_bin == end_bin) {
        # Region is within a single bin
        split_regions[[length(split_regions) + 1]] <- data.frame(
          start_pos = reg$start_pos,
          end_pos = reg$end_pos,
          group = reg$group,
          bin = start_bin
        )
      } else {
        # Region spans multiple bins - split it
        for (b in start_bin:end_bin) {
          bin_start_global <- (b - 1) * bin_width + 1
          bin_end_global <- b * bin_width

          seg_start <- max(reg$start_pos, bin_start_global)
          seg_end <- min(reg$end_pos, bin_end_global)

          if (seg_start <= seg_end) {
            split_regions[[length(split_regions) + 1]] <- data.frame(
              start_pos = seg_start,
              end_pos = seg_end,
              group = reg$group,
              bin = b
            )
          }
        }
      }
    }

    if (length(split_regions) > 0) {
      sig_regions_split <- dplyr::bind_rows(split_regions)

      # Convert to schematic coordinates
      sig_regions_split <- sig_regions_split %>%
        dplyr::mutate(
          schematic_start = global_to_schematic(start_pos),
          schematic_end = global_to_schematic(end_pos)
        )

      # Join with bin max y values and set bar position just above data
      sig_regions_split <- sig_regions_split %>%
        dplyr::left_join(bin_max_y, by = "bin") %>%
        dplyr::mutate(
          bar_y = max_y + y_range * dplyr::case_when(
            group == "Retained" ~ 0.06,
            group == "Excluded" ~ 0.03,
            TRUE ~ 0.045
          )
        )

      # Add significance bars as segments (thinner lines)
      plot <- plot +
        ggplot2::geom_segment(
          data = sig_regions_split,
          ggplot2::aes(x = schematic_start, xend = schematic_end,
                       y = bar_y, yend = bar_y,
                       color = group),
          linewidth = 1.2,
          lineend = "butt",
          inherit.aes = FALSE
        )
    }
  }

  return(plot)
}
