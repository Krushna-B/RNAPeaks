# Shared internal pipeline for plot_gene() and plot_region().
#
# Both public functions resolve a GTF and select transcript(s), then call the pipeline
#   1. resolves the filter / clip window (caller-supplied for region; derived
#      from transcripts for gene)
#   2. validates and prepares the BED into per-track peaks
#   3. builds the gene or region structure (region structure is clipped to
#      the window so the plot matches what the user requested)
#   4. renders (draw_plot validates the region/peaks contracts itself)
#
#' @noRd
#' @family plot
plot_peaks_pipeline <- function(transcripts,
                                bed,
                                is_region,
                                window     = NULL,
                                bam_files  = NULL,
                                peaks_opts = peaks_options(),
                                style      = peaks_plot_style()) {

  # 1. Resolve the active window. Gene mode window is found through transcript
  # Region mode window is passed from input params
  if (is.null(window)) {
    window <- list(start = min(transcripts$start),
                   end   = max(transcripts$end))
  }

  filter <- list(
    chr      = as.character(transcripts$seqnames[1]),
    start    = window$start,
    end      = window$end,
    strand   = as.character(transcripts$strand[1]),
    collapse = peaks_opts$collapse
  )

  # 2. Validate BED (restricting to `include` tracks), then prepare into peaks
  cli::cli_progress_step("Validating BED")
  bed      <- check_bed(bed, split_col = peaks_opts$split_by,
                        include = peaks_opts$include)
  cli::cli_progress_step("Filtering peaks ({nrow(bed)} rows)")
  peaks_df <- prepare_bed(
    bed,
    filter       = filter,
    order        = list(by           = peaks_opts$order_by,
                        in_          = peaks_opts$order_in,
                        max_proteins = peaks_opts$max_proteins),
    track_height = style$peak_height
  )

  # No peaks in window
  if (is.null(peaks_df)) {
    cli::cli_progress_done()
    return(invisible(NULL))
  }

  # 3. Place gene baseline above the peak stack
  center <- max(peaks_df$y_end) + style$gene_offset + style$exon_height / 2
  layout <- list(
    center           = center,
    exon_height      = style$exon_height,
    utr_height       = style$utr_height,
    gene_lane_height = style$gene_lane_height
  )

  # 4. Build the gene or region structure
  region <- if (is_region) {
    build_region_structure(transcripts, layout, window = window)
  } else {
    build_gene_structure(transcripts, layout)
  }

  # 5. BAM coverage tracks stack above the gene structure
  if (length(bam_files) > 0L) {
    cli::cli_progress_step(
      "Loading BAM coverage ({length(bam_files)} track{?s})"
    )
  }
  bam_tracks <- prepare_bam_tracks(
    bam_files = bam_files,
    chr       = filter$chr,
    start     = window$start,
    end       = window$end,
    base_y    = max(region$y_end) + style$bam_gap,
    style     = style
  )

  # 6. Render
  cli::cli_progress_step("Rendering plot")
  out <- draw_plot(region = region, peaks = peaks_df, bam_tracks = bam_tracks,
                    is_region = is_region, style = style)
  cli::cli_progress_done()
  out
}
