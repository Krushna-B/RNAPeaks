# Shared internal pipeline for plot_gene() and plot_region().
#
# Both public functions resolve a GTF and select transcript(s), then call the pipeline
#   1. derives the filter window from `transcripts`
#   2. validates and prepares the BED into per-track peaks
#   3. builds the gene or region structure
#   4. renders (draw_plot validates the region/peaks contracts itself)
#
#' @noRd
plot_peaks_pipeline <- function(transcripts,
                                bed,
                                is_region,
                                peaks_opts = peaks_options(),
                                style      = peaks_plot_style()) {

  # 1. Derive filter window from transcripts
  filter <- list(
    chr      = as.character(transcripts$seqnames[1]),
    start    = min(transcripts$start),
    end      = max(transcripts$end),
    strand   = as.character(transcripts$strand[1]),
    omit     = peaks_opts$omit,
    collapse = peaks_opts$collapse
  )

  # 2. Validate BED, then prepare into peaks
  bed      <- check_bed(bed, split_col = peaks_opts$split_by)
  peaks_df <- prepare_bed(
    bed,
    filter       = filter,
    order        = list(by           = peaks_opts$order_by,
                        in_          = peaks_opts$order_in,
                        max_proteins = peaks_opts$max_proteins),
    track_height = style$peak_height
  )

  # No peaks in window
  if (is.null(peaks_df)) return(invisible(NULL))

  # 3. Place gene baseline above the peak stack
  center <- max(peaks_df$y_end) + style$gene_offset + style$exon_height / 2
  layout <- list(
    center      = center,
    exon_height = style$exon_height,
    utr_height  = style$utr_height
  )

  # 4. Build the gene or region structure
  region <- if (is_region) {
    build_region_structure(transcripts, layout)
  } else {
    build_gene_structure(transcripts, layout)
  }

  # 5. Render (draw_plot validates the region + peaks contracts)
  draw_plot(region = region, peaks = peaks_df,
            is_region = is_region, style = style)
}
