#' Skipped Exon Splicing Map
#'
#' Computes per-position peak binding frequency around the four splice
#' boundaries of skipped-exon events and renders the standard four-region
#' splicing map.
#'
#' @param events Data frame of rMATS SE.MATS events.
#' @param bed_file BED file path or BED data frame of peaks.
#' @param opts Result of [splicing_options()].
#' @param style Result of [splicing_style()].
#' @param title Plot title.
#'
#' @return A list with `plot` (ggplot) and `data` (list of `frequency` and
#'   `significance` tables).
#'
#' @family splicing_maps
#' @export
skipped_exon_splicing_map <- function(events, bed_file,
                                       opts  = splicing_options(),
                                       style = splicing_style(),
                                       title = "") {
  .run_sm_entry("skipped exon splicing map", {
    .assert_sm_entry_args(title, bed_file = bed_file, has_bed = TRUE)

    bed_gr <- .peaks_to_granges(bed_file)
    scorer <- function(regions_gr, n_events, n_regions, region_width) {
      peaks_scorer(regions_gr, bed_gr, n_events, n_regions, region_width)
    }
    event_map_pipeline(
      events  = events,
      schema  = event_schema_se,
      scorer  = scorer,
      opts    = opts,
      style   = style,
      plot_fn = plot_event_map,
      title   = title
    )
  })
}
