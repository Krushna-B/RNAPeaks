#' Skipped Exon Splicing Map
#'
#' Computes per-position peak binding frequency around the four splice
#' boundaries of skipped-exon events and renders the standard four-region
#' splicing map. Events are filtered into Negative / Positive / Control groups
#' by `opts$psi_cutoff`, scored against the supplied peaks, and the Control
#' distribution is bootstrapped for the SD ribbon. Negative and Positive
#' groups are tested per-position against Control.
#'
#' @param events Data frame of rMATS SE.MATS events. Required columns are
#'   listed by [event_schema_se]; `IncLevelDifference`, `FDR`, and `PValue`
#'   are used for group assignment.
#' @param bed_file BED file path or BED data frame of peaks (chr, start, end,
#'   name, score, strand). Validated via [check_bed()].
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
}
