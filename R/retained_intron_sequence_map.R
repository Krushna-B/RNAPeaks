#' Retained Intron Sequence Map
#'
#' Computes per-position motif frequency around the two splice boundaries of
#' retained-intron events and renders the two-region sequence map.
#'
#' @param events Data frame of rMATS RI.MATS events. Required columns are
#'   listed by [event_schema_ri].
#' @param sequence Character vector of motifs (IUPAC ambiguity codes
#'   supported).
#' @param genome `NULL` (hg38 default), one of `"hg38"` / `"mm10"` / `"mm39"`,
#'   or a `BSgenome` instance.
#' @param opts Result of [splicing_options()].
#' @param style Result of [splicing_style()].
#' @param title Plot title.
#'
#' @return A list with `plot` (ggplot) and `data` (list of `frequency` and
#'   `significance` tables).
#'
#' @family sequence_maps
#' @export
retained_intron_sequence_map <- function(events, sequence,
                                          genome = NULL,
                                          opts   = splicing_options(),
                                          style  = splicing_style(),
                                          title  = "") {
  .run_sm_entry("retained intron sequence map", {
    .assert_sm_entry_args(title)

    motifs <- .normalize_motifs(sequence)
    bsg    <- .resolve_genome(genome)
    scorer <- function(regions_gr, n_events, n_regions, region_width) {
      motif_scorer(regions_gr, bsg, motifs, n_events, n_regions, region_width)
    }
    event_map_pipeline(
      events  = events,
      schema  = event_schema_ri,
      scorer  = scorer,
      opts    = opts,
      style   = style,
      plot_fn = plot_event_map,
      title   = title
    )
  })
}
