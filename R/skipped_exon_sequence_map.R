#' Skipped Exon Sequence Map
#'
#' Computes per-position motif frequency around the four splice boundaries of
#' skipped-exon events and renders the standard four-region sequence map.
#'
#' @param events Data frame of rMATS SE.MATS events.
#' @param sequence Character vector of motifs (IUPAC ambiguity codes
#'   supported).
#' @param genome `NULL` (hg38 default), one of `"hg38"` / `"mm10"` / `"mm39"`,
#'   or a `BSgenome` instance.
#' @param opts Result of [splicing_options()].
#' @param style Result of [splicing_style()].
#' @param title Plot title.
#' @param motif_mode One of `"combined"` (default all motifs pooled into
#'   a single map) or `"individual"` (one map per motif).
#'
#' @return For `motif_mode = "combined"`: a list with `plot` (ggplot) and
#'   `data` (list of `frequency` and `significance` tables). For
#'   `motif_mode = "individual"`: a named list of those results, one entry
#'   per motif.
#'
#' @family sequence_maps
#' @export
skipped_exon_sequence_map <- function(events, sequence,
                                       genome     = NULL,
                                       opts       = splicing_options(),
                                       style      = splicing_style(),
                                       title      = "",
                                       motif_mode = "combined") {
  wrap_sm_errors("skipped exon sequence map", {
    validate_sm_inputs(events, opts, style, title,
                       sequence = sequence, genome = genome,
                       motif_mode = motif_mode)

    motifs <- .normalize_motifs(sequence)
    bsg    <- .resolve_genome(genome)
    prep   <- .prepare_sequence_map_prep(events, event_schema_se, opts, bsg, motifs)

    .run_sequence_map(motifs, motif_mode, prep,
                      event_schema_se, opts, style, title)
  })
}
