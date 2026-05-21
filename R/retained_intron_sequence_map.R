#' Retained Intron Sequence Map
#'
#' Computes per-position motif frequency around the two splice boundaries of
#' retained-intron events and renders the two-region sequence map.
#'
#' @param events Data frame of rMATS RI.MATS events.
#' @param sequence Character vector of motifs (IUPAC ambiguity codes
#'   supported).
#' @param genome `NULL` (hg38 default), one of `"hg38"` / `"mm10"` / `"mm39"`,
#'   or a `BSgenome` instance.
#' @param opts Result of [splicing_options()].
#' @param style Result of [splicing_style()].
#' @param title Plot title.
#' @param motif_mode One of `"combined"` (default — all motifs pooled into
#'   a single map) or `"individual"` (one map per motif).
#'
#' @return For `motif_mode = "combined"`: a list with `plot` (ggplot) and
#'   `data` (list of `frequency` and `significance` tables). For
#'   `motif_mode = "individual"`: a named list of those results, one entry
#'   per motif.
#'
#' @family sequence_maps
#' @export
retained_intron_sequence_map <- function(events, sequence,
                                          genome     = NULL,
                                          opts       = splicing_options(),
                                          style      = splicing_style(),
                                          title      = "",
                                          motif_mode = "combined") {
  #Wraps Error's thrown
  wrap_sm_errors("retained intron sequence map", {
    #Validate Input Params
    validate_sm_inputs(events, opts, style, title,
                       sequence = sequence, genome = genome,
                       motif_mode = motif_mode)
    #Normalize all motifs and switch U to T
    motifs <- .normalize_motifs(sequence)

    #Select Genome based on user genome input (Default is "hg38" if user doesn't provide)
    bsg    <- .resolve_genome(genome)

    #Build Preparation (Filter out events, build regions into GRanges, Extract Sequences for those Regions)
    prep   <- .prepare_sequence_map_prep(events, event_schema_ri, opts, bsg, motifs)

    #Finish rest of pipeline (Finding hits, calculating significance, and plotting)
    .run_sequence_map(motifs, motif_mode, prep,
                      event_schema_ri, opts, style, title)
  })
}
