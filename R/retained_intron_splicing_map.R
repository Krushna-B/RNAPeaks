#' Retained Intron Splicing Map
#'
#' Computes per-position peak binding frequency around the two splice
#' boundaries of retained-intron events and renders the two-region splicing
#' map.
#'
#' @param events Data frame of rMATS RI.MATS events.
#' @param bed_file BED file path or BED data frame of peaks (chr, start, end,
#'   name, score, strand). Validated via `check_bed()`.
#' @param opts Result of [splicing_options()].
#' @param style Result of [splicing_style()].
#' @param title Plot title.
#'
#' @return A list with `plot` (ggplot) and `data` (list of `frequency` and
#'   `significance` tables).
#'
#' @family splicing_maps
#' @export
retained_intron_splicing_map <- function(events, bed_file,
                                          opts  = splicing_options(),
                                          style = splicing_style(),
                                          title = "") {
  #Wraps Error's thrown
  wrap_sm_errors("retained intron splicing map", {
    # Required args have no default; report a forgotten one as a clear
    # validation error instead of a base-R "missing argument" error.
    if (missing(events)) {
      abort_invalid_arg(c(
        "{.arg events} is required.",
        "i" = "Supply a data frame of rMATS RI.MATS events."
      ))
    }
    if (missing(bed_file)) {
      abort_invalid_arg(c(
        "{.arg bed_file} is required.",
        "i" = "Supply a BED file path or a BED data frame of peaks."
      ))
    }
    #Validate Input Params
    validate_sm_inputs(events, opts, style, title, bed_file = bed_file)

    #Validate bed file and turn into GRanges
    bed_gr <- .peaks_to_granges(bed_file)
    #Choose Scorer function (peaks_scorer for splicing)
    scorer <- function(regions_gr, n_regions, region_width, group_name) {
      peaks_scorer(regions_gr, bed_gr, n_regions, region_width)
    }
    #Run pipeline
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
