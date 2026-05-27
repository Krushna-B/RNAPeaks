#' Plot RBP peak binding density over 5' and 3' UTRs
#'
#' For every protein-coding gene in `gtf` (or the subset named by
#' `transcripts`), picks the transcript with the longest exonic 5' UTR and
#' the one with the longest exonic 3' UTR. Peaks from `bed` that fall in
#' the UTR are binarized per-bp on the spliced transcript, resampled to
#' `n_bins` per side, and averaged across genes. UTRs shorter than
#' `n_bins` keep each bp at its true position and pad the rest with 0.
#' The result is two density curves (one per BED track) shown as
#' `[5' UTR] // [3' UTR]`.
#'
#' @param bed Validated BED data frame, BED file path, or named list of
#'   either; one density curve is drawn per element.
#' @param gtf Optional GTF file path. If `NULL`, the bundled
#'   annotation for `species` is used.
#' @param transcripts Optional character vector of transcript ids to
#'   restrict to. Default uses every protein-coding transcript.
#' @param species One of `"hg38"`, `"mm10"`, `"mm39"`. Ignored when
#'   `gtf` is supplied.
#' @param moving_average Window size for moving-average smoothing of the
#'   per-bin curve. `0` or `NULL` disables smoothing.
#' @param style Output of [utr_style()].
#' @param track_colors Optional named character vector of colors keyed
#'   by BED track name. When `NULL`, a single track is drawn black and
#'   multiple tracks use a `Set1` palette.
#' @param title Plot title.
#'
#' @return A list with `plot` (ggplot) and `data` (per-bin density frame).
#' @export
#' @family plot
plot_utr_binding <- function(bed,
                             gtf            = NULL,
                             transcripts    = NULL,
                             species        = "hg38",
                             moving_average = 5L,
                             style          = utr_style(),
                             track_colors   = NULL,
                             title          = "") {
  tryCatch(
    {
      check_string(title, "title")
      species <- normalize_str(species)
      check_string(species, "species", choices = c("hg38", "mm10", "mm39"))
      if (!is.null(moving_average)) {
        check_scalar_int(moving_average, "moving_average", min = 0)
      }

      cli::cli_progress_step("Loading GTF")
      gtf <- get_GTF(species = species, file = gtf)

      cli::cli_progress_step("Building UTR event table")
      built <- build_utr_events(gtf, transcripts = transcripts)
      events <- built$events
      if (nrow(events) == 0L) {
        abort_not_found("No genes with usable UTRs were found.")
      }
      n5 <- sum(events$utr5_len > 0)
      n3 <- sum(events$utr3_len > 0)
      cli::cli_alert_info(
        "Kept {nrow(events)} gene{?s}: {n5} with 5' UTR, {n3} with 3' UTR."
      )

      cli::cli_progress_step("Preparing BED tracks")
      tracks <- .prep_utr_bed_tracks(bed)

      cli::cli_progress_step("Scoring utr region")
      n_bins <- event_schema_utr$n_bins
      freq_df <- .score_all_tracks(tracks, built, n_bins,
                                   n5 = n5, n3 = n3)
      cli::cli_progress_done()

      # Per-track smoothing across the two concatenated panels.
      freq_df$moving_avg <- .smooth_utr(freq_df, moving_average, n_bins)

      cli::cli_progress_step("Rendering plot")
      p <- plot_utr_map(
        data         = freq_df,
        schema       = event_schema_utr,
        style        = style,
        track_colors = track_colors,
        title        = title
      )
      cli::cli_progress_done()

      list(plot = p, data = freq_df)
    },
    error = function(cnd) {
      if (inherits(cnd, "rnapeaks_error")) {
        cli::cli_abort("Failed to generate UTR binding plot.", parent = cnd)
      } else {
        cli::cli_abort(
          c("Failed to generate UTR binding plot.",
            "x" = "An unexpected error occurred."),
          parent = cnd
        )
      }
    }
  )
}

# Normalize the `bed` argument to a named list of reduced GRanges.
.prep_utr_bed_tracks <- function(bed) {
  items <- if (is.data.frame(bed)) {
    list(bed)
  } else if (is.character(bed)) {
    as.list(bed)
  } else if (is.list(bed)) {
    bed
  } else {
    abort_invalid_arg(c(
      "{.arg bed} must be a BED data frame, file path(s), or a named list of either.",
      "x" = "Got {.cls {class(bed)[1]}}."
    ))
  }
  if (length(items) == 0L) {
    abort_invalid_arg("{.arg bed} must contain at least one track.")
  }
  ok <- vapply(items, function(b) {
    is.data.frame(b) || (is.character(b) && length(b) == 1L && !is.na(b))
  }, logical(1L))
  if (!all(ok)) {
    bad <- which(!ok)[1L]
    abort_invalid_arg(c(
      "Each element of {.arg bed} must be a BED data frame or a single file path.",
      "x" = "Element {bad} is {.cls {class(items[[bad]])[1]}}.",
      "i" = "If you have multiple BEDs, use {.code list(bed1, bed2)}."
    ))
  }
  nm <- names(items)
  if (is.null(nm)) nm <- rep("", length(items))
  nm[nm == ""] <- paste0("bed", which(nm == ""))
  if (anyDuplicated(nm)) {
    abort_invalid_arg(c(
      "{.arg bed} has duplicate track names.",
      "x" = "Duplicates: {.val {unique(nm[duplicated(nm)])}}."
    ))
  }

  lapply(stats::setNames(items, nm), function(b) {
    if (is.character(b)) {
      if (!file.exists(b)) {
        abort_not_found(c("BED file does not exist.",
                          "x" = "Path: {.path {b}}."))
      }
      b <- utils::read.table(b)
    }
    df <- check_bed(b)
    gr <- GenomicRanges::makeGRangesFromDataFrame(
      df,
      seqnames.field     = "chr",
      start.field        = "start",
      end.field          = "end",
      strand.field       = "strand",
      keep.extra.columns = FALSE
    )
    GenomicRanges::reduce(gr)
  })
}

# Build the per-bin frequency frame across all tracks and both sides.
.score_all_tracks <- function(tracks, built, n_bins, n5, n3) {
  rows <- list()
  for (g in names(tracks)) {
    s5 <- score_utr_metagene(built$utr5_pieces, tracks[[g]],
                             n_events = n5, n_bins = n_bins)
    s3 <- score_utr_metagene(built$utr3_pieces, tracks[[g]],
                             n_events = n3, n_bins = n_bins)
    rows[[g]] <- data.frame(
      global_position    = seq_len(2L * n_bins),
      region_idx         = rep(c(1L, 2L), each = n_bins),
      position_in_region = rep(seq_len(n_bins), times = 2L),
      frequency          = c(s5$density, s3$density),
      group              = g,
      n_events           = c(rep(s5$n, n_bins), rep(s3$n, n_bins)),
      stringsAsFactors   = FALSE
    )
  }
  do.call(rbind, rows)
}

# Per-track moving average. Smoothing stays within each utr pannel.
.smooth_utr <- function(df, window, n_bins) {
  if (is.null(window) || window <= 1L) return(df$frequency)
  out <- df$frequency
  for (g in unique(df$group)) {
    idx <- df$group == g
    out[idx] <- apply_moving_average(df$frequency[idx], window,
                                     n_regions = 2L, region_width = n_bins)
  }
  out
}


#' Internal: UTR plotter
#'
#' @keywords internal
#' @noRd
plot_utr_map <- function(data, schema, style, track_colors = NULL,
                          title = "") {
  layout <- schema$plot_layout(n_bins = schema$n_bins)

  data$schematic_position <-
    layout$region_starts[data$region_idx] + data$position_in_region
  data$group      <- factor(data$group, levels = unique(data$group))
  data$plot_group <- paste(data$region_idx, data$group, sep = ":")

  y <- .y_geometry(data)

  # Resolve per-track colors. Single track defaults to black
  track_levels <- levels(data$group)
  colors <- if (!is.null(track_colors)) {
    missing <- setdiff(track_levels, names(track_colors))
    if (length(missing) > 0L) {
      abort_invalid_arg(c(
        "{.arg track_colors} is missing entries for some tracks.",
        "x" = "Missing: {.val {missing}}."
      ))
    }
    track_colors[track_levels]
  } else if (length(track_levels) == 1L) {
    stats::setNames(style$single_track_color, track_levels)
  } else {
    stats::setNames(rep(style$palette, length.out = length(track_levels)),
                    track_levels)
  }

  legend_labels <- vapply(track_levels, function(g) {
    n5 <- data$n_events[data$group == g & data$region_idx == 1L][1L]
    n3 <- data$n_events[data$group == g & data$region_idx == 2L][1L]
    sprintf("%s [n = %s / %s]", g,
            format(n5, big.mark = ","),
            format(n3, big.mark = ","))
  }, character(1L))

  ggplot2::ggplot(
    data,
    ggplot2::aes(x = schematic_position, y = moving_avg, color = group,
                  group = plot_group)
  ) +
    ggplot2::geom_line(linewidth = style$line_width,
                       alpha     = style$line_alpha) +
    schema$build_schematic_layers(layout, style, y$y_min, y$exon_height) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_manual(values = colors, labels = legend_labels,
                                 name = "BED track") +
    ggplot2::scale_y_continuous(
      limits = c(y$y_min - y$exon_height * 2,
                 y$y_max * 1.05 + y$y_range * 0.15)
    ) +
    ggplot2::labs(x = NULL, y = style$ylab, title = title) +
    .plot_theme(style)
}
