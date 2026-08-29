#' Plot RBP peak binding density over 5' and 3' UTRs
#'
#' For every protein-coding gene in `gtf` (or the subset named by
#' `transcripts`), picks the transcript with the longest exonic 5' UTR and
#' the one with the longest exonic 3' UTR. Peaks from `bed` that fall in
#' the UTR are binarized per-bp on the spliced transcript, resampled to
#' `n_bins` per side, and averaged across genes. UTRs shorter than
#' `n_bins` keep each bp at its true position and pad the rest with 0.
#' The 5' and 3' UTRs are returned as two independent plots, each with its
#' own schematic (a thin UTR box abutting the CDS) and percentage markers
#' running 5' -> 3' along the window.
#'
#' @param bed Validated BED data frame, BED file path, or named list of
#'   either; one density curve is drawn per element.
#' @param gtf Optional GTF file path. If `NULL`, the bundled
#'   annotation for `species` is used.
#' @param transcripts Optional character vector restricting which transcripts
#'   contribute. Accepts Ensembl transcript ids (`ENST…`), gene ids (`ENSG…`),
#'   or gene symbols (e.g. `"CXCR4"`); gene-level ids expand to all of that
#'   gene's transcripts. Default uses every protein-coding transcript. Draws a
#'   single pooled curve per BED track; use `gene_groups` for split curves.
#'   Mutually exclusive with `gene_groups`.
#' @param gene_groups Optional *named* list of id vectors (same id forms as
#'   `transcripts`), one curve per group drawn with its own line type and
#'   labelled in the legend. A gene may appear in more than one group. Ids
#'   that match nothing are skipped; a group with no usable UTRs is dropped
#'   with a warning. Mutually exclusive with `transcripts`.
#' @param species One of `"hg38"`, `"mm10"`, `"mm39"`. Ignored when
#'   `gtf` is supplied.
#' @param moving_average Window size for moving-average smoothing of the
#'   per-bin curve. `0` or `NULL` disables smoothing.
#' @param style Output of [utr_style()]. Track colors come from the style:
#'   `single_track_color` for one track, and `palette` for several. A
#'   *named* `palette` (keyed by BED track name) assigns colors per track
#'   regardless of order; an unnamed one is applied positionally.
#' @param title Optional plot title. When supplied it prefixes each
#'   panel's heading (e.g. `"<title> - 5' UTR"`); otherwise each plot is
#'   titled `"5' UTR"` / `"3' UTR"`.
#'
#' @return A list with two elements, `utr5` and `utr3`. Each is itself a
#'   list with `plot` (ggplot) and `data` (per-bin density frame for that
#'   side).
#' @export
#' @family plot
plot_utr_binding <- function(bed,
                             gtf            = NULL,
                             transcripts    = NULL,
                             gene_groups    = NULL,
                             species        = "hg38",
                             moving_average = 5L,
                             style          = utr_style(),
                             title          = "") {
  tryCatch(
    {
      # `bed` is required and has no default; catch it here so a forgotten
      # argument reports as a clear validation error instead of falling
      # through to the generic "unexpected error" branch below.
      if (missing(bed)) {
        abort_invalid_arg(c(
          "{.arg bed} is required.",
          "i" = "Supply a BED data frame, a file path, or a named list of either."
        ))
      }

      check_string(title, "title")
      species <- normalize_str(species)
      check_string(species, "species", choices = c("hg38", "mm10", "mm39"))
      if (!is.null(moving_average)) {
        check_scalar_int(moving_average, "moving_average", min = 0)
      }
      groups_spec <- .resolve_utr_group_spec(gene_groups, transcripts)

      cli::cli_progress_step("Loading GTF")
      gtf <- get_GTF(species = species, file = gtf)

      cli::cli_progress_step("Building UTR event table")
      build_ids <- if (is.null(groups_spec)) NULL
                   else unique(unlist(groups_spec, use.names = FALSE))
      built <- build_utr_events(gtf, transcripts = build_ids)
      events <- built$events
      if (nrow(events) == 0L) {
        abort_not_found("No genes with usable UTRs were found.")
      }

      group_events      <- .assign_group_events(groups_spec, events, gtf)
      show_group_legend <- length(group_events) > 1L

      n5 <- sum(events$utr5_len > 0)
      n3 <- sum(events$utr3_len > 0)
      cli::cli_alert_info(
        "Kept {nrow(events)} gene{?s}: {n5} with 5' UTR, {n3} with 3' UTR."
      )

      # Capture the symbol the caller passed (e.g. `K562_bed`) so a single
      # unnamed data frame is labelled with its variable name rather than a
      # generic "bed1".
      bed_sym  <- substitute(bed)
      bed_name <- if (is.symbol(bed_sym)) as.character(bed_sym) else NULL

      cli::cli_progress_step("Preparing BED tracks")
      tracks <- .prep_utr_bed_tracks(bed, default_name = bed_name)

      cli::cli_progress_step("Scoring utr region")
      n_bins <- event_schema_utr$n_bins
      side5 <- .score_one_side(tracks, built$utr5_pieces, events,
                               group_events, "utr5_len", n_bins)
      side3 <- .score_one_side(tracks, built$utr3_pieces, events,
                               group_events, "utr3_len", n_bins)
      cli::cli_progress_done()

      # Per-track moving average within each (single-region) panel.
      side5$moving_avg <- .smooth_side(side5, moving_average, n_bins)
      side3$moving_avg <- .smooth_side(side3, moving_average, n_bins)

      cli::cli_progress_step("Rendering plots")
      ttl5 <- if (nzchar(title)) paste0(title, " \u2014 5' UTR") else "5' UTR"
      ttl3 <- if (nzchar(title)) paste0(title, " \u2014 3' UTR") else "3' UTR"
      p5 <- plot_utr_side_map(side5, event_schema_utr, style, "utr5",
                              title = ttl5, show_group_legend = show_group_legend)
      p3 <- plot_utr_side_map(side3, event_schema_utr, style, "utr3",
                              title = ttl3, show_group_legend = show_group_legend)
      cli::cli_progress_done()

      list(
        utr5 = list(plot = p5, data = side5),
        utr3 = list(plot = p3, data = side3)
      )
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
# `default_name` is the caller's symbol for a single unnamed data frame.
.prep_utr_bed_tracks <- function(bed, default_name = NULL) {
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
  for (i in seq_along(items)) {
    if (nm[i] != "") next
    it <- items[[i]]
    nm[i] <- if (is.character(it)) {
      tools::file_path_sans_ext(basename(it))
    } else if (length(items) == 1L && !is.null(default_name)) {
      default_name
    } else {
      paste0("bed", i)
    }
  }
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

# Normalize the grouping args to a named list of id vectors, or NULL when no
# grouping is requested (all protein-coding genes as one pooled curve).
# transcripts is the single-group shorthand; the two args are exclusive.
.resolve_utr_group_spec <- function(gene_groups, transcripts) {
  if (!is.null(gene_groups) && !is.null(transcripts)) {
    abort_invalid_arg(c(
      "Pass either {.arg transcripts} or {.arg gene_groups}, not both.",
      "i" = "{.arg transcripts} is the single-group case of {.arg gene_groups}."
    ))
  }
  if (!is.null(gene_groups)) {
    nm <- names(gene_groups)
    if (!is.list(gene_groups) || length(gene_groups) == 0L ||
        is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
      abort_invalid_arg(
        "{.arg gene_groups} must be a non-empty, uniquely named list of id vectors."
      )
    }
    ok <- vapply(gene_groups, function(v) {
      is.character(v) && length(v) > 0L && !anyNA(v)
    }, logical(1L))
    if (!all(ok)) {
      abort_invalid_arg(c(
        "Each element of {.arg gene_groups} must be a non-empty character vector.",
        "x" = "Bad group{?s}: {.val {nm[!ok]}}."
      ))
    }
    return(gene_groups)
  }
  if (!is.null(transcripts)) return(list(`All genes` = transcripts))
  NULL
}

# Map each group to the event rows it covers. NULL spec -> one group over all
# events. Groups that resolve to no usable UTR are warned and dropped; an
# all-empty result aborts.
.assign_group_events <- function(groups_spec, events, gtf) {
  if (is.null(groups_spec)) {
    return(list(`All genes` = seq_len(nrow(events))))
  }
  out <- list()
  for (g in names(groups_spec)) {
    res <- resolve_gene_ids(gtf, groups_spec[[g]])
    if (length(res$unmatched)) {
      cli::cli_alert_info("Group {.val {g}}: no match for {.val {res$unmatched}}.")
    }
    idx <- which(events$gene_id %in% res$gene_ids)
    if (length(idx) == 0L) {
      cli::cli_warn("Dropping group {.val {g}}: no genes with usable UTRs.")
      next
    }
    out[[g]] <- idx
  }
  if (length(out) == 0L) {
    abort_not_found("No gene groups had genes with usable UTRs.")
  }
  out
}

# Build the per-bin frequency frame for one UTR side across every
# (BED track, gene group) pair. Each group is averaged over its own events.
.score_one_side <- function(tracks, pieces, events, group_events, len_col, n_bins) {
  rows <- list()
  for (g in names(group_events)) {
    evt_idx <- group_events[[g]]
    n_ev    <- sum(events[[len_col]][evt_idx] > 0)
    sub     <- pieces[pieces$event_idx %in% evt_idx, , drop = FALSE]
    for (t in names(tracks)) {
      s <- score_utr_side(sub, tracks[[t]], n_events = n_ev, n_bins = n_bins)
      rows[[paste(g, t, sep = "\r")]] <- data.frame(
        position_in_region = seq_len(n_bins),
        frequency          = s$density,
        track              = t,
        gene_group         = g,
        n_events           = s$n,
        stringsAsFactors   = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# Moving average within each (track, group) curve of a single-region panel.
.smooth_side <- function(df, window, n_bins) {
  if (is.null(window) || window <= 1L) return(df$frequency)
  out   <- df$frequency
  combo <- interaction(df$track, df$gene_group, drop = TRUE)
  for (g in levels(combo)) {
    idx <- combo == g
    out[idx] <- apply_moving_average(df$frequency[idx], window,
                                     n_regions = 1L, region_width = n_bins)
  }
  out
}


#' Internal: single-side UTR plotter
#'
#' Draws one UTR side (`"utr5"` or `"utr3"`) as a single-region density
#' plot with a gene-anatomy schematic (UTR box + adjacent CDS) and
#' percentage markers beneath it.
#'
#' @keywords internal
#' @noRd
plot_utr_side_map <- function(data, schema, style, side, title = "",
                              show_group_legend = FALSE) {
  layout <- schema$plot_layout(n_bins = schema$n_bins)

  data$schematic_position <- layout$region_start + data$position_in_region
  data$track      <- factor(data$track,      levels = unique(data$track))
  data$gene_group <- factor(data$gene_group, levels = unique(data$gene_group))

  # Size the schematic to the displayed height. UTR curves sit on an
  # elevated baseline, so scaling boxes to the data range (as the splicing
  # maps do) would collapse them to slivers.
  vals  <- data$moving_avg[!is.na(data$moving_avg)]
  y_max <- if (length(vals)) max(vals) else 0.01
  if (y_max <= 0) y_max <- 0.01
  exon_height <- y_max * (style$schematic_height %||% 0.06)
  y_top       <- 0

  track_levels <- levels(data$track)
  group_levels <- levels(data$gene_group)

  # Per-track colors: named palette maps by name, unnamed is positional, a
  # lone track falls back to single_track_color.
  pal <- style$palette
  colors <- if (!is.null(names(pal)) && any(nzchar(names(pal)))) {
    missing <- setdiff(track_levels, names(pal))
    if (length(missing) > 0L) {
      abort_invalid_arg(c(
        "{.arg palette} is missing entries for some tracks.",
        "x" = "Missing: {.val {missing}}."
      ))
    }
    pal[track_levels]
  } else if (length(track_levels) == 1L) {
    stats::setNames(style$single_track_color, track_levels)
  } else {
    stats::setNames(rep(pal, length.out = length(track_levels)), track_levels)
  }

  # Per-group line types, same named / positional rule as the palette.
  lty <- style$linetypes
  linetypes <- if (!is.null(names(lty)) && any(nzchar(names(lty)))) {
    missing <- setdiff(group_levels, names(lty))
    if (length(missing) > 0L) {
      abort_invalid_arg(c(
        "{.arg linetypes} is missing entries for some gene groups.",
        "x" = "Missing: {.val {missing}}."
      ))
    }
    lty[group_levels]
  } else {
    stats::setNames(rep(lty, length.out = length(group_levels)), group_levels)
  }

  # The event count varies by gene group but not by track, so it annotates
  # whichever legend distinguishes the groups. With one group it stays on
  # the track legend, matching the ungrouped plot.
  with_n <- function(x, col) {
    vapply(x, function(v) {
      n <- data$n_events[data[[col]] == v][1L]
      sprintf("%s [n = %s]", v, format(n, big.mark = ","))
    }, character(1L))
  }
  if (show_group_legend) {
    track_labels <- track_levels
    group_labels <- with_n(group_levels, "gene_group")
  } else {
    track_labels <- with_n(track_levels, "track")
    group_labels <- group_levels
  }

  # x range includes the CDS block on the side-appropriate end.
  if (identical(side, "utr5")) {
    x_lim <- c(layout$region_start, layout$bin_width + layout$cds_width)
  } else {
    x_lim <- c(-layout$cds_width, layout$bin_width)
  }

  ggplot2::ggplot(
    data,
    ggplot2::aes(x = schematic_position, y = moving_avg, color = track,
                  linetype = gene_group,
                  group = interaction(track, gene_group))
  ) +
    ggplot2::geom_line(linewidth = style$line_width,
                       alpha     = style$line_alpha) +
    schema$build_schematic_layers(layout, style, y_top, exon_height, side) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_manual(values = colors, labels = track_labels,
                                 name = "BED track") +
    ggplot2::scale_linetype_manual(
      values = linetypes, labels = group_labels,
      name   = style$linetype_legend_name %||% "Gene group",
      guide  = if (show_group_legend) "legend" else "none"
    ) +
    ggplot2::scale_x_continuous(limits = x_lim,
                                expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(
      limits = c(-exon_height * 3.2, y_max * 1.1)
    ) +
    ggplot2::labs(x = NULL, y = style$ylab, title = title) +
    .plot_theme(style)
}
