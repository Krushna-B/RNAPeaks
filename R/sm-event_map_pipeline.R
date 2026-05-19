#' Internal: splicing / sequence map pipeline
#'
#' Stages:
#'   1. filter events into groups.
#'   2. score each group via `scorer` -> 0/1 matrix of
#'      `[n_events x (n_regions * region_width)]`.
#'   3. bootstrap Control mean + SD on the raw matrix.
#'   4. per-position significance (Negative / Positive vs. Control).
#'   5. assemble frequency frame.
#'   6. apply moving-average smoothing.
#'   7. render via `plot_fn`.
#'
#' @param events Data frame with `schema$required_cols`.
#' @param schema Event-type schema list.
#' @param scorer `function(regions_gr, n_events, n_regions, region_width)`
#'   returning an `[n_events x (n_regions * region_width)]` 0/1 matrix.
#' @param opts Result of [splicing_options()].
#' @param style Result of [splicing_style()].
#' @param plot_fn `function(data, schema, style, opts, significance, title)`
#'   returning a ggplot.
#' @param title Plot title.
#'
#' @return `list(plot, data)` with `data = list(frequency, significance)`.
#'   `frequency` carries raw `frequency` / `frequency_sd` plus smoothed
#'   `moving_avg` / `moving_avg_sd` columns.
#'
#' @keywords internal
#' @noRd
event_map_pipeline <- function(events, schema, scorer, opts, style,
                                plot_fn, title = "") {
  #Validate Params
  .validate_pipeline_inputs(schema, scorer, opts, style, plot_fn)

  #Pull options from style and opts constructors
  region_width    <- schema$region_width(opts$width_exon, opts$width_intron)
  n_regions       <- schema$n_regions
  total_positions <- as.integer(n_regions * region_width)

  #Filter groups
  groups         <- filter_events(events, schema, opts)
  .report_group_sizes(groups, opts) #Print for verbose

  #Score all events in one pass, then slice rows back by group.
  non_empty <- groups[vapply(groups, nrow, 0L) > 0L]
  if (length(non_empty) == 0L) {
    score_matrices <- lapply(groups, function(.) matrix(0L, 0L, total_positions))
  } else {
    all_events <- do.call(rbind, non_empty)
    regions_gr <- schema$build_regions(all_events, opts$width_exon, opts$width_intron)
    M_all      <- scorer(
      regions_gr,
      n_events     = nrow(all_events),
      n_regions    = n_regions,
      region_width = region_width
    )
    .check_score_matrix_shape(M_all, nrow(all_events), total_positions) #Validate shape

    score_matrices <- lapply(names(groups), function(g) {
      if (g %in% names(non_empty)) M_all[all_events$group == g, , drop = FALSE]
      else                         matrix(0L, 0L, total_positions)
    })
    names(score_matrices) <- names(groups)
  }


  #Bootstrap Control mean + SD (skip if Control wasn't requested)
  control_stats <- if ("Control" %in% names(score_matrices)) {
    bootstrap_control(
      control_score_matrix = score_matrices$Control,
      n_pos                = NROW(groups$Positive),
      n_neg                = NROW(groups$Negative),
      opts                 = opts
    )
  } else NULL

  #Per-position significance (Negative / Positive vs. Control)
  sig_df <- .significance_table(score_matrices, opts, style)

  #Build Frequency Frame
  freq_df <- .assemble_frequency_frame(
    score_matrices, control_stats, n_regions, region_width
  )

  #Moving average for display only — does not feed stats
  freq_df$moving_avg    <- .smooth_by_group(
    freq_df$frequency,    freq_df$group, opts$moving_average,
    n_regions, region_width
  )
  freq_df$moving_avg_sd <- .smooth_by_group(
    freq_df$frequency_sd, freq_df$group, opts$moving_average,
    n_regions, region_width
  )

  #Create Plot
  plot <- plot_fn(
    data         = freq_df,
    schema       = schema,
    style        = style,
    opts         = opts,
    significance = sig_df,
    title        = title
  )

  #Return plot, and data
  list(plot = plot, data = list(frequency = freq_df, significance = sig_df))
}


#Helpers

#Input Validation
.validate_pipeline_inputs <- function(schema, scorer, opts, style, plot_fn) {
  if (!is.list(schema)) {
    abort_invalid_arg("{.arg schema} must be a list.")
  }
  required_fields <- c("required_cols", "n_regions", "region_width",
                       "build_regions", "group_set")
  missing_fields  <- setdiff(required_fields, names(schema))
  if (length(missing_fields) > 0L) {
    abort_invalid_arg(c(
      "{.arg schema} is missing required field{?s}.",
      "x" = "Missing: {.val {missing_fields}}."
    ))
  }
  if (!is.function(schema$region_width) || !is.function(schema$build_regions)) {
    abort_invalid_arg(
      "{.arg schema$region_width} and {.arg schema$build_regions} must be functions."
    )
  }
  if (!is.character(schema$required_cols) || length(schema$required_cols) == 0L) {
    abort_invalid_arg("{.arg schema$required_cols} must be a non-empty character vector.")
  }
  if (!is.numeric(schema$n_regions) || length(schema$n_regions) != 1L ||
      schema$n_regions < 1L) {
    abort_invalid_arg("{.arg schema$n_regions} must be a positive integer.")
  }
  if (!is.function(scorer))  abort_invalid_arg("{.arg scorer} must be a function.")
  if (!is.function(plot_fn)) abort_invalid_arg("{.arg plot_fn} must be a function.")
  if (!is.list(opts) || is.null(opts$width_exon) || is.null(opts$moving_average)) {
    abort_invalid_arg("{.arg opts} must be a {.fn splicing_options} result.")
  }
  if (!is.list(style) || is.null(style$show_significance)) {
    abort_invalid_arg("{.arg style} must be a {.fn splicing_style} result.")
  }
}

#CLI information for verbose
.report_group_sizes <- function(groups, opts) {
  if (!isTRUE(opts$verbose)) return(invisible())
  sizes <- vapply(groups, nrow, 0L)
  cli::cli_inform(
    "Group sizes: {paste(names(sizes), sizes, sep = '=', collapse = ', ')}"
  )
}

#Matrix shape validation
.check_score_matrix_shape <- function(M, expected_rows, expected_cols) {
  if (!is.matrix(M)) {
    abort_invalid_arg(c(
      "{.arg scorer} must return a matrix.",
      "x" = "Got {.cls {class(M)[1]}}."
    ))
  }
  if (nrow(M) != expected_rows || ncol(M) != expected_cols) {
    abort_invalid_arg(c(
      "{.arg scorer} returned the wrong shape.",
      "i" = "Expected {.val {expected_rows}}x{.val {expected_cols}}.",
      "x" = "Got {nrow(M)}x{ncol(M)}."
    ))
  }
}


#Per-position significance vs. Control
.significance_table <- function(score_matrices, opts, style) {
  if (!isTRUE(style$show_significance))       return(NULL)
  if (!"Control" %in% names(score_matrices))  return(NULL)
  if (nrow(score_matrices$Control) == 0L) {
    cli::cli_warn("Skipping significance: Control group has no events.")
    return(NULL)
  }

  rows <- list()
  for (g in intersect(c("Negative", "Positive"), names(score_matrices))) {
    M <- score_matrices[[g]]
    if (nrow(M) == 0L) next
    out       <- test_per_position(M, score_matrices$Control, opts)
    out$group <- g
    rows[[g]] <- out
  }
  if (length(rows) == 0L) return(NULL)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#frequency frame, one row per group x position
.assemble_frequency_frame <- function(score_matrices, control_stats,
                                       n_regions, region_width) {
  total_positions <- as.integer(n_regions * region_width)
  group_names     <- names(score_matrices)
  n_groups        <- length(group_names)

  freqs <- lapply(score_matrices, function(M) {
    if (nrow(M) > 0L) colMeans(M) else rep(0, total_positions)
  })
  sds <- rep(list(rep(0, total_positions)), n_groups)
  names(sds) <- group_names
  if (!is.null(control_stats)) {
    freqs$Control <- control_stats$mean_per_position
    sds$Control   <- control_stats$sd_per_position
  }

  data.frame(
    global_position    = rep(seq_len(total_positions), times = n_groups),
    region_idx         = rep(rep(seq_len(n_regions),   each = region_width), times = n_groups),
    position_in_region = rep(rep(seq_len(region_width), times = n_regions),  times = n_groups),
    frequency          = unlist(freqs, use.names = FALSE),
    frequency_sd       = unlist(sds,   use.names = FALSE),
    group              = rep(group_names, each = total_positions),
    n_events           = rep(vapply(score_matrices, nrow, 0L), each = total_positions),
    stringsAsFactors   = FALSE
  )
}

#Per-group moving average
.smooth_by_group <- function(values, groups, window, n_regions, region_width) {
  out <- values
  for (g in unique(groups)) {
    idx      <- groups == g
    out[idx] <- apply_moving_average(values[idx], window, n_regions, region_width)
  }
  out
}
