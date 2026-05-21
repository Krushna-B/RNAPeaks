#' Internal: splicing / sequence map pipeline
#'
#' Stages:
#'   1. partition events into group index sets.
#'   2. score each non-empty group via `scorer` -> hit list.
#'   3. tabulate per-position counts per group.
#'   4. bootstrap Control mean + SD from a sparse hit matrix.
#'   5. per-position significance (Negative / Positive vs. Control).
#'   6. assemble frequency frame.
#'   7. moving-average smoothing.
#'   8. render via `plot_fn`.
#'
#' @param events Data frame with `schema$required_cols`. Ignored if `prep`
#'   is supplied.
#' @param schema Event-type schema list.
#' @param scorer `function(regions_gr, n_regions, region_width, group_name)`
#'   returning `list(event_id, col_idx)` for hits in that group.
#' @param opts Result of [splicing_options()].
#' @param style Result of [splicing_style()].
#' @param plot_fn `function(data, schema, style, opts, significance, title)`
#'   returning a ggplot.
#' @param title Plot title.
#' @param prep Optional precomputed result of `prepare_event_map()`. When
#'   supplied, `events` is not used.
#'
#' @return `list(plot, data)` with `data = list(frequency, significance)`.
#'
#' @keywords internal
#' @noRd
event_map_pipeline <- function(events = NULL, schema, scorer, opts, style,
                                plot_fn, title = "", prep = NULL) {
  #Validate input params
  .validate_pipeline_inputs(schema, scorer, opts, style, plot_fn)

  #Build out preparation if prep is not provided (Splicing Map, no prep, Sequence Map prep
  # due to possibility of individual or combined mode)
  if (is.null(prep)) {
    if (is.null(events)) {
      abort_invalid_arg("Either {.arg events} or {.arg prep} must be supplied.")
    }
    prep <- prepare_event_map(events, schema, opts)
  }

  groups_idx       <- prep$groups_idx
  regions_by_group <- prep$regions_by_group
  n_regions        <- prep$n_regions
  region_width     <- prep$region_width
  n_positions      <- prep$total_positions

  #Verbose mode print size to cli
  .report_group_sizes(groups_idx, opts)

  #Per each group find hits
  per_group <- vector("list", length(groups_idx))
  names(per_group) <- names(groups_idx)
  for (g in names(groups_idx)) {
    n_g <- length(groups_idx[[g]])
    if (n_g == 0L) {
      per_group[[g]] <- list(counts = integer(n_positions), n = 0L, hits = NULL)
      next
    }
    cli::cli_progress_step("Scoring {.val {g}} ({n_g} events)")
    #Pass inputs to specific scorer for hits
    regions_g <- regions_by_group[[g]]
    hits_g    <- scorer(regions_g, n_regions, region_width, group_name = g)
    #Check shape
    .check_hits_shape(hits_g, n_positions)

    #Tabulate hits per position
    counts_g  <- tabulate(hits_g$col_idx, nbins = n_positions)
    per_group[[g]] <- list(
      counts = counts_g,
      n      = n_g,
      hits   = if (identical(g, "Control")) hits_g else NULL
    )
  }
  cli::cli_progress_done()

  #Bootstrap Controls
  control_stats <- if ("Control" %in% names(per_group)) {
    n_pos <- if (!is.null(per_group$Positive)) per_group$Positive$n else 0L
    n_neg <- if (!is.null(per_group$Negative)) per_group$Negative$n else 0L
    bootstrap_control(
      ctrl_hits   = per_group$Control$hits,
      n_control   = per_group$Control$n,
      n_positions = n_positions,
      n_pos       = n_pos,
      n_neg       = n_neg,
      opts        = opts
    )
  } else NULL

  #Build significance df
  cli::cli_progress_step("Computing significance")
  sig_df <- .significance_table(per_group, opts, style)

  #Build frequency df + smoothing
  cli::cli_progress_step("Assembling frequency table")
  freq_df <- .assemble_frequency_frame(per_group, control_stats,
                                        n_regions, region_width)

  #Moving averages on mean and std
  freq_df$moving_avg    <- .smooth_by_group(
    freq_df$frequency,    freq_df$group, opts$moving_average,
    n_regions, region_width
  )
  freq_df$moving_avg_sd <- .smooth_by_group(
    freq_df$frequency_sd, freq_df$group, opts$moving_average,
    n_regions, region_width
  )

  #Return plot + data
  cli::cli_progress_step("Rendering plot")
  plot <- plot_fn(
    data         = freq_df,
    schema       = schema,
    style        = style,
    opts         = opts,
    significance = sig_df,
    title        = title
  )
  cli::cli_progress_done()

  list(plot = plot, data = list(frequency = freq_df, significance = sig_df))
}


#' Internal: prepare shared event-map inputs
#'
#' Normalizes chromosomes, applies the read-count filter, partitions events
#' into group index sets, and pre-builds per-group region `GRanges`.
#' @noRd
prepare_event_map <- function(events, schema, opts) {
  cli::cli_progress_step("Preparing events ({nrow(events)} rows)")
  #Validate events frame based on scheme type
  check_data_frame(events, "events", required = schema$required_cols)
  #Validate columns types
  .check_event_column_types(events, schema)
  #Normalize chr column in events
  events$chr <- normalize_chr(events$chr)

  #Precompute per-position axis geometry
  region_width    <- schema$region_width(opts$width_exon, opts$width_intron)
  n_regions       <- schema$n_regions
  total_positions <- as.integer(n_regions * region_width)

  #Filter low-coverage events and partition into Negative / Positive / Control
  filtered   <- filter_events(events, schema, opts)
  events     <- filtered$events
  groups_idx <- filtered$groups_idx

  #Build per-group region GRanges (mcols carry event_id + region_idx)
  regions_by_group <- lapply(groups_idx, function(idx) {
    if (length(idx) == 0L) return(GenomicRanges::GRanges())
    schema$build_regions(events[idx, , drop = FALSE],
                         opts$width_exon, opts$width_intron)
  })

  cli::cli_progress_done()
  list(
    groups_idx       = groups_idx,
    events           = events,
    regions_by_group = regions_by_group,
    n_regions        = n_regions,
    region_width     = region_width,
    total_positions  = total_positions
  )
}


#Every required event column except the text ones (chr, strand, GeneID)
#must be numeric
.check_event_column_types <- function(events, schema) {
  text_cols    <- c("chr", "strand", "GeneID")
  numeric_cols <- setdiff(schema$required_cols, text_cols)
  is_num       <- vapply(numeric_cols,
                          function(col) is.numeric(events[[col]]),
                          logical(1L))
  not_numeric  <- numeric_cols[!is_num]
  if (length(not_numeric) > 0L) {
    abort_invalid_arg(c(
      "{.arg events} has non-numeric column{?s}: {.field {not_numeric}}.",
      "i" = "Coordinates and statistics must be numeric. Did your reader (e.g. {.fn read.csv} with default {.code stringsAsFactors}) keep them as character?"
    ))
  }
}

#Validate param inputs schema params, and ensure scorer,and plots are functions
#(Dev validation for new schemas) not public exposed API
.validate_pipeline_inputs <- function(schema, scorer, opts, style, plot_fn) {
  if (!is.list(schema)) {
    abort_invalid_arg("{.arg schema} must be a list.")
  }
  required_fields <- c("required_cols", "n_regions", "region_width",
                       "build_regions", "group_set")
  missing_fields  <- setdiff(required_fields, names(schema))
  if (length(missing_fields) > 0L) {
    abort_invalid_arg(
      "{.arg schema} is missing required field{?s}: {.val {missing_fields}}."
    )
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
}

#Print out sizes of Groups to cli for verbose mode
.report_group_sizes <- function(groups_idx, opts) {
  if (!isTRUE(opts$verbose)) return(invisible())
  sizes <- vapply(groups_idx, length, 0L)
  cli::cli_inform(
    "Group sizes: {paste(names(sizes), sizes, sep = '=', collapse = ', ')}"
  )
}

#Validate hits shape
#(Dev validation for new scorers) not public exposed API
.check_hits_shape <- function(h, n_positions) {
  if (!is.list(h) || is.null(h$event_id) || is.null(h$col_idx)) {
    abort_invalid_arg(c(
      "{.arg scorer} must return a list with {.field event_id} and {.field col_idx}.",
      "x" = "Got {.cls {class(h)[1]}}."
    ))
  }
  if (length(h$event_id) != length(h$col_idx)) {
    abort_invalid_arg(c(
      "{.arg scorer} returned mismatched {.field event_id} / {.field col_idx} lengths.",
      "x" = "{length(h$event_id)} vs {length(h$col_idx)}."
    ))
  }
  if (length(h$col_idx) > 0L &&
      (min(h$col_idx) < 1L || max(h$col_idx) > n_positions)) {
    abort_invalid_arg(c(
      "{.arg scorer} returned {.field col_idx} outside {.val 1..{n_positions}}.",
      "x" = "Range: {.val {range(h$col_idx)}}."
    ))
  }
}

#Inputs: Named list for groups of counts for each position
#And input params for options and style
#Return dataframe with signficance per position
.significance_table <- function(per_group, opts, style) {
  #Skips if signifance shoudl be skipped or empty
  if (!isTRUE(style$show_significance))      return(NULL)
  if (!"Control" %in% names(per_group))      return(NULL)
  if (per_group$Control$n == 0L) {
    cli::cli_warn("Skipping significance: Control group has no events.")
    return(NULL)
  }


  rows <- list()
  #Run signifance test per position on each group
  for (g in intersect(c("Negative", "Positive"), names(per_group))) {
    if (per_group[[g]]$n == 0L) next
    out <- test_per_position(
      group_counts   = per_group[[g]]$counts,
      n_group        = per_group[[g]]$n,
      control_counts = per_group$Control$counts,
      n_control      = per_group$Control$n,
      opts           = opts
    )
    out$group <- g
    rows[[g]] <- out
  }
  if (length(rows) == 0L) return(NULL)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#Build's data frame for plotting
.assemble_frequency_frame <- function(per_group, control_stats,
                                       n_regions, region_width) {
  total_positions <- as.integer(n_regions * region_width)
  group_names     <- names(per_group)
  n_groups        <- length(group_names)

  freqs <- vector("list", n_groups); names(freqs) <- group_names
  sds   <- vector("list", n_groups); names(sds)   <- group_names
  for (g in group_names) {
    n_g <- per_group[[g]]$n
    freqs[[g]] <- if (n_g > 0L) per_group[[g]]$counts / n_g
                  else          rep(0, total_positions)
    sds[[g]]   <- rep(0, total_positions)
  }
  if (!is.null(control_stats)) {
    freqs$Control <- control_stats$mean_per_position
    sds$Control   <- control_stats$sd_per_position
  }

  n_by_group <- vapply(per_group, `[[`, integer(1L), "n")

  data.frame(
    global_position    = rep(seq_len(total_positions), times = n_groups),
    region_idx         = rep(rep(seq_len(n_regions),   each = region_width), times = n_groups),
    position_in_region = rep(rep(seq_len(region_width), times = n_regions),  times = n_groups),
    frequency          = unlist(freqs, use.names = FALSE),
    frequency_sd       = unlist(sds,   use.names = FALSE),
    group              = rep(group_names, each = total_positions),
    n_events           = rep(n_by_group, each = total_positions),
    stringsAsFactors   = FALSE
  )
}

#Moving average on each group
.smooth_by_group <- function(values, groups, window, n_regions, region_width) {
  out <- values
  for (g in unique(groups)) {
    idx      <- groups == g
    out[idx] <- apply_moving_average(values[idx], window, n_regions, region_width)
  }
  out
}
