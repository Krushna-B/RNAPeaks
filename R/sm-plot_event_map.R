#' Internal: unified splicing / sequence map plotter
#'
#
#' @param data Frequency frame from `event_map_pipeline()`.
#' @param schema Event-type schema list.
#' @param style Result of [splicing_style()].
#' @param opts Result of [splicing_options()] (needs `width_exon`,
#'   `width_intron`, `psi_cutoff` for legend labels).
#' @param significance Per-position p-value table or `NULL`.
#' @param title Plot title.
#'
#' @return A `ggplot` object.
#'
#' @keywords internal
#' @noRd
plot_event_map <- function(data, schema, style, opts,
                            significance = NULL, title = "") {
  layout <- schema$plot_layout(opts$width_exon, opts$width_intron)
  data   <- .add_schematic_position(data, schema, layout, opts)
  y      <- .y_geometry(data)

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x     = schematic_position,
      y     = moving_avg,
      color = group,
      group = plot_group
    )
  ) +
    .ribbon_layer(data, style) +
    ggplot2::geom_line(linewidth = style$line_width, alpha = style$line_alpha) +
    .boundary_lines_layer(layout, style) +
    schema$build_schematic_layers(layout, style, y$y_min, y$exon_height) +
    ggplot2::scale_fill_identity() +
    .significance_bars(significance, schema, layout, opts, y, style) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    .group_color_scale(data, style, opts) +
    ggplot2::scale_y_continuous(
      limits = c(y$y_min - y$exon_height * 1.5,
                 y$y_max * 1.05 + y$y_range * 0.15)
    ) +
    ggplot2::labs(x = NULL, y = style$ylab, title = title) +
    .plot_theme(style)
}


#Helpers

#Add per-row x position and the (region x group) plotting group key.
.add_schematic_position <- function(data, schema, layout, opts) {
  data$schematic_position <-
    layout$region_starts[data$region_idx] + data$position_in_region

  present <- intersect(c("Negative", "Positive", "Control"), unique(data$group))
  data$group      <- factor(data$group, levels = present)
  data$plot_group <- paste(data$region_idx, data$group, sep = ":")
  data
}

#Y-axis geometry: data range plus space below for the schematic.
.y_geometry <- function(data) {
  valid <- data$moving_avg[!is.na(data$moving_avg)]
  if (length(valid) == 0L) {
    y_max <- 0.01; y_min_data <- 0
  } else {
    y_max <- max(valid); y_min_data <- min(valid)
  }
  y_range     <- y_max - y_min_data
  if (y_range == 0) y_range <- max(y_max * 0.1, .Machine$double.eps)
  exon_height <- y_range * 0.08
  y_min       <- min(0, y_min_data) - y_range * 0.05
  list(y_min = y_min, y_max = y_max,
       y_range = y_range, exon_height = exon_height)
}

#Control SD ribbon.
.ribbon_layer <- function(data, style) {
  ribbon <- data[data$moving_avg_sd > 0, , drop = FALSE]
  if (nrow(ribbon) == 0L) return(NULL)
  ribbon$ymin        <- ribbon$moving_avg - ribbon$moving_avg_sd
  ribbon$ymax        <- ribbon$moving_avg + ribbon$moving_avg_sd
  ribbon$ribbon_fill <- unlist(style$group_colors[as.character(ribbon$group)],
                                use.names = FALSE)
  ggplot2::geom_ribbon(
    data    = ribbon,
    mapping = ggplot2::aes(x = schematic_position,
                            ymin = ymin, ymax = ymax,
                            group = plot_group,
                            fill  = ribbon_fill),
    alpha = style$ribbon_alpha, color = NA, inherit.aes = FALSE
  )
}

#Dashed vertical lines at each exon/intron junction inside each region.
.boundary_lines_layer <- function(layout, style) {
  ggplot2::geom_vline(
    xintercept = layout$region_starts + layout$boundary_offsets,
    linetype   = "dashed",
    color      = style$boundary_col,
    linewidth  = 0.5
  )
}

#Significance bars above the data.
.significance_bars <- function(significance, schema, layout, opts, y, style) {
  if (is.null(significance) || nrow(significance) == 0L) return(NULL)
  sig <- significance[which(significance$significant), , drop = FALSE]
  if (nrow(sig) == 0L) return(NULL)

  region_width <- schema$region_width(opts$width_exon, opts$width_intron)
  sig$region_idx         <- as.integer(ceiling(sig$position / region_width))
  sig$position_in_region <- sig$position - (sig$region_idx - 1L) * region_width
  sig$x <- layout$region_starts[sig$region_idx] + sig$position_in_region

  sig <- sig[order(sig$group, sig$region_idx, sig$position), , drop = FALSE]

  key      <- paste(sig$group, sig$region_idx, sep = "::")
  new_run  <- c(TRUE, key[-1L] != key[-nrow(sig)] |
                       diff(sig$position) != 1L)
  run_id   <- cumsum(new_run)

  splits <- split(sig, run_id)
  bars   <- do.call(rbind, lapply(splits, function(d) {
    data.frame(x = min(d$x), xend = max(d$x),
               group = d$group[1L], stringsAsFactors = FALSE)
  }))
  bars$y    <- y$y_max + y$y_range * 0.05
  bars$yend <- bars$y

  ggplot2::geom_segment(
    data    = bars,
    mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = group),
    linewidth = 2.5, inherit.aes = FALSE, show.legend = FALSE
  )
}

#Legend labels.
.group_color_scale <- function(data, style, opts) {
  present <- levels(data$group)
  n_by_group <- vapply(present, function(g) {
    rows <- data$n_events[data$group == g]
    if (length(rows) == 0L) NA_integer_ else as.integer(rows[[1L]])
  }, integer(1L))

  fmt_n <- function(n) {
    if (is.na(n)) "?" else format(n, big.mark = ",")
  }

  labels <- vapply(present, function(g) {
    switch(g,
      Negative = sprintf("\u0394\u03a8 < %g [n = %s]",
                         opts$psi_cutoff[1L], fmt_n(n_by_group[[g]])),
      Positive = sprintf("\u0394\u03a8 > %g [n = %s]",
                         opts$psi_cutoff[2L], fmt_n(n_by_group[[g]])),
      Control  = sprintf("Control [n = %s]", fmt_n(n_by_group[[g]])),
      g
    )
  }, character(1L))

  values <- unlist(style$group_colors[present])

  ggplot2::scale_color_manual(values = values, labels = labels,
                               name = "Event group")
}

#Theme.
.plot_theme <- function(style) {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_blank(),
      axis.ticks.x     = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_text(size  = style$axis_text_size,
                                                color = "black"),
      axis.ticks.y     = ggplot2::element_line(color = "black"),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.title       = ggplot2::element_text(hjust = 0.5,
                                                size  = style$title_size,
                                                color = style$title_color,
                                                face  = "bold.italic"),
      legend.position  = style$legend_position
    )
}
