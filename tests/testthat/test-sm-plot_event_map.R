# Tests for the splicing / sequence map plotter in R/sm-plot_event_map.R
#
# Following the house style for plots: assert the data and geometry the helpers
# compute (positions, y-geometry, significance-bar coordinates, legend labels)
# rather than rendered pixels. Layer-adding helpers are checked for the layer
# they emit (or NULL), and plot_event_map is checked to assemble a ggplot.
#
# Layout constants use width_exon = 10, width_intron = 20 with the RI schema
# (2 regions): region_width = 31, region_starts = c(0, 110), offsets = c(10, 20).

opts_fixture <- function() splicing_options(width_exon = 10, width_intron = 20)

# --- .add_schematic_position ----------------------------------------------

test_that(".add_schematic_position places x and builds the region:group key", {
  layout <- event_schema_ri$plot_layout(10, 20)
  data <- data.frame(
    region_idx         = c(1L, 1L, 2L),
    position_in_region = c(1L, 2L, 1L),
    group              = c("Positive", "Positive", "Control"),
    stringsAsFactors   = FALSE
  )
  out <- .add_schematic_position(data, event_schema_ri, layout, opts_fixture())
  expect_equal(out$schematic_position, c(1, 2, 111))     # 0+1, 0+2, 110+1
  expect_equal(out$plot_group, c("1:Positive", "1:Positive", "2:Control"))
  # group becomes a factor in canonical order, keeping only present levels
  expect_equal(levels(out$group), c("Positive", "Control"))
})

# --- .y_geometry ----------------------------------------------------------

test_that(".y_geometry derives the range and schematic band from the curve", {
  y <- .y_geometry(data.frame(moving_avg = c(0, 0.5, 1, NA)))
  expect_equal(y$y_max, 1)
  expect_equal(y$y_range, 1)
  expect_equal(y$exon_height, 0.08)             # 8% of the range
  expect_equal(y$y_min, -0.05)                  # min(0, data_min) - 5% range
})

test_that(".y_geometry handles a flat curve and an all-NA curve", {
  flat <- .y_geometry(data.frame(moving_avg = c(0.5, 0.5)))
  expect_equal(flat$y_range, 0.05)              # zero range -> 10% of y_max
  empty <- .y_geometry(data.frame(moving_avg = NA_real_))
  expect_equal(empty$y_max, 0.01)               # fallback when nothing is valid
})

# --- .ribbon_layer --------------------------------------------------------

test_that(".ribbon_layer is drawn only when a group carries SD", {
  base <- data.frame(
    schematic_position = 1:2, moving_avg = c(0.5, 0.5),
    plot_group = "1:Control", group = "Control"
  )
  expect_null(.ribbon_layer(cbind(base, moving_avg_sd = c(0, 0)), splicing_style()))
  lyr <- .ribbon_layer(cbind(base, moving_avg_sd = c(0.1, 0.1)), splicing_style())
  expect_s3_class(lyr, "ggproto")
})

# --- .boundary_lines_layer ------------------------------------------------

test_that(".boundary_lines_layer marks each junction at region_start + offset", {
  layout <- event_schema_ri$plot_layout(10, 20)
  lyr <- .boundary_lines_layer(layout, splicing_style())
  expect_s3_class(lyr, "ggproto")
  expect_equal(lyr$data$xintercept, c(10, 130))   # c(0,110) + c(10,20)
})

# --- .significance_bars ---------------------------------------------------

test_that(".significance_bars collapses contiguous significant positions per run", {
  # Positive significant at positions 1,2,3 (contiguous, region 1), 5 (region 1
  # but broken by a gap), and 40 (region 2). Each run becomes one segment.
  sig <- data.frame(
    position    = c(1, 2, 3, 5, 40),
    significant = TRUE,
    group       = "Positive",
    stringsAsFactors = FALSE
  )
  layout <- event_schema_ri$plot_layout(10, 20)
  y   <- list(y_max = 1, y_range = 1)
  lyr <- .significance_bars(sig, event_schema_ri, layout, opts_fixture(),
                            y, splicing_style())
  expect_s3_class(lyr, "ggproto")
  # region_starts[region] + position_in_region: run1 x1..3, run2 x5, run3 x119
  expect_equal(lyr$data$x,    c(1, 5, 119))
  expect_equal(lyr$data$xend, c(3, 5, 119))
  expect_equal(unique(lyr$data$y), 1.05)          # y_max + 5% of range
})

test_that(".significance_bars returns NULL when nothing is significant", {
  expect_null(.significance_bars(NULL, event_schema_ri,
                                 event_schema_ri$plot_layout(10, 20),
                                 opts_fixture(), list(y_max = 1, y_range = 1),
                                 splicing_style()))
  none <- data.frame(position = 1:3, significant = FALSE, group = "Positive")
  expect_null(.significance_bars(none, event_schema_ri,
                                 event_schema_ri$plot_layout(10, 20),
                                 opts_fixture(), list(y_max = 1, y_range = 1),
                                 splicing_style()))
})

# --- .group_color_scale ---------------------------------------------------

test_that(".group_color_scale labels each group with its cutoff and event count", {
  data <- data.frame(
    group    = factor(c("Negative", "Positive", "Control"),
                      levels = c("Negative", "Positive", "Control")),
    n_events = c(5L, 4L, 3L)
  )
  sc <- .group_color_scale(data, splicing_style(), splicing_options())
  expect_s3_class(sc, "ScaleDiscrete")
  expect_equal(sc$name, "Event group")
  expect_equal(sc$labels, c(
    Negative = sprintf("ΔΨ < %g [n = %s]", -0.1, "5"),
    Positive = sprintf("ΔΨ > %g [n = %s]",  0.1, "4"),
    Control  = "Control [n = 3]"
  ))
})

# --- .plot_theme ----------------------------------------------------------

test_that(".plot_theme returns a ggplot theme", {
  expect_s3_class(.plot_theme(splicing_style()), "theme")
})

# --- plot_event_map -------------------------------------------------------

test_that("plot_event_map assembles a ggplot from a frequency frame", {
  # A minimal but complete frequency frame: two groups over region 1, carrying
  # the columns event_map_pipeline would hand the plotter.
  data <- data.frame(
    region_idx         = 1L,
    position_in_region = rep(1:3, times = 2),
    group              = rep(c("Positive", "Control"), each = 3),
    moving_avg         = c(0.2, 0.4, 0.2, 0.1, 0.1, 0.1),
    moving_avg_sd      = c(0, 0, 0, 0.05, 0.05, 0.05),
    n_events           = rep(c(4L, 3L), each = 3),
    stringsAsFactors   = FALSE
  )
  p <- plot_event_map(data, event_schema_ri, splicing_style(), opts_fixture(),
                      significance = NULL, title = "RI")
  expect_s3_class(p, "ggplot")
})
