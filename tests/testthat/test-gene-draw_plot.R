# Tests for R/gene-draw_plot.R
#
# The big draw_plot() assembler is exercised structurally (returns a ggplot,
# validates inputs). Its pure helpers are tested for exact behaviour, and the
# layer builders are checked by inspecting the ggproto layer / theme objects.

# --- gene_title -----------------------------------------------------------

test_that("gene_title joins gene names and appends the first transcript id", {
  region <- data.frame(gene_name = c("SRSF1", "SRSF1"), gene_id = "G1",
                       transcript_id = "ENST1")
  expect_equal(gene_title(region), "SRSF1 (ENST1)")
})

test_that("gene_title falls back to gene_id and omits an absent transcript", {
  r1 <- data.frame(gene_name = NA_character_, gene_id = "G1", transcript_id = "ENST1")
  expect_equal(gene_title(r1), "G1 (ENST1)")
  r2 <- data.frame(gene_name = "A", gene_id = "G1", transcript_id = NA_character_)
  expect_equal(gene_title(r2), "A")
})

# --- gene_label_text ------------------------------------------------------

test_that("gene_label_text uses gene_name unless it is NA/empty, then gene_id", {
  region <- data.frame(gene_name = c("A", NA, ""), gene_id = c("G1", "G2", "G3"))
  expect_equal(gene_label_text(region), c("A", "G2", "G3"))
})

# --- background_table -----------------------------------------------------

test_that("background_table returns one band per unique y-range spanning xlim", {
  peaks <- data.frame(y_start = c(0, 0.3, 0), y_end = c(0.3, 0.6, 0.3))
  bg <- background_table(peaks, c(100, 1000))
  expect_equal(nrow(bg), 2)               # two distinct (y_start, y_end) bands
  expect_equal(bg$y_start, c(0, 0.3))     # sorted ascending
  expect_true(all(bg$x_start == 100))
  expect_true(all(bg$x_end == 1000))
})

# --- find_left_margin -----------------------------------------------------

test_that("find_left_margin returns the pad for no labels and grows with label length", {
  expect_equal(find_left_margin(character(0), 4), 8)  # pad_pt default
  wide   <- find_left_margin("WWWWWWWWWWWW", 4)
  narrow <- find_left_margin("i", 4)
  expect_gt(wide, narrow)
  expect_gt(narrow, 8)
})

# --- add_highlight_band ---------------------------------------------------

test_that("add_highlight_band adds a layer only when a visible highlight is set", {
  g0 <- ggplot2::ggplot()
  n0 <- length(g0$layers)

  none <- add_highlight_band(g0, c(0, 1000), peaks_plot_style(highlight = NULL))
  expect_equal(length(none$layers), n0)

  inside <- add_highlight_band(g0, c(0, 1000), peaks_plot_style(highlight = c(100, 200)))
  expect_equal(length(inside$layers), n0 + 1)

  # highlight entirely outside the padded range -> no layer
  outside <- add_highlight_band(g0, c(0, 1000), peaks_plot_style(highlight = c(2000, 3000)))
  expect_equal(length(outside$layers), n0)
})

# --- region_gene_labels ---------------------------------------------------

test_that("region_gene_labels builds one label per gene at the visual-left edge", {
  region <- data.frame(gene_id = c("G1", "G1", "G2"),
                       gene_name = c("A", "A", "B"),
                       y_start = c(0, 0.5, 1),
                       y_end   = c(0.3, 0.8, 1.3))
  layer <- region_gene_labels(region, xlim = c(100, 1000), style = peaks_plot_style())
  expect_s3_class(layer, "LayerInstance")
  d <- layer$data
  expect_setequal(d$label, c("A", "B"))
  # center = (min y_start(0) + max y_end(0.8)) / 2 for gene G1
  expect_equal(d$y[d$label == "A"], 0.4)
  # x = xlim[1] - axis_pad_bp(500) * gene_label_x_offset(0.25) = 100 - 125
  expect_true(all(d$x == -25))
})

# --- plot_theme -----------------------------------------------------------

test_that("plot_theme returns a ggplot theme object", {
  th <- plot_theme(peaks_plot_style(), left_margin = 40)
  expect_s3_class(th, "theme")
})

# --- draw_plot (structural) -----------------------------------------------

make_region <- function() {
  build_gene_structure(
    make_gtf()[make_gtf()$transcript_id == "ENST00000001", ],
    list(center = 0, exon_height = 0.5, utr_height = 0.3)
  )
}
make_peaks <- function() {
  data.frame(start = 200, end = 250, y_start = -0.3, y_end = 0,
             group_name = "SRSF1", stringsAsFactors = FALSE)
}

test_that("draw_plot returns a ggplot for gene and region modes", {
  expect_s3_class(draw_plot(make_region(), make_peaks()), "ggplot")
  expect_s3_class(draw_plot(make_region(), make_peaks(), is_region = TRUE), "ggplot")
})

test_that("draw_plot validates its region and peaks frames", {
  region <- make_region(); peaks <- make_peaks()
  expect_error(draw_plot(region[, setdiff(names(region), "y_start")], peaks),
               class = "rnapeaks_error_invalid_arg")
  expect_error(draw_plot(region, peaks[, setdiff(names(peaks), "group_name")]),
               class = "rnapeaks_error_invalid_arg")
})

test_that("draw_plot rejects an invalid strand", {
  region <- make_region(); region$strand <- "."
  expect_error(draw_plot(region, make_peaks()), class = "rnapeaks_error_invalid_arg")
})

test_that("draw_plot in region mode requires a gene_id column", {
  region <- make_region(); region$gene_id <- NULL
  expect_error(draw_plot(region, make_peaks(), is_region = TRUE),
               class = "rnapeaks_error_invalid_arg")
})
