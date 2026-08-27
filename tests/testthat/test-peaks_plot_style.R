# Tests for peaks_plot_style() in R/params-peaks_plot_style.R
#
# ~60 visual knobs, each delegating to a check_* validator (the validators
# themselves are exhaustively tested in test-utils-conditions.R). Here we pin:
#   (1) the list contract,
#   (2) the special / cross-argument logic unique to this constructor,
#   (3) that EVERY argument's validation is actually wired up, via a table of
#       one known-bad value per argument.

# One rejecting value per argument. The choice reflects that argument's
# documented constraint: colors -> a non-color, non-negative numbers -> -1,
# unit intervals -> 1.5, integer floors -> below the floor, choices -> a
# non-member, flags -> NA, ranges -> descending.
bad_values <- list(
  exon_color                = "notacolor",
  exon_height               = -1,
  utr_color                 = "notacolor",
  utr_height                = -1,
  gene_offset               = -1,
  gene_lane_height          = -1,
  intron_color              = "notacolor",
  intron_linewidth          = -1,
  intron_arrow_len_in       = -1,
  total_arrows              = -1,
  max_per_intron            = -1,
  min_intron_frac           = 1.2,
  min_intron_bp             = -1,
  shaft_frac                = 1.2,
  peak_color                = "notacolor",
  peak_height               = -1,
  peak_alpha                = 1.5,
  peak_border_color         = "notacolor",
  peak_border_linewidth     = -1,
  band_even_fill            = "notacolor",
  band_odd_fill             = "notacolor",
  band_sep_color            = "notacolor",
  band_sep_linewidth        = -1,
  protein_label_size        = -1,
  protein_label_color       = "notacolor",
  protein_label_x_offset_bp = -1,
  strand_label_size         = -1,
  strand_label_color        = "notacolor",
  strand_label_x_offset     = -1,
  strand_label_y_offset     = "x",       # free numeric -> reject non-numeric
  gene_label_size           = -1,
  gene_label_color          = "notacolor",
  gene_label_x_offset       = -1,
  title_size                = -1,
  title_color               = "notacolor",
  subtitle_size             = -1,
  subtitle_color            = "notacolor",
  subtitle_sep              = 1,         # must be a string
  axis_title_size           = -1,
  axis_text_size            = -1,
  axis_pad_bp               = -1,
  axis_breaks_n             = 1,         # min = 2
  five_to_three             = NA,
  plot_top_margin           = -1,
  plot_right_margin         = -1,
  plot_bottom_margin        = -1,
  plot_left_margin          = -1,        # NULL allowed, but -1 rejected
  bam_fill_color            = "notacolor",
  bam_fill_alpha            = 1.5,
  bam_label_size            = -1,
  bam_axis_text_size        = -1,
  bam_ylim                  = c(100, 0), # descending range
  bam_track_height          = -1,
  bam_gap                   = -1,
  highlight                 = c(200, 100),
  highlight_color           = "notacolor",
  highlight_opacity         = 1.5,
  show_junctions            = NA,
  junction_color            = "notacolor",
  junction_linetype         = "squiggle",
  junction_linewidth        = -1,
  junction_alpha            = 1.5
)

test_that("defaults produce a valid list with documented values", {
  st <- peaks_plot_style()
  expect_type(st, "list")
  expect_equal(st$exon_color, "navy")
  expect_false(st$five_to_three)
  expect_null(st$plot_left_margin)
})

test_that("the returned list contains exactly the documented parameters", {
  expect_setequal(names(peaks_plot_style()), names(formals(peaks_plot_style)))
})

test_that("the bad-value table covers every argument of peaks_plot_style", {
  expect_setequal(names(bad_values), names(formals(peaks_plot_style)))
})

test_that("every argument rejects a known-bad value with an invalid_arg error", {
  for (arg in names(bad_values)) {
    expect_error(
      do.call(peaks_plot_style, stats::setNames(list(bad_values[[arg]]), arg)),
      class = "rnapeaks_error_invalid_arg",
      info  = arg
    )
  }
})

# --- special: gene_lane_height clamp --------------------------------------

test_that("gene_lane_height below exon_height warns and clamps up to exon_height", {
  expect_warning(
    st <- peaks_plot_style(gene_lane_height = 0.1, exon_height = 0.5),
    regexp = "clamping"
  )
  expect_equal(st$gene_lane_height, 0.5)
})

test_that("gene_lane_height at or above exon_height is left untouched", {
  expect_silent(st <- peaks_plot_style(gene_lane_height = 0.9, exon_height = 0.5))
  expect_equal(st$gene_lane_height, 0.9)
})

# --- special: highlight via normalize_coord + range check -----------------

test_that("highlight accepts NULL, numeric, or comma-grouped numeric strings", {
  expect_null(peaks_plot_style(highlight = NULL)$highlight)
  expect_equal(peaks_plot_style(highlight = c(100, 200))$highlight, c(100, 200))
  expect_equal(peaks_plot_style(highlight = c("1,000", "2,000"))$highlight, c(1000, 2000))
})

# --- special: NULL-allowed numeric args -----------------------------------

test_that("plot_left_margin and bam_ylim accept NULL (auto) or valid values", {
  expect_null(peaks_plot_style(plot_left_margin = NULL)$plot_left_margin)
  expect_equal(peaks_plot_style(plot_left_margin = 40)$plot_left_margin, 40)
  expect_null(peaks_plot_style(bam_ylim = NULL)$bam_ylim)
  expect_equal(peaks_plot_style(bam_ylim = c(0, 100))$bam_ylim, c(0, 100))
})

# --- special: NA-allowed color (peak border) ------------------------------
# Unlike splicing_style / utr_style, this constructor keeps allow_na = TRUE on
# its colors, so NA is a meaningful value meaning "no border".

test_that("peak_border_color may be NA (no border)", {
  expect_true(is.na(peaks_plot_style(peak_border_color = NA)$peak_border_color))
})

# --- choices --------------------------------------------------------------

test_that("junction_linetype accepts a valid ggplot linetype", {
  expect_equal(peaks_plot_style(junction_linetype = "dotted")$junction_linetype, "dotted")
})
