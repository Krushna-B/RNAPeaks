# Tests for prepare_bed() in R/bed-prepare.R
#
# Contract:
#   Filter peaks, choose a display order, and assign each protein track a y band.
#   - order$by "Count" ranks proteins by peak count; "Alphabetical" by name.
#   - order$in_ overrides ordering with an explicit vector (warns on absent names).
#   - order$max_proteins caps how many proteins are shown (warns when it truncates).
#   - The FIRST protein in the resolved order is drawn at the TOP (highest y).
#   - Returns a frame carrying rank / y_start / y_end; NULL when nothing survives
#     filtering.
#
# Fixtures use targets AAA (1 peak) and ZZZ (2 peaks) so Count and Alphabetical
# orderings are distinguishable.

full_window <- list(chr = "1", start = 0, end = 1000,
                    strand = "+", collapse = 0)

bed_two_targets <- function() {
  rbind(
    make_checked_bed(start = c(100, 300), end = c(150, 350), target = "ZZZ"),
    make_checked_bed(start = 700, end = 750, target = "AAA")
  )
}

# --- basic layout ---------------------------------------------------------

test_that("returns a frame with rank / y_start / y_end", {
  out <- prepare_bed(
    bed_two_targets(), filter = full_window,
    order = list(by = "Count", in_ = NULL, max_proteins = 40),
    track_height = 0.3
  )
  expect_true(all(c("rank", "y_start", "y_end") %in% colnames(out)))
  expect_equal(out$y_end, out$y_start + 0.3)
})

test_that("Count order puts the most-bound protein on top (highest rank / y)", {
  out <- prepare_bed(
    bed_two_targets(), filter = full_window,
    order = list(by = "Count", in_ = NULL, max_proteins = 40),
    track_height = 0.3
  )
  # ZZZ has more peaks -> ranked first -> drawn on top (larger rank value).
  expect_true(all(out$rank[out$group_name == "ZZZ"] >
                    out$rank[out$group_name == "AAA"]))
})

test_that("Alphabetical order ranks by name, not by count", {
  out <- prepare_bed(
    bed_two_targets(), filter = full_window,
    order = list(by = "Alphabetical", in_ = NULL, max_proteins = 40),
    track_height = 0.3
  )
  # AAA sorts first -> on top, despite having fewer peaks.
  expect_true(all(out$rank[out$group_name == "AAA"] >
                    out$rank[out$group_name == "ZZZ"]))
})

test_that("an unknown order$by aborts", {
  expect_error(
    prepare_bed(bed_two_targets(), filter = full_window,
                order = list(by = "Nonsense", in_ = NULL, max_proteins = 40)),
    class = "rnapeaks_error_invalid_arg"
  )
})

# --- explicit ordering ----------------------------------------------------

test_that("order$in_ imposes an explicit order and puts its first entry on top", {
  out <- prepare_bed(
    bed_two_targets(), filter = full_window,
    order = list(by = "Count", in_ = c("AAA", "ZZZ"), max_proteins = 40),
    track_height = 0.3
  )
  expect_true(all(out$rank[out$group_name == "AAA"] >
                    out$rank[out$group_name == "ZZZ"]))
})

test_that("order$in_ naming an absent target warns but still succeeds", {
  expect_warning(
    prepare_bed(bed_two_targets(), filter = full_window,
                order = list(by = "Count", in_ = c("AAA", "ZZZ", "NOPE"),
                             max_proteins = 40)),
    regexp = "not present"
  )
})

# --- max_proteins ---------------------------------------------------------

test_that("max_proteins truncates to the top-ranked proteins with a warning", {
  expect_warning(
    out <- prepare_bed(bed_two_targets(), filter = full_window,
                       order = list(by = "Count", in_ = NULL, max_proteins = 1)),
    regexp = "top"
  )
  # Only ZZZ (highest count) survives.
  expect_equal(unique(out$group_name), "ZZZ")
})

# --- pass-through of empty filter -----------------------------------------

test_that("prepare_bed returns NULL when filtering removes every peak", {
  empty_window <- modifyList(full_window, list(start = 0, end = 10))
  out <- suppressMessages(
    prepare_bed(bed_two_targets(), filter = empty_window,
                order = list(by = "Count", in_ = NULL, max_proteins = 40))
  )
  expect_null(out)
})

# --- argument validation --------------------------------------------------

test_that("an empty or non-frame bed aborts", {
  expect_error(
    prepare_bed(bed_two_targets()[0, ], filter = full_window,
                order = list(by = "Count", in_ = NULL, max_proteins = 40)),
    class = "rnapeaks_error_invalid_bed"
  )
})

test_that("filter and order must be lists", {
  expect_error(
    prepare_bed(bed_two_targets(), filter = "x",
                order = list(by = "Count", in_ = NULL, max_proteins = 40)),
    class = "rnapeaks_error_invalid_arg"
  )
})

test_that("max_proteins < 1 and non-positive track_height abort", {
  expect_error(
    prepare_bed(bed_two_targets(), filter = full_window,
                order = list(by = "Count", in_ = NULL, max_proteins = 0)),
    class = "rnapeaks_error_invalid_arg"
  )
  expect_error(
    prepare_bed(bed_two_targets(), filter = full_window,
                order = list(by = "Count", in_ = NULL, max_proteins = 40),
                track_height = 0),
    class = "rnapeaks_error_invalid_arg"
  )
})
