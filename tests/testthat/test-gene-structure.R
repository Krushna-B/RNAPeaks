# Tests for R/gene-structure.R
#
# These verify the actual coordinate math, not just that the functions run:
#   - feature y-bands = center +/- (height / 2)
#   - introns span the gap between consecutive features (prev_end+1 .. next_start-1)
#     with y = center +/- exon_height/8
#   - region lanes stack at center + (i-1)*gene_lane_height
#   - clip_to_window trims to [start,end] and recomputes width

# --- build_gene_structure: exon-only (non-coding) geometry ----------------

test_that("build_gene_structure places exon bands at center +/- exon_height/2 and fills introns", {
  enst1 <- make_gtf()[make_gtf()$transcript_id == "ENST00000001", ]  # exons 100-300, 900-1100
  layout <- list(center = 0, exon_height = 0.5, utr_height = 0.3)
  out <- build_gene_structure(enst1, layout)

  exons <- out[out$type == "exon", ]
  expect_equal(nrow(exons), 2)
  expect_equal(exons$y_start, c(-0.25, -0.25))
  expect_equal(exons$y_end,   c(0.25, 0.25))

  intron <- out[out$type == "intron", ]
  expect_equal(nrow(intron), 1)
  expect_equal(intron$start, 301)      # prev_end (300) + 1
  expect_equal(intron$end,   899)      # next_start (900) - 1
  expect_equal(intron$y_start, -0.0625) # center - exon_height/8
  expect_equal(intron$y_end,    0.0625)
  expect_equal(intron$width, 599)       # 899 - 301 + 1
})

# --- build_gene_structure: protein-coding uses CDS + UTR with own heights --

test_that("protein-coding transcripts use CDS/UTR features with per-type heights", {
  tx <- data.frame(
    seqnames = "1", start = c(100, 300, 900), end = c(150, 500, 1100),
    strand = "+", type = c("UTR", "CDS", "CDS"),
    gene_biotype = "protein_coding", width = c(51, 201, 201),
    stringsAsFactors = FALSE
  )
  out <- build_gene_structure(tx, list(center = 0, exon_height = 0.5, utr_height = 0.3))

  utr <- out[out$type == "UTR", ]
  cds <- out[out$type == "CDS", ]
  expect_equal(utr$y_start, -0.15)   # utr_height/2
  expect_equal(unique(cds$y_start), -0.25)  # exon_height/2

  introns <- out[out$type == "intron", ]
  expect_equal(introns$start, c(151, 501))
  expect_equal(introns$end,   c(299, 899))
})

test_that("build_gene_structure returns NULL when no plottable features exist", {
  tx <- data.frame(seqnames = "1", start = 100, end = 1100, strand = "+",
                   type = "transcript", stringsAsFactors = FALSE)
  expect_null(build_gene_structure(tx, list(center = 0, exon_height = 0.5, utr_height = 0.3)))
})

# --- intron_rows ----------------------------------------------------------

test_that("intron_rows returns nothing for <= 1 feature or book-ended features", {
  one <- data.frame(start = 100, end = 300, type = "exon")
  expect_equal(nrow(intron_rows(one, 0, 0.5)), 0)

  # gap of exactly 0 (next_start == prev_end + 1) is not an intron
  adjacent <- data.frame(start = c(100, 301), end = c(300, 500), type = "exon")
  expect_equal(nrow(intron_rows(adjacent, 0, 0.5)), 0)
})

test_that("intron_rows sorts features and spans each real gap", {
  feats <- data.frame(start = c(900, 100), end = c(1100, 300), type = "exon")  # unsorted
  out <- intron_rows(feats, center = 0, exon_height = 0.8)
  expect_equal(out$start, 301)
  expect_equal(out$end,   899)
  expect_equal(out$y_start, -0.1)  # center - 0.8/8
  expect_equal(out$y_end,    0.1)
})

# --- clip_to_window -------------------------------------------------------

test_that("clip_to_window trims overlapping features and recomputes width", {
  df <- data.frame(start = c(100, 900), end = c(300, 1100), width = c(201, 201))
  out <- clip_to_window(df, 200, 1000)
  expect_equal(out$start, c(200, 900))
  expect_equal(out$end,   c(300, 1000))
  expect_equal(out$width, c(101, 101))  # end - start + 1
})

test_that("clip_to_window drops features that fall entirely outside the window", {
  df <- data.frame(start = c(100, 1500), end = c(300, 1600))
  out <- clip_to_window(df, 0, 1000)
  expect_equal(nrow(out), 1)
  expect_equal(out$start, 100)
})

test_that("clip_to_window returns an empty frame unchanged", {
  df <- data.frame(start = numeric(0), end = numeric(0))
  expect_equal(nrow(clip_to_window(df, 0, 100)), 0)
})

# --- build_region_structure: lane stacking --------------------------------

test_that("build_region_structure stacks lanes at center + (i-1)*gene_lane_height", {
  g <- make_gtf()
  tx <- g[g$transcript_id %in% c("ENST00000001", "ENST00000003"), ]
  layout <- list(center = 0, exon_height = 0.5, utr_height = 0.3, gene_lane_height = 1)
  out <- build_region_structure(tx, layout, list(start = 0, end = 6000))

  # ENST1 starts left (100) -> lane 1 at center 0; ENST3 (5000) -> lane 2 at center 1
  e1 <- out[out$transcript_id == "ENST00000001" & out$type == "exon", ]
  e3 <- out[out$transcript_id == "ENST00000003" & out$type == "exon", ]
  expect_equal(unique(e1$y_start), -0.25)
  expect_equal(unique(e3$y_start), 0.75)   # 1 - 0.25
  expect_equal(unique(e3$y_end),   1.25)
})

test_that("build_region_structure clips to the window, dropping out-of-view transcripts", {
  g <- make_gtf()
  tx <- g[g$transcript_id %in% c("ENST00000001", "ENST00000003"), ]
  layout <- list(center = 0, exon_height = 0.5, utr_height = 0.3, gene_lane_height = 1)
  out <- build_region_structure(tx, layout, list(start = 0, end = 1000))
  expect_false("ENST00000003" %in% out$transcript_id)  # 5000-5600 is outside
  expect_equal(max(out$end), 1000)                      # ENST1 exon2 clipped to window
})

test_that("build_region_structure requires gene_lane_height in the layout", {
  tx <- make_gtf()[make_gtf()$transcript_id == "ENST00000001", ]
  expect_error(
    build_region_structure(tx, list(center = 0, exon_height = 0.5, utr_height = 0.3),
                           list(start = 0, end = 2000)),
    class = "rnapeaks_error_invalid_arg"
  )
})

# --- make_strand_labels ---------------------------------------------------

test_that("make_strand_labels places 5'/3' tags offset from the gene ends", {
  tx <- data.frame(start = c(100, 900), end = c(300, 1100), strand = "+",
                   y_start = c(-0.25, -0.25))
  labs <- make_strand_labels(tx, axis_pad_bp = 500, x_offset = 100)
  expect_equal(labs$left$Label, "5'")
  expect_equal(labs$left$X, 0)       # gene_min(100) - 100
  expect_equal(labs$right$Label, "3'")
  expect_equal(labs$right$X, 1200)   # gene_max(1100) + 100
  expect_equal(labs$left$Y, -0.25)
})

test_that("make_strand_labels caps x_offset at axis_pad_bp and flips labels on - strand", {
  tx <- data.frame(start = c(100, 900), end = c(300, 1100), strand = "-",
                   y_start = c(-0.25, -0.25))
  labs <- make_strand_labels(tx, axis_pad_bp = 500, x_offset = 1000, y_offset = 0.5)
  expect_equal(labs$left$Label, "3'")     # minus strand reverses tags
  expect_equal(labs$left$X, -400)         # gene_min(100) - min(1000, 500)
  expect_equal(labs$left$Y, 0.25)         # -0.25 + y_offset
})

test_that("make_strand_labels rejects negative offsets", {
  tx <- data.frame(start = 100, end = 300, strand = "+", y_start = 0)
  expect_error(make_strand_labels(tx, axis_pad_bp = -1), class = "rnapeaks_error_invalid_arg")
  expect_error(make_strand_labels(tx, axis_pad_bp = 500, x_offset = -1),
               class = "rnapeaks_error_invalid_arg")
})

# --- validators -----------------------------------------------------------

test_that("require_layout enforces presence, type, and positivity of heights", {
  expect_error(require_layout("x"), class = "rnapeaks_error_invalid_arg")
  expect_error(require_layout(list(center = 0, exon_height = 0.5)), class = "rnapeaks_error_invalid_arg")
  expect_error(require_layout(list(center = 0, exon_height = 0, utr_height = 0.3)),
               class = "rnapeaks_error_invalid_arg")
  expect_silent(require_layout(list(center = 0, exon_height = 0.5, utr_height = 0.3)))
})

test_that("require_transcript_cols rejects non-frames, empty frames, and missing columns", {
  expect_error(require_transcript_cols(list(), c("start")), class = "rnapeaks_error_invalid_arg")
  expect_error(require_transcript_cols(data.frame(start = numeric(0)), "start"),
               class = "rnapeaks_error_invalid_arg")
  expect_error(require_transcript_cols(data.frame(a = 1), c("start", "end")),
               class = "rnapeaks_error_invalid_arg")
})

test_that("require_window validates the region window", {
  expect_error(require_window("x"), class = "rnapeaks_error_invalid_arg")
  expect_error(require_window(list(start = 0)), class = "rnapeaks_error_invalid_arg")
  expect_error(require_window(list(start = "a", end = 1)), class = "rnapeaks_error_invalid_arg")
  expect_error(require_window(list(start = 100, end = 0)), class = "rnapeaks_error_invalid_arg")
  expect_silent(require_window(list(start = 0, end = 100)))
})

test_that("require_positive_scalar rejects negatives and non-scalars", {
  expect_error(require_positive_scalar(-1, "v"), class = "rnapeaks_error_invalid_arg")
  expect_error(require_positive_scalar(c(1, 2), "v"), class = "rnapeaks_error_invalid_arg")
  expect_silent(require_positive_scalar(0, "v"))
})
