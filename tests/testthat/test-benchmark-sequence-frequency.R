# Benchmark tests for calculate_sequence_frequency
# Run with: testthat::test_file("tests/testthat/test-benchmark-sequence-frequency.R")
#
# These tests compare performance of parallel vs sequential processing.
# They require BSgenome.Hsapiens.UCSC.hg38 to be installed.

# Skip benchmarks in automated testing (they take time)
skip_on_cran()
skip_on_ci()

# Check if required packages are available
skip_if_not_installed("bench")
skip_if_not_installed("future")
skip_if_not_installed("future.apply")
skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

test_that("benchmark: calculate_sequence_frequency parallel vs sequential", {

  # Load required packages

  library(BSgenome.Hsapiens.UCSC.hg38)
  library(bench)

  # Create mock SEMATS-like data for testing
  # Simulating ~1000 events (4000 bins) for reasonable benchmark time
  set.seed(42)
  n_events <- 1000

  # Generate random genomic coordinates on chr1
  chr1_length <- 248956422
  random_starts <- sort(sample(1e6:(chr1_length - 1e6), n_events))

  mock_semats <- data.frame(
    chr = rep("chr1", n_events),
    strand = sample(c("+", "-"), n_events, replace = TRUE),
    GeneID = paste0("GENE", seq_len(n_events)),
    upstreamES = random_starts,
    upstreamEE = random_starts + 100,
    exonStart_0base = random_starts + 400,
    exonEnd = random_starts + 500,
    downstreamES = random_starts + 800,
    downstreamEE = random_starts + 900,
    PValue = runif(n_events, 0, 0.1),
    FDR = runif(n_events, 0, 0.1),
    IncLevelDifference = runif(n_events, -0.5, 0.5),
    stringsAsFactors = FALSE
  )

  # Build bins matrix
  bins_gr <- make_bins_matrix(mock_semats,
                               WidthIntoExon = 50,
                               WidthIntoIntron = 250)

  genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
  bin_width <- 300
  sequence <- "YGCY"

  # Ensure we start with sequential plan
  future::plan(future::sequential)

  # Run benchmark
  message("\n=== Benchmarking calculate_sequence_frequency ===")
  message(sprintf("Events: %d, Bins: %d, Sequence: %s\n", n_events, length(bins_gr), sequence))

  bm <- bench::mark(
    sequential = {
      calculate_sequence_frequency(bins_gr, sequence, genome, bin_width, cores = 1)
    },
    parallel_2cores = {
      calculate_sequence_frequency(bins_gr, sequence, genome, bin_width, cores = 2)
    },

    parallel_4cores = {
      calculate_sequence_frequency(bins_gr, sequence, genome, bin_width, cores = 4)
    },
    iterations = 3,
    check = TRUE,  # Verify all methods return same result
    filter_gc = FALSE
  )

  # Print results
  message("\n=== Benchmark Results ===")
  print(bm[, c("expression", "min", "median", "mem_alloc", "n_itr")])

  # Basic assertions - all should complete without error

  expect_true(all(bm$n_itr >= 1))

  # Clean up
 future::plan(future::sequential)
})


test_that("benchmark: scaling with number of events", {

  skip("Long running benchmark - run manually")

  library(BSgenome.Hsapiens.UCSC.hg38)
  library(bench)

  genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
  bin_width <- 300
  sequence <- "YGCY"

  # Test different sizes
  event_sizes <- c(500, 1000, 2000, 8000)

  results <- list()

  for (n_events in event_sizes) {
    message(sprintf("\n--- Testing with %d events ---", n_events))

    set.seed(42)
    chr1_length <- 248956422
    random_starts <- sort(sample(1e6:(chr1_length - 1e6), n_events))

    mock_semats <- data.frame(
      chr = rep("chr1", n_events),
      strand = sample(c("+", "-"), n_events, replace = TRUE),
      GeneID = paste0("GENE", seq_len(n_events)),
      upstreamES = random_starts,
      upstreamEE = random_starts + 100,
      exonStart_0base = random_starts + 400,
      exonEnd = random_starts + 500,
      downstreamES = random_starts + 800,
      downstreamEE = random_starts + 900,
      PValue = runif(n_events, 0, 0.1),
      FDR = runif(n_events, 0, 0.1),
      IncLevelDifference = runif(n_events, -0.5, 0.5),
      stringsAsFactors = FALSE
    )

    bins_gr <- make_bins_matrix(mock_semats, WidthIntoExon = 50, WidthIntoIntron = 250)

    future::plan(future::sequential)

    bm <- bench::mark(
      sequential = {
        calculate_sequence_frequency(bins_gr, sequence, genome, bin_width, cores = 1)
      },
      parallel_4cores = {
        calculate_sequence_frequency(bins_gr, sequence, genome, bin_width, cores = 4)
      },
      iterations = 2,
      check = TRUE,
      filter_gc = FALSE
    )

    bm$n_events <- n_events
    results[[as.character(n_events)]] <- bm
  }

  # Combine and print results
  all_results <- do.call(rbind, results)
  message("\n=== Scaling Results ===")
  print(all_results[, c("expression", "n_events", "median", "mem_alloc")])

  expect_true(TRUE)  # Just check it runs

  future::plan(future::sequential)
})


test_that("correctness: parallel and sequential return identical results", {

  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

  library(BSgenome.Hsapiens.UCSC.hg38)

  # Small dataset for quick correctness check
  set.seed(123)
  n_events <- 100

  chr1_length <- 248956422
  random_starts <- sort(sample(1e6:(chr1_length - 1e6), n_events))

  mock_semats <- data.frame(
    chr = rep("chr1", n_events),
    strand = sample(c("+", "-"), n_events, replace = TRUE),
    GeneID = paste0("GENE", seq_len(n_events)),
    upstreamES = random_starts,
    upstreamEE = random_starts + 100,
    exonStart_0base = random_starts + 400,
    exonEnd = random_starts + 500,
    downstreamES = random_starts + 800,
    downstreamEE = random_starts + 900,
    PValue = runif(n_events, 0, 0.1),
    FDR = runif(n_events, 0, 0.1),
    IncLevelDifference = runif(n_events, -0.5, 0.5),
    stringsAsFactors = FALSE
  )

  bins_gr <- make_bins_matrix(mock_semats, WidthIntoExon = 50, WidthIntoIntron = 250)
  genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
  bin_width <- 300

  # Run both methods
  future::plan(future::sequential)

  result_seq <- calculate_sequence_frequency(bins_gr, "YGCY", genome, bin_width, cores = 1)
  result_par <- calculate_sequence_frequency(bins_gr, "YGCY", genome, bin_width, cores = 2)

  # Check they're identical

  expect_equal(result_seq$global_position, result_par$global_position)
  expect_equal(result_seq$match_count, result_par$match_count)

  future::plan(future::sequential)
})
