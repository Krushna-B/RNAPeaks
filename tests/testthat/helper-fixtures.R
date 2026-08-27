# Shared fixtures for RNAPeaks tests. testthat sources helper-*.R before any
# test file, so these constructors are available everywhere.

# Raw BED exactly as accepted by check_bed(): six positional columns
# (chr, start, end, name, score, strand). `protein` fills the name column (pos 4)
# so it can double as a split_col target
make_raw_bed <- function(n = 3,
                         chr = "chr1",
                         start = seq(100, by = 200, length.out = n),
                         end = start + 50,
                         protein = "SRSF1",
                         score = 100,
                         strand = "+") {
  data.frame(
    chr = rep(chr, length.out = n),
    start = start,
    end = end,
    name = rep(protein, length.out = n),
    score = rep(score, length.out = n),
    strand = rep(strand, length.out = n),
    stringsAsFactors = FALSE
  )
}

# Validated BED shaped like check_bed()'s output: canonical coordinate columns
# plus a `target` column
make_checked_bed <- function(chr = "1",
                            start = c(100, 300, 500),
                            end = start + 50,
                            strand = "+",
                            target = "SRSF1") {
  n <- length(start)
  data.frame(
    chr = rep(chr, length.out = n),
    start = as.numeric(start),
    end = as.numeric(end),
    strand = rep(strand, length.out = n),
    target = rep(target, length.out = n),
    stringsAsFactors = FALSE
  )
}
