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

# Canonical (already chr-normalized) GTF data frame for gene / selection tests.
# Two + strand genes on chr "1": SRSF1 (a long transcript ENST...1 and a shorter
# ENST...2) and U2AF2 (ENST...3) in a separate window; plus a - strand gene
# MINUS (ENST...4). Carries transcript_name so name-based lookup can be tested.
make_gtf <- function() {
  row <- function(type, start, end, strand, gid, gname, tid, tname, width) {
    data.frame(seqnames = "1", start = start, end = end, strand = strand,
               type = type, gene_id = gid, gene_name = gname,
               transcript_id = tid, transcript_name = tname, width = width,
               stringsAsFactors = FALSE)
  }
  rbind(
    # SRSF1, + strand: long transcript (width 1000) and short transcript (width 400)
    row("transcript", 100, 1100, "+", "ENSG00000136450", "SRSF1", "ENST00000001", "SRSF1-201", 1000),
    row("exon",       100,  300, "+", "ENSG00000136450", "SRSF1", "ENST00000001", "SRSF1-201", 200),
    row("exon",       900, 1100, "+", "ENSG00000136450", "SRSF1", "ENST00000001", "SRSF1-201", 200),
    row("UTR",        100,  150, "+", "ENSG00000136450", "SRSF1", "ENST00000001", "SRSF1-201", 50),
    row("transcript", 100,  500, "+", "ENSG00000136450", "SRSF1", "ENST00000002", "SRSF1-202", 400),
    row("exon",       100,  500, "+", "ENSG00000136450", "SRSF1", "ENST00000002", "SRSF1-202", 400),
    # U2AF2, + strand, separate window
    row("transcript", 5000, 5600, "+", "ENSG00000063244", "U2AF2", "ENST00000003", "U2AF2-201", 600),
    row("exon",       5000, 5600, "+", "ENSG00000063244", "U2AF2", "ENST00000003", "U2AF2-201", 600),
    # MINUS, - strand
    row("transcript", 2000, 2400, "-", "ENSG00000000009", "MINUS", "ENST00000004", "MINUS-201", 400),
    row("exon",       2000, 2400, "-", "ENSG00000000009", "MINUS", "ENST00000004", "MINUS-201", 400)
  )
}

# Minimal on-disk GTF for the plot entry points: gene AAA (ENSG1) with one + strand
# transcript ENST1 on chr 1, exons 100-400 and 800-1100. Returns the file path.
# Import requires rtracklayer, so guard callers with skip_if_not_installed().
write_min_gtf <- function() {
  path <- tempfile(fileext = ".gtf")
  writeLines(c(
    '1\tsrc\ttranscript\t100\t1100\t.\t+\t.\tgene_id "ENSG1"; gene_name "AAA"; transcript_id "ENST1";',
    '1\tsrc\texon\t100\t400\t.\t+\t.\tgene_id "ENSG1"; gene_name "AAA"; transcript_id "ENST1";',
    '1\tsrc\texon\t800\t1100\t.\t+\t.\tgene_id "ENSG1"; gene_name "AAA"; transcript_id "ENST1";'
  ), path)
  path
}
