# Internal function to pair computed intron intervals with plotting Y coordinates
Compute_Intron_Positions <- function(Gene_s, Intron_s) {

  # Compute genomic intron intervals from the ordered gene features
  Introns_Positions <- get_Intron_Positions(Gene_s)

  if (is.null(Introns_Positions) || nrow(Introns_Positions) == 0L) {
    return(NULL)
  }

  # Copy y_start / y_end from provided Intron_s from Gene Structure Data Frame
  Introns_Positions$y_start <- Intron_s$y_start
  Introns_Positions$y_end <- Intron_s$y_end

  # Midline in y for drawing intron baselines/arrows
  Intron_s$mid_y <- (Intron_s$y_start + Intron_s$y_end) / 2

  list(Introns_Positions = Introns_Positions, Intron_s = Intron_s)
}


# ----Helper Functions-------
# Get Intron positions from the gaps between coding features
get_Intron_Positions <- function(gene_df) {
  # Keep only rows relevant to exon/intron structure
  types <- c("intron", "exon", "CDS", "five_prime_utr", "three_prime_utr",
             "five_prime_UTR", "three_prime_UTR")
  gene_df <- gene_df[gene_df$type %in% types, ]

  # Order features by genomic position
  gene_df <- gene_df[order(gene_df$start, gene_df$end), ]

  # Blank row
  empty_introns <- function() data.frame(
    seqnames = character(0),
    strand = character(0),
    intron_id = integer(0),
    start = integer(0),
    end = integer(0),
    length_bp = integer(0),
    stringsAsFactors = FALSE
  )

  # Need at least two features to have a gap (intron)
  n <- nrow(gene_df)
  if (n < 2) {
    return(NULL)
  }

  # Gap between end of feature n and start of n+1
  start <- gene_df$end[-n]
  end_prev <- gene_df$start[-1]

  # Keep only proper gaps
  keep <- end_prev > start

  # Build intron intervals as closed ranges: [end_n, start_n+1]
  data.frame(
    seqnames = gene_df$seqnames[1],
    strand = gene_df$strand[1],
    intron_id = seq_len(n - 1)[keep],
    start = start[keep] + 1,
    end = end_prev[keep] - 1,
    length_bp = (end_prev[keep] - start[keep]),
    row.names = NULL,
    check.names = FALSE
  )
}
