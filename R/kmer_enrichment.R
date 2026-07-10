#' Compare k-mer enrichment between two sets of sequences
#'
#' Counts every length-`k` k-mer across two sets of sequences, normalizes each
#' set to a relative frequency, and reports the per-k-mer difference
#' `freq_a - freq_b`. Each set may be gene ids, transcript ids, BED peaks, or a
#' mix of gene / transcript ids.
#'
#' @param set_a,set_b The two sets to compare. Each is a BED data frame, a list
#'   of BED data frames, a `.bed` path, or a character vector of gene /
#'   transcript ids.
#' @param k K-mer length.
#' @param gtf Optional path to a local GTF file. If `NULL`, the bundled
#'   annotation for `species` is used
#' @param species One of `"hg38"`, `"mm10"`, or `"mm39"`. Selects both the
#'   bundled GTF and the default genome. Ignored where `gtf` / `genome` are
#'   supplied.
#' @param genome A BSgenome object or a genome key (`"hg38"` / `"mm10"` /
#'   `"mm39"`) used for sequence extraction. Defaults to `species`.
#' @param label_a,label_b Display names for the two sets, used in the table
#'   context and plot axes.
#' @param top_n Number of most strongly enriched k-mers (by absolute
#'   difference) to label in the scatter plot
#' @param title Plot title.
#'
#' @return A list with `table` (a data frame sorted by `difference`, one row per
#'   k-mer with `kmer`, `freq_a`, `freq_b`, `count_a`, `count_b`, `difference`,
#'   and `rank`) and `plots` (a list of two ggplots: `scatter` and `rank`).
#'
#' @export
kmer_enrichment <- function(set_a, set_b, k,
                            gtf     = NULL,
                            species = "hg38",
                            genome  = NULL,
                            label_a = "Set A",
                            label_b = "Set B",
                            top_n   = 20,
                            title   = "") {
  wrap_sm_errors("k-mer enrichment", {
    #Ensure required args
    if (missing(set_a) || missing(set_b)) {
      absent <- c(if (missing(set_a)) "set_a", if (missing(set_b)) "set_b")
      abort_invalid_arg(c(
        "Required argument{?s} missing: {.arg {absent}}.",
        "i" = "Each set is a BED, a {.file .bed} path, or a vector of gene / transcript ids."
      ))
    }
    if (missing(k)) {
      abort_invalid_arg(c(
        "{.arg k} is required.",
        "i" = "Supply a k-mer length, e.g. {.val 4}."
      ))
    }

    #Validate scalars, Cap k at 12, 4^12.
    check_scalar_int(k, "k", min = 1L, max = 12L)
    check_scalar_int(top_n, "top_n", min = 0L)
    check_string(label_a, "label_a")
    check_string(label_b, "label_b")
    check_string(title, "title")
    species <- normalize_str(species)

    #Resolve the genome once
    genome <- .resolve_genome(genome %||% species)

    #Only touch the GTF if a set if id's are provided
    needs_gtf <- !.is_bed_set(set_a) || !.is_bed_set(set_b)
    if (needs_gtf) {
      cli::cli_progress_step("Loading GTF")
      gtf <- get_GTF(species = species, file = gtf)
    }

    #Turn each set into sequences, then into normalized k-mer frequencies
    cli::cli_progress_step("Resolving sequences for {.val {label_a}}")
    seqs_a <- .resolve_kmer_set(set_a, gtf, genome, "set_a")
    cli::cli_progress_step("Resolving sequences for {.val {label_b}}")
    seqs_b <- .resolve_kmer_set(set_b, gtf, genome, "set_b")

    cli::cli_progress_step("Counting {k}-mers")
    kf_a <- .kmer_freq(seqs_a, k)
    kf_b <- .kmer_freq(seqs_b, k)
    cli::cli_progress_done()

    #Assemble the sorted enrichment table and the two plots.
    tbl <- .build_kmer_table(kf_a, kf_b)
    plots <- list(
      scatter = .plot_kmer_scatter(tbl, label_a, label_b, top_n, title),
      rank    = .plot_kmer_rank(tbl, label_a, label_b, title)
    )

    list(table = tbl, plots = plots)
  })
}


#Check's if set is a BED input
#Criteria: If a df, or list of df, or file path
.is_bed_set <- function(set) {
  if (is.data.frame(set)) return(TRUE)
  if (is.list(set) && !is.data.frame(set) && length(set) > 0L &&
      all(vapply(set, is.data.frame, logical(1L)))) {
    return(TRUE)
  }
  if (is.character(set) && length(set) == 1L && !is.na(set) &&
      file.exists(set)) {
    return(TRUE)
  }
  FALSE
}


#Resolve one set to a DNAStringSet of sequences
.resolve_kmer_set <- function(set, gtf, genome, arg) {
  if (.is_bed_set(set)) {
    gr <- .bed_set_to_granges(set)
  } else if (is.character(set) && length(set) > 0L && !anyNA(set) &&
             all(nzchar(set))) {
    gr <- .ids_to_span_granges(set, gtf)
  } else {
    abort_invalid_arg(c(
      "{.arg {arg}} must be a BED, a {.file .bed} path, or a non-empty vector of gene / transcript ids.",
      "x" = "Got {.cls {class(set)[1]}}."
    ))
  }

  #Pull the sequences for the GRanges

  gr <- .align_seqnames_to_genome(gr, genome) #Normalize function
  if (length(gr) == 0L) {
    abort_not_found(c(
      "No sequences could be extracted for {.arg {arg}}.",
      "i" = "Check that {.arg genome} matches the chromosome naming of the input."
    ))
  }
  #Return to Seqs
  BSgenome::getSeq(genome, gr)
}


#Combine one or more BEDs into a single reduced GRanges of peaks
.bed_set_to_granges <- function(set) {
  #One data frame or path goes straight through the shared BED reader
  if (is.data.frame(set) || is.character(set)) {
    return(.peaks_to_granges(set))
  }
  #For list,read each, concatenate, then collapse overlaps
  grs <- lapply(set, .peaks_to_granges)
  GenomicRanges::reduce(do.call(c, grs))
}


#Turn a vector of gene / transcript ids into one pre-mRNA span per id
.ids_to_span_granges <- function(ids, gtf) {
  #Collect the region bound by each id into parallel vectors, then build a single GRanges
  seqn <- character(length(ids))
  strd <- character(length(ids))
  lo   <- integer(length(ids))
  hi   <- integer(length(ids))

  for (i in seq_along(ids)) {
    id <- ids[i]
    #Ensembl transcript ids are used directly; everything else resolves
    # to the gene's longest transcript inside select_transcript
    is_tx <- grepl("^(ENST|ENSMUST)\\d", toupper(id))
    rows  <- if (is_tx) {
      select_transcript(gtf, TxID = id)
    } else {
      select_transcript(gtf, geneID = id)
    }
    seqn[i] <- as.character(rows$seqnames[1])
    strd[i] <- as.character(rows$strand[1])
    lo[i]   <- min(rows$start)
    hi[i]   <- max(rows$end)
  }

  #Create GRanges
  GenomicRanges::GRanges(
    seqnames = seqn,
    ranges   = IRanges::IRanges(start = lo, end = hi),
    strand   = strd
  )
}


#Count every k-mer over the sequences and normalize to a relative frequency.
.kmer_freq <- function(seqs, k) {
  if (length(seqs) == 0L) {
    abort_not_found("No sequences to count k-mers over.")
  }
  total <- Biostrings::oligonucleotideFrequency(seqs, width = k,
                                                simplify.as = "collapsed")
  denom <- sum(total)
  if (denom == 0L) {
    abort_not_found(c(
      "No valid {k}-mers found in the sequences.",
      "i" = "The extracted sequences may be shorter than {.arg k} or all ambiguous."
    ))
  }
  list(count = total, freq = total / denom)
}


#Join the two frequency vectors into the sorted enrichment table
.build_kmer_table <- function(kf_a, kf_b) {
  kmers <- names(kf_a$freq)
  tbl <- data.frame(
    kmer       = kmers,
    freq_a     = as.numeric(kf_a$freq),
    freq_b     = as.numeric(kf_b$freq[kmers]),
    count_a    = as.numeric(kf_a$count),
    count_b    = as.numeric(kf_b$count[kmers]),
    difference = as.numeric(kf_a$freq) - as.numeric(kf_b$freq[kmers]),
    stringsAsFactors = FALSE
  )
  #Most enriched in set A - B at the top
  tbl <- tbl[order(-tbl$difference), , drop = FALSE]
  tbl$rank <- seq_len(nrow(tbl))
  rownames(tbl) <- NULL
  tbl
}


#Scatter of per-k-mer frequencies with a y = x reference line
#Points above the diagonal are enriched in set A, below in set B
.plot_kmer_scatter <- function(tbl, label_a, label_b, top_n, title) {
  labels <- if (top_n > 0L) {
    tbl[order(-abs(tbl$difference)), , drop = FALSE][seq_len(min(top_n, nrow(tbl))), , drop = FALSE]
  } else {
    tbl[0, , drop = FALSE]
  }

  ggplot2::ggplot(tbl, ggplot2::aes(x = freq_b, y = freq_a)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey60") +
    ggplot2::geom_point(color = "#3B6FB6", alpha = 0.6, size = 1.4) +
    ggplot2::geom_text(
      data    = labels,
      mapping = ggplot2::aes(label = kmer),
      size    = 3, vjust = -0.6, check_overlap = TRUE
    ) +
    ggplot2::labs(
      x     = sprintf("%s k-mer frequency", label_b),
      y     = sprintf("%s k-mer frequency", label_a),
      title = title
    ) +
    ggplot2::theme_classic()
}


#Rank curve: k-mers ordered by enrichment against their difference
#Shows the spread from the most A-enriched down to the most B-enriched.
.plot_kmer_rank <- function(tbl, label_a, label_b, title) {
  ggplot2::ggplot(tbl, ggplot2::aes(x = rank, y = difference)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey60") +
    ggplot2::geom_point(color = "#3B6FB6", alpha = 0.6, size = 1.2) +
    ggplot2::labs(
      x     = "k-mer rank (by enrichment)",
      y     = sprintf("frequency difference (%s - %s)", label_a, label_b),
      title = title
    ) +
    ggplot2::theme_classic()
}
