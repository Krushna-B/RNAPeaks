#!/usr/bin/env Rscript

# Fast R rewrite of the Python Genomic Peak Analysis Tool, using IRanges
# for all interval algebra.
#
# Public API:
#   generate_control_peaks(
#     input  = "peaks.bed",         # file path OR data.frame/data.table
#     anno   = "annotation.bed",    # file path OR data.frame/data.table
#     gene   = "genes.bed",         # file path OR data.frame/data.table
#     pool   = 4,
#     output = "./result/",
#     seed   = 1234
#   )
#
# Output columns (no header):
#   chr, peak_start, peak_end, name, strand, control_start, control_end
#
# Coordinate convention
# ---------------------
# BED is 0-based half-open, so position count of [s, e) is e - s. The Python
# script uses portion's closed intervals but only ever reads
# upper - lower as the "length", which matches BED widths.
#
# Internally we store every interval as IRanges(s, e - 1L). Then:
#   * width(ir) == e - s        (matches Python lengths)
#   * intersect / setdiff / reduce match the BED interpretation
#   * To recover BED coordinates: bed_start = start(ir),
#     bed_end = end(ir) + 1L
#
# Random control choices may differ from Python because R and Python use
# different RNG internals; the algorithm and distribution are preserved.



REGIONS <- c("UTR3", "UTR5", "CDS", "exon")
MIN_REGION_LENGTH <- 50L
DEFAULT_EXTENSION_LENGTH <- 100L
SLIDING_WINDOW_SIZE <- 200L
SLIDING_WINDOW_STEP <- 50L

# -----------------------------
# IRanges helpers
# -----------------------------
# All `_bed`-suffixed helpers take or return BED-style (start, end) pairs.

empty_ir <- function() IRanges::IRanges()

ir_bed <- function(start, end) {
  # Build an IRanges from BED-style coordinates.
  start <- as.integer(start)
  end <- as.integer(end)
  keep <- !is.na(start) & !is.na(end) & (end > start)
  if (!any(keep)) return(empty_ir())
  IRanges::IRanges(start = start[keep], end = end[keep] - 1L)
}

ir_bed_starts <- function(ir) IRanges::start(ir)
ir_bed_ends <- function(ir) IRanges::end(ir) + 1L
ir_widths <- function(ir) IRanges::width(ir)

ir_intersect <- function(a, b) {
  if (length(a) == 0L || length(b) == 0L) return(empty_ir())
  IRanges::intersect(a, b)
}

ir_setdiff <- function(a, b) {
  if (length(a) == 0L) return(empty_ir())
  if (length(b) == 0L) return(IRanges::reduce(a))
  IRanges::setdiff(a, b)
}

ir_union <- function(a, b) {
  if (length(a) == 0L) return(IRanges::reduce(b))
  if (length(b) == 0L) return(IRanges::reduce(a))
  IRanges::union(a, b)
}

ir_contains_single <- function(start, end, container_ir) {
  # TRUE iff the single BED interval [start, end) lies entirely inside
  # container_ir.
  if (length(container_ir) == 0L) return(FALSE)
  q <- ir_bed(start, end)
  if (length(q) == 0L) return(FALSE)
  inter <- ir_intersect(q, container_ir)
  length(inter) == 1L &&
    IRanges::start(inter) == IRanges::start(q) &&
    IRanges::end(inter) == IRanges::end(q)
}

ir_overlaps_any <- function(start, end, ir) {
  if (length(ir) == 0L) return(FALSE)
  q <- ir_bed(start, end)
  if (length(q) == 0L) return(FALSE)
  IRanges::countOverlaps(q, ir) > 0L
}

# Flat-vector interval helpers. starts/ends are BED-style half-open and assumed
# sorted, non-overlapping (i.e. produced from a reduced IRanges).

flat_extract <- function(ir) {
  list(starts = IRanges::start(ir), ends = IRanges::end(ir) + 1L)
}

# Single (qs, qe) "any overlap" against sorted non-overlapping starts/ends.
flat_overlaps_any_scalar <- function(qs, qe, starts, ends) {
  n <- length(starts)
  if (n == 0L) return(FALSE)
  i <- findInterval(qs, starts)
  if (i >= 1L && ends[i] > qs) return(TRUE)
  if (i + 1L <= n && starts[i + 1L] < qe) return(TRUE)
  FALSE
}

# Vectorized "fully contained in some interval".
# Returns logical vector parallel to qs/qe.
flat_contained_batch <- function(qs, qe, starts, ends) {
  n <- length(starts)
  out <- logical(length(qs))
  if (n == 0L) return(out)
  i <- findInterval(qs, starts)
  hit <- i >= 1L
  if (!any(hit)) return(out)
  out[hit] <- qe[hit] <= ends[i[hit]]
  out
}

splice_slice <- function(splice_sites, python_start_index_zero_based) {
  if (length(splice_sites) == 0L) return(integer())
  first <- as.integer(python_start_index_zero_based) + 1L
  if (first > length(splice_sites)) return(integer())
  splice_sites[seq.int(first, length(splice_sites), by = 2L)]
}

# -----------------------------
# Input normalization helpers
# -----------------------------

is_file_like_input <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

coerce_to_dt <- function(x, label = "input") {
  if (is_file_like_input(x)) {
    return(data.table::fread(x, header = FALSE, sep = "\t", fill = TRUE, quote = ""))
  }
  if (is.data.frame(x) || data.table::is.data.table(x)) {
    return(data.table::as.data.table(data.table::copy(x)))
  }
  stop(label, " must be either a file path, data.frame, or data.table.")
}

# -----------------------------
# Input readers
# -----------------------------

read_peaks <- function(x) {
  raw <- coerce_to_dt(x, label = "input peaks")
  raw[, line_id := .I]

  named_required <- c("chr", "start", "end", "strand", "FC", "transcript")

  if (all(named_required %in% names(raw))) {
    peaks <- raw[
      strand %in% c("+", "-"),
      .(
        line_id = line_id,
        chr = as.character(chr),
        start = as.integer(start),
        end = as.integer(end),
        strand = as.character(strand),
        FC = as.character(FC),
        transcript = as.character(transcript)
      )
    ]
  } else {
    if (ncol(raw) < 12L) {
      stop(
        "Peak input must either have named columns ",
        paste(named_required, collapse = ", "),
        " or at least 12 columns matching the original peak-file format."
      )
    }

    peaks <- data.table::data.table(
      line_id = raw$line_id,
      chr = as.character(raw[[1L]]),
      start = as.integer(raw[[2L]]),
      end = as.integer(raw[[3L]]),
      strand = as.character(raw[[6L]]),
      FC = as.character(raw[[7L]]),
      transcript = as.character(raw[[12L]])
    )

    peaks <- peaks[strand %in% c("+", "-")]
  }

  peaks[, peak_range := paste0(start, "-", end)]

  # Match Python build_peak_dict: one record per (chr, peak_range), first FC
  # retained, all transcripts appended in file order.
  grouped <- peaks[
    ,
    .(
      first_line = min(line_id),
      start = start[1L],
      end = end[1L],
      FC = FC[1L],
      transcripts = list(transcript)
    ),
    by = .(chr, peak_range)
  ]

  data.table::setorder(grouped, first_line)
  grouped[]
}

read_annotation <- function(x) {
  raw <- coerce_to_dt(x, label = "annotation")
  raw[, line_id := .I]

  named_required <- c("chr", "start", "end", "name", "strand", "gene")

  if (all(named_required %in% names(raw))) {
    anno <- raw[
      ,
      .(
        line_id = line_id,
        chr = as.character(chr),
        start = as.integer(start),
        end = as.integer(end),
        name = as.character(name),
        strand = as.character(strand),
        gene = as.character(gene)
      )
    ]

    if ("transcript" %in% names(raw)) {
      anno[, transcript := as.character(raw$transcript)]
    }
    if ("region" %in% names(raw)) {
      anno[, region := as.character(raw$region)]
    }
  } else {
    if (ncol(raw) < 7L) {
      stop(
        "Annotation input must either have named columns ",
        paste(named_required, collapse = ", "),
        " or at least 7 columns matching the original annotation-file format."
      )
    }

    anno <- data.table::data.table(
      line_id = raw$line_id,
      chr = as.character(raw[[1L]]),
      start = as.integer(raw[[2L]]),
      end = as.integer(raw[[3L]]),
      name = as.character(raw[[4L]]),
      strand = as.character(raw[[6L]]),
      gene = as.character(raw[[7L]])
    )
  }

  if (!"transcript" %in% names(anno) || !"region" %in% names(anno)) {
    name_parts <- data.table::tstrsplit(anno$name, "_", fixed = TRUE, keep = 1:2)
    if (!"transcript" %in% names(anno)) {
      anno[, transcript := as.character(name_parts[[1L]])]
    }
    if (!"region" %in% names(anno)) {
      anno[, region := as.character(name_parts[[2L]])]
    }
  }

  anno[]
}

read_genes <- function(x) {
  raw <- coerce_to_dt(x, label = "gene annotation")
  raw[, line_id := .I]

  named_required <- c("chr", "start", "end", "gene")

  if (all(named_required %in% names(raw))) {
    genes <- raw[
      ,
      .(
        line_id = line_id,
        chr = as.character(chr),
        start = as.integer(start),
        end = as.integer(end),
        gene = as.character(gene)
      )
    ]
  } else {
    if (ncol(raw) < 4L) {
      stop(
        "Gene input must either have named columns ",
        paste(named_required, collapse = ", "),
        " or at least 4 columns: chr, start, end, gene."
      )
    }

    genes <- data.table::data.table(
      line_id = raw$line_id,
      chr = as.character(raw[[1L]]),
      start = as.integer(raw[[2L]]),
      end = as.integer(raw[[3L]]),
      gene = as.character(raw[[4L]])
    )
  }

  genes[]
}

# -----------------------------
# Annotation preprocessing
# -----------------------------

build_transcript_annotations <- function(anno_chr) {
  if (nrow(anno_chr) == 0L) return(list())

  # Capture file-order-last gene per transcript BEFORE sorting, to match
  # Python's `anno_dict[chrom][tx]["gene"] = fields[6]` overwrite-on-every-line
  # semantics.
  data.table::setorder(anno_chr, line_id)
  file_last <- anno_chr[!duplicated(transcript, fromLast = TRUE),
                       .(transcript, gene)]
  gene_by_tx_file <- setNames(file_last$gene, file_last$transcript)

  # One sort. Then everything is vectorized base-R from here on — no
  # data.table group operations, no per-transcript R loops calling [.data.table.
  data.table::setorder(anno_chr, transcript, region, start)

  n <- nrow(anno_chr)
  tx <- anno_chr$transcript
  rg <- anno_chr$region
  st <- anno_chr$start
  en <- anno_chr$end
  strands <- anno_chr$strand

  # Per-transcript first/last row indices (transcript-sorted, so contiguous).
  tx_change <- if (n == 1L) TRUE else c(TRUE, tx[-1L] != tx[-n])
  tx_first_idx <- which(tx_change)
  tx_last_idx <- c(tx_first_idx[-1L] - 1L, n)
  tx_names <- tx[tx_first_idx]
  n_tx <- length(tx_names)

  strand_by_tx <- setNames(strands[tx_first_idx], tx_names)
  gene_by_tx <- gene_by_tx_file[tx_names]

  # Reduce (start, end) per (transcript, region) using a flat vectorized sweep.
  in_regions <- rg %in% REGIONS
  if (!any(in_regions)) {
    red_tx_idx <- integer(0L)
    red_reg <- character(0L)
    red_rs <- integer(0L)
    red_re <- integer(0L)
  } else {
    sub_idx <- which(in_regions)
    tx_sub <- tx[sub_idx]
    rg_sub <- rg[sub_idx]
    s_sub <- st[sub_idx]
    e_sub <- en[sub_idx]
    m <- length(sub_idx)

    group_change <- if (m == 1L) TRUE else
      c(TRUE, tx_sub[-1L] != tx_sub[-m] | rg_sub[-1L] != rg_sub[-m])
    group_id <- cumsum(group_change)

    # Vectorized "cummax with reset at each group boundary" without a grouped
    # data.table call. Pattern: cummax(e) overall - boundary correction. Easier:
    # use a single R loop over groups via ave-style trick, or compute per-group.
    # The simplest reliable approach is one tight C-like R for-loop here.
    ce <- integer(m)
    ce[1L] <- e_sub[1L]
    if (m > 1L) {
      g_prev <- group_id[1L]
      for (i in 2L:m) {
        gi <- group_id[i]
        if (gi != g_prev) {
          ce[i] <- e_sub[i]
          g_prev <- gi
        } else {
          prev_ce <- ce[i - 1L]
          ce[i] <- if (e_sub[i] > prev_ce) e_sub[i] else prev_ce
        }
      }
    }

    ce_lag <- c(0L, ce[-m])
    new_block <- group_change | (s_sub > ce_lag)
    bs_idx <- which(new_block)
    nb <- length(bs_idx)
    be_idx <- if (nb == 1L) m else c(bs_idx[-1L] - 1L, m)

    red_tx <- tx_sub[bs_idx]
    red_reg <- rg_sub[bs_idx]
    red_rs <- s_sub[bs_idx]
    red_re <- ce[be_idx]
    red_tx_idx <- match(red_tx, tx_names)
  }

  # Bucket reduced intervals into per-transcript per-region lists.
  out <- vector("list", n_tx)
  names(out) <- tx_names

  # Pre-init each transcript with empty regions_flat.
  for (k in seq_len(n_tx)) {
    out[[k]] <- list(
      regions_flat = list(),
      region_order = character(0L),
      strand = strand_by_tx[[k]],
      gene = gene_by_tx[[k]],
      splice_sites = integer(0L),
      tx_start = NA_integer_,
      tx_end = NA_integer_
    )
  }

  # Walk reduced intervals; since red_tx_idx and red_reg are grouped by
  # (tx_idx, region), find boundaries and slice.
  if (length(red_tx_idx) > 0L) {
    keyA <- red_tx_idx
    keyB <- red_reg
    bn <- length(keyA)
    rb_change <- if (bn == 1L) TRUE else
      c(TRUE, keyA[-1L] != keyA[-bn] | keyB[-1L] != keyB[-bn])
    rb_first <- which(rb_change)
    rb_last <- c(rb_first[-1L] - 1L, bn)

    for (j in seq_along(rb_first)) {
      a <- rb_first[j]; b <- rb_last[j]
      tx_i <- keyA[a]
      reg <- keyB[a]
      out[[tx_i]]$regions_flat[[reg]] <- list(
        starts = red_rs[a:b],
        ends = red_re[a:b]
      )
    }
  }

  # Per-transcript splice sites from raw exon rows (anno_chr already sorted by tx).
  exon_mask <- rg == "exon"
  if (any(exon_mask)) {
    ex_tx <- tx[exon_mask]
    ex_s <- st[exon_mask]
    ex_e <- en[exon_mask]
    # Each tx contributes 2 sites per exon row; combine starts and ends.
    sites_combined <- c(ex_s, ex_e)
    tx_combined <- c(ex_tx, ex_tx)
    sites_split <- split(sites_combined, tx_combined)
    # Sort and attach.
    for (txn in names(sites_split)) {
      tx_i <- match(txn, tx_names)
      if (is.na(tx_i)) next
      sl <- sort(sites_split[[txn]])
      out[[tx_i]]$splice_sites <- sl
      out[[tx_i]]$tx_start <- sl[1L]
      out[[tx_i]]$tx_end <- sl[length(sl)]
      if (length(sl) > 2L) {
        mid <- sl[2L:(length(sl) - 1L)]
        i_s <- mid[seq.int(1L, length(mid), by = 2L)]
        i_e <- mid[seq.int(2L, length(mid), by = 2L)]
        out[[tx_i]]$regions_flat[["intron"]] <- flat_reduce(i_s, i_e)
      } else {
        out[[tx_i]]$regions_flat[["intron"]] <- list(starts = integer(0L),
                                                    ends = integer(0L))
      }
    }
  }

  # Finalize region_order per transcript.
  for (k in seq_len(n_tx)) {
    rn <- names(out[[k]]$regions_flat)
    out[[k]]$region_order <- c(REGIONS[REGIONS %in% rn], "intron"[("intron" %in% rn)])
  }

  out
}

build_gene_map <- function(genes_chr) {
  if (nrow(genes_chr) == 0L) return(list())

  data.table::setorder(genes_chr, line_id)

  # Match Python's last-write-wins: keep only the last row per gene name.
  # R's `[[` on a named list with duplicates returns the FIRST match, so we
  # have to drop the earlier duplicates explicitly.
  genes_chr <- genes_chr[!duplicated(gene, fromLast = TRUE)]

  out <- vector("list", nrow(genes_chr))
  for (i in seq_len(nrow(genes_chr))) {
    out[[i]] <- list(start = genes_chr$start[i], end = genes_chr$end[i])
  }
  setNames(out, genes_chr$gene)
}

# Intersect a single interval (a_s, a_e) with sorted-reduced (b_s, b_e).
# Returns flat list(starts, ends).
flat_intersect_one <- function(a_s, a_e, b_s, b_e) {
  if (length(b_s) == 0L) return(list(starts = integer(0L), ends = integer(0L)))
  ov <- b_s < a_e & b_e > a_s
  if (!any(ov)) return(list(starts = integer(0L), ends = integer(0L)))
  list(
    starts = pmax(b_s[ov], a_s),
    ends = pmin(b_e[ov], a_e)
  )
}

# Set difference: a minus b. Both inputs are sorted-reduced flat lists.
flat_setdiff <- function(a_s, a_e, b_s, b_e) {
  na <- length(a_s)
  if (na == 0L) return(list(starts = integer(0L), ends = integer(0L)))
  if (length(b_s) == 0L) return(list(starts = a_s, ends = a_e))

  # Result is at most na + length(b_s) intervals.
  cap <- na + length(b_s) + 1L
  out_s <- integer(cap)
  out_e <- integer(cap)
  n <- 0L

  for (i in seq_len(na)) {
    cur <- a_s[i]
    a_end <- a_e[i]
    # b intervals overlapping [cur, a_end)
    ov <- which(b_s < a_end & b_e > cur)
    if (length(ov) == 0L) {
      n <- n + 1L; out_s[n] <- cur; out_e[n] <- a_end
      next
    }
    for (jk in ov) {
      js <- b_s[jk]; je <- b_e[jk]
      if (js > cur) {
        n <- n + 1L; out_s[n] <- cur; out_e[n] <- js
      }
      if (je > cur) cur <- je
      if (cur >= a_end) break
    }
    if (cur < a_end) {
      n <- n + 1L; out_s[n] <- cur; out_e[n] <- a_end
    }
  }
  list(starts = out_s[seq_len(n)], ends = out_e[seq_len(n)])
}

# -----------------------------
# Overlap classification (pass 1)
# -----------------------------

find_overlap_for_peak <- function(peak_start, peak_end, tx) {
  # Returns a character vector of qualifying region names for the peak.
  # Non-intron region qualifies iff some single piece of the region fully
  # contains the peak (since pieces are non-overlapping and reduced, this is
  # equivalent to "fully contained in any merged region piece"). Intron
  # qualifies iff any overlap exists.
  out <- character()

  for (region in tx$region_order) {
    rf <- tx$regions_flat[[region]]
    if (is.null(rf) || length(rf$starts) == 0L) next

    if (region == "intron") {
      if (flat_overlaps_any_scalar(peak_start, peak_end, rf$starts, rf$ends)) {
        out <- c(out, region)
      }
      next
    }

    # Full coverage: peak fits entirely inside one piece.
    if (flat_contained_batch(peak_start, peak_end, rf$starts, rf$ends)[1L]) {
      out <- c(out, region)
    }
  }

  out
}

choose_peak_region <- function(region_names) {
  if ("CDS" %in% region_names) return("CDS")
  if ("UTR3" %in% region_names) return("UTR3")
  if ("UTR5" %in% region_names) return("UTR5")
  if ("exon" %in% region_names) return("exon")
  "intron"
}

# -----------------------------
# Control generation
# -----------------------------

sample_simple_control <- function(peak_start, peak_end, peak_length, cl_starts, cl_ends, final_region, rng) {
  # Python materializes every integer start in
  #   range(control_region[0], control_region[1] - peak_length)
  # We compute the same uniform-over-positions sample without enumerating.
  if (length(cl_starts) == 0L) return(NULL)

  starts <- cl_starts
  ends <- cl_ends
  n_starts <- ends - peak_length - starts

  keep <- n_starts > 0L
  if (!any(keep)) return(NULL)

  starts <- starts[keep]
  n_starts <- n_starts[keep]

  total <- sum(n_starts)
  chosen_global <- rng$pick(total)

  cumulative <- cumsum(n_starts)
  chosen_row <- which(cumulative >= chosen_global)[1L]
  previous_total <- if (chosen_row == 1L) 0L else cumulative[chosen_row - 1L]
  offset <- chosen_global - previous_total - 1L

  control_start <- starts[chosen_row] + offset
  control_end <- control_start + peak_length

  list(
    peak_start = as.integer(peak_start),
    peak_end = as.integer(peak_end),
    control_start = as.integer(control_start),
    control_end = as.integer(control_end),
    region_type = final_region
  )
}

sample_uniform_flat <- function(starts, ends, rng) {
  if (length(starts) == 0L) return(NULL)
  i <- rng$pick(length(starts))
  list(start = starts[i], end = ends[i])
}

# Merge a small unsorted set of candidate intervals (starts/ends) by sorting
# and walking. Returns a list(starts, ends) of disjoint, sorted intervals.
flat_reduce <- function(s, e) {
  if (length(s) == 0L) return(list(starts = integer(0L), ends = integer(0L)))
  ord <- order(s, e)
  s <- s[ord]; e <- e[ord]
  out_s <- integer(length(s))
  out_e <- integer(length(s))
  n <- 0L
  cur_s <- s[1L]; cur_e <- e[1L]
  if (length(s) > 1L) {
    for (k in 2L:length(s)) {
      if (s[k] <= cur_e) {
        if (e[k] > cur_e) cur_e <- e[k]
      } else {
        n <- n + 1L
        out_s[n] <- cur_s; out_e[n] <- cur_e
        cur_s <- s[k]; cur_e <- e[k]
      }
    }
  }
  n <- n + 1L
  out_s[n] <- cur_s; out_e[n] <- cur_e
  list(starts = out_s[seq_len(n)], ends = out_e[seq_len(n)])
}

get_cds_control <- function(peak_start, peak_end, peak_length, cl_starts, cl_ends, rf, tx, rng) {
  # rf = tx$regions_flat[["CDS"]] (sorted, reduced).
  # Python's _get_cds_control iterates region_list and overwrites exon_coords
  # on each hit, so the LAST overlapping CDS interval wins.
  hit <- which(rf$starts < peak_end & rf$ends > peak_start)
  if (length(hit) == 0L) return(NULL)
  k <- hit[length(hit)]
  exon_coords <- c(rf$starts[k], rf$ends[k])

  dist_to_start <- as.integer(peak_start) - exon_coords[1L]
  dist_to_end <- exon_coords[2L] - as.integer(peak_end)

  if (dist_to_start < dist_to_end) {
    splice_sites <- splice_slice(tx$splice_sites, 1L)
    distance <- dist_to_start
  } else {
    splice_sites <- splice_slice(tx$splice_sites, 0L)
    distance <- dist_to_end
  }

  if (length(splice_sites) == 0L) return(NULL)

  control_pos <- as.integer(splice_sites) + as.integer(distance)
  cand_starts <- control_pos
  cand_ends <- control_pos + peak_length

  ok <- flat_contained_batch(cand_starts, cand_ends, cl_starts, cl_ends)
  if (!any(ok)) return(NULL)

  picked <- flat_reduce(cand_starts[ok], cand_ends[ok])
  chosen <- sample_uniform_flat(picked$starts, picked$ends, rng)
  if (is.null(chosen)) return(NULL)

  list(
    peak_start = as.integer(peak_start),
    peak_end = as.integer(peak_end),
    control_start = as.integer(chosen$start),
    control_end = as.integer(chosen$end),
    region_type = "CDS"
  )
}

get_intron_control <- function(peak_start, peak_end, peak_length, cl_starts, cl_ends, rf, tx, rng) {
  # rf = tx$regions_flat[["intron"]].
  region_coords <- NULL
  intersect_length <- 0L
  is_exon_overlap <- FALSE
  is_5_prime_ss <- FALSE

  # Last-overlapping intron wins.
  if (length(rf$starts) > 0L) {
    hit <- which(rf$starts < peak_end & rf$ends > peak_start)
    if (length(hit) > 0L) {
      k <- hit[length(hit)]
      region_coords <- c(rf$starts[k], rf$ends[k])
      intersect_length <- min(peak_end, rf$ends[k]) - max(peak_start, rf$starts[k])
    }
  }

  if (intersect_length == 0L) {
    erf <- tx$regions_flat[["exon"]]
    if (!is.null(erf) && length(erf$starts) > 0L) {
      hit <- which(erf$starts < peak_end & erf$ends > peak_start)
      if (length(hit) > 0L) {
        k <- hit[length(hit)]
        region_coords <- c(erf$starts[k], erf$ends[k])
        intersect_length <- min(peak_end, erf$ends[k]) - max(peak_start, erf$starts[k])
        is_exon_overlap <- TRUE
      }
    }
  }

  if (intersect_length == 0L || is.null(region_coords)) return(NULL)

  dist_to_start <- as.integer(peak_start) - region_coords[1L]
  dist_to_end <- region_coords[2L] - as.integer(peak_end)

  final_region <- NULL
  splice_sites <- integer()
  distance <- NA_integer_

  if (peak_length == intersect_length) {
    if (dist_to_start < dist_to_end) {
      if (is_exon_overlap) {
        splice_sites <- splice_slice(tx$splice_sites, 0L)
        final_region <- if (tx$strand == "+") "3pexon" else "5pexon"
        distance <- dist_to_start
      } else {
        splice_sites <- splice_slice(tx$splice_sites, 1L)
        final_region <- if (tx$strand == "+") "5pintron" else "3pintron"
        distance <- dist_to_start
      }
    } else {
      if (is_exon_overlap) {
        splice_sites <- splice_slice(tx$splice_sites, 1L)
        final_region <- if (tx$strand == "+") "5pexon" else "3pexon"
        distance <- dist_to_end
        is_5_prime_ss <- TRUE
      } else {
        splice_sites <- splice_slice(tx$splice_sites, 2L)
        final_region <- if (tx$strand == "+") "3pintron" else "5pintron"
        distance <- dist_to_end
        is_5_prime_ss <- TRUE
      }
    }
  } else {
    # Python uses two independent ifs here; both can fire.
    if (dist_to_start < 0L) {
      splice_sites <- splice_slice(tx$splice_sites, 1L)
      final_region <- if (tx$strand == "+") "5pss" else "3pss"
      distance <- dist_to_start
    }
    if (dist_to_end < 0L) {
      splice_sites <- splice_slice(tx$splice_sites, 2L)
      final_region <- if (tx$strand == "+") "3pss" else "5pss"
      distance <- dist_to_end
      is_5_prime_ss <- TRUE
    }
  }

  if (is.null(final_region) || length(splice_sites) == 0L || is.na(distance)) {
    return(NULL)
  }

  if (is_5_prime_ss) {
    control_pos <- as.integer(splice_sites) - as.integer(distance)
    cand_starts <- control_pos - peak_length
    cand_ends <- control_pos
  } else {
    control_pos <- as.integer(splice_sites) + as.integer(distance)
    cand_starts <- control_pos
    cand_ends <- control_pos + peak_length
  }

  ok <- flat_contained_batch(cand_starts, cand_ends, cl_starts, cl_ends)
  if (!any(ok)) return(NULL)

  picked <- flat_reduce(cand_starts[ok], cand_ends[ok])
  chosen <- sample_uniform_flat(picked$starts, picked$ends, rng)
  if (is.null(chosen)) return(NULL)

  list(
    peak_start = as.integer(peak_start),
    peak_end = as.integer(peak_end),
    control_start = as.integer(chosen$start),
    control_end = as.integer(chosen$end),
    region_type = final_region
  )
}

get_control <- function(peak_start, peak_end, region, cl_starts, cl_ends, tx, rng) {
  peak_start <- as.integer(peak_start)
  peak_end <- as.integer(peak_end)
  peak_length <- peak_end - peak_start

  if (peak_length < MIN_REGION_LENGTH) return(NULL)
  if (peak_length <= 0L) return(NULL)
  if (length(cl_starts) == 0L) return(NULL)

  rf <- tx$regions_flat[[region]]
  if (is.null(rf)) rf <- list(starts = integer(0L), ends = integer(0L))

  if (region %in% c("UTR3", "UTR5", "exon")) {
    return(sample_simple_control(
      peak_start, peak_end, peak_length, cl_starts, cl_ends, region, rng
    ))
  }
  if (region == "CDS") {
    return(get_cds_control(
      peak_start, peak_end, peak_length, cl_starts, cl_ends, rf, tx, rng
    ))
  }
  if (region == "intron") {
    return(get_intron_control(
      peak_start, peak_end, peak_length, cl_starts, cl_ends, rf, tx, rng
    ))
  }

  NULL
}

# Sliding window matching Python sliding_window():
# step == window_size; final block emitted only if remaining length > min_len;
# otherwise the final block re-emits the previous window_str (a Python bug we
# faithfully preserve).
sliding_window_python <- function(start, end,
                                  window_size = SLIDING_WINDOW_SIZE,
                                  min_len = SLIDING_WINDOW_STEP) {
  start <- as.integer(start)
  end <- as.integer(end)

  seq_len_total <- end - start + 1L
  if (seq_len_total <= 0L) {
    return(list(window_start = integer(), window_end = integer()))
  }

  out_start <- integer()
  out_end <- integer()

  last_start <- NA_integer_
  last_end <- NA_integer_

  i <- 0L
  while (i < seq_len_total) {
    if (i + window_size > seq_len_total) {
      if (seq_len_total - i > min_len) {
        last_start <- start + i
        last_end <- end
      }
    } else {
      last_start <- start + i
      last_end <- start + i + window_size
    }
    out_start <- c(out_start, last_start)
    out_end <- c(out_end, last_end)
    i <- i + window_size
  }

  list(window_start = out_start, window_end = out_end)
}

# -----------------------------
# Per-chromosome processing
# -----------------------------

process_chromosome <- function(chromosome, chrom_index, peaks, anno, genes, seed, rng) {
  # In R-native mode we reseed per chromosome so parallel mclapply workers
  # don't share state. In Python-RNG mode the RNG is a single global stream
  # owned by the embedded Python interpreter, so we skip per-chrom reseeding
  # to mirror Python's behavior.
  if (identical(rng$mode, "R")) {
    set.seed(as.integer(seed) + as.integer(chrom_index))
  }

  peaks_chr <- peaks[chr == chromosome]
  anno_chr <- anno[chr == chromosome]
  genes_chr <- genes[chr == chromosome]

  if (nrow(peaks_chr) == 0L || nrow(anno_chr) == 0L) {
    return(data.table::data.table())
  }

  tx_annotations <- build_transcript_annotations(data.table::copy(anno_chr))
  gene_map <- build_gene_map(data.table::copy(genes_chr))

  # Pass 1: build peak_region_dict and total peak IRanges per strand.
  n_peaks <- nrow(peaks_chr)
  peak_starts <- peaks_chr$start
  peak_ends <- peaks_chr$end
  peak_ranges <- peaks_chr$peak_range
  peak_tx_lists <- peaks_chr$transcripts

  peak_region_dict <- vector("list", n_peaks)
  names(peak_region_dict) <- peak_ranges
  peak_order <- character(0L)
  peak_order_cap <- n_peaks
  peak_order_n <- 0L
  peak_order_buf <- character(peak_order_cap)

  plus_starts <- integer(0L); plus_ends <- integer(0L)
  minus_starts <- integer(0L); minus_ends <- integer(0L)
  plus_buf_starts <- integer(n_peaks); plus_buf_ends <- integer(n_peaks)
  minus_buf_starts <- integer(n_peaks); minus_buf_ends <- integer(n_peaks)
  plus_n <- 0L; minus_n <- 0L

  for (p_i in seq_len(n_peaks)) {
    ps <- peak_starts[p_i]
    pe <- peak_ends[p_i]
    pr <- peak_ranges[p_i]
    transcripts <- peak_tx_lists[[p_i]]

    if (is.na(ps) || is.na(pe) || pe <= ps) next
    peak_length <- pe - ps

    for (tx_name in transcripts) {
      tx <- tx_annotations[[tx_name]]
      if (is.null(tx)) next

      if (tx$strand == "+") {
        plus_n <- plus_n + 1L
        plus_buf_starts[plus_n] <- ps
        plus_buf_ends[plus_n] <- pe
      } else {
        minus_n <- minus_n + 1L
        minus_buf_starts[minus_n] <- ps
        minus_buf_ends[minus_n] <- pe
      }

      qualifying <- find_overlap_for_peak(ps, pe, tx)
      if (length(qualifying) == 0L) next

      bucket <- peak_region_dict[[pr]]
      if (is.null(bucket)) {
        bucket <- list()
        peak_order_n <- peak_order_n + 1L
        peak_order_buf[peak_order_n] <- pr
      }
      for (rt in qualifying) {
        bucket[[rt]] <- c(bucket[[rt]], tx_name)
      }
      peak_region_dict[[pr]] <- bucket
    }
  }

  if (plus_n > 0L) {
    plus_starts <- plus_buf_starts[seq_len(plus_n)]
    plus_ends <- plus_buf_ends[seq_len(plus_n)]
  }
  if (minus_n > 0L) {
    minus_starts <- minus_buf_starts[seq_len(minus_n)]
    minus_ends <- minus_buf_ends[seq_len(minus_n)]
  }

  tpp <- if (plus_n > 0L) flat_reduce(plus_starts, plus_ends) else list(starts = integer(0L), ends = integer(0L))
  tpm <- if (minus_n > 0L) flat_reduce(minus_starts, minus_ends) else list(starts = integer(0L), ends = integer(0L))

  if (peak_order_n == 0L) {
    return(data.table::data.table())
  }
  peak_order <- peak_order_buf[seq_len(peak_order_n)]

  tpp_starts <- tpp$starts; tpp_ends <- tpp$ends
  tpm_starts <- tpm$starts; tpm_ends <- tpm$ends

  # Running total_control per strand maintained as sorted flat vectors.
  # Inserts are O(N) via append; checks are O(log N) via findInterval.
  tcp_starts <- integer(0L); tcp_ends <- integer(0L)
  tcm_starts <- integer(0L); tcm_ends <- integer(0L)

  peak_start_map <- setNames(peaks_chr$start, peaks_chr$peak_range)
  peak_end_map <- setNames(peaks_chr$end, peaks_chr$peak_range)
  peak_fc_map <- setNames(peaks_chr$FC, peaks_chr$peak_range)

  # Output buffers, pre-allocated to peak_order length (one row per accepted peak).
  res_n <- 0L
  res_chr <- character(peak_order_n)
  res_pstart <- integer(peak_order_n)
  res_pend <- integer(peak_order_n)
  res_name <- character(peak_order_n)
  res_strand <- character(peak_order_n)
  res_cstart <- integer(peak_order_n)
  res_cend <- integer(peak_order_n)

  processed_windows <- new.env(parent = emptyenv(), hash = TRUE, size = 1024L)

  for (peak_range in peak_order) {
    bucket <- peak_region_dict[[peak_range]]
    peak_region <- choose_peak_region(names(bucket))

    peak_start <- as.integer(peak_start_map[[peak_range]])
    peak_end <- as.integer(peak_end_map[[peak_range]])
    peak_fc <- as.character(peak_fc_map[[peak_range]])
    peak_length <- peak_end - peak_start

    peak_transcripts <- bucket[[peak_region]]

    for (tx_name in peak_transcripts) {
      tx <- tx_annotations[[tx_name]]
      if (is.null(tx)) next

      strand_key <- if (tx$strand == "+") "plus" else "minus"

      if (strand_key == "plus") {
        tp_starts <- tpp_starts; tp_ends <- tpp_ends
      } else {
        tp_starts <- tpm_starts; tp_ends <- tpm_ends
      }

      gene_name <- tx$gene
      gene_rec <- gene_map[[gene_name]]
      if (!is.null(gene_rec)) {
        bl <- flat_intersect_one(gene_rec$start, gene_rec$end, tp_starts, tp_ends)
        bl_starts <- bl$starts; bl_ends <- bl$ends
      } else {
        bl_starts <- tp_starts; bl_ends <- tp_ends
      }

      if (peak_region == "intron") {
        if (is.na(tx$tx_start) || is.na(tx$tx_end)) next
        ri_starts <- tx$tx_start; ri_ends <- tx$tx_end
      } else {
        rf <- tx$regions_flat[[peak_region]]
        if (is.null(rf)) { ri_starts <- integer(0L); ri_ends <- integer(0L) }
        else { ri_starts <- rf$starts; ri_ends <- rf$ends }
      }

      cl <- flat_setdiff(ri_starts, ri_ends, bl_starts, bl_ends)
      cl_starts <- cl$starts; cl_ends <- cl$ends

      if (peak_length <= 400L) {
        win_starts <- peak_start
        win_ends <- peak_end
      } else {
        ws <- sliding_window_python(peak_start, peak_end,
                                    SLIDING_WINDOW_SIZE, SLIDING_WINDOW_STEP)
        win_starts <- ws$window_start
        win_ends <- ws$window_end
      }

      peak_processed <- FALSE

      for (w_i in seq_along(win_starts)) {
        window_start <- win_starts[w_i]
        window_end <- win_ends[w_i]
        window_key <- paste0(window_start, "-", window_end)

        if (exists(window_key, envir = processed_windows, inherits = FALSE)) next

        control <- get_control(window_start, window_end, peak_region, cl_starts, cl_ends, tx, rng)
        if (is.null(control)) next
        if (is.na(control$control_start) || control$control_start == 0L) next

        cs <- control$control_start
        ce <- control$control_end

        if (flat_overlaps_any_scalar(cs, ce, tp_starts, tp_ends)) next
        if (strand_key == "plus") {
          if (flat_overlaps_any_scalar(cs, ce, tcp_starts, tcp_ends)) next
        } else {
          if (flat_overlaps_any_scalar(cs, ce, tcm_starts, tcm_ends)) next
        }

        # Accept: insert into the sorted flat vectors.
        if (strand_key == "plus") {
          pos <- findInterval(cs, tcp_starts)
          tcp_starts <- append(tcp_starts, cs, after = pos)
          tcp_ends <- append(tcp_ends, ce, after = pos)
        } else {
          pos <- findInterval(cs, tcm_starts)
          tcm_starts <- append(tcm_starts, cs, after = pos)
          tcm_ends <- append(tcm_ends, ce, after = pos)
        }

        assign(window_key, cs, envir = processed_windows)

        res_n <- res_n + 1L
        res_chr[res_n] <- chromosome
        res_pstart[res_n] <- control$peak_start
        res_pend[res_n] <- control$peak_end
        res_name[res_n] <- paste(gene_name, control$region_type, peak_fc, sep = "_")
        res_strand[res_n] <- tx$strand
        res_cstart[res_n] <- cs
        res_cend[res_n] <- ce

        peak_processed <- TRUE
      }

      if (peak_processed) break
    }
  }

  if (res_n == 0L) return(data.table::data.table())

  data.table::data.table(
    chr = res_chr[seq_len(res_n)],
    peak_start = as.character(res_pstart[seq_len(res_n)]),
    peak_end = as.character(res_pend[seq_len(res_n)]),
    name = res_name[seq_len(res_n)],
    strand = res_strand[seq_len(res_n)],
    control_start = as.character(res_cstart[seq_len(res_n)]),
    control_end = as.character(res_cend[seq_len(res_n)])
  )
}

# -----------------------------
# Main user-facing function
# -----------------------------

make_rng <- function(mode, seed) {
  mode <- match.arg(mode, c("R", "python"))
  if (mode == "python") {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("rng = 'python' requires the 'reticulate' package. ",
           "Install with install.packages('reticulate').")
    }
    py_random <- reticulate::import("random", convert = TRUE)
    py_random$seed(as.integer(seed))
    list(
      mode = "python",
      # Returns an R 1-based index in {1..n}. random.choice(seq) is
      # equivalent to seq[random.randrange(len(seq))].
      pick = function(n) {
        py_random$randrange(as.integer(n)) + 1L
      }
    )
  } else {
    set.seed(as.integer(seed))
    list(
      mode = "R",
      pick = function(n) sample.int(as.integer(n), 1L)
    )
  }
}

generate_control_peaks <- function(
    input,
    anno,
    gene,
    pool = 1L,
    output = "./result/",
    output_file = NULL,
    output_basename = NULL,
    seed = 1234L,
    rng = c("R", "python"),
    return_results = TRUE,
    write_output = TRUE,
    verbose = TRUE) {
  start_time <- Sys.time()

  if (missing(input) || is.null(input)) {
    stop("You must provide an input peak file path or peak data frame.")
  }
  if (missing(anno) || is.null(anno)) {
    stop("You must provide an annotation file path or annotation data frame.")
  }
  if (missing(gene) || is.null(gene)) {
    stop("You must provide a gene annotation file path or gene data frame.")
  }

  rng_mode <- match.arg(rng)
  pool <- as.integer(pool)
  seed <- as.integer(seed)
  if (is.na(pool) || pool < 1L) pool <- 1L
  if (is.na(seed)) seed <- 1234L

  # Python RNG is a single global stream held by the embedded interpreter.
  # mclapply forks would each get a private copy of that interpreter state
  # and diverge from Python's reference behavior, so force sequential mode.
  if (rng_mode == "python" && pool > 1L) {
    if (verbose) {
      message("rng = 'python' forces pool = 1 to preserve a single RNG stream.")
    }
    pool <- 1L
  }

  rng_state <- make_rng(rng_mode, seed)

  if (verbose) message("Building annotation dict...")
  anno_dt <- read_annotation(anno)
  if (verbose) {
    message(sprintf("Annotation table loaded: %.2f s",
                    as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  }

  genes_dt <- read_genes(gene)
  if (verbose) {
    message(sprintf("Gene table loaded: %.2f s",
                    as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  }

  peaks_dt <- read_peaks(input)
  if (verbose) {
    message(sprintf("Peak table loaded: %.2f s",
                    as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  }

  chromosomes <- unique(peaks_dt$chr)
  chromosomes <- chromosomes[chromosomes %in% unique(anno_dt$chr)]

  if (length(chromosomes) == 0L) {
    warning("No chromosomes from the peak input were found in the annotation input.")
    total_results <- data.table::data.table()
  } else {
    cores <- pool
    if (.Platform$OS.type == "windows" && cores > 1L) {
      if (verbose) message("parallel::mclapply is not multicore on Windows; using 1 core.")
      cores <- 1L
    }

    if (verbose) message("Processing peaks...")

    worker <- function(k) {
      process_chromosome(
        chromosome = chromosomes[k],
        chrom_index = k,
        peaks = peaks_dt,
        anno = anno_dt,
        genes = genes_dt,
        seed = seed,
        rng = rng_state
      )
    }

    if (cores > 1L) {
      result_list <- parallel::mclapply(seq_along(chromosomes), worker,
                                        mc.cores = cores, mc.preschedule = TRUE)
    } else {
      result_list <- lapply(seq_along(chromosomes), worker)
    }

    total_results <- data.table::rbindlist(result_list, fill = TRUE)
  }

  if (verbose) {
    message(sprintf("Peak processing completed: %.2f s",
                    as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  }

  if (write_output) {
    dir.create(output, showWarnings = FALSE, recursive = TRUE)
    if (is.null(output_file)) {
      if (is.null(output_basename)) {
        if (is_file_like_input(input)) {
          output_basename <- tools::file_path_sans_ext(basename(input))
        } else {
          output_basename <- "control_peaks"
        }
      }
      output_file <- file.path(output, paste0(output_basename, "_control.bed"))
    }

    data.table::fwrite(total_results, file = output_file,
                       sep = "\t", col.names = FALSE, quote = FALSE, na = "")

    if (verbose) {
      message("Analysis completed successfully!")
      message("Output written to: ", output_file)
    }
  } else {
    output_file <- NULL
    if (verbose) {
      message("Analysis completed successfully!")
      message("No output file written because write_output = FALSE.")
    }
  }

  if (verbose) {
    message(sprintf("Total time: %.2f s",
                    as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  }

  if (return_results) return(total_results)
  invisible(output_file)
}
