#' Per-base coverage for a BAM file over a genomic region.
#' @family gene
#' @noRd
compute_bam_coverage <- function(bam_path, chr, start, end) {

  #check if the bam file is valid and can be opened
  if (!is.character(bam_path) || length(bam_path) != 1L || is.na(bam_path) ||
      !nzchar(bam_path)) {
    abort_invalid_arg("{.arg bam_path} must be a single non-empty file path.")
  }
  if (!file.exists(bam_path)) {
    abort_not_found(c(
      "BAM file does not exist.",
      "x" = "Path: {.path {bam_path}}."
    ))
  }
  bai_path <- paste0(bam_path, ".bai")
  if (!file.exists(bai_path)) {
    abort_not_found(c(
      "BAM index (.bai) not found.",
      "x" = "Expected: {.path {bai_path}}.",
      "i" = "Index the BAM first, e.g. with {.code Rsamtools::indexBam()} or {.code samtools index}."
    ))
  }

  bam_file <- tryCatch(
    Rsamtools::BamFile(bam_path),
    error = function(e) abort_invalid_arg(c(
      "Failed to open {.arg bam_path} as a BAM file.",
      "x" = "{conditionMessage(e)}",
      "i" = "Path: {.path {bam_path}}."
    ))
  )

  bam_chroms <- GenomeInfoDb::seqnames(
    Rsamtools::seqinfo(Rsamtools::BamFile(bam_path))
  )

  #check chromosome name convention
  chr_resolved <- if (chr %in% bam_chroms) {
    chr
  } else if (paste0("chr", chr) %in% bam_chroms) {
    paste0("chr", chr)
  } else if (sub("^chr", "", chr) %in% bam_chroms) {
    sub("^chr", "", chr)
  } else {
    abort_not_found(
      "Chromosome {.val {chr}} not in BAM header {.path {bam_path}}."
    )
  }

  #ranges o fregion of interest
  roi <- GenomicRanges::GRanges(
    seqnames = chr_resolved,
    ranges   = IRanges::IRanges(start = start, end = end)
  )

  #loading genomic alignments
  ga  <- GenomicAlignments::readGAlignments(
    bam_path, param = Rsamtools::ScanBamParam(which = roi)
  )
  cov <- GenomicAlignments::coverage(ga)


  if (!chr_resolved %in% names(cov)) {
    return(data.frame(pos = seq(start, end), coverage = 0L))
  }

  #if the end of the region is less than the region
  data.frame(
    pos      = seq(start, end),
    coverage = as.numeric(cov[[chr_resolved]][start:end])
  )
}


#' Build BAM ribbon / label / scale data frames for draw_plot().]
#' @family gene
#' @noRd
prepare_bam_tracks <- function(bam_files, chr, start, end, base_y, style) {
  if (is.null(bam_files) || !length(bam_files)) return(NULL)
  bam_files <- resolve_bam_names(bam_files)

  covs <- lapply(bam_files, compute_bam_coverage,
                 chr = chr, start = start, end = end)

  if (!is.null(style$bam_ylim)) {
    y_lo <- style$bam_ylim[1]
    y_hi <- style$bam_ylim[2]
  } else {
    y_lo <- 0
    y_hi <- max(vapply(covs, function(d) max(d$coverage, na.rm = TRUE), numeric(1)))
    if (y_hi <= y_lo) y_hi <- y_lo + 1
  }

  band_h <- style$bam_track_height
  n      <- length(bam_files)
  nms    <- names(bam_files)

  rib <- lab <- sc <- vector("list", n)
  for (i in seq_len(n)) {
    floor_y   <- base_y + (i - 1L) * band_h
    ceiling_y <- floor_y + band_h
    d         <- covs[[i]]
    clipped   <- pmin(pmax(d$coverage, y_lo), y_hi)

    d$track <- nms[i]
    d$y_bot <- floor_y
    d$y_top <- floor_y + (clipped - y_lo) / (y_hi - y_lo) * band_h

    rib[[i]] <- d[, c("pos", "track", "y_bot", "y_top")]
    lab[[i]] <- data.frame(track = nms[i], y = (floor_y + ceiling_y) / 2,
                           stringsAsFactors = FALSE)
    sc[[i]]  <- data.frame(
      track = nms[i],
      y     = c(floor_y, ceiling_y),
      text  = c(as.character(y_lo), scales::comma(y_hi)),
      stringsAsFactors = FALSE
    )
  }

  list(
    ribbons      = do.call(rbind, rib),
    labels       = do.call(rbind, lab),
    scales       = do.call(rbind, sc),
    total_height = n * band_h
  )
}


# Fill in missing names from the filename.
resolve_bam_names <- function(bam_files) {
  nm <- names(bam_files)
  if (is.null(nm)) nm <- rep("", length(bam_files))
  blank <- is.na(nm) | !nzchar(nm)
  nm[blank] <- tools::file_path_sans_ext(basename(bam_files[blank]))
  names(bam_files) <- nm
  bam_files
}
