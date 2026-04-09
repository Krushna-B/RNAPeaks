
#' Compute per-base read coverage for a BAM file over a genomic region.
#'
#' @param bam_path  Character. Path to a sorted, indexed BAM file (.bai must exist).
#' @param chr       Character. Chromosome name as it appears in the BAM header (e.g. "chr12" or "12").
#' @param start     Integer. Region start position (1-based).
#' @param end       Integer. Region end position (1-based, inclusive).
#
#' @return A data.frame with columns:
#'   \describe{
#'    \item{pos}{Integer genomic position.}
#'     \item{coverage}{Integer read depth at that position.}
#'   }
Compute_BAM_Coverage <- function(bam_path, chr, start, end) {

  # Resolve chromosome name against the BAM header — handles "19" vs "chr19"
  bam_seqlevels <- Rsamtools::seqinfo(Rsamtools::BamFile(bam_path))
  bam_chroms    <- GenomeInfoDb::seqnames(bam_seqlevels)

  chr_resolved <- if (chr %in% bam_chroms) {
    chr
  } else if (paste0("chr", chr) %in% bam_chroms) {
    paste0("chr", chr)
  } else if (sub("^chr", "", chr) %in% bam_chroms) {
    sub("^chr", "", chr)
  } else {
    stop(sprintf(
      "Chromosome '%s' not found in BAM header. Available: %s",
      chr, paste(head(bam_chroms, 10), collapse = ", ")
    ))
  }

  roi <- GenomicRanges::GRanges(
    seqnames = chr_resolved,
    ranges   = IRanges::IRanges(start = start, end = end)
  )

  param <- Rsamtools::ScanBamParam(which = roi)
  ga    <- GenomicAlignments::readGAlignments(bam_path, param = param)

  cov     <- GenomicAlignments::coverage(ga)
  chr_key <- chr_resolved

  if (!chr_key %in% names(cov)) {
    # Return zero coverage if chromosome not found in BAM
    return(data.frame(pos = seq(start, end), coverage = 0L))
  }

  cov_vec <- as.numeric(cov[[chr_key]][start:end])

  data.frame(
    pos      = seq(start, end),
    coverage = cov_vec
  )
}


# Derive a display label from a BAM file path (filename without extension).
BAM_Label_From_Path <- function(bam_path) {
  tools::file_path_sans_ext(basename(bam_path))
}
