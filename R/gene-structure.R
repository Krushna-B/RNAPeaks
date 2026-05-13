# Feature types for protein coding transcripts
UTR_TYPES                    <- c("five_prime_utr", "three_prime_utr")
PROTEIN_CODING_FEATURE_TYPES <- c("CDS", UTR_TYPES)

# Columns the transcript df wil carry
TRANSCRIPT_COLS <- c("seqnames", "start", "end", "width", "strand", "type",
                     "gene_id", "gene_name", "gene_biotype",
                     "transcript_id", "transcript_name",
                     "y_start", "y_end")





build_gene_structure <- function(transcript, layout) {
  #Validate transcript df cols
  require_transcript_cols(transcript,
                          c("seqnames", "start", "end", "strand", "type"))

  #Validate layout list params
  require_layout(layout)

  #Identify biotype of transcript either protein coding or not
  biotype <- stats::na.omit(transcript$gene_biotype)[1]
  if (length(biotype) != 1 || is.na(biotype)) biotype <- "other"

  #Select which cols to filter out based on biotype
  feature_types <- switch(biotype,
    protein_coding = PROTEIN_CODING_FEATURE_TYPES,
    "exon"
  )

  #Select for features
  features <- transcript[transcript$type %in% feature_types, , drop = FALSE]
  if (!nrow(features)) return(NULL)

  #Build y-axis coordinates for UTRs and Exons
  half_h <- ifelse(features$type %in% UTR_TYPES,
                   layout$utr_height  / 2,
                   layout$exon_height / 2)
  features$y_start <- layout$center - half_h
  features$y_end   <- layout$center + half_h

  #Keep only required rows
  features <- features[, intersect(TRANSCRIPT_COLS, colnames(features)), drop = FALSE]

  #Add intron positions which are computed
  rbind(features, intron_rows(features, layout$center, layout$exon_height))
}


#Compute intron positions between exon and utr features
intron_rows <- function(features, center, exon_height) {
  #Return out if only 1 feature
  if (nrow(features) <= 1) return(features[0, , drop = FALSE])

  #Order features and ensure start > end
  features <- features[order(features$start, features$end), , drop = FALSE]
  n          <- nrow(features)
  prev_end   <- features$end[-n]
  next_start <- features$start[-1]
  keep       <- next_start > prev_end

  #If no features return
  if (!any(keep)) return(features[0, , drop = FALSE])

  #Computer intron start and end positions
  half    <- exon_height / 8
  introns <- features[rep(1, sum(keep)), , drop = FALSE]
  introns$type    <- "intron"
  introns$start   <- prev_end[keep] + 1
  introns$end     <- next_start[keep] - 1

  #Compute intron y-coordinates
  introns$y_start <- center - half
  introns$y_end   <- center + half
  if ("width" %in% colnames(introns)) {
    introns$width <- introns$end - introns$start + 1
  }
  row.names(introns) <- NULL
  return(introns)
}




# Place 5'/3' tags just outside the gene.
make_strand_labels <- function(transcript,
                               offset_frac = 0.05,
                               offset_min  = 100,
                               offset_max  = 1000) {
  require_transcript_cols(transcript,
                          c("start", "end", "strand", "y_start"))
  require_positive_scalar(offset_frac, "offset_frac")
  require_positive_scalar(offset_min,  "offset_min")
  require_positive_scalar(offset_max,  "offset_max")
  if (offset_max < offset_min) {
    abort_invalid_arg(c(
      "{.arg offset_max} must be >= {.arg offset_min}.",
      "x" = "Got {.code offset_min} = {.val {offset_min}}, {.code offset_max} = {.val {offset_max}}."
    ))
  }

  gene_min <- min(transcript$start)
  gene_max <- max(transcript$end)
  y_pos    <- min(transcript$y_start)
  offset   <- clamp((gene_max - gene_min + 1) * offset_frac,
                    offset_min, offset_max)

  labs <- if (identical(transcript$strand[1], "+")) c("5'", "3'") else c("3'", "5'")
  list(
    left  = data.frame(Label = labs[1], X = gene_min - offset, Y = y_pos),
    right = data.frame(Label = labs[2], X = gene_max + offset, Y = y_pos)
  )
}

#Helper for Strand labels
clamp <- function(x, lo, hi) min(max(x, lo), hi)



#Validate transcript df cols
require_transcript_cols <- function(transcript, cols) {
  if (!is.data.frame(transcript)) {
    abort_invalid_arg(c(
      "{.arg transcript} must be a data frame.",
      "x" = "You supplied {.cls {class(transcript)[1]}}."
    ))
  }
  if (nrow(transcript) == 0L) {
    abort_invalid_arg("{.arg transcript} has no rows.")
  }
  missing <- setdiff(cols, colnames(transcript))
  if (length(missing)) {
    abort_invalid_arg(c(
      "{.arg transcript} is missing required column{?s}: {.field {missing}}.",
      "i" = "Required: {.field {cols}}."
    ))
  }
}

#Validate layout list params
require_layout <- function(layout) {
  if (!is.list(layout)) {
    abort_invalid_arg(c(
      "{.arg layout} must be a list.",
      "x" = "You supplied {.cls {class(layout)[1]}}."
    ))
  }
  required <- c("center", "exon_height", "utr_height")
  missing  <- setdiff(required, names(layout))
  if (length(missing)) {
    abort_invalid_arg(c(
      "{.arg layout} is missing key{?s}: {.field {missing}}.",
      "i" = "Required: {.field {required}}."
    ))
  }
  for (k in required) {
    v <- layout[[k]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      key_label <- paste0("layout$", k)
      abort_invalid_arg(c(
        "{.code {key_label}} must be a single non-NA numeric value.",
        "x" = "Got {.cls {class(v)[1]}} of length {length(v)}."
      ))
    }
  }
  if (layout$exon_height <= 0 || layout$utr_height <= 0) {
    abort_invalid_arg(c(
      "{.code layout$exon_height} and {.code layout$utr_height} must be positive.",
      "x" = "Got {.val {layout$exon_height}} and {.val {layout$utr_height}}."
    ))
  }
}

# Verify a value is a single non-NA non-negative numeric.
require_positive_scalar <- function(value, arg) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < 0) {
    abort_invalid_arg(c(
      "{.arg {arg}} must be a single non-negative number.",
      "x" = "Got {.cls {class(value)[1]}} = {.val {value}}."
    ))
  }
}



