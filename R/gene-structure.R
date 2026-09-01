# Feature types for protein coding transcripts
#UTR Type after gtf normalization
UTR_TYPES                    <- "UTR"
PROTEIN_CODING_FEATURE_TYPES <- c("CDS", UTR_TYPES)

# Columns the transcript df wil carry
TRANSCRIPT_COLS <- c("seqnames", "start", "end", "width", "strand", "type",
                     "gene_id", "gene_name", "gene_biotype",
                     "transcript_id", "transcript_name",
                     "y_start", "y_end")



#' Build a single-transcript structure data frame.
#'
#' @param transcript GTF rows for one transcript with `seqnames`, `start`,
#'   `end`, `strand`, `type`, and (for protein coding) `gene_biotype`.
#' @param layout Named list with `center`, `exon_height`, `utr_height`.
#' @return Data frame of exon/UTR/intron rows with `y_start` and `y_end`
#'   coordinates, or `NULL` if no plottable features are present.
#' @noRd
#' @family gene
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


#Build a multi-transcript region structure stacked vertically
#  - One lane per transcript_id
build_region_structure <- function(transcripts, layout, window) {
  require_transcript_cols(transcripts,
                          c("seqnames", "start", "end", "strand", "type",
                            "gene_id", "transcript_id"))
  require_layout(layout)
  require_window(window)
  if (is.null(layout$gene_lane_height)) {
    abort_invalid_arg("{.code layout$gene_lane_height} is required for region structures.")
  }

  # Order lanes left-to-right by each transcript's leftmost feature start
  transcript_ids <- unique(stats::na.omit(transcripts$transcript_id))
  if (!length(transcript_ids)) {
    abort_invalid_arg("{.arg transcripts} has no usable {.field transcript_id} values.")
  }
  tx_starts <- vapply(transcript_ids, function(t) {
    min(transcripts$start[transcripts$transcript_id == t])
  }, numeric(1))
  transcript_ids <- transcript_ids[order(tx_starts)]

  # Build each transcript at its stacked center; lane i sits at
  # center + i*lane_height
  lane_height <- layout$gene_lane_height
  parts <- lapply(seq_along(transcript_ids), function(i) {
    rows <- transcripts[!is.na(transcripts$transcript_id) &
                        transcripts$transcript_id == transcript_ids[i],
                        , drop = FALSE]
    per_layout <- layout
    per_layout$center <- layout$center + (i - 1L) * lane_height
    build_gene_structure(rows, per_layout)
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(transcripts[0, , drop = FALSE])

  region <- do.call(rbind, parts)

  # Clip features to the user window so the plot shows exactly what was asked
  clip_to_window(region, window$start, window$end)
}


# Trim feature start/end to [start, end] and drop rows that fall fully outside.
# Width is recomputed.
clip_to_window <- function(df, start, end) {
  if (!nrow(df)) return(df)
  df$start <- pmax(df$start, start)
  df$end   <- pmin(df$end,   end)
  df <- df[df$end >= df$start, , drop = FALSE]
  if ("width" %in% names(df)) df$width <- df$end - df$start + 1L
  df
}


# Validate window list passed to build_region_structure()
require_window <- function(window) {
  if (!is.list(window)) {
    abort_invalid_arg(c(
      "{.arg window} must be a list with {.field start} and {.field end}.",
      "x" = "Got {.cls {class(window)[1]}}."
    ))
  }
  missing <- setdiff(c("start", "end"), names(window))
  if (length(missing)) {
    abort_invalid_arg(c(
      "{.arg window} is missing required field{?s}: {.field {missing}}."
    ))
  }
  for (k in c("start", "end")) {
    v <- window[[k]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
      abort_invalid_arg(c(
        "{.code window${k}} must be a single non-NA numeric value.",
        "x" = "Got {.cls {class(v)[1]}} of length {length(v)}."
      ))
    }
  }
  if (window$start > window$end) {
    abort_invalid_arg(c(
      "{.code window$start} must be <= {.code window$end}.",
      "x" = "Got start = {.val {window$start}}, end = {.val {window$end}}."
    ))
  }
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
  keep       <- next_start > prev_end +1L

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


# Place 5'/3' tags just outside the gene. x_offset is measured from the gene
# ends and capped at axis_pad_bp so labels stay inside the padded panel.
make_strand_labels <- function(transcript, axis_pad_bp, x_offset = axis_pad_bp, y_offset = 0) {
  require_transcript_cols(transcript,
                          c("start", "end", "strand", "y_start"))
  require_positive_scalar(axis_pad_bp, "axis_pad_bp")
  require_positive_scalar(x_offset, "x_offset")

  transcript$strand <- as.character(transcript$strand)

  gene_min <- min(transcript$start)
  gene_max <- max(transcript$end)
  # Anchor at the box bottom so the tags sit slightly below the (centered)
  # gene name rather than overlapping it. Nudge with strand_label_y_offset.
  y_pos    <- min(transcript$y_start) + y_offset
  off      <- min(x_offset, axis_pad_bp)

  labs <- if (identical(transcript$strand[1], "+")) c("5'", "3'") else c("3'", "5'")
  list(
    left  = data.frame(Label = labs[1], X = gene_min - off, Y = y_pos),
    right = data.frame(Label = labs[2], X = gene_max + off, Y = y_pos)
  )
}



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




