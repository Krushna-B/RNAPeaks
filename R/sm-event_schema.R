#' Internal: event-type schemas

#' Schema slots:
#' \describe{
#'   \item{name}{Human-readable label.}
#'   \item{required_cols}{Character vector of columns the input table must
#'     contain. Universal columns plus event-type anchors.}
#'   \item{n_regions}{Number of regions per event.}
#'   \item{region_width}{`function(width_exon, width_intron)` returning the
#'     nominal bp width of each region. Same for every region within a
#'     schema.}
#'   \item{build_regions}{`function(events, width_exon, width_intron)`
#'     returning a `GRanges` of `nrow(events) * n_regions` ranges in
#'     event-then-region order, carrying `event_id` and `region_idx` in
#'     `mcols`.
#'   \item{region_labels}{Length `n_regions` character vector of facet labels
#'     (transcript order, 5' -> 3').}
#'   \item{schematic}{String key for the plotter's bottom-strip
#'     drawing}
#'   \item{group_set}{`c("Negative","Positive","Control")` for trio-mode
#'     event types, or `"Single"` for single-distribution modes (e.g. UTR).}
#' }
#'
#' @keywords internal
#' @name event_schema



# Columns every event-table schema must provide.
.event_required_base <- c(
  "chr", "strand", "GeneID",
  "PValue", "FDR", "IncLevelDifference"
)


#' @rdname event_schema
#' @noRd
event_schema_se <- list(
  name          = "Skipped Exon",

  required_cols = c(
    .event_required_base,
    "exonStart_0base", "exonEnd",
    "upstreamES", "upstreamEE",
    "downstreamES", "downstreamEE"
  ),

  n_regions     = 4L,

  region_width  = function(width_exon, width_intron) {
    as.integer(width_exon + width_intron + 1L)
  },

  build_regions = function(events, width_exon, width_intron) {
    n <- nrow(events)
    if (n == 0L) {
      return(GenomicRanges::GRanges())
    }

    upS <- events$upstreamES
    upE <- events$upstreamEE
    eS  <- events$exonStart_0base
    eE  <- events$exonEnd
    dnS <- events$downstreamES
    dnE <- events$downstreamEE

    # R1: upstream exon end into 1st intron
    s1 <- pmax(upE - width_exon,   upS); e1 <- pmin(upE + width_intron, eS)
    # R2: exon start into intron
    s2 <- pmax(eS  - width_intron, upE); e2 <- pmin(eS  + width_exon,   eE)
    # R3: exon End into exon
    s3 <- pmax(eE  - width_exon,   eS);  e3 <- pmin(eE  + width_intron, dnS)
    # R4: downstream exon start into intron
    s4 <- pmax(dnS - width_intron, eE);  e4 <- pmin(dnS + width_exon,   dnE)

    starts <- cbind(s1, s2, s3, s4)
    ends   <- cbind(e1, e2, e3, e4)

    GenomicRanges::GRanges(
      seqnames   = rep(events$chr, each = 4L),
      ranges     = IRanges::IRanges(start = as.vector(t(starts)),
                                    end   = as.vector(t(ends))),
      strand     = rep(events$strand, each = 4L),
      event_id   = rep(seq_len(n), each = 4L),
      region_idx = rep(seq_len(4L), times = n)
    )
  },

  region_labels = c(
    "Upstream exon | intron",
    "Intron | cassette exon",
    "Cassette exon | intron",
    "Intron | downstream exon"
  ),

  schematic     = "se",

  group_set     = c("Negative", "Positive", "Control")
)


#' @rdname event_schema
#' @noRd
event_schema_ri <- list(
  name          = "Retained Intron",

  required_cols = c(
    .event_required_base,
    "upstreamES", "upstreamEE",
    "downstreamES", "downstreamEE"
  ),

  n_regions     = 2L,

  region_width  = function(width_exon, width_intron) {
    as.integer(width_exon + width_intron + 1L)
  },

  build_regions = function(events, width_exon, width_intron) {
    n <- nrow(events)
    if (n == 0L) {
      return(GenomicRanges::GRanges())
    }

    upS <- events$upstreamES
    upE <- events$upstreamEE
    dnS <- events$downstreamES
    dnE <- events$downstreamEE

    # R1: upstream exon end into retained intron
    s1 <- pmax(upE - width_exon,   upS); e1 <- pmin(upE + width_intron, dnS)
    # R2: retained intron end into downstream exon
    s2 <- pmax(dnS - width_intron, upE); e2 <- pmin(dnS + width_exon,   dnE)

    starts <- cbind(s1, s2)
    ends   <- cbind(e1, e2)

    GenomicRanges::GRanges(
      seqnames   = rep(events$chr, each = 2L),
      ranges     = IRanges::IRanges(start = as.vector(t(starts)),
                                    end   = as.vector(t(ends))),
      strand     = rep(events$strand, each = 2L),
      event_id   = rep(seq_len(n), each = 2L),
      region_idx = rep(seq_len(2L), times = n)
    )
  },

  region_labels = c(
    "Upstream exon | intron",
    "Intron | downstream exon"
  ),

  schematic     = "ri",

  group_set     = c("Negative", "Positive", "Control")
)
