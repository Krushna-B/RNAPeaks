# Columns build_intron_arrows() always returns
ARROW_COLS <- c("x", "xend", "y", "yend")



build_intron_arrows <- function(introns,
                                xlim,
                                arrows_per_view = 12,
                                min_intron_frac = 0.02,
                                min_intron_bp   = 0,
                                max_per_intron  = 6L,
                                shaft_frac      = 0.4) {
  #Validate inputs
  require_intron_df(introns)
  require_xlim(xlim)
  require_positive_scalar(arrows_per_view, "arrows_per_view")
  require_unit_fraction(min_intron_frac,   "min_intron_frac")
  require_positive_scalar(min_intron_bp,   "min_intron_bp")
  require_positive_scalar(max_per_intron,  "max_per_intron")
  require_unit_fraction(shaft_frac,        "shaft_frac")

  #Empty input or zero-width view -> empty result
  if (!nrow(introns)) return(empty_arrows())
  view_span <- xlim[2] - xlim[1]
  if (view_span <= 0) return(empty_arrows())

  #Clip each intron to the visible window
  vis_start <- pmax(introns$start, xlim[1])
  vis_end   <- pmin(introns$end,   xlim[2])
  vis_len   <- pmax(0, vis_end - vis_start)
  vis_frac  <- vis_len / view_span

  #Arrows per intron: density proportional to visible fraction
  #Length smaller than min are ignored
  n <- ifelse(vis_frac < min_intron_frac | vis_len < min_intron_bp,
              0L,
              pmin(max_per_intron,
                   pmax(1L, as.integer(round(vis_frac * arrows_per_view)))))

  if (!any(n > 0L)) return(empty_arrows())

  #y-line per intron (supports several upstream column shapes)
  y_line <- resolve_intron_y(introns)

  #Direction per intron from strand
  dir <- ifelse(introns$strand == "-", -1, 1)

  #Emit segments
  rows <- lapply(seq_len(nrow(introns)), function(i) {
    if (n[i] == 0L) return(NULL)
    spacing <- vis_len[i] / (n[i] + 1L)
    centers <- vis_start[i] + spacing * seq_len(n[i])
    shaft   <- shaft_frac * spacing
    data.frame(
      x    = centers - dir[i] * shaft / 2,
      xend = centers + dir[i] * shaft / 2,
      y    = y_line[i],
      yend = y_line[i],
      row.names = NULL
    )
  })

  out <- do.call(rbind, rows)
  if (is.null(out)) empty_arrows() else out
}




#Empty result with the canonical columns
empty_arrows <- function() {
  data.frame(x = numeric(0), xend = numeric(0), y = numeric(0), yend = numeric(0))
}

#Pick a y-coordinate per intron from whichever column shape is present
resolve_intron_y <- function(introns) {
  if ("y" %in% names(introns)) return(introns$y)
  if ("y_mid" %in% names(introns)) return(introns$y_mid)
  if (all(c("y_start", "y_end") %in% names(introns))) {
    return((introns$y_start + introns$y_end) / 2)
  }
  rep(0, nrow(introns))
}




#Validate introns data frame has the columns we need
require_intron_df <- function(introns) {
  if (!is.data.frame(introns)) {
    abort_invalid_arg(c(
      "{.arg introns} must be a data frame.",
      "x" = "You supplied {.cls {class(introns)[1]}}."
    ))
  }
  required <- c("start", "end")
  missing  <- setdiff(required, colnames(introns))
  if (length(missing)) {
    abort_invalid_arg(c(
      "{.arg introns} is missing required column{?s}: {.field {missing}}.",
      "i" = "Required: {.field {required}}."
    ))
  }
}

#Validate xlim is a length-2 numeric in ascending order
require_xlim <- function(xlim) {
  if (!is.numeric(xlim) || length(xlim) != 2L || anyNA(xlim)) {
    abort_invalid_arg(c(
      "{.arg xlim} must be a length-2 numeric vector with no NAs.",
      "x" = "Got {.cls {class(xlim)[1]}} of length {length(xlim)}."
    ))
  }
  if (xlim[2] < xlim[1]) {
    abort_invalid_arg(c(
      "{.arg xlim} must be in ascending order.",
      "x" = "Got {.code xlim[1]} = {.val {xlim[1]}}, {.code xlim[2]} = {.val {xlim[2]}}."
    ))
  }
}

#Validate a value is in [0, 1]
require_unit_fraction <- function(value, arg) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value < 0 || value > 1) {
    abort_invalid_arg(c(
      "{.arg {arg}} must be a single number in {.val {c(0, 1)}}.",
      "x" = "Got {.cls {class(value)[1]}} = {.val {value}}."
    ))
  }
}
