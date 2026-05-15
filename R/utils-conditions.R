#'Internal: classed-conditions helpers for RNAPeak's Errors
#' @keywords internal

#1. Any BED Related Issue
abort_invalid_bed <- function(message, ..., call = parent.frame()){
  cli::cli_abort(
    message,
    ...,
    class = c("rnapeaks_error_invalid_bed","rnapeaks_error"),
    call = call,
    .envir = call
  )
}

#2. Any Argument Related Issue
abort_invalid_arg <- function(message, ..., call = parent.frame()){
  cli::cli_abort(
    message,
    ...,
    class = c("rnapeaks_error_invalid_arg","rnapeaks_error"),
    call = call,
    .envir = call
  )
}

#3. Any lookup not found
abort_not_found <- function(message, ..., call = parent.frame()){
  cli::cli_abort(
    message,
    ...,
    class = c("rnapeaks_error_not_found","rnapeaks_error"),
    call = call,
    .envir = call
  )
}

#4. Invalid GTF
abort_invalid_gtf <- function(message, ..., call = parent.frame()){
  cli::cli_abort(
    message,
    ...,
    class = c("rnapeaks_error_invalid_gtf","rnapeaks_error"),
    call = call,
    .envir = call
  )
}

#5. Cannot determine GTF Species
abort_species_unknown <- function(..., call = parent.frame()){
  cli::cli_abort(
    c("Cannot infer species from {.arg gtf}.",
      "i" = "Expected gene_id starting with {.val ENSG} (human) or {.val ENSMUSG} (mouse)."),
    ...,
    class = c("rnapeaks_error_species_unknown",
              "rnapeaks_error_invalid_gtf",
              "rnapeaks_error"),
    call = call,
    .envir = call
  )
}


# More Internal validation used by peaks_plot_style()
.abort_check <- function(message, call) {
  cli::cli_abort(
    message,
    class = c("rnapeaks_error_invalid_arg", "rnapeaks_error"),
    call  = call,
    .envir = parent.frame()
  )
}

# Single non-NA numeric in [min, max].
check_scalar_number <- function(x, arg, min = -Inf, max = Inf,
                                call = parent.frame()) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    .abort_check(c(
      "{.arg {arg}} must be a single non-NA number.",
      "x" = "Got {.cls {class(x)[1]}} of length {length(x)}."
    ), call = call)
  }
  if (x < min || x > max) {
    .abort_check(c(
      "{.arg {arg}} must lie in [{min}, {max}].",
      "x" = "Got {.val {x}}."
    ), call = call)
  }
  invisible(x)
}

# Single non-NA whole number in [min, max].
check_scalar_int <- function(x, arg, min = -Inf, max = Inf,
                             call = parent.frame()) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != trunc(x)) {
    .abort_check(c(
      "{.arg {arg}} must be a single non-NA whole number.",
      "x" = "Got {.cls {class(x)[1]}} {.val {x}}."
    ), call = call)
  }
  if (x < min || x > max) {
    .abort_check(c(
      "{.arg {arg}} must lie in [{min}, {max}].",
      "x" = "Got {.val {x}}."
    ), call = call)
  }
  invisible(x)
}

# Single non-NA number in [0, 1].
check_unit_interval <- function(x, arg, call = parent.frame()) {
  check_scalar_number(x, arg, min = 0, max = 1, call = call)
}

# NULL or a single non-NA number in [min, max].
check_scalar_number_or_null <- function(x, arg, min = -Inf, max = Inf,
                                        call = parent.frame()) {
  if (is.null(x)) return(invisible(x))
  check_scalar_number(x, arg, min = min, max = max, call = call)
}

# Single non-NA TRUE/FALSE.
check_flag <- function(x, arg, call = parent.frame()) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .abort_check(c(
      "{.arg {arg}} must be a single {.val TRUE} or {.val FALSE}.",
      "x" = "Got {.cls {class(x)[1]}} of length {length(x)}."
    ), call = call)
  }
  invisible(x)
}

# Single non-NA character. If `choices` is supplied, value must be in it.
check_string <- function(x, arg, choices = NULL, call = parent.frame()) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) {
    .abort_check(c(
      "{.arg {arg}} must be a single non-NA string.",
      "x" = "Got {.cls {class(x)[1]}} of length {length(x)}."
    ), call = call)
  }
  if (!is.null(choices) && !x %in% choices) {
    .abort_check(c(
      "{.arg {arg}} must be one of {.or {.val {choices}}}.",
      "x" = "Got {.val {x}}."
    ), call = call)
  }
  invisible(x)
}

# Single color: a named color, hex string, or (when allow_na) NA.
check_color <- function(x, arg, allow_na = TRUE, call = parent.frame()) {
  if (length(x) != 1L) {
    .abort_check(c(
      "{.arg {arg}} must be a single color.",
      "x" = "Got {length(x)} values."
    ), call = call)
  }
  if (is.na(x)) {
    if (allow_na) return(invisible(x))
    .abort_check("{.arg {arg}} cannot be {.val NA}.", call = call)
  }
  ok <- tryCatch({grDevices::col2rgb(x); TRUE}, error = function(e) FALSE)
  if (!ok) {
    .abort_check(c(
      "{.arg {arg}} is not a recognised color.",
      "x" = "Got {.val {x}}.",
      "i" = "Use a named color (e.g. {.val red}) or a hex string (e.g. {.val \"#FF0000\"})."
    ), call = call)
  }
  invisible(x)
}

# NULL or a length-2 numeric c(min, max) with min < max and no NAs.
check_range_or_null <- function(x, arg, call = parent.frame()) {
  if (is.null(x)) return(invisible(x))
  if (!is.numeric(x) || length(x) != 2L || anyNA(x)) {
    .abort_check(c(
      "{.arg {arg}} must be {.val NULL} or a length-2 numeric.",
      "x" = "Got {.cls {class(x)[1]}} of length {length(x)}."
    ), call = call)
  }
  if (x[1] >= x[2]) {
    .abort_check(c(
      "{.arg {arg}} must satisfy {.code min < max}.",
      "x" = "Got {.code c({x[1]}, {x[2]})}."
    ), call = call)
  }
  invisible(x)
}

# Trim whitespace on a single string input. NULL passes through.
normalize_str <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1L && !is.na(x)) return(trimws(x))
  x
}

# Canonicalise a chromosome name: trim, uppercase, strip a leading "CHR".
# So "chr19", "Chr19", " CHR19" all become "19"; "chrX" -> "X". NULL passes
# through. The same canonical form is used by the BED and GTF normalizers,
# so values from any source compare directly.
normalize_chr <- function(x) {
  if (is.null(x)) return(NULL)
  sub("^CHR", "", toupper(trimws(as.character(x))))
}

# Accept numeric, or character with commas / surrounding whitespace
# (e.g. "43,190,342"), and return numeric. NULL passes through. Aborts with
# a clear message if a character element cannot be parsed.
normalize_coord <- function(x, arg = "x") {
  if (is.null(x)) return(NULL)
  if (is.character(x)) {
    cleaned <- trimws(gsub(",", "", x, fixed = TRUE))
    cleaned[cleaned == ""] <- NA_character_
    coerced <- suppressWarnings(as.numeric(cleaned))
    bad <- is.na(coerced) & !is.na(x)
    if (any(bad)) {
      abort_invalid_arg(c(
        "{.arg {arg}} must be numeric or a numeric string (commas allowed).",
        "x" = "Could not parse: {.val {x[bad]}}."
      ))
    }
    return(coerced)
  }
  x
}

# Non-empty data frame containing every column in `required`.
check_data_frame <- function(df, arg, required = character(),
                             call = parent.frame()) {
  if (!is.data.frame(df)) {
    .abort_check(c(
      "{.arg {arg}} must be a data frame.",
      "x" = "Got {.cls {class(df)[1]}}."
    ), call = call)
  }
  if (nrow(df) == 0L) {
    .abort_check("{.arg {arg}} has no rows.", call = call)
  }
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols)) {
    .abort_check(c(
      "{.arg {arg}} is missing required column{?s}: {.field {missing_cols}}.",
      "i" = "Required: {.field {required}}."
    ), call = call)
  }
  invisible(df)
}


