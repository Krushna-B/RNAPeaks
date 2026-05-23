#' Validate and normalize BED frame input
#'
#' Checks BED data for required columns and valid coordinates, then combines
#' the input into a single data frame with standardized column names.
#'
#' @param bed A data frame or list of data frames with at least 6 columns.
#'   Columns are mapped by position: 1=chr, 2=start, 3=end, 6=strand.
#'
#' @param split_col Optional positive integer giving the column index used to
#'   group rows. Every bed must have at least `split_col` columns, and the
#'   index cannot point at a canonical position (1=chr, 2=start, 3=end,
#'   6=strand).
#'
#' @return A single combined, validated bed data frame
#' @noRd
#' @family bed
check_bed <- function(bed, split_col = NULL) {
  #Check if input is a list c(Bed1,...Bedn)
  beds <- if (is.data.frame(bed)) list(bed) else bed

  #Check input length is not 0 and assert all beds are dataframes
  if (!is.list(beds) || length(beds) == 0L ||
      !all(vapply(beds, is.data.frame, logical(1)))){
    abort_invalid_bed(c(
      "{.arg bed} must be a non-empty list of bed files.",
      "x" = "You supplied {.cls {class(bed)[1]}}."
    ))
  }

  #Check if each bed file has 6 cols
  ncols <- vapply(beds, ncol, integer(1))
  if (any(ncols < 6L)){
    bad_bed <- which(ncols < 6L)
    abort_invalid_bed(c(
      "{.arg bed} must have at least 6 columns.",
      "x" = "Bed {bad_bed[1]} has {ncols[bad_bed[1]]} column{?s}.",
      "i" = "Expected order: chr, start, end, name, score, strand."
    ))
  }

  #Check if all beds have same number of cols
  if (length(unique(ncols)) > 1L) {
    abort_invalid_bed(c(
      "All beds in {.arg bed} must have the same number of columns.",
      "x" = "Got column counts: {.val {ncols}}.",
      "i" = "Cannot combine beds with mismatched shapes."
    ))
  }

  #Rename columns 1,2,3,6 to (chr,start,end,strand)
  CANONICAL <- c("chr", "start", "end", "strand")
  beds <- lapply(beds, function(b){
    colnames(b)[c(1L, 2L, 3L, 6L)] <- CANONICAL
    return(b)
  })

  #Check all colnames are the same for each bed
  ref_names  <- colnames(beds[[1]])
  same_names <- vapply(beds, function(b) identical(colnames(b), ref_names),
                       logical(1))
  if (!all(same_names)) {
    bad <- which(!same_names)[1]
    abort_invalid_bed(c(
      "All beds in {.arg bed} must share identical column names.",
      "x" = "Bed {bad} has different column names from bed 1."
    ))
  }

  if (!is.null(split_col)) {
    if (!is.numeric(split_col) || length(split_col) != 1L || is.na(split_col) ||
        split_col != trunc(split_col) || split_col < 1L) {
      abort_invalid_arg("{.arg split_col} must be a single positive integer or {.code NULL}.")
    }
    split_col <- as.integer(split_col)
    if (split_col %in% c(1L, 2L, 3L, 6L)) {
      abort_invalid_arg(c(
        "{.arg split_col} cannot point at a canonical column.",
        "x" = "You supplied {.val {split_col}}.",
        "i" = "Canonical positions: 1=chr, 2=start, 3=end, 6=strand."
      ))
    }
    if (any(ncols < split_col)) {
      bad <- which(ncols < split_col)
      abort_invalid_arg(c(
        "{.arg split_col} = {split_col} exceeds available columns.",
        "x" = "Bed {bad[1]} has {ncols[bad[1]]} column{?s}."
      ))
    }
  }

  #Adding target column
  if (!is.null(split_col)) {
    beds <- lapply(beds, function(b){ b$target <- b[[split_col]]; return(b) })
  } else if (length(beds) >= 1L){
    nm <- names(beds)
    if (is.null(nm)){ nm <- rep("", length(beds)) }
      nm[nm == ""] <- paste0("bed", which(nm == ""))
      beds <- Map(function(b, lab) { b$target <- lab; b }, beds, nm)
  }

  #Combine into one bed and normalize
  combined <- do.call(rbind, beds)
  combined$chr <- normalize_chr(combined$chr) #All upper case, '1', 'X', 'Y'
  combined$start <- suppressWarnings(as.numeric(combined$start))
  combined$end   <- suppressWarnings(as.numeric(combined$end))    #may produce NA which will be ignored
  combined$strand <- as.character(combined$strand)

  #Remove any rows with Start and End coordinates are NA (Warning)
  bad_rows <- is.na(combined$start) | is.na(combined$end)
  if (any(bad_rows)) {
    cli::cli_warn(c(
      "Dropped {sum(bad_rows)} row{?s} with non-numeric {.field start} or {.field end}.",
      "i" = "First dropped row: {which(bad_rows)[1]}."
    ))
    combined <- combined[!bad_rows, , drop = FALSE]
    if (nrow(combined) == 0L) {
      abort_invalid_bed("All rows were dropped; {.arg bed} has no valid coordinates.")
    }
  }

  #Check all start coordinates are positive
  if (any(combined$start < 0)) {
    bad <- which(combined$start < 0)
    abort_invalid_bed(c(
      "{.arg bed} has negative {.field start} coordinates.",
      "x" = "{length(bad)} bad row{?s} (first: row {bad[1]}, value {combined$start[bad[1]]})."
    ))
  }

  #Check all end coordinates are positive
  if (any(combined$end < 0)) {
    bad <- which(combined$end < 0)
    abort_invalid_bed(c(
      "{.arg bed} has negative {.field end} coordinates.",
      "x" = "{length(bad)} bad row{?s} (first: row {bad[1]}, value {combined$end[bad[1]]})."
    ))
  }

  #Check for all coordinates start < end
  if (any(combined$end < combined$start)) {
    bad <- which(combined$end < combined$start)
    abort_invalid_bed(c(
      "{.arg bed} has rows where {.field end} < {.field start}.",
      "x" = "{length(bad)} bad row{?s} (first: row {bad[1]})."
    ))
  }

  #Check that strand is either (+) or (-)
  if (!all(combined$strand %in% c("+", "-"))) {
    bad <- which(!combined$strand %in% c("+", "-"))
    abort_invalid_bed(c(
      "{.arg bed} {.field strand} must be {.val +} or {.val -}.",
      "x" = "{length(bad)} bad row{?s} (first: row {bad[1]}, value {.val {combined$strand[bad[1]]}})."
    ))
  }

  return(combined)
}
