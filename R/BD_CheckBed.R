#' Validate and normalize BED frame input
#'
#' Checks BED data for required columns and valid coordinates, then combines
#' the input into a single data frame with standardized column names.
#'
#' @param bed A data frame or list of data frames with at least 6 columns.
#'   Columns are mapped by position: 1=chr, 2=start, 3=end, 6=strand.
#'
#' @param split_col Optional. Name of column used to group
#'   rows. Must exist in every bed and cannot be one of the
#'   canonical columns (`chr`, `start`, `end`, `strand`).
#'
#' @return A single combined, validated bed data frame
#'
#' @export
#' @family bed
checkBed <- function(bed, split_col = NULL) {
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

  #Check if split_col is provided and if it is valid
  if (!is.null(split_col)) {
    if (!is.character(split_col) || length(split_col) != 1L || is.na(split_col)) {
      abort_invalid_arg("{.arg split_col} must be a single column name or {.code NULL}.")
    }
    if (split_col %in% CANONICAL) {
      abort_invalid_arg(c(
        "{.arg split_col} cannot be one of the canonical columns.",
        "x" = "You supplied {.val {split_col}}.",
        "i" = "Canonical (positional): {.field {CANONICAL}}."
      ))
    }
    has_col <- vapply(beds, function(b) split_col %in% colnames(b), logical(1))
    if (!all(has_col)) {
      bad <- which(!has_col)
      abort_invalid_arg(c(
        "Column {.field {split_col}} not found in bed {bad[1]}.",
        "i" = "Available: {.field {colnames(beds[[bad[1]]])}}."
      ))
    }
  }

  #Adding target column
  if (!is.null(split_col)) {
    beds <- lapply(beds, function(b){ b$target <- b[[split_col]]; return(b) })
  } else if (length(beds) >= 2L){
    nm <- names(beds)
    if (is.null(nm)){ nm <- rep("", length(beds)) }
      nm[nm == ""] <- paste0("bed", which(nm == ""))
      beds <- Map(function(b, lab) { b$target <- lab; b }, beds, nm)
  }

  #Combine into one bed and normalize
  combined <- do.call(rbind, beds)
  combined$chr <- sub("^CHR", "", toupper(as.character(combined$chr))) #All upper case, '1', 'X', 'Y'
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

  cli::cli_alert_success("BED files validated successfully")
  return(combined)
}
