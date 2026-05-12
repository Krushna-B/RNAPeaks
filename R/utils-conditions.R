#'Internal: classed-conditions helpers for RNAPeek's Errors
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


