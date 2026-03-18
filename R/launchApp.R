#' Launch RNAPeaks Shiny Application
#'
#' Starts an interactive web application for visualizing RNA-binding protein
#' peaks on gene structures. The app provides a user-friendly interface for
#' all RNAPeaks functionality without requiring R programming knowledge.
#'
#' @param launch.browser Logical. If TRUE (default), the app opens in the
#'   default web browser. If FALSE, returns the app URL for manual access.
#' @param port The TCP port for the Shiny server. Default is NULL (random port).
#' @param host The IPv4 address for the server. Default is "127.0.0.1" (localhost).
#'   Use "0.0.0.0" to allow external connections.
#'
#' @return Invisibly returns the Shiny app object. The function runs until
#'   the user closes the browser window or stops the R session.
#'
#' @details
#' The Shiny app provides:
#' \itemize{
#'   \item File upload for BED files or use of sample data
#'   \item Single gene visualization (PlotGene)
#'   \item Genomic region visualization (PlotRegion)
#'   \item Splicing map analysis (createSplicingMap)
#'   \item Sequence motif mapping (createSequenceMap)
#'   \item Customizable styling options
#'   \item Export to PDF, or CSV
#' }
#'
#' @section Requirements:
#' The shiny package must be installed. Install with:
#' \code{install.packages("shiny")}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Launch the app in your default browser
#'   launchApp()
#'
#'   # Launch on a specific port
#'   launchApp(port = 3838)
#'
#'   # Allow external connections (for sharing on local network)
#'   launchApp(host = "0.0.0.0", port = 3838)
#' }
launchApp <- function(launch.browser = TRUE, port = NULL, host = "127.0.0.1") {

  # Check if shiny is installed
if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the app. ",
         "Install it with: install.packages('shiny')",
         call. = FALSE)
  }

  # Set max file upload size to 500 MB
  options(shiny.maxRequestSize = 500 * 1024^2)

  # Launch the app from GitHub
  message("Starting RNAPeaks Shiny app...")
  message("Close the browser window or press Ctrl+C to stop.")

  shiny::runGitHub(
    repo      = "RNAPeaks",
    username  = "Krushna-B",
    subdir    = "deployment",
    launch.browser = launch.browser,
    port      = port,
    host      = host
  )
}
