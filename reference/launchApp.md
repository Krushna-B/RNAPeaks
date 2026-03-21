# Launch RNAPeaks Shiny Application

Starts an interactive web application for visualizing RNA-binding
protein peaks on gene structures. The app provides a user-friendly
interface for all RNAPeaks functionality without requiring R programming
knowledge.

## Usage

``` r
launchApp(launch.browser = TRUE, port = NULL, host = "127.0.0.1")
```

## Arguments

- launch.browser:

  Logical. If TRUE (default), the app opens in the default web browser.
  If FALSE, returns the app URL for manual access.

- port:

  The TCP port for the Shiny server. Default is NULL (random port).

- host:

  The IPv4 address for the server. Default is "127.0.0.1" (localhost).
  Use "0.0.0.0" to allow external connections.

## Value

Invisibly returns the Shiny app object. The function runs until the user
closes the browser window or stops the R session.

## Details

The Shiny app provides:

- File upload for BED files or use of sample data

- Single gene visualization (PlotGene)

- Genomic region visualization (PlotRegion)

- Splicing map analysis (createSplicingMap)

- Sequence motif mapping (createSequenceMap)

- Customizable styling options

- Export to PDF, or CSV

## Requirements

The shiny package must be installed. Install with:
`install.packages("shiny")`

## Examples

``` r
if (FALSE) { # \dontrun{
  # Launch the app in your default browser
  launchApp()

  # Launch on a specific port
  launchApp(port = 3838)

  # Allow external connections (for sharing on local network)
  launchApp(host = "0.0.0.0", port = 3838)
} # }
```
