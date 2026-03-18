# global.R - Setup for RNAPeaks Shiny App

# Increase max file upload size to 500 MB
options(shiny.maxRequestSize = 500 * 1024^2)

# Set BioConductor repositories
options(repos = BiocManager::repositories())

# Install RNAPeaks from GitHub if not available
if (!requireNamespace("RNAPeaks", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("Krushna-B/RNAPeaks")
}

# Load RNAPeaks and its dependencies
library(RNAPeaks)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)
library(BiocGenerics)
library(GenomeInfoDb)
library(Biostrings)
library(scales)
library(grid)
library(slider)
library(magrittr)




