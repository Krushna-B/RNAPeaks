# global.R - Setup for RNAPeaks Shiny App

# Increase max file upload size to 500 MB
options(shiny.maxRequestSize = 500 * 1024^2)

# Set BioConductor repositories
options(repos = BiocManager::repositories())

# Load required packages
library(shiny)
library(bslib)
library(bsicons)
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
library(progress)
library(magrittr)

# Get the app directory (works on shinyapps.io)
app_dir <- getwd()
message("App directory: ", app_dir)
message("Files in app directory: ", paste(list.files(app_dir), collapse = ", "))

# Always source bundled R files and data for shinyapps.io compatibility
message("Sourcing bundled R files...")

r_dir <- file.path(app_dir, "R")
data_dir <- file.path(app_dir, "data")

message("R directory: ", r_dir)
message("Data directory: ", data_dir)
message("R files: ", paste(list.files(r_dir), collapse = ", "))
message("Data files: ", paste(list.files(data_dir), collapse = ", "))

# BED file operations
source(file.path(r_dir, "BD_CheckBed.R"))
source(file.path(r_dir, "BD_FilterBed.R"))
source(file.path(r_dir, "BD_OrderPeak.R"))
source(file.path(r_dir, "BD_Prepare_Bed.R"))
source(file.path(r_dir, "BD_AnotBed.R"))

# Gene plotting
source(file.path(r_dir, "GN_GetGene.R"))
source(file.path(r_dir, "GN_Build_Gene_Structure_Helper_Functions.R"))
source(file.path(r_dir, "GN_Compute_Intron_Positions.R"))
source(file.path(r_dir, "GN_Make_Intron_Arrows.R"))
source(file.path(r_dir, "GN_Get_Plot_by_Gene.R"))
source(file.path(r_dir, "GN_draw_Gene_Plot.R"))
source(file.path(r_dir, "GN_PlotGene.R"))

# Region plotting
source(file.path(r_dir, "RG_GetRegion.R"))
source(file.path(r_dir, "RG_Build_Region_Structure.R"))
source(file.path(r_dir, "RG_Make_Intron_Arrows.R"))
source(file.path(r_dir, "RG_Plot_Region_MultiGene.R"))
source(file.path(r_dir, "RG_PlotRegion.R"))

# Splicing and Sequence maps
source(file.path(r_dir, "SM_HelperFunctions.R"))
source(file.path(r_dir, "SM_CreateSequenceMap.R"))
source(file.path(r_dir, "Splicing_Map.R"))

# Data is loaded directly in app.R for simplicity
message("Data loading handled in app.R")

message("RNAPeaks Shiny app initialized successfully")
