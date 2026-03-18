# global.R - Setup for RNAPeaks Shiny App

# Increase max file upload size to 500 MB
options(shiny.maxRequestSize = 500 * 1024^2)

# Set BioConductor repositories
options(repos = BiocManager::repositories())

# Core packages
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

# Source RNAPeaks functions directly from the package R/ directory
r_dir <- file.path(dirname(getwd()), "R")

source(file.path(r_dir, "BD_CheckBed.R"))
source(file.path(r_dir, "BD_FilterBed.R"))
source(file.path(r_dir, "BD_OrderPeak.R"))
source(file.path(r_dir, "BD_Prepare_Bed.R"))
source(file.path(r_dir, "BD_AnotBed.R"))

source(file.path(r_dir, "GN_GetGene.R"))
source(file.path(r_dir, "GN_Build_Gene_Structure_Helper_Functions.R"))
source(file.path(r_dir, "GN_Compute_Intron_Positions.R"))
source(file.path(r_dir, "GN_Make_Intron_Arrows.R"))
source(file.path(r_dir, "GN_Get_Plot_by_Gene.R"))
source(file.path(r_dir, "GN_draw_Gene_Plot.R"))
source(file.path(r_dir, "GN_PlotGene.R"))

source(file.path(r_dir, "RG_GetRegion.R"))
source(file.path(r_dir, "RG_Build_Region_Structure.R"))
source(file.path(r_dir, "RG_Make_Intron_Arrows.R"))
source(file.path(r_dir, "RG_Plot_Region_MultiGene.R"))
source(file.path(r_dir, "RG_PlotRegion.R"))

source(file.path(r_dir, "SM_HelperFunctions.R"))
source(file.path(r_dir, "SM_CreateSequenceMap.R"))
source(file.path(r_dir, "Splicing_Map.R"))




