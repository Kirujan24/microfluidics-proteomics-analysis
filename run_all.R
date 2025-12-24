# run_all.R

# Master script for Paper-1:
# Microfluidic enrichment and multi-omics analysis of urinary extracellular vesicles

# This script sources the individual analysis modules in the order used in the manuscript.

# Clear workspace for clean reproducibility
rm(list = ls())

# Set a consistent seed (if any random steps are added later)
set.seed(123)

# Run NTA analysis (Figure 3)
source("scripts/nta_analysis.R")

# Run proteomics analysis (Figures 5–7)
source("scripts/proteomics_analysis.R")



