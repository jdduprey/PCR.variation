# ==========================================================
# getting calcofi data into format to plot amp eff. and microscopy counts
# Joe Duprey, Helen Casendino, Kai  
# Last edited 01/06/2021
# ==========================================================

library(tidyverse)
library(dplyr)
library(here)

# load data
micro_counts <- readRDS(here("special_issue_manuscript", "data", "microscopy_tech_nReads.RDS"))
ASV_reads <- readRDS(here("special_issue_manuscript", "data", "mifish_tech_nReads.RDS"))

counts_reads <- left_join(ASV_reads, micro_counts)
