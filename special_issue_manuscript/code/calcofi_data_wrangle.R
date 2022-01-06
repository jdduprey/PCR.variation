# ==========================================================
# getting calcofi data into format to plot amp eff. and microscopy counts
# Joe Duprey, Helen Casendino, Kai  
# Last edited 01/06/2021
# ==========================================================

library(tidyverse)
library(dplyr)
library(here)
library(ggplot2)

# load data
micro_counts <- readRDS(here("special_issue_manuscript", "data", "microscopy_tech_nReads.RDS"))
mifish_reads <- readRDS(here("special_issue_manuscript", "data", "mifish_tech_nReads.RDS"))

counts_plus_reads <- left_join(mifish_reads, micro_counts)

print(unique(counts_plus_reads$tech_rep))

# lets get the data into wide form 
counts_plus_reads <- counts_plus_reads %>%
  mutate(bio_rep = paste(station_id, ext_rep, sep = "_")) %>%
  select(-c(Sample, ID_sebastes, Unique_ID)) %>%
  distinct() 

write_csv(counts_plus_reads, here("special_issue_manuscript", "data", "calcofi_reads_counts_long.csv"))

counts_plus_reads_wide <- counts_plus_reads %>%
  pivot_wider(names_from = tech_rep, values_from = mifish_reads) 
  
write_csv(counts_plus_reads_wide, here("special_issue_manuscript", "data", "calcofi_reads_counts_wide.csv"))

# ggplot(counts_plus_reads, aes(x=mifish_reads)) + 
#   geom_histogram()
# 
# ggplot(counts_plus_reads, aes(x=larval_counts)) + 
#   geom_histogram()

