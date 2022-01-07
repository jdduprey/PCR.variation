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

# calculating proportion reads 
counts_plus_reads %>% 
  unite(., Sample_ID, c("station_id","ext_rep","tech_rep"), sep=":" ) %>% 
  dplyr::select(-larval_counts) %>% 
  filter(., !is.na(mifish_reads)) %>% 
  group_by(Sample_ID) %>% 
  dplyr::summarise(tot_reads=sum(mifish_reads)) -> mifish_tot_reads

counts_plus_reads %>% 
  unite(., Sample_ID, c("station_id","ext_rep","tech_rep"), sep=":") %>% 
  dplyr::select(-mifish_reads) %>% 
  filter(., !is.na(larval_counts)) %>% 
  group_by(Sample_ID) %>% 
  dplyr::summarise(tot_counts=sum(larval_counts)) -> larvae_tot_counts

cal_props <- counts_plus_reads %>% 
  unite(., Sample_ID, c("station_id","ext_rep","tech_rep"), sep=":",remove="F") %>% 
  left_join(mifish_tot_reads) %>% 
  left_join(larvae_tot_counts) %>% 
  mutate(., prop_reads= mifish_reads/tot_reads,
         prop_counts=larval_counts/tot_counts)

calcofi_amps <- read.csv(here("special_issue_manuscript", "data", "alphas_oceanic.csv"))
calcofi_amps <- calcofi_amps %>% rename("ID_mifish" = "Species")
complete_cal_props_long <- left_join(cal_props, calcofi_amps, by = "ID_mifish" )

write_csv(complete_cal_props_long, here("special_issue_manuscript", "data", "complete_cal_props_long.csv"))


# cal_props %>%
#   dplyr::select(-prop_counts) %>% 
#   drop_na(prop_reads) %>% 
#   group_by(Sample_ID) %>%
#   dplyr::summarise(prop_check = sum(prop_reads)) -> sanity_check

# cal_props_csv <- cal_props %>%
#   dplyr::select(-prop_reads) %>%
#   drop_na(prop_reads)

# cal_prop_reads_wide <- cal_props %>%
#   drop_na(prop_reads) %>% 
#   dplyr::select(-c(Sample_ID, prop_counts, tot_reads, tot_counts)) %>% 
#   pivot_wider(names_from = tech_rep, values_from = prop_reads)
