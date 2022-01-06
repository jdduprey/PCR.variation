# anchovy, sardine data explore


library(tidyverse)
library(dplyr)
library(here)
library(ggplot2)

counts_plus_reads <- read_csv(here("special_issue_manuscript", "data", "calcofi_reads_counts_long.csv"))

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

counts_plus_reads %>% 
  unite(., Sample_ID, c("station_id","ext_rep","tech_rep"), sep=":",remove="F") %>% 
  left_join(mifish_tot_reads) %>% 
  left_join(larvae_tot_counts) %>% 
  mutate(., prop_reads= mifish_reads/tot_reads,
         prop_counts=larval_counts/tot_counts)

just_anchovies <- counts_plus_reads %>% 
  unite(., Sample_ID, c("station_id","ext_rep","tech_rep"), sep=":",remove="F") %>% 
  left_join(mifish_tot_reads) %>% 
  left_join(larvae_tot_counts) %>% 
  mutate(., prop_reads= mifish_reads/tot_reads,
         prop_counts=larval_counts/tot_counts) %>% 
  filter(., ID_mifish=="Engraulis mordax") %>% 
  filter(., larval_counts >0) %>% 
  dplyr::select(station_id,larval_counts,prop_counts) %>% distinct() %>% 
  arrange(desc(prop_counts))

