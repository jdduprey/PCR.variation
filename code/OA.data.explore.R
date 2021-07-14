library(dplyr)
library(tidyverse)
library(vegan)
library(ggplot2)

ASV.table <- read.csv('../data/OA.ASV_table_all_together.csv')

# # filter to just miseq run 1 
# miseq1 <- ASV.table %>%
#    filter(Miseq_run %in% 1)


tech1 <- ASV.table %>%
  separate(col=sample, into=c('bio','tech'), sep='[.]', remove=FALSE) %>%
  filter(tech %in% 1) %>%
  rename(tech1reads=nReads)# %>%
  #select(-c(1:4))

reads.avg <- ASV.table %>%
  separate(col=sample, into=c('bio','tech'), sep='[.]', remove=FALSE) %>%
  select(bio, tech, Hash, nReads) %>%
  pivot_wider(names_from=tech, values_from=nReads, values_fill=0)
  
test <- left_join(reads.avg, tech1)
# distribution of log(nReads) for 
ggplot(miseq1, aes(x=log_reads)) +
  geom_histogram(color='black', fill='lightblue', bins=50)

# summary statistics? 
summary(miseq1)

# ======================================================
# NMDS/Bray Curtis - interested in comparing single techreplicate community 
# to 3 technical replicate community 
# just try to make NMDS with all techs for single miseq

# ======================================================
# make community matrix - extract columns with abundance information
# do we need to flip the table in order to do this?

#read tide paper again
#ask Ramon/Zach about code to analyse variation in data 
#zack and ramon say permanova + companion betadisper test

