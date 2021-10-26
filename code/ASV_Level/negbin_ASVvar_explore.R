# Exploring Ryan's Negbin ASV variation code 
# 10/25/21

library(tidyverse)
library(vegan)
library(stringi)
library(fitdistrplus)


# Read in data and clean it a bit (not removing low read outliers)
raw.ASVs <- read.csv('/Users/helencasendino/Desktop/PCR.variation/input/raw_ASVs.csv')

reads_long <- raw.ASVs %>%
  filter(str_detect(sample, "Ostrich", TRUE)) %>% # remove ostrich samples
  filter(str_detect(sample, "Kangaroo", TRUE)) %>% # remove kangaroo samples
  filter(str_detect(sample, "K+", TRUE)) %>%  # remove negative control (??) samples
  filter(str_detect(sample, "k+", TRUE)) %>% 
  separate(col = sample, into = c('bio', 'tech'), sep = '[.]', remove = FALSE) 

reads_long <- reads_long %>% 
  mutate(bio = gsub("_", "", reads_long$bio)) %>% 
  mutate(sample = gsub("_", "", reads_long$sample))


# first couple bits are creating two functions: 
#first to get the variance of a negbin distribution, given the parameters mu and phi 
# Second to use fitdist() to find the most likely parameters mu and phi (which R calls "size" in this case), given some data vector. 
# Then the tidyr code takes a dataset (here, called dat.mifish, but you can use your dataset), groups by species (ID_mifish) and biological sample (station_id),
# and then fits a negative binomial distribution  to each vector of 3 observations for each species at each site. 
# Then it plots the result.


varNB <- function(mu, phi){
  mu + ((mu^2)/phi)
}
getNB <- function(x){
  require(fitdistrplus)
  
  return(
    (fitdist(unlist(x), "nbinom"))$estimate
  )
}

a <- reads_long %>% 
  filter(bio %in% unique(reads_long$bio)) %>% 
  group_by(bio, Hash) %>% 
  filter(sum(nReads)>0) %>%
  filter(n() > 1) %>%
  dplyr::select(nReads) %>% 
  nest() %>%
  mutate(mlmodels = map(data,getNB))

a$mlmodels %>% 
  as.data.frame() %>% 
  t() %>% as_tibble() %>% 
  rename(phi = size) %>% 
  mutate(estVariance = varNB(mu, phi)) %>% 
  pivot_longer(-mu) %>% 
  ggplot(aes(x = log(mu), y = log(value))) +
  geom_point() +
  facet_wrap(~name)