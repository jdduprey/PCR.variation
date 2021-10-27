#####=====ASV variation using Negative Binomial Parameters=======
# Uses Ryan's code to fit PCR reads for each ASV (usually 3) to a negbin distribution,
# and estimating mu and phi
# 10/25/21
#####===================================================

###==Dependencies====
library(tidyverse)
library(stringi)
library(fitdistrplus)
library(ggpubr)
library(here)
#===============

#####=======READ IN DATA====
raw.ASVs <- read.csv(here("input/raw_ASVs.csv"))

reads_long <- raw.ASVs %>%
  filter(str_detect(sample, "Ostrich", TRUE)) %>% # remove ostrich samples
  filter(str_detect(sample, "Kangaroo", TRUE)) %>% # remove kangaroo samples
  filter(str_detect(sample, "K+", TRUE)) %>%  # remove negative control (??) samples
  filter(str_detect(sample, "k+", TRUE)) %>% 
  separate(col = sample, into = c('bio', 'tech'), sep = '[.]', remove = FALSE) 

# Standardizing 
reads_long <- reads_long %>% 
  mutate(bio = gsub("_", "", reads_long$bio)) %>% 
  mutate(sample = gsub("_", "", reads_long$sample))


###======Ryan's code, updated for Hood Canal Data===========

# "first couple bits are creating two functions: 
#first to get the variance of a negbin distribution, given the parameters mu and phi 
# Second to use fitdist() to find the most likely parameters mu and phi (which R calls "size" in this case), given some data vector. 
# Then the tidyr code takes a dataset (here, called dat.mifish, but you can use your dataset), groups by species (ID_mifish) and biological sample (station_id),
# and then fits a negative binomial distribution  to each vector of 3 observations for each species at each site. 
# Then it plots the result."

varNB <- function(mu, phi){ # calculates variance
  mu + ((mu^2)/phi)
} # calculates variance
getNB <- function(x){
  require(fitdistrplus)
  
  return(
    (fitdist(unlist(x), "nbinom"))$estimate
  )
} # estimates mu and phi

# Below: estimates parameters for each Hash (2-3 PCR reads/Hash)
# nested_df needs cols: Miseq_run, bio (bottle, ie PO2017A), Hash, and nReads
# spits out list of: 
      # 1) df of mu, estimated variance and phi, 
#   and 2) plot of mu vs. estimated variance and mu vs. phi

estimates_pars_byHash <- function(df){
  
  list_NBpars <- df %>% 
    filter(bio %in% unique(df$bio)) %>% 
    group_by(bio, Hash) %>% 
    filter(sum(nReads)>0) %>%
    filter(n() > 1) %>%
    dplyr::select(nReads) %>% 
    nest() %>%
    mutate(mlmodels = map(data,getNB))
  
  df_NBpars <- list_NBpars$mlmodels %>%  
    as.data.frame() %>% 
    t() %>% as_tibble() %>% 
    rename(phi = size) %>% 
    mutate(estVariance = varNB(mu, phi)) %>% 
    pivot_longer(-mu)

  NBpars_plot <- df_NBpars %>% 
    ggplot(aes(x = log(mu), y = log(value))) +
    geom_point() +
    facet_wrap(~name)
  
  ### Just exploratory (adds list of original 3 reads to df_NBpars) - Helen 
  # df_NBpars_reads <- df_NBpars %>% 
  #  pivot_wider(names_from = name) %>% 
  #  mutate(Reads = list_NBpars$data)
  
  return(list(df_NBpars, df_NBpars_reads, NBpars_plot))
}

NB_ouput <- estimates_pars_byHash(reads_long)



# Helen: If the curve of mu vs. variance is dependent on phi, we also want to 
# keep track of the slope of mu vs. variance, to compare to the slope of other kinds of replicates?

# a <- NB_ouput[[1]] %>% pivot_wider(names_from = name )
# summary(lm(data = a, estVariance ~ mu))
  
