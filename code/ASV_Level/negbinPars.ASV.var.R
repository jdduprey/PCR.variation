##### =====ASV variation using Negative Binomial Parameters=======
# Uses Ryan's code to fit PCR reads for each ASV (usually 3) to a negbin distribution,
# and estimating mu and phi
# 10/28/21
##### ===================================================

### ==Dependencies====
library(tidyverse)
library(stringi)
library(fitdistrplus)
library(ggpubr)
library(here)
# ===============

##### =======READ IN DATA====
all_data <- read_csv(here("data/all_data.csv"))


### ======Ryan's code, updated for Hood Canal Data===========

# "first couple bits are creating two functions:
# first to get the variance of a negbin distribution, given the parameters mu and phi
# Second to use fitdist() to find the most likely parameters mu and phi (which R calls "size" in this case), given some data vector.
# Then the tidyr code takes a dataset (here, called dat.mifish, but you can use your dataset), groups by species (ID_mifish) and biological sample (station_id),
# and then fits a negative binomial distribution  to each vector of 3 observations for each species at each site.
# Then it plots the result."

varNB <- function(mu, phi) {
  mu + ((mu^2) / phi)
} # calculates variance

getNB <- function(x) {
  require(fitdistrplus)

  return(
    (fitdist(unlist(x), "nbinom"))$estimate
  )
} # estimates mu and phi

# Below: estimates parameters for each Hash (2-3 PCR reads/Hash)
# df needs cols: Miseq_run, bio (bottle, ie PO2017A), Hash, and nReads
# rep_type is either "technical" for PCR variation or "biological" for bottle var
# spits out list of:
# 1) df of mu, estimated variance and phi, and original reads for each ASV
#   and 2) plot of mu vs. estimated variance and mu vs. phi


estimates_pars_byHash <- function(df, rep_type) {
  if (rep_type == "technical") {
    grouping_df <- df %>%
      mutate(bottle = paste(site, bio, sep = "_"))
      filter(bottle %in% unique(df$bottle)) %>%
      group_by(bottle, hash)
  }

  if (rep_type == "biological") {
    grouping_df <- df %>%
      filter(site %in% unique(df$site)) %>% # WHAT IS THIS LINE DOING?
      group_by(site, hash)
  }

  list_NBpars <- grouping_df %>%
    filter(sum(reads) > 0) %>%
    filter(n() > 1) %>%
    dplyr::select(reads) %>%
    nest() %>%
    mutate(mlmodels = map(data, getNB))

  df_NBpars <- list_NBpars$mlmodels %>%
    as.data.frame() %>%
    t() %>%
    as_tibble() %>%
    rename(phi = size) %>%
    mutate(estVariance = varNB(mu, phi)) %>%
    mutate(Reads = list_NBpars$data) %>%
    pivot_longer(-c(mu, Reads))

  NBpars_plot <- df_NBpars %>%
    ggplot(aes(x = log(mu), y = log(value))) +
    geom_point() +
    facet_wrap(~name)

  return(list(df_NBpars, NBpars_plot))
}

NB_ouput <- estimates_pars_byHash(all_data, rep_type = "biological")


# Helen: If the curve of mu vs. variance is dependent on phi, we also want to
# keep track of the slope of mu vs. variance, to compare to the slope of other kinds of replicates?
