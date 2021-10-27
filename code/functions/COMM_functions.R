# ===================================
# saving functions so they can be called in other scripts
# 10/26/2021
# find_PCR_BCDs()

# Computes pairwise Bray-Curtis dissimilarities among technical replicates from 
# each biological sample. Returns a vector of those values.
# Args:
#   df: A data frame with columns for bio, tech, hash, prop, and seq_run.
#   seq_runs: A vector containing the numbers of the seq runs to be analyzed.
find_PCR_BCDs <- function(df, seq_runs = unique(df$seq_run)) {  # Includes all seq runs by default
  df <- df %>% 
    ungroup() %>%  # Just in case
    filter(seq_run %in% seq_runs) %>% 
    select(bio, tech, hash, prop)
  
  bottles <- unique(df$bio)
  bcds_upper_bound <- choose(max(df$tech), 2) * length(bottles)
  bcds <- rep(NA, bcds_upper_bound)  # Vector to store BCDs
  i <- 1
  
  # Iterate over each bottle and calculate BCDs among PCRs
  for (bottle in bottles) {
    # Data wrangling in preparation for the vegdist function
    bottle_data <- df %>% 
      filter(bio == bottle) %>% 
      pivot_wider(names_from = tech, values_from = prop, values_fill = 0) %>% 
      select(matches("\\d{1, }"))  # Columns that have numbers as names (i.e., represent PCRs)
    flipped_bottle_data <- t(bottle_data)
    print(flipped_bottle_data)
    num_PCR_pairs <- choose(nrow(flipped_bottle_data), 2)
    dis <- vegdist(flipped_bottle_data)
    bcds[i:(i + num_PCR_pairs - 1)] <- as.vector(dis)
    i <- i + num_PCR_pairs
  }
  
  bcds <- as.numeric(na.omit(bcds))
  
  # Sanity check
  print("Bray-Curtis Dissimilarities:")
  print(bcds)
  print(paste("Mean BCD value:", mean(bcds), sep = " "))
  
  # Plot
  hist(
    bcds,
    col = viridis::plasma(3, 0.4, 0.7),
    main = "PCR Variation",
    xlab = "Pairwise Bray-Curtis Dissimilarities"
  )
  
  bcds
}
