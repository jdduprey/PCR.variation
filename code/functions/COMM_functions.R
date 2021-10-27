# ==========================================================================
# saving functions so they can be called in other scripts
# 10/26/2021
# find_PCR_BCDs()
# find_bottle_BCDs()
# ==========================================================================

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

# ==========================================================================
# ==========================================================================

# Computes pairwise Bray-Curtis dissimilarities across different biological
# replicates (but not among PCRs from the same bottle). Returns a vector of
# those values and a histogram.
# Args:
#   df: A data frame with columns for bio, tech, hash, prop, and seq_run.
find_bottle_BCDs <- function(df) {
  # Find a length for bcds that is certain to contain all calculated values. We
  # do this by finding the maximum number of replicates from any individual site.
  temp <- df %>%
    mutate(site = gsub(".{1}$", "", bio)) %>%  # Separate site and bottle info from each other
    mutate(bottle = stri_sub(bio, -1)) %>% 
    select(site, bottle, tech) %>% 
    group_by(site) %>% 
    unique() %>%  # Get rid of duplicate rows bc we don't care about the number of hashes
    mutate(n = n())  # Number of replicates for that site
  
  max_replicates <- max(temp$n)
  df <- nest_data(df)  # Data wrangling using helper function
  
  # Notes:
  # a) max_replicates ^ 2 is the max size of any BCD matrix generated from a
  # single site
  # b) We multiply by 0.5 bc at least half the values from each dis_mat matrix
  # will be removed (to eliminate duplicates and intra-bottle comparisons)
  bcds_upper_bound <- ceiling((max_replicates ^ 2) * 0.5 * length(df$site))
  bcds <- rep(NA, bcds_upper_bound)
  i <- 1
  
  for (j in 1:length(df$site)) {
    event <- df$data[[j]]  # A data frame corresponding to one sampling event
    sub_tib <- event %>% select(!c(hash))
    flip_tib <- t(sub_tib)
    dis <- vegdist(flip_tib)
    dis_mat <- as.matrix(dis)
    
    # Removes all cells corresponding to BCDs within the same bottle by using
    # row and column names. Also removes all the duplicates below the diagonal.
    for (col_num in 1:ncol(dis_mat)) {
      # Extracting the first letter (A, B, C, etc) from the column's name
      col_letter <- substr(colnames(dis_mat)[col_num], 1, 1)
      for (row_num in 1:nrow(dis_mat)) {
        # Extracting the first letter (A, B, C, etc) from the row's name
        row_letter <- substr(rownames(dis_mat)[row_num], 1, 1)
        if (! col_letter > row_letter) {
          dis_mat[row_num, col_num] <- NA
        }
      }
    }
    
    bcd_vector <- na.omit(as.vector(dis_mat))  # All BCDs for that sampling event as vector
    len <- length(bcd_vector)
    if (len > 0) {
      bcds[i:(i + len - 1)] <- bcd_vector  # Add the BCD values to our long vector
      i <- i + len  # Update index
    }
  }
  
  bcds <- as.numeric(na.omit(bcds))
  
  # Sanity check
  print("Bray-Curtis Dissimilarities:")
  print(bcds)
  print(paste("Mean BCD value:", mean(bcds), sep = " "))
  
  # Plot
  hist(
    as.numeric(bcds[-1]),
    col = viridis::viridis(3, 0.4, 0.7),
    main = "Bottle Variation",
    xlab = "Bray-Curtis Pairwise Distances"
  )
  
  bcds
}

bottle_output <- find_bottle_BCDs(long_props)
