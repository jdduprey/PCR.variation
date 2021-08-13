### eDNA Index function
# Helen Casendino
# Aug 13 2021

# make a function to use on 1 dataframe across all PCR columns
# ===================================================
# df: DATAFRAME with 1 column for each PCR run (data = nReads)
# PCRcols: vector of PCR column indices (ie, 3:5)
# ===================================================

calc_eDNA_index <- function(df, PCRcols) { 
  
  for(i in 1:length(PCRcols)){
    
    nRead.col <- colnames(df)[PCRcols[i]] # col name of PCR column
    
    df[,nRead.col] <- df %>%
      group_by(bio) %>% # "convert read counts to proportion with biological sample" 
      mutate(totalReads = sum(.data[[nRead.col]]),
             propReads = .data[[nRead.col]]/totalReads) %>% 
      
      group_by(Hash) %>% # "scale the resulting proportion to the largest observed proportion across samples"
      mutate(index = propReads/(max(propReads)), .keep = "none") %>% ungroup() %>% select(-Hash)
  }
  colnames(df)[PCRcols] <- c("PCR1_index", "PCR2_index", "PCR3_index") # rename columns to "index"
  return(df)
} 