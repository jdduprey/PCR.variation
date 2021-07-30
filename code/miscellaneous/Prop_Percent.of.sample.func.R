# Prop_Percent.of.sample function
# Written by Helen Casendino 7/30/21

# Rationale: So using flat proportion cutoffs (0.02, 0.04, etc.), it seems like Hashes with average proportions from 0 to 0.02 are most responsible for causing PCR variability.
# I want to see what rough percentage the average Hash proportions < 0.02 make up of each bio sample's total Hashes. 

# Function calls: df = dataframe with a column ("bio") for all biological samples, 
#                         and three columns with Hash proportions for each PCR replicate (PCR1_prop, PCR2_prop,PCR3_prop)
#                 bioVector = vector of all unique bio sample IDs
#                 value = the proportion value that you want to look for in each bio sample
# Function output: a vector the same length as bioVector, that gives the percent of the avg Hash proportions that are at or below the value.
#                 So, if the first element of my output is is 0.95 and my input value was 0.02, I know that for my first bio sample, 
#                 Hashes with average proportions at or below 0.02 make up 95% of the data

library(tidyverse)

Prop_Percent.of.sample <- function(df, bioVector, value){
 
  loc_vector <- rep(NA, length = length(bioVector) )
 
   for(i in 1:length(bioVector)){
    meanprops <- df[df$bio == bioVector[i],] %>% group_by(Hash) %>% summarise(meanProp.PCR = mean(c(PCR1_prop, PCR2_prop,PCR3_prop))) %>% arrange(desc(meanProp.PCR))
    first_instance <- (which(meanprops$meanProp.PCR <= value))[1]
    total_hash_n <- length(meanprops$meanProp.PCR)
    loc_vector[i] <- 1-(first_instance/total_hash_n) 
    print(i)
  }
return(loc_vector)
}