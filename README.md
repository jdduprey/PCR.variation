# Technical & Biological Variation

Project using existing eDNA metabarcoding data to quantify variation arising from technical and biological processes. Uses metabarcoding data from Hood Canal and San Juan Island (Gallego et al. 2020)
 * [Hood Canal edna repo](https://github.com/ramongallego/eDNA.and.Ocean.Acidification.Gallego.et.al.2020)

## Progress
1. Starting with ASV hash table scripts, filters out seriously low read PCR replicates and convert read data to proportional data (Hash.Proportions.Rmd)
2. Finds Bray-Curtis dissimilarities among technical replicates (fn: find_PCR_BCDs) and biological replicates (fn: find_bottle_BCDs) (braycurtis.comm.var.Rmd)
3. IN PROGRESS: Quantifies ASV-level variation and the impact of rare species by fitting the proportions of each ASV across replicates to a negbin distribution and observing how variation changes with relative Hash rarity. 
