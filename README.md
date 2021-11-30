# Technical & Biological Variation in eDNA Metabarcoding Processes 

This project analyzes existing eDNA metabarcoding data to quantify variation arising from technical and biological processes. We use metabarcoding data from: 

 * [Hood Canal eDNA data](https://github.com/ramongallego/eDNA.and.Ocean.Acidification.Gallego.et.al.2020) - Gallego 2020
 * [Eel Grass Halo data](https://github.com/invertdna/EelgrassHalo) - Jacobs-Palmer 2020 
 

## Progress
1. Starting with ASV hash table scripts, filters out seriously low read PCR replicates and convert read data to proportional data `Hash.Proportions.Rmd`
2. Finds Bray-Curtis dissimilarities among technical replicates `find_PCR_BCDs()` and biological replicates `find_bottle_BCDs()` in `braycurtis.comm.var.Rmd`
3. IN PROGRESS: Quantifies ASV-level variation and the impact of rare species by fitting the proportions of each ASV across replicates to a negbin distribution and observing how variation changes with relative Hash rarity. 

## Additonal Datasets to Input into Community Level Functions
* [Puget Sound Urbanization Gradient](https://datadryad.org/stash/dataset/doi:10.5061/dryad.04tq4)
