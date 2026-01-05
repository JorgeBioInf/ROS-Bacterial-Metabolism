# Analysis of structural features susceptible to ROS
Once structural characteristics in proteins had been retrieved, manual curation of the dataset and statistical analyses were performed. All scripts used for these purposes can be found here. 

Details of each script workflow is included in their source code. 

## What does each script do?
1. `BRENDA_ROS.ipynb`: BRENDA database exploration in search of experimental evidence, plots and dataset management.
2. `Enrichment_Analysis.ipynb`: performs different enrichment analysis concerning structural features.
3. `EC_number_retrieval`: creates dictionary [UniProt_EC.json](../../data/UniProt_EC.json) for subsequent analyses.
4. `metrics_distributions.R`: plots distance-related metrics distributions. 
