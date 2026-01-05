# Some data retrieved during analysis
Intermediate datasets, various results of enrichment analysis and other results can be found here.

### Datasets with ROS susceptibility information for each protein
- `ROS_summary_full.tsv`: with the structural information gathered from [monomeric processing](../scripts/Proteome_processing_monomers).
- `ROS_summary_BRENDA_full.txt`: expanded dataset with aditional information. Obtained from the analysis through [BRENDA_ROS.py](scripts/Stats_analysis/BRENDA_ROS.py) script.
- `Pputida_evidenced.txt`: subset of fully annotated enzymes. Comprises all proteins experimentally evidenced to be susceptible to ROS.
  
### Dictionaries used in the workflow
- `ID_relationships.json`: contains the relationship between KEGG IDs and their corresponding UniProt ID.
- `UniProt_EC.json`: stores the relationship between UniProt IDs and their associated Enzyme Commission number (EC).

### Enrichment analysis of structural features
- `SF_byClades_enrichment.xlsx`: enrichment analysis of structural features between a specific enzyme class and all other classes and enrichment of each structural feature between experimentally evidenced and no evidenced data, divided by protein class.
- `EC.tsv`: enrichment analysis of enzyme classes within the experimentally evidenced subset of proteins.

### Metabolic analysis results
- `fba_comparison.tsv`: flux changes across all reactions when the original model is constrained using ROS susceptibility information.
- `fva_comparison.tsv`: changes in flux ranges across ROS-inhibitable reactions when the original model is constrained using ROS susceptibility information.
- `metab_enrichment.tsv`: enrichment of metabolic subsystems in peripherical metabolism. 
