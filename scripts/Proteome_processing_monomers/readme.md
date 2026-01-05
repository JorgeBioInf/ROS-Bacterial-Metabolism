# Find structural features susceptible to ROS in monomeric proteins of *P. putida*
Here are stored all the scripts used to extract potentially **ROS-inhibitable structural features** present in **monomeric proteins**. Nevertheless, this workflow is adapted specifically to the proteome available in the **genomic-scale model iJN1480** of *P. putida*.

## Followed workflow
Each script's details can be found inside their source code

### 1. Proteome extraction and identification of proteins

1. `PP_metab_proc.py`: Extracts KEGG identifiers from GEM `iJN1480.txt' and processess them.

  **INPUT**

  -`iJN1480.txt`: GEM model
  -`(ID_realtionships.json): optional dictionary with KEGG-UniProt identifiers relationship. Accelerates the process. An example can be found here.
