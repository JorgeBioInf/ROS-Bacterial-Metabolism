# Find structural features susceptible to ROS in monomeric proteins of *P. putida*
Here are stored all the scripts used to extract potentially **ROS-inhibitable structural features** present in **monomeric proteins**. Nevertheless, this workflow is adapted specifically to the proteome available in the **genomic-scale model iJN1480** of *P. putida*.

## Followed workflow
Each script's details can be found inside their source code.

### 1. Proteome extraction and identification of proteins

1. `PP_metab_proc.py`: Extracts KEGG identifiers from GEM `iJN1480.txt` and processes them.
   1. **INPUT**
      - `iJN1480.txt`: GEM model.
      - `(ID_relationships.json)`: optional dictionary with KEGG–UniProt identifiers relationship, accelerates the process. An example can be found [here](../../data/ID_relationships.json).
   2. **OUTPUT**
      - `Monomers_to_model/`: folder with FASTA sequences of proteins without available predicted structure.
      - `Monomers_predictions/`: folder with the monomeric `.pdb` structures.
      - `Complex_fastas/`: folder with FASTA sequences of extracted complexes.
      - Optional: `ID_relationships.json` dictionary.


2. `monomeric_subfolders.py`: Creates subfolders for each protein stored in `Monomers_predictions/` folder.

3. `AF.sh`: Recursively feeds AlphaFold2 with sequences stored in `Monomers_to_model` or `Complex_fastas` (with inner subfolders for each protein/complex). Requires `AF_BestRanked_and_ipTM.py` to extract the best prediction in each case. Currently NOT available. 
   1. **INPUT**
      - `Monomers_to_model/` or `Complex_fastas`.
      - Type of analysis should be specified (monomeric or multimeric).
   2. **OUTPUT**
      - `ranked_0.pdb` file for each protein/complex.

3. `Name_transformer.py`: Changes `ranked_0.pdb` predictions file names to `{PP_XXXX|pWW0_XXXX}.pdb`.
   1. **INPUT**
      - `Monomers_predictions`.
   2. **OUTPUT**
      - `{{PP_XXXX|pWW0_XXXX}.pdb}`.

### 2. UniProt information extraction

1. `Uniprot_Entries_Retrieval.py`: Retrieves UniProt entries of all monomers in `.json` format.
   1. **INPUT**
      - `Monomers_predictions/`.
      - `ID_relationships.json`.
      - Type of analysis should be specified (monomeric or multimeric).
   2. **OUTPUT**
      - `{KEGG_ID}_UniProt_Features.json` for each monomer (KEGG_ID = PP_XXXX or pWW0_XXXX). 
   

### 3. Search of ROS-susceptible structural features

1. `ROS_SR_finder.py`: Extracts susceptibility information of functional sites in the protein. 
   1. **INPUT**
      - `Monomers_predictions/`.
      - Type of analysis should be specified (monomeric or multimeric).
   2. **OUTPUT**
      - `{KEGG_ID}_Susceptibility_Scores.json` for each monomer.
     
2. `Disulfide_Bonds.py`: Calculates potential disulfide bonds in proteins.  
   1. **INPUT**
      - `Monomers_predictions/`.
   2. **OUTPUT**
      - `{KEGG_ID}_Disulfide_Bonds.json` for each monomer.
   
2. `Cofactors_and_CheBI(_and_AUX).py`: Obtains cofactors with transition metals within each proteic structure. Version ending in `"_and_AUX"` uses an auxiliary dataset for remaining CheBI identifiers (aux_chemical_data.txt).   
   1. **INPUT**
      - `Monomers_predictions/`.
   2. **OUTPUT**
      - `cofactors.json` dictionary with the cofactors of each protein.
     
### 4. Classification summary
Its necessary to execute all previous steps in the exposed order. However, this can be done directly through the script used in this last step.

2. `ROS_summary.py`: Generates table with each monomer's characteristics. 
   1. **INPUT**
      - `Monomers_predictions/`.
   2. **OUTPUT**
      - [`ROS_summary_full.tsv`](../../data/ROS_summary_full.tsv).
     
This table is latter curated and statistically analysed. Subsequently, metabolic modelling can be confidently performed integrating ROS susceptibility information. 
