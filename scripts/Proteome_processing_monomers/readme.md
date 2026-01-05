# Find structural features susceptible to ROS in monomeric proteins of *P. putida*
Here are stored all the scripts used to extract potentially **ROS-inhibitable structural features** present in **monomeric proteins**. Nevertheless, this workflow is adapted specifically to the proteome available in the **genomic-scale model iJN1480** of *P. putida*.

## Followed workflow
Each script's details can be found inside their source code.

### 1. Proteome extraction and identification of proteins

1. `PP_metab_proc.py`: Extracts KEGG identifiers from GEM `iJN1480.txt` and processes them.
   1. **INPUT**
      - `iJN1480.txt`: GEM model.
      - `(ID_relationships.json)`: optional dictionary with KEGG–UniProt identifiers relationship, accelerates the process. An example can be found here.
   2. **OUTPUT**
      - `Monomers_to_model/`: folder with FASTA sequences of proteins without available predicted structure.
      - `Monomers_predictions/`: folder with the monomeric `.pdb` structures.
      - `Complex_fastas/`: folder with FASTA sequences of extracted complexes.
      - Optional: `ID_relationships.json` dictionary.

2. `monomeric_subfolders.py`: Creates subfolders for each protein stored in `Monomers_predictions/` folder.


### 2. UniProt information extraction

1. `Uniprot_Entries_Retrieval.py`: Retrieves UniProt entries of all monomers in `.json` format.
   1. **INPUT**
      - `Monomers_predictions/`.
      - `ID_relationships.json`.
      - Type of analysis should be specified (monomeric or multimeric).
   2. **OUTPUT**
      - `{KEGG_ID}_UniProt_Features.json` for each monomer.
   
