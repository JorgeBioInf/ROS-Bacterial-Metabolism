
"""
Monomer_name_transformer.py
----------------------
Changes the name of .pdb files for monomers ({UniProt_ID}.pdb --> {Protein_ID}.pdb)

Author: Jorge Marcos Fernández
Date: 2025-09-20
Version: 1.0

Usage:
    python Monomer_name_transformer.py folder_with_predictions ID_relationships.json

Dependencies:
    - Python 3.10+
    - pathlib
    - sys
    - os
    - json

Notes:
    - Requires file ID_relationships.json

"""


# PACKAGES
from pathlib import Path
import sys
import os
import json


# CHECK ARGUMENTS
if len(sys.argv) != 3:
    print('Use: python3 Monomer_name_transformer.py folder_with_predictions ID_relationships.json')
    sys.exit(1)


# MAIN PROGRAM

with open(sys.argv[2], 'r') as dic:
    ID_dic = json.load(dic)

Uni_to_KEGG = {v: k for k, v in ID_dic.items()}


folder = Path(sys.argv[1])

# Transform the name of each subfolder in folder_with_predictions
for file_path in folder.iterdir():
    if file_path.is_file():
        
        Uni_id = file_path.stem
        
        if not Uni_id in Uni_to_KEGG.keys():
            print(f'Error when processing file {file_path.name}')
            print(f'No correspondence found for UniProt ID {Uni_id}')
            continue
            
        # Search for associated KEGG id
        KEGG = Uni_to_KEGG[Uni_id]
        
        new_name = f"{KEGG}.pdb"
        new_path = folder / new_name
        
        # Avoid writing already existing files
        if new_path.exists():
            print(f"Warning: {new_name} already exists. Skipping renaming of {file_path.name}")
            continue
        
        # Rename
        file_path.rename(new_path)
        print(f"{file_path.name} -> {new_name}")

