"""
AF_BestRanked_and_ipTM.py
----------------------
Retrieves UniProt entries for each protein in folder_with_predictions

Author: Jorge Marcos Fernández
Date: 2025-09-26
Version: 1.0

Usage:
    python Uniprot_Entries_Retrieval.py folder_with_predictions monomer|multimer ID_relationships.json

Output:
    {protein_ID}_UniProt_Features.json

Dependencies:
    - requests
    . pathlib
    - sys
    - json

Notes:
    - Requires folder_with_predictions and ID_relationships.json ready
"""


import requests
import sys
from pathlib import Path
import json


if len(sys.argv) != 4:
    print("Use: python Uniprot_Entries_Retrieval.py folder_with_predictions type_of_prediction ID_relationship_dictionary")
    print("Type of prediction = monomer if the directory contains only monomeric predictions")
    print("Type of prediction = multimer if the directory contains only multimeric predictinos")
    sys.exit(1)

if sys.argv[2] != 'multimer' and sys.argv[2] != 'monomer':
    print(f"Error: unknown prediction type {sys.argv[2]}")
    print("Type of prediction = monomer if the directory contains only monomeric predictions")
    print("Type of prediction = multimer if the directory contains only multimeric predictions")
    sys.exit(1) 


# Choose whether monomers or multimers are being processed
mon = True
mult = False

if sys.argv[2] == 'multimer':
    mon = False
    mult = True



# Load KEGG_ID - UniProt-ID relationship
dict_path = Path(sys.argv[3])

with open(dict_path, 'r') as dict_file:
    dic = json.load(dict_file)


# MAIN PROGRAM
folder = Path(sys.argv[1])

if not folder.is_dir():
    print(f'Error: no folder named {folder.name}')
    sys.exit(1)

# Monomers
if mon:
    for subfolder in folder.iterdir():
        ID = subfolder.name
        
        out_filename = f'{subfolder.name}_UniProt_Features.json'
        out_path = subfolder / out_filename
        
        if out_path.exists():
            print(f'File {out_filename} already exists. Skipping')
            continue
        
        print(f'Processing protein {ID} ...')
        
        if ID not in dic.keys():
            print(f'Error with protein {ID}. No PP_XXXX ID found.')
            continue
        
        Uni_ID = dic[ID]
        
        if not Uni_ID:
            print(f'Error with protein {ID}. No UniProt ID found.')
            continue
        
        url = f'https://rest.uniprot.org/uniprotkb/{Uni_ID}.json'
        response = requests.get(url)
        
        if not response.ok:
            print(f'There was an error when searching protein {Uni_ID} in UniProt')
            continue
        
        data = response.json()
        
        with open(out_path, 'w') as out_file:
            json.dump(data, out_file, indent=4)       


# Multimers            
if mult:
    for subfolder in folder.iterdir():
        
        IDs = subfolder.name.split('-')
        print(f'Processing {subfolder.name} ...')
        
        for ID in IDs: 
            if ID not in dic.keys():
                print(f'Error with protein {ID}. No PP_XXXX ID found.')
                continue
            
            Uni_ID = dic[ID]

            if not Uni_ID:
                print(f'Error with protein {ID}. No UniProt ID found.')
                continue

            url = f'https://rest.uniprot.org/uniprotkb/{Uni_ID}.json'
            response = requests.get(url)

            if not response.ok:
                print(f'There was an error when searching protein {Uni_ID} in UniProt')
                continue
            
            data = response.json()
        
            out_filename = f'{ID}_UniProt_Features.json'
            out_path = subfolder / out_filename
        
            with open(out_path, 'w') as out_file:
                json.dump(data, out_file, indent=4)

    


