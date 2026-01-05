"""
Cofactors_and_CheBI_and_AUX.py
----------------------
Retrieves cofactors containing transition metals for each monomer in the proteome

Author: Jorge Marcos Fernández
Date: 2025-09-20
Version: 1.0

Usage:
    python Cofactors_and_CheBI_and_AUX.py folder_with_predictions chemical_data.tsv.gz

Output:
    - cofactors.json

Dependencies:
    - json
    - pandas
    - sys
    - pathlib

Notes:
    - Requires proteome folder Monomer_predictions/ with files {protein_ID}_UniProt_Features.json in the same directory
    - Requires file aux_chemical_data.txt in the same directory
"""

import json
import sys
from pathlib import Path
import pandas as pd


# FUNCTIONS
def transition_metal(chemical_info, aux_chemical_info, CheBI_id, transition_metals):
    """Returns chemical formula (str) of a given CheBI_id if it contains a transition metal"""
    TM = False
    
    try:
        formula = chemical_info.loc[chemical_info["compound_id"].astype(str) == str(CheBI_id), "formula"].iloc[0]
        for tm in transition_metals:
            if tm in formula:
                TM = True
                break
    except Exception as e:
        print(f'Error: CheBI_id {CheBI_id} not found in original chemical dataset!')
        print('Retrieving from auxiliary data ...')
        formula = aux_chemical_info.loc[aux_chemical_info["compound_id"].astype(str) == str(CheBI_id), "formula"].iloc[0]
        for tm in transition_metals:
            if tm in formula:
                TM = True
                break
    
    if TM:
        return formula
    else:
        return None
    

def get_cofactors(data, chemical_info, aux_chemical_info, transition_metals):
    """Returns a dictionary with all cofactors and their CheBI IDs for a given protein
    providing that they contain a transition metal"""

    out_list = []
    used_cof = {}

    # Look for cofactors and CheBI IDs
    if 'comments' in data:
        comments = data['comments']
        c = False
        id_ = False

        for i in range(0, len(comments)):
            if 'cofactors' in comments[i]:
                more_info = comments[i]['cofactors']

                # Explore all anotated cofactors
                for cof in range(0, len(more_info)):
                    
                    # Cofactor name
                    if 'name' in more_info[cof]:
                        cofactor = more_info[cof]['name']
                        c = True

                    # CheBI id
                    if 'cofactorCrossReference' in more_info[cof]:
                        crossref = more_info[cof]['cofactorCrossReference']
                        if crossref['database'] == 'ChEBI': 
                            chebi_id = crossref['id']
                            id_ = True

                    # Retrieve formula                     
                    if id_ and c:  
                        
                        if cofactor in used_cof and used_cof[cofactor]: # Cofactor and formula already retrieved
                            formula = used_cof[cofactor]

                        elif cofactor in used_cof and not used_cof[cofactor]:  # Already analysed cofactor without formula
                            formula = None
            
                        else:  # New cofactor
                            chebi_id_num = chebi_id.replace('CHEBI:','') 
                            formula = transition_metal(chemical_info,aux_chemical_info, chebi_id_num, transition_metals)
                            used_cof[cofactor] = formula

                    # Store information
                    if formula:
                        out_dict = {'Cofactor name':cofactor,
                                    'CheBI ID': chebi_id, 
                                    'Formula': formula}
                        
                        out_list.append(out_dict)
        
    return out_list
                    

# ARGUMENT CHECK

if len(sys.argv) != 3:
    print('Use: python3 ROS_summary.py target_directory chemical_data.tsv.gz')
    sys.exit(1)

# Check input directory
try:
    target_dir = Path(sys.argv[1])
except:
    print(f'Error: no folder named {sys.argv[1]}!')
    sys.exit(1)

# Read chemical information
try:
    df = pd.read_csv(sys.argv[2], compression='gzip', sep = '\t', header=0)
    chemical_info = df[["compound_id", "formula"]]
except Exception as e:
    print(f'Error when loading {sys.argv[2]}: {e}')
    sys.exit(1)


# Read auxiliar chemical information
try: 
    df = pd.read_csv("aux_chemical_data.txt", sep = '\t', header=0)
    aux_chemical_info = df[["compound_id", "formula"]]
except Exception as e:
    print(f'Error when loading aux_chemical_data.txt: {e}')
    sys.exit(1)


# Transition metals list
TM = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', \
      'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', \
        'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', \
            'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn']


# Main program
folder = Path('Monomers_predictions/')
output = {}

for subfolder in folder.iterdir():

    prot = subfolder.name

    target_file = f'{prot}_UniProt_Features.json'
    target_path = subfolder / target_file

    try:
        with open(target_path, 'r') as file:
            data = json.load(file)
    except:
        print(f'Error: no file {target_file} in {subfolder}')
        continue

    # Get cofactors
    cofactors = get_cofactors(data, chemical_info, aux_chemical_info, TM)

    # Store
    if cofactors:
        output[prot] = cofactors
          
with open('cofactors.json', 'w') as file:
    json.dump(output, file, indent=4)
    
    
    

