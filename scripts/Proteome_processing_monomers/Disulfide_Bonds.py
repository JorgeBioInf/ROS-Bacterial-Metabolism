
"""
Disulfide_Bonds.py
----------------------
Computes and retrieves disulfide bonds for each protein in the proteome

Author: Jorge Marcos Fernández
Date: 2025-10-01
Version: 1.0

Usage:
    python Disulfide_Bonds.py folder_with_predictions

Output:
    {protein_ID}_Disulfide_Bonds.txt

Dependencies:
    - Bio.PDB
    - numpy
    - sys
    - pandas
    - pathlib
    - json
    - os

Notes:
    - Requires structures (ranked_0.pdb) and uniprot features ({protein_ID}_UniProt_Features.json) for each protein
"""


# PACKAGES
from Bio.PDB import PDBParser # type: ignore
from Bio.PDB.vectors import calc_dihedral  # type: ignore
import numpy as np # type: ignore
import argparse
import sys
import pandas as pd  # type: ignore
from pathlib import Path  # type: ignore
import json
import os


## FUNCTIONS
def anotated_diS(uniprot_info):
    """Extracts anotated disulfide bonds which might have not been detected"""
    domains = []

    # Get region info
    if 'features' not in uniprot_info.keys():
        return None
    
    regions = uniprot_info['features']

    # Search for catalytic domains: 
    for region in regions:
        reg = region['type']  
        if reg == 'Disulfide bond':

            # Start and end positions
            start = region['location']['start']['value']
            end = region['location']['end']['value']
            pos = f'{start}_{end}'
            if pos not in domains:
                domains.append(pos)
    
    return domains


# CHECK ARGUMENTS

# A target folder should be given as input
if len(sys.argv) != 2:
    print("Use: python Disulfide_Bonds.py folder_with_predictions")
    sys.exit(1)

if not os.path.isdir(sys.argv[1]):
    print(f'Error: unable to find {sys.argv[1]} directory')
    sys.exit(1)


# Define thresholds
PLDDT_LIMIT = 50
B_FACTOR_LIMIT = 30
QMEAN_LIMIT = 0.7

D_INF = 1.5
D_SUP = 2.5
A_INF = 84
A_SUP = 96


### MAIN PROGRAM

# Iterate through the target folder 
folder = Path(sys.argv[1])
count = 0 

for subfolder in folder.iterdir():

    data = []

    # Check if .pdb file is found
    pdb_name = f'{subfolder.name}.pdb'
    pdb_path = subfolder / pdb_name

    if not os.path.exists(pdb_path):
        print(f'Error: no {pdb_name} file found for {subfolder.name}')
        continue

    # Check if folder has already been processed
    out_filename = f'{subfolder.name}_Disulfide_Bonds.txt'
    out_path = subfolder / out_filename

    if os.path.exists(out_path):
        print(f'Skipping folder {subfolder.name}: results already available in {out_filename}')
        count += 1
        continue

    # Check if UniProt info file in subfolder
    uniprot_filename = f'{subfolder.name}_UniProt_Features.json'
    uniprot_filepath = subfolder / uniprot_filename

    if not os.path.exists(uniprot_filepath):
        print(f'Error: no {uniprot_filename} found in folder {subfolder.name}')
        print('Please, try executing Uniprot_Entries_Retrieval.py')
        continue

    # Create PDB parser from Bio.PDB
    parser = PDBParser(QUIET=True)
    protein = parser.get_structure("candidate", pdb_path)

    # Read UniProt info file
    with open(uniprot_filepath, 'r') as uf:
            uniprot_info = json.load(uf)
        
    # Get domains from UniProt info
    anotated = anotated_diS(uniprot_info)
    if anotated:
        data = anotated
        count += 1

    # Type of protein (experimental or AlphaFold)
    mAF = True
    mSM = False

    if "SWISS-MODEL" in protein.header["name"].upper():
        mSM = True
        mAF = False

    ## 1. DISULFIDE BONDS   # Based on previous code by Víctor Sanz Salvador

    # CYS detection according to their quality score (depends on the type of model)
    cys_res = []
    for residue in protein.get_residues():

        # AlphaFold prediction
        if mAF:
            if residue.resname == "CYS":  # Check if cysteine
                # Quality control
                if residue["C"].get_bfactor() > PLDDT_LIMIT:
                    cys_res.append(residue)

        # SWISS-MODEL prediction
        elif mSM:
            if residue.resname == "CYS":
                # Quality control
                if residue["C"].get_bfactor() > QMEAN_LIMIT:
                    cys_res.append(residue)

        # Experimental structure
        else:
            if residue.resname == "CYS":
                if residue["C"].get_bfactor() < B_FACTOR_LIMIT:
                    cys_res.append(residue)

    # Compute all disulfide bonds between every possible pair of cysteines
    diS = []

    for i, cys1 in enumerate(cys_res):
        for j, cys2 in enumerate(cys_res):
            if i < j:
                # Extract atom coordinates
                sulf1 = cys1["SG"]
                sulf2 = cys2["SG"]

                # Compute distance
                d = sulf1 - sulf2

                if D_INF <= d <= D_SUP:
                    c1 = cys1["CB"].get_vector()
                    c2 = cys2["CB"].get_vector()
                    s1 = sulf1.get_vector()
                    s2 = sulf2.get_vector()

                    dihedral = calc_dihedral(c1, s1, s2, c2) * (180 / np.pi)

                    if A_INF <= abs(dihedral) <= A_SUP:
                        diS_bond = [cys1, cys2, d, dihedral]
                        diS.append(diS_bond)

    # Retrieving the results
    if diS:
        count += 1
        for bond in diS:
            res1 = f"{bond[0].get_id()[1]}"
            res2 = f"{bond[1].get_id()[1]}"
            bond = f"{res1}_{res2}"
            if bond not in data:
                data.append(bond)

    # Write output file
    if data:
        with open(out_path, 'w') as file:
            file.write(str(data))
            print(f'Successfully written {out_path}')

print(f'{count} proteins with disulfide bonds were found!')
