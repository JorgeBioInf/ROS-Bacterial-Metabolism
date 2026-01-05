"""
ROS_summary.py
----------------------
Merges ROS susceptibility information retrieved from previous programs to generate a summary table.
Computes disulfide bonds mean minimum distance to functional sites of the protein (when applicable).

Author: Jorge Marcos Fernández
Date: 2025-11-02
Version: 1.0

Usage:
    python ROS_summary.py target_directory

Output:
    - ROS_summary.tsv

Dependencies:
    - json
    - Bio
    - sys
    - numpy
    - pathlib
    - subprocess
    - pandas
    - ast

Notes:
    - target_directory must contain all information gathered from previous programs in the workflow. Hence,
    all these need to be previously executed correctly. However, they can be directly executed through this program.
    - For executing previous codes, these must be found in the same directory.
"""


## MODULES
import json
from Bio.PDB import PDBParser
import sys
import numpy as np 
from pathlib import Path
import subprocess
import pandas as pd
import ast

## FUNCTIONS
def calculate_distance(res1, res2):
    "Returns C-alpha distance between two residues"
    diff  = res1["SG"].coord - res2["CA"].coord
    return np.sqrt(np.sum(diff * diff))

def get_residue_by_number(protein, res_num):
    "Returns residue name from its number"
    for res in protein.get_residues():
        if res.id[1] == res_num:
            return res
    return None

def DiS_AS_distance(DiS_info, ROS_info, protein):
    """Calculates the mean minimum distance between S carbons participating in
    disulfide bonds and the Active Sites or Binding Sites residues"""
    
    distances = np.zeros((len(DiS_info), len(ROS_info)-1), np.float64)
    
    for row, DiS in enumerate(DiS_info):
        S_pos = DiS.split('_')
        S_res = []
        
        for S in S_pos: 
            S_res.append(get_residue_by_number(protein, int(S)))
            
        for col, (desc, info) in enumerate(ROS_info.items()):  
            if desc == 'Total scores':
                continue
                
            AS_res = []  
            start_pos = int(info['start'])
            end_pos = int(info['end'])
            positions = list(range(start_pos, end_pos+1))
            
            for pos in positions:
                AS_res.append(get_residue_by_number(protein, pos))
            
            # Calculate distances between S atoms and AS residues
            d = []
            for res in AS_res:
                d1 = calculate_distance(S_res[0], res)
                d2 = calculate_distance(S_res[1], res)
                # Mean distance between S atoms and a specific residue
                mean_dist = (d1 + d2) / 2    
                d.append(mean_dist)
            
            # Mean distance between S atoms and all Active Site's residues
            if d:
                distances[row, col] = np.mean(d)
            else:
                distances[row, col] = None
    
    # Minimum distance between each S-S and all Active Sites / Binding Sites
    if distances.any():
        min_dist = np.min(distances, axis=1)
    else:
        min_dist = 0

    # Mean minimum distance
    mmd = np.mean(min_dist)
    
    return mmd

## ARGUMENT CHECK

# Check if input directory exists
if len(sys.argv) != 2:
    print('Use: python3 ROS_summary.py target_directory')
    sys.exit(1)

try:
    target_dir = Path(sys.argv[1])
except:
    print(f'Error: no folder named {sys.argv[1]}!')
    sys.exit(1)

# Read ID_relationships.json
try:
    with open('ID_relationships.json', 'r') as file:
        KEGG2UniProt = json.load(file)
except Exception as e:
    print(f'Error when loading ID_relationships.json: {e}')
    sys.exit(1)


## INTERFACE
print('Welcome! This program processes a list of proteins and retrieves a ROS susceptibility summary')
print('The following python scripts should be executed for the target folder:')
print('[1] Uniprot_Entries_Retrieval.py')
print('[2] ROS_SR.py')
print('[3] Disulfide_Bonds.py')
print('[4] Cofactors_and_CheBI.py')
print('')
print('You can choose to execute them throughout this program if you have not previously done so')
print('')
print('Remember that ID_relationship.json file is needed')

get_uni_info = False
get_ROS_info = False
get_DiS_info = False
get_cof_info = False

while True:
    i0 = input('\nPlease, tell which kind of analysis do you want [monomer | multimer] ').lower()
    if i0 in ('monomer', 'multimer'):
        analysis = i0
        break
    print(f'Invalid option "{i0}". Please select [monomer | multimer]')

while True:
    i1 = input('Do you have PP_XXXX_UniProt_Features.json files ready? [y|n] ').lower()
    if i1 in ('y', 'n'):
        get_uni_info = (i1 == 'n')
        break
    print(f'Invalid option "{i1}". Please select [y|n]')

while True:
    i2 = input('Do you have PP_XXXX_Susceptibiliy_Scores.json files ready? [y|n] ').lower()
    if i2 in ('y', 'n'):
        get_ROS_info = (i2 == 'n')
        break
    print(f'Invalid option "{i2}". Please select [y|n]')

while True:
    i3 = input('Do you have PP_XXXX_Disulfide_Bonds.txt files ready? [y|n] ').lower()
    if i3 in ('y', 'n'):
        get_DiS_info = (i3 == 'n')
        break
    print(f'Invalid option "{i3}". Please select [y|n]')

while True:
    i4 = input('Do you have cofactors.txt file ready? [y|n] ').lower()
    if i4 in ('y', 'n'):
        get_cof_info = (i4 == 'n')
        break
    print(f'Invalid option "{i4}". Please select [y|n]')


## EXTRACT RESULTS

# If UniProt info is needed
if get_uni_info:
    print('\nRetrieving Uniprot entries for each protein ...')
    result = subprocess.run(
        ["python3", "Uniprot_Entries_Retrieval.py", target_dir, i0, "ID_relationships.json"])

# If ROS susceptibility info is needed
if get_ROS_info:
    opt = f'{i0}s'
    print('\nRetrieving ROS susceptibility scores ...')
    result = subprocess.run(
    ["python3", "ROS_SR_DEFINITIVO.py", target_dir, opt])

# If disulfide bonds info is needed
if get_DiS_info:
    print('\nCalculating disulfide bonds ...')
    result = subprocess.run(
    ["python3", "Disulfide_Bonds.py", target_dir])
    
# Get cofactors information
if get_cof_info:
    print('\nRetrieving susceptible cofactors ...')
    result = subprocess.run(
    ["python3", "Cofactors_and_CheBI.py", target_dir, "chemical_data.tsv.gz"])


# Read cofactors.json
try:
    with open('cofactors.json', 'r') as file:
        cofactors = json.load(file)
except Exception as e:
    print(f'Error when loading cofactors.json: {e}')
    sys.exit(1)


## MAIN PROGRAM
KEGG_IDs = []
UniProt_IDs = []
gene_names = []
protein_names = []
AS_residue_number = []
AS_mean_distance = []
AS_cys = []
BS_residue_number = []
BS_mean_distance = []
BS_cys = []
DiS_number = []
DiS_mean_distance = []
cofactor = []

for subfolder in sorted(target_dir.iterdir()):
    
    prot = subfolder.name
    print(f'Processing protein {prot} ...')
    
    # Read files
    uniprot_filename = f'{prot}_UniProt_Features.json'
    uniprot_path = subfolder / uniprot_filename
    ROS_filename = f'{prot}_Susceptibility_Scores.json'
    ROS_path = subfolder / ROS_filename
    DiS_filename = f'{prot}_Disulfide_Bonds.txt'
    DiS_path = subfolder / DiS_filename
    pdb_filename = f'{prot}.pdb'
    pdb_path = subfolder / pdb_filename
    
    try:
        with open(uniprot_path, 'r') as uf:
            uniprot_info = json.load(uf)
    except Exception as e:
        print(f'Error reading {uniprot_filename} file')
        print('Please, check your files')
        continue
    
    try:
        with open(ROS_path, 'r') as rf:
            ROS_info = json.load(rf)
    except Exception as e:
        print(f'Error reading {ROS_filename} file')
        print('Please, check your files')
        continue
    
    try:
        with open(DiS_path, 'r') as df:
            DiS_info = ast.literal_eval(df.read())

    except Exception as e:
        # print(f'Error reading {DiS_filename} file', file = sys.stderr)
        DiS_info = []

    try:
        parser = PDBParser(QUIET=True)
        protein = parser.get_structure("candidate", pdb_path)
    except Exception as e:
        print(f'Error reading {pdb_filename} file')
        continue

    # Get gene and protein names
    try:
        gene_name = uniprot_info['genes'][0]['geneName']['value']
    except:
        gene_name = None
        
    try:
        prot_name = uniprot_info['proteinDescription']['recommendedName']['fullName']['value']
    except:
        try:
            prot_name = uniprot_info['proteinDescription']['submissionNames'][0]['fullName']['value']
        except:
            prot_name = None

    # Calculate disulfide bonds' mean minimum distance to Active Sites or Binding Sites
    if DiS_info:
        DiS_score = DiS_AS_distance(DiS_info, ROS_info, protein)
    else:
        DiS_score = 0
    
    # Get cofactors
    all_cof = []
    if prot in cofactors:
        cof_info = cofactors[prot]
        for i in cof_info:
            cof_name = i['Cofactor name']
            cof_formula = i['Formula']
            cof = f'{cof_name} ({cof_formula})'
            all_cof.append(cof)

    # Store results  
    KEGG_IDs.append(prot)
    UniProt_IDs.append(KEGG2UniProt[prot])
    gene_names.append(gene_name)
    protein_names.append(prot_name)
    AS_residue_number.append(ROS_info['Total scores']['s1'])
    AS_mean_distance.append(round(ROS_info['Total scores']['s2'],2))
    AS_cys.append(ROS_info['Total scores']['s3'])
    BS_residue_number.append(ROS_info['Total scores']['s4'])
    BS_mean_distance.append(round(ROS_info['Total scores']['s5'],2))
    BS_cys.append(ROS_info['Total scores']['s6'])
    DiS_number.append(len(DiS_info))
    DiS_mean_distance.append(round(DiS_score,2))
    cofactor.append(' AND '.join(all_cof))


## SUMARY TABLE
data = {
    "KEGG ID": KEGG_IDs,
    "UniProt ID": UniProt_IDs,
    "Gene": gene_names,
    "Protein": protein_names,
    "# SR in AS": AS_residue_number, 
    "Mean SR - AS minimum distance": AS_mean_distance,
    "cys in AS": AS_cys,
    "# SR in BS": BS_residue_number, 
    "Mean SR - BS minimum distance": BS_mean_distance,
    "cys in BS": BS_cys,
    "# Disulfide bonds": DiS_number,
    "Mean S-S - AS/BS minimum distance": DiS_mean_distance,
    "Cofactor": cofactor,
}

df = pd.DataFrame(data)
df.to_csv("ROS_summary.tsv", sep="\t", index=False)

print('\nFile ROS_summary.tsv succesfully written!')
