"""
ROS_SR_finder.py
----------------------
Searches ROS-susceptible residues within functional sites of the protein
Returns a .json dictionary for each protein with the following information:
    1) Number of susceptible residues in Active Sites 
    2) Mean minimum distance between susceptible residues outside functional sites and suceptible residues in Active Sites
    3) Number of cysteine residues in Active Sites
    4) Number of susceptible residues in Binding Sites 
    5) Mean minimum distance between susceptible residues outside functional sites and suceptible residues in Binding Sites
    6) Number of cysteine residues in Binding Sites

Author: Jorge Marcos Fernández
Date: 2025-10-10
Version: 1.0

Usage:
    python ROS_SR_finder.py folder_with_predictions monomers|multimers

Output:
    {protein_ID}_Susceptibility_Scores.json

Dependencies:
    - Bio
    - numpy
    - sys
    - pathlib
    - json
    - os

Notes:
    - Requires structures ({protein_ID}.pdb) and uniprot features ({protein_ID}_UniProt_Features.json) for each protein
"""

## LIBRARIES
import json
from Bio.PDB import PDBParser # type: ignore
import sys
import os
from Bio.SeqUtils import seq1 # type: ignore
import numpy as np # type: ignore
from pathlib import Path


## FUNCTIONS
def get_domains(uniprot_info):
    "Returns each domain"
    domains = []

    # Get region info
    if 'features' not in uniprot_info.keys():
        print('Error: No features found')
        return None
    
    regions = uniprot_info['features']

    # Search for catalytic domains: 
    for region in regions:
        reg = region['type']  
        if reg == 'Active site' or reg == 'Site' or reg == 'Binding site':

            # Description
            desc = region['description'] 

            # Start and end positions
            start = region['location']['start']['value']
            end = region['location']['end']['value']

            # Store domain
            dom = {'Type': reg, 'Desc': desc, 'start': start, 'end': end}
            domains.append(dom)
    
    if domains:
        return domains
    
    else:
        return None
            
def calculate_distance(res1, res2):
    "Returns C-alpha distance between two residues"
    diff  = res1["CA"].coord - res2["CA"].coord
    return np.sqrt(np.sum(diff * diff))

def get_residue_by_number(structure, res_num):
    "Returns redidue name from its number"
    for res in structure.get_residues():
        if res.id[1] == res_num:
            return res
    return None

def total_scores(scores):
    """Returns total susceptibility score for each domain in each chain"""
    s1 = 0
    s2 = []
    mean_s2 = 0
    s3 = 0
    s4 = 0
    s5 = []
    mean_s5 = 0
    s6 = 0 

    for (desc, info) in scores.items():
        if info['type'] == 'Active site':
            s1 += info['score1']
            s2.append(info['score2'])
            s3 += info['score3']
        else:
            s4 += info['score1']
            s5.append(info['score2'])
            s6 += info['score3']

    if s2:
        mean_s2 = sum(s2)/len(s2)
    if s4:
        mean_s5 = sum(s5)/len(s5)
    tot_scores = {'s1':s1,
                  's2':mean_s2,
                  's3':s3,
                  's4':s4,
                  's5':mean_s5,
                  's6':s6
                  }
    
    return tot_scores

def susceptibility_score(uniprot_info, protein, susceptible, domains, chains, chain):
    """Returns susceptibility score for each domain in each chain, considering susceptible
     amino acids in the Active Site and others in the structure"""
     
    scores = {}
    
    # Get full sequence
    seq = uniprot_info['sequence']['value']

    # Explore each domain
    for domain in domains: 

        # Get domain features
        CA = []
        Others = []
        start = domain['start']
        end = domain['end']
        typ = domain['Type']

        subseq = seq[start-1:end] 

        # Check whether UniProt and .pdb domains have the same sequence
        pdb_subseq = chains[chain][start-1:end]

        if subseq != pdb_subseq:
            print(f'Error: PDB and UniProt sequences of {typ} do not match!')
            print(f'Error: PDB and UniProt sequences of {typ} do not match!', file=sys.stderr)
            continue

        # Get susceptible amino acids in CS
        count = 0 
        cys_count = 0

        for idx, aa in enumerate(subseq): 
            CA.append(f'{aa}_{start + idx}')
            if aa in susceptible:
                count += 1
            if aa == 'C':
                cys_count += 1

        # Get other susceptible amino acids
        for idx, aa in enumerate(seq):
            aa_name = f'{aa}_{idx + 1}'
            if aa in susceptible and aa_name not in CA:
                Others.append(f'{aa}_{idx + 1}')

        # Check if any list is empty
        if len(Others) == 0:
            print(f"No susceptible residues outside domain {start}-{end} (chain {chain})")
            info = {'type':typ,
                    'start':start,
                    'end':end,
                    'score1':count,
                    'score2':0,
                    'score3':cys_count
                    }
            
            key = f"{domain['Type']}_{domain['start']}_{domain['end']}"
            scores[key] = info
            continue

        # Compute distance between CS residues and other susceptible residues
        distances = np.zeros((len(CA), len(Others)), np.float64)
        for row, res1 in enumerate(CA):
            for col, res2 in enumerate(Others):
                resnum1 = int(res1.split('_')[1])
                resnum2 = int(res2.split('_')[1])
                r1 = get_residue_by_number(protein, resnum1)
                r2 = get_residue_by_number(protein, resnum2)

                if r1 is None or r2 is None:
                    print(f'Error calculating distance between {res1} and {res2}')
                    print(f'Error calculating distance between {res1} and {res2}', file=sys.stderr)
                    continue

                distances[row, col] = calculate_distance(r1,r2) 

        # Get distance between each CA residue and its closest susceptible residue
        min_dist = np.min(distances, axis=1)
        mean_dist = np.mean(min_dist)

        # Final score
        info = {'type':typ,
                'start':start,
                'end':end,
                'score1':count,
                'score2':mean_dist,
                'score3':cys_count
                }
        
        key = f"{domain['Type']}_{domain['start']}_{domain['end']}"
        scores[key] = info
    
    scores['Total scores'] = total_scores(scores)

    return scores

def write_empty(subfolder):

    scores = {}
    info = {'s1':0,
            's2':0,
            's3':0,
            's4':0, 
            's5':0,
            's6':0
            }

    scores['Total scores'] = info

    name = subfolder.name
    out_name = f'{name}_Susceptibility_Scores.json'
    out_pth = subfolder / out_name

    with open(out_pth, 'w') as file:
        json.dump(scores, file, indent=4)


## ARGUMENT CHECK
if len(sys.argv) != 3:
    print("Use: python ROS_SR.py folder_with_predictions monomers|multimers")
    sys.exit(1)
    
if sys.argv[2] != 'multimers' and sys.argv[2] != 'monomers':
    print(f"Error: unknown prediction type {sys.argv[2]}")
    print("Type of prediction = monomers if the directory contains only monomeric predictions")
    print("Type of prediction = multimers if the directory contains only multimeric predictions")
    sys.exit(1)


# Check whether monomers or multimers are processed
mon = True
mult = False

if sys.argv[2] == 'multimers':
    mon = False
    mult = True

# Check if input directory exists
try:
    folder = Path(sys.argv[1])
except:
    print(f'Error: no folder named {sys.argv[1]}!')

# Susceptible aminoacids
susceptible = ['C', 'M', 'Y', 'W', 'H', 'L', 'R', 'P', 'T']


## MAIN PROGRAM
# Monomers
if mon:
    print('Running program for monomers...\n')
    AS_number = 0
    FILES_number = 0
    for subfolder in folder.iterdir():

        FILES_number += 1

        name = subfolder.name
        pdb_filename = f'{name}.pdb'
        pdb_path = subfolder / pdb_filename
        uniprot_filename = f'{name}_UniProt_Features.json'
        json_path = subfolder / uniprot_filename

        print(f'Protein {name}', file = sys.stderr)

        # Check if files exist in folder
        if not os.path.exists(pdb_path):
            print(f'\nNo {pdb_filename} file found in folder {name}')
            print(f'No {pdb_filename} file found in folder {name}', file = sys.stderr)
            continue

        if not os.path.exists(json_path):
            print(f'\nNo {uniprot_filename} file found in folder {name}')
            print(f'No {uniprot_filename} file found in folder {name}', file = sys.stderr)
            print('Try executing Uniprot_Entries_Retrieval.py')
            continue

        # Read UniProt info
        print(f'\nProcessing {name} ...')
        
        with open(json_path, 'r') as uf:
            uniprot_info = json.load(uf)
        
        # Get domains from UniProt info
        domains = get_domains(uniprot_info)
        if not domains:
            print(f'No Sites found for protein {name}')
            print(f'No Sites found for protein {name}', file=sys.stderr)
            write_empty(subfolder)
            continue

        AS_number += 1

        # Read .pdb
        parser = PDBParser(QUIET=True)
        protein = parser.get_structure("candidate", pdb_path)

        # Get full sequence from .pdb
        chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in protein.get_chains()}
               
        # Retrieve scores dictionary for each domain
        scores = susceptibility_score(uniprot_info, protein, susceptible, domains, chains, 'A')
        
        # Save in a new file
        if any(scores.values()):
            out_filename = f'{name}_Susceptibility_Scores.json'
            out_path = subfolder / out_filename
            
            with open(out_path, 'w') as out_file:
                json.dump(scores, out_file, indent=4)

        else:
            print('No score for protein {name}')

# Multimers
else:
    print('Running program for multimers...\n')
    AS_number = 0
    FILES_number = 0

    for subfolder in folder.iterdir():
        
        FILES_number += 1
        results = {}
        
        ## Process .pdb file
        name = subfolder.name
        pdb_filename = f'{name}.pdb'
        pdb_path = subfolder / pdb_filename
        print(f'\nProcessing {name} ...')
        
        # Check if files exist in folder
        if not os.path.exists(pdb_path):
            print(f'\nNo {pdb_filename} file found in folder {name}')
            print(f'No {pdb_filename} file found in folder {name}', file=sys.stderr)
            continue
        
        # Get complex structure
        parser = PDBParser(QUIET=True)
        protein = parser.get_structure("candidate", pdb_path)

        # Get full sequence from .pdb
        chains = {chain.id:seq1(''.join(residue.resname for residue in chain)) for chain in protein.get_chains()}

        ## Process .json files
        # One uniprot filename for each individual protein
        prots = name.split('-')
        uniprot_filenames = [f'{prot}_UniProt_Features.json' for prot in prots]
        ACT_site = False
        
        for uniprot_filename, prot, (chain_id, seq) in zip(uniprot_filenames, prots, chains.items()):
            json_path = subfolder / uniprot_filename
        
            if not os.path.exists(json_path):
                print(f'\nNo {uniprot_filename} file found in folder {name}')
                print(f'nNo {uniprot_filename} file found in folder {name}', file = sys.stderr)
                print('Try executing Uniprot_Entries_Retrieval.py')
                break

            # Read UniProt file
            with open(json_path, 'r') as uf:
                uniprot_info = json.load(uf)

            # Get domains from UniProt info
            domains = get_domains(uniprot_info)
            if not domains:
                print(f'No Active Site found for protein {prot}')
                print(f'No Active Site found for protein {prot}', file=sys.stderr)
                continue
        
            ACT_site = True
            # Retrieve scores dictionary for each domain
            scores = susceptibility_score(uniprot_info, protein, susceptible, domains, chains, chain_id)
            
            # Store the results for that protein
            results[prot] = scores
        
        if ACT_site:
            AS_number += 1

        # Save results in file
        if any(results.values()):
            out_filename = f'{name}_Susceptibility_Scores.json'
            out_path = subfolder / out_filename
        
            with open(out_path, 'w') as out_file:
                json.dump(results, out_file, indent=4)
                
        else:
            print(f'No results for complex {name}\n')

percentage = AS_number / FILES_number * 100
print(f'\n{AS_number} susceptible sites found in {FILES_number} files ({percentage}%)')
