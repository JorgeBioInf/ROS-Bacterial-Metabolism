"""
AF_BestRanked_and_ipTM.py
----------------------
Retrieves the best structural model from a given AlphaFold2 output and deletes unnecessary files

Author: Jorge Marcos Fernández
Date: 2025-09-26
Version: 1.0

Usage:
    python ipTM.py folder_with_predictions monomers|multimers 

Output:
    ranked_0.pdb 

Dependencies:
    - pathlib
    . shutil
    - sys
    - json
    - os
    - pickle

Notes:
    - ipTM score is also stored for complexes
"""


# MODULES
import json
import pickle
import os
from pathlib import Path
import sys
import shutil


### MAIN PROGRAM

# Check arguments
# A target folder, output folder and type of prediction should be entered

if len(sys.argv) < 2:
    print("Use: python ipTM.py folder_with_predictions type_of_prediction")
    print("Type of prediction = monomers if the directory contains only monomeric predictions")
    print("Type of prediction = multimers if the directory contains only multimeric predictinos")
    sys.exit(1)

if sys.argv[2] != 'multimer' and sys.argv[2] != 'monomer':
    print(f"Error: unknown prediction type {sys.argv[2]}")
    print("Type of prediction = monomers if the directory contains only monomeric predictions")
    print("Type of prediction = multimers if the directory contains only multimeric predictinos")
    sys.exit(1)


# Control variables to distinguish between monomers and multimers
mult = False
if sys.argv[2] == 'multimer':
    mult = True

    
# Initialize variables 
folder = Path(sys.argv[1])
target_filename = 'ranking_debug.json'

   
## MONOMERS & COMPLEXES
# Restart control variables
ipTM = False
r0 = False
keep_files = [f'Error_{folder.name}.txt'] # Keep error log file

print(f'Processing prediction {folder.name}')

# Retrieve best model if exists
ranked0 = 'ranked_0.pdb'
ranked0_path = folder / ranked0

if os.path.exists(ranked0_path):
    r0 = True
    keep_files.append(ranked0)

## COMPLEXES  
if mult == True:
    
    # Get the name of the best model's .pkl file 
    file_path = folder / target_filename

    if not os.path.exists(file_path):
        print(f'Error: unable to find file {target_filename} for protein {folder.name}')
        print('Please, check Error file\n')

        # Delete all files except error log
        keep_files.pop(ranked0)
        for item in folder.iterdir():
            if item.name not in keep_files:
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
        exit
    
    with open(file_path, 'r') as json_data:
        d = json.load(json_data)
        best_model_name = d["order"][0]

    bm_filename = f'result_{best_model_name}.pkl'
    file_path2 = folder / bm_filename

    # Get ipTM score 
    with open(file_path2, 'rb') as f:
        data = pickle.load(f)
        pred_ipTM = data['iptm']

        if pred_ipTM is not None: 
            ipTM = True
            filename_ipTM = f'{folder.name}_ipTM.txt'
            filepath_ipTM = folder / filename_ipTM

            with open(filepath_ipTM, 'w') as file:
                file.write(str(pred_ipTM))
    
    if ipTM == True: 
        keep_files.append(filename_ipTM)

# If ranked_0 found for monomers or ranked_0 and ipTM found for complexes, delete rest of the files
if mult == True:
    if r0 == True and ipTM == True:
        print(f'All information found for model {folder.name}')
        print(f'Deleting unnecessary files')
        print('')
        
        # Delete files
        for item in folder.iterdir():
            if item.name not in keep_files:
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
        
    elif r0 == False:
        print(f'Error when searching for {ranked0} file in {folder.name} prediction')
        print(f'Deleting files')
        print('')

        for item in folder.iterdir():
            if item.name not in keep_files:
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
    
    elif ipTM == False:
        print(f'Error when searching for ipTM value in {folder.name} prediction')
        print('Deleting files')
        print('')

        for item in folder.iterdir():
            if item.name not in keep_files:
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
                
else:
    if r0 == True:
        print(f'All information found for model {folder.name}')
        print(f'Deleting unnecessary files')
        print('')
            
        # Delete files
        for item in folder.iterdir():
            if item.name not in keep_files:
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
                
    else:
        print(f'Error when searching for {ranked0} file in {folder.name} prediction')
        print('Deleting files')
        print('')

        for item in folder.iterdir():
            if item.name not in keep_files:
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)



