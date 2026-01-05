"""
Name_transformer.py
----------------------
Changes ranked_0.pdb predictions filenames to {PP_XXXX}.pdb 

Author: Jorge Marcos Fernández
Date: 2025-10-03
Version: 1.0

Usage:
    python Name_transformer target_directory

Output:
    - ranked_0.pdb --> PP_XXXX.pdb or pWW0_XXXX.pdb

Dependencies:
    - sys
    - os
    - pathlib

Notes:
    - Requires target_folder to be divided into subfolders for each predicion
    - PP_XXXX or pWW0_XXXX identifiers are stored in each subfolder's name
"""

# Change ranked_0.pdb filenames to PP_XXXX identifiers stored in folder's name
# ranked_0.pdb --> PP_XXXX.pdb
#

# PACKAGES
from pathlib import Path
import sys
import os

# Check arguments
if len(sys.argv) != 2:
    print("Use: python Name_transformer.py folder_with_ranked_0_predictions")
    sys.exit(1)

# Process folder
folder = Path(sys.argv[1])
target_filename = 'ranked_0.pdb'

for subfolder in folder.iterdir():
    
    name = subfolder.name
    pdb_path = subfolder / target_filename
    
    # Check if ranked_0 file exists
    if not os.path.exists(pdb_path):
        print(f'No {target_filename} file found in folder {subfolder.name}')
        continue

    # New name is retrieved from the subfolder name   
    new_name = f'{subfolder.name}.pdb' 
    new_path = subfolder / new_name
    
    # Check whether the ranked_0 file has already been processed
    if new_path.exists():
        print(f'File {new_name} already exists in {subfolder.name}, skipping...')
        continue
        
    # Rename file
    pdb_path.rename(new_path)
    print(f"{pdb_path.name} -> {new_name}")

