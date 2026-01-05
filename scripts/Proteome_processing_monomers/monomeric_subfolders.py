"""
monomeric_subfolders.py
----------------------
Given a folder with monomeric predictions, creates a subfolder {pdb_name}/ for each model.

Author: Jorge Marcos Fernández
Date: 2025-09-30
Version: 1.0

Usage:
    python monomeric_subfolders.py 

Output:
    - {pdb_name}/ subfolders for each prediction

Dependencies:
    - shutil
    - pathlib

Notes:
    - Requires folder with monomeric predictions. Name can be changed directly in the code.
"""

import shutil
from pathlib import Path

folder = Path("./Monomers_predictions/")  # <-- Change to folder where monomeric predictions can be found

for pdb_file in folder.glob("*.pdb"):

    new_folder = folder / pdb_file.stem  
    new_folder.mkdir(exist_ok=True)
    new_path = new_folder / pdb_file.name  
    
    shutil.move(str(pdb_file), str(new_path))
    
    print(f"{pdb_file.name} -> {new_folder}/")

