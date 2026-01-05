"""
PP_metab_proc.py
----------------------
Processes a metabolic model of Pseudomonas putida and identifies all monomers and complexes.
Monomers --> Returns UniProt ID and structural prediction from Uniprot when available. If not, returns FASTA sequence for further modelling.
Complexes --> Returns multi-FASTA sequence for structure prediction

Author: Jorge Marcos Fernández
Date: 2025-09-15
Version: 1.0

Usage:
    python PP_metab_proc.py 

Output:
    Monomers_predictions/ : folder with all monomeric .pdb structures (if available)
    Monomers_to_model/ : folder with FASTA files for remaining monomers
    Complex_fastas/ : folder with multi-FASTA sequences of identified complexes

Dependencies:
    - pandas
    - re
    - numpy
    - sys
    - itertools
    - json
    - os
    - requests

Notes:
    - Requires model iJN1480.txt in the same directory. Model name might be changed in code. 
    - ID_relationships.json dictionary in the same directory accelerates analysis. If not available, it is created.
    - Might trigger errors due to API server saturation. Consider increasing the waiting time between queries. 
"""


# MODULES
import pandas as pd # type: ignore
import numpy as np # type: ignore
import re
import os
import requests # type: ignore
from itertools import product
import time
import json

# ======================================================================

# SCRIPT PREPARATION

# Subdirectories names and paths
direct = 'Complex_fastas'
path = os.path.join(os.getcwd(), direct)

direct2 = 'Monomers_predictions'
path2 = os.path.join(os.getcwd(), direct2)

direct3 = 'Monomers_to_model'
path3 = os.path.join(os.getcwd(), direct3)

# Create subdirectories if they do not exist
if not os.path.exists(path): 
            os.makedirs(path)

if not os.path.exists(path2): 
            os.makedirs(path2)

if not os.path.exists(path3): 
    os.makedirs(path3)

# Load dictionary with relation PP_XXXX - UniProt ID if exists. If not, create it
id_rel_filename = 'ID_relationships.json'

if os.path.exists(id_rel_filename):
    print(f'KEGG_ID - UniProt ID relationship retrieved from {id_rel_filename}')
    with open('ID_relationships.json', 'r') as file:   
        KEGG2UniProt = json.load(file)
else:
    print(f'No {id_rel_filename} found, generating new dictionary')
    KEGG2UniProt = {}

# Empty data structures 
KEGG2Seq = {}    # PP_XXXX - sequence
Not_found = []  # Will contain not found IDs

# =========================================================================

### PART 1. Access to protein's PP_XXXX and pWW0_XXX identifiers

# Load dataset
model_name = 'iJN1480.txt' # Change if needed
data = pd.read_csv(model_name, sep="\t")

# Get Gen-reaction associations
genes = data['GPR'].dropna().unique()

# Store IDs for each monomer (PP_XXXX or pWW0_XXX)
monomers = []

for ID in genes: 
    
    if not 'and' in ID:  # Monomer condition
        splitted = ID.split(' or ')
        r = re.compile('PP_\d{4}|pWW0_\d+') # Find protein's Uniprot ID
        newlist = list(filter(r.match, splitted))
        
        for mon in newlist:
            mon = mon.strip().replace('(','').replace(')','')
            if mon not in monomers: # Avoid repetitions
                monomers.append(mon)


# Store IDs for each complex protein (PP_XXXX or pWW0_XXX)
complexes = []

for ID in genes:    
    
    if 'and' in ID: # Complex condition
        
        # Opt 1: Only 1 possible complex (A and B)
        if not 'or' in ID:
            r = re.findall('PP_\d{4}|pWW0_\d+', ID)
            complexes.append(r)
        
        # Opt 2: Different possible monomers for the same complex ((A or B) and C)
        # All possible combinations must be stored
        elif ' or ' in ID and not ') or' in ID and not 'or (' in ID:
            perma_prot = []
            variable_prot = []
            
            split1 = ID.split(' and ')
            
            for element in split1:
                
                # Fixed proteins (followed or preceded by 'and')
                if not 'or' in element: 
                    r = re.findall('PP_\d{4}|pWW0_\d+', element)
                    perma_prot.append(r)
                
                # Variable proteins (followed or preceded by 'or')
                else:  
                    r = re.findall('PP_\d{4}|pWW0_\d+', element) # Store a list for each option
                    variable_prot.append(r)
                
                # Flatten perma_prot (list of lists to list) 
                perma_flat = [item for sublist in perma_prot for item in sublist]
            
            # Merge all possible combinations (cartesian product of lists)
            complexes.extend([list(t) + perma_flat for t in product(*variable_prot)])
            
            
        # Opt 3: Distincition between different possible complexes ((A and B) or (C and D))
        elif ') or (' in ID:
            split1 = ID.split(') or (')
            
            # If no more possible combinations found (A and B) or (C and D))
            if all("or" not in x for x in split1):
                for element in split1:
                    r = re.findall('PP_\d{4}|pWW0_\d+', element) 
                    complexes.append(r)
            
            # If more possible combinations ((A and B) or (C and (E or D))
            else:           
                for element in split1:
                    
                    if "or" in element and 'and' not in element:  # (A or B)
                        r = re.findall('PP_\d{4}|pWW0_\d+', element)
                        variable_prot.append(r)

                    elif "or" in element and 'and' in element:   # (A or B) and C
                        split2 = element.split('and')
                        perma_prot = []
                        variable_prot = []
                        
                        for element2 in split2:
                
                            # Fixed proteins (followed or preceded by 'and')
                            if not 'or' in element2: 
                                r = re.findall('PP_\d{4}|pWW0_\d+', element2)
                                perma_prot.append(r)

                            # Variable proteins (followed or preceded by 'or')
                            else:  
                                r = re.findall('PP_\d{4}|pWW0_\d+', element2) # Store a list for each option
                                variable_prot.append(r)
                        
                        # Flatten perma_prot (list of lists to list) 
                        perma_flat = [item for sublist in perma_prot for item in sublist]
                        
                        # Cartesian product
                        complexes.extend([list(t) + perma_flat for t in product(*variable_prot)])
            
                    else:  # A and B
                        r = re.findall('PP_\d{4}|pWW0_\d+', element)
                        complexes.append(r)

            
        # Opt 4: Many distinctions ((A and B) or C or (D and (E or F)) / (A and B) or C)
        else:
            perma_prot = []
            variable_prot = []
            
            # Split when ') or' OR 'or (' are found --> Subdivide the problem
            split1 = re.split(r'\)\s*or\s*|\s*or\s*\(', ID)   
            
            for element in split1:
                
                # If element is 'A' add it directly into monomers list'
                if not 'or' in element and not 'and' in element:
                    r = re.findall('PP_\d{4}|pWW0_\d+', element)
                    if r not in monomers:
                        monomers.extend(r)
                
                # If element is a complex (A and B and C ...) add it into complexes list
                elif not 'or' in element:
                    r = re.findall('PP_\d{4}|pWW0_\d+', element)
                    complexes.append(r)
                
                # If (A or B)
                elif 'or' in element and not 'and' in element:
                     r = re.findall('PP_\d{4}|pWW0_\d+', element)
                     monomers.extend([mon for mon in r if mon not in monomers])
                
                # If more distinctions are included (e.g. (A and B) or C) --> Same procedure as before
                else:
                    split2 = element.split('or')
                    
                    for element2 in split2:
                        perma_prot = []
                        variable_prot = []
                        
                        # Fixed proteins (followed or preceded by 'and')
                        if not 'or' in element2: 
                            r = re.findall('PP_\d{4}|pWW0_\d+', element2)
                            perma_prot.append(r)

                        # Variable proteins (followed or preceded by 'or')
                        else:  
                            r = re.findall('PP_\d{4}|pWW0_\d+', element2) # Store a list for each option
                            variable_prot.append(r)
                        
                        # Flatten perma_prot (list of lists to list) 
                        perma_flat = [item for sublist in perma_prot for item in sublist]
                        
                        # Cartesian product
                        complexes.extend([list(t) + perma_flat for t in product(*variable_prot)])
                        

## Data purification

# 1. Delete empty lists (re delivers an empty list [] if no match is found)
complexes = [c for c in complexes if c]

# 2. Repeated and permuted complexes are eliminated
unique_complexes = []
for c in complexes:
    s = tuple(sorted(c))  # Sort complex so permutation variability is avoided
    if s not in unique_complexes:
        unique_complexes.append(s)


## Data flux control
print('Total entries:', len(genes))
print('Total monomers:', len(monomers))
print('Total possible complexes', len(unique_complexes))


# ===========================================================================================

### PART 2. Uniprot ID and sequences retrieval

# Generate plasmid dictionary with the relation pWW0_pXXX identifier - RefSeq identifier
print('\nExtracting plasmid dictionary...')

url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
params = {
    "db": "nuccore",
    "id": "NC_003350.1",  # accession of pWW0 plasmid in genbank
    "rettype": "gb",
    "retmode": "text"
}

resp = requests.get(url, params=params)
gb_record = resp.text

# Dictionary with relationship pWWO identifier - RefSeq identifier
plasmid_dic = {}

# Divide text for each CDS
gb_list = gb_record.split('CDS')

# For each CDS
for element in gb_list:
    m_locus = re.search(r'/old_locus_tag="([^"]+)"', element)  # pWWO identifier
    m_prot  = re.search(r'/protein_id="([^"]+)"', element) # RefSeq identifier

    if m_locus and m_prot: # Check if not empty
        plasmid_dic[m_locus.group(1)] = m_prot.group(1)   # Avoid duplicates

print('Done\n')


## 2.1. MONOMERS

# Already known structures: monomers
Uni_ids = []
used_m = []

print('Working on monomers...')

# Sequence and UniProt ID retrieval
for ID in list(monomers):

    # Check whether the monomer has already been processed
    filename = f"{ID}.fasta"
    filepath = os.path.join(direct3, filename)

    if os.path.exists(filepath):
        print(f"File {filename} already exists!")
        monomers.remove(ID)
        continue
    
    # PLASMID GENES (pWW0_XXX)
    if ID.startswith('pWW0'):

        # Change format so it can be used to parse the databases (pWW0_XXX --> pWWO_pXXX)
        ID_org = ID
        ID = ID.replace('0_', 'O_p')

        # Check if ID exists in NCBI plasmid pWW0 database:
        if ID not in plasmid_dic.keys():
            print(f'Error: protein {ID} not found in NCBI plasmid pWW0 database. Skipping...')
            Not_found.append(f'Error: protein {ID} not found in NCBI plasmid pWW0 database')
            continue

        # Retrieve UniProt information (UniProt ID and Sequence)
        RefSeq = plasmid_dic[ID]
        url = f'https://rest.uniprot.org/uniprotkb/search?query={RefSeq}&format=json'
        response = requests.get(url)
        text = response.json()

        if "results" not in text or len(text["results"]) == 0:  # Check access
            print(f'There was an error searching for protein {ID} in UniProt')
            Not_found.append(f'There was an error searching for protein {ID} in UniProt')
            continue

        s = text['results'][0]['sequence']['value']
        u = text['results'][0]['primaryAccession']

        KEGG2UniProt[ID_org] = u
        KEGG2Seq[ID_org] = s


    # NORMAL GENES
    else:
        # Check if not found in KEGG2UniProt dictionary or if empty entry
        if ID not in KEGG2UniProt or not KEGG2UniProt[ID]:  

            url = f'https://rest.kegg.jp/get/ppu:{ID}'
            response = requests.get(url)

            if not response.ok:  # Check access
                print(f'There was an error searching for protein {ID} in KEGG database')
                Not_found.append(f'There was an error searching for protein {ID} in KEGG database')
                continue
        
            info = response.text

            # Retrieve UniProt ID 
            r = re.findall(r'UniProt:\s+(\S+)', info)
            
            # Check if UniProt ID is found
            if not r:
                print(f'No UniProt ID was found for protein {ID}')
                Not_found.append(f'No UniProt ID was found for protein {ID}')
                continue
            
            Uni_id = r[0]  
            Uni_ids.append(Uni_id)

            KEGG2UniProt[ID] = Uni_id
            used_m.append(ID)

            # Retrieve sequence
            subseqs = re.findall('\\s{8}[A-Z]+$', info, re.M) 

            if not subseqs:
                print(f'No sequence was found for protein {ID} in Uniprot')
                Not_found.append(f'No sequence was found for protein {ID} in Uniprot')
                continue

            s = ''.join(subseqs).replace(' ', '')  # Avoid white spaces
            KEGG2Seq[ID] = s
        

        else:  # Retrieve sequence directly
            Uni_id = KEGG2UniProt[ID]
            url = f"https://rest.uniprot.org/uniprotkb/{Uni_id}.fasta"
            response = requests.get(url)

            if not response.ok:
                print(f'No UniProt fasta was found for protein {ID}')
                Not_found.append(f'No UniProt ID was found for protein {ID}')
                continue

            fasta = response.text
            s = "".join(fasta.splitlines()[1:]) # Retrieves sequence only
            KEGG2Seq[ID] = s

print('Completed\n')


## 2.2. COMPLEXES

complexes = unique_complexes
new_complexes = []

print('Working on complexes...')

for com in complexes:

    # Check whether the complex has already been processed
    name = '-'.join(com)
    filename = f"{name}.fasta"
    filepath = os.path.join(direct, filename)

    if os.path.exists(filepath):
        print(f"File {filename} already exists!")  
        continue
    
    # New complex is being processed
    if not os.path.exists(filepath):
        new_complexes.append(com)

    # For each protein in each complex, get its UniProt ID and sequence
    for ID in com: 

        # NEW PLASMID GENES (pWW0_XXX)
        if ID.startswith('pWW0'):
            # Check if already processed (Uniprot ID + Sequence)
            if ID in KEGG2Seq and ID in KEGG2UniProt and KEGG2Seq[ID] and KEGG2UniProt[ID]:
                continue

            # If no Uniprot Accession
            if ID not in KEGG2UniProt or not KEGG2UniProt[ID]:

            # Change format so it can be used to parse the databases (pWW0_XXX --> pWWO_pXXX)
                ID_org = ID
                ID = ID.replace('0_', 'O_p')

            # Check if ID exists in NCBI plasmid pWW0 database:
                if ID not in plasmid_dic.keys():
                    print(f'Error: protein {ID} not found in NCBI plasmid pWW0 database. Skipping...')
                    Not_found.append(f'Error: protein {ID} not found in NCBI plasmid pWW0 database')
                    continue

                # Retrieve UniProt information (UniProt ID and Sequence)
                RefSeq = plasmid_dic[ID]
                url = f'https://rest.uniprot.org/uniprotkb/search?query={RefSeq}&format=json'
                response = requests.get(url)
                text = response.json()

                if "results" not in text or len(text["results"]) == 0:  # Check access
                    print(f'There was an error searching for protein {ID} in UniProt')
                    Not_found.append(f'There was an error searching for protein {ID} in UniProt')
                    continue

                s = text['results'][0]['sequence']['value']
                u = text['results'][0]['primaryAccession']

                KEGG2UniProt[ID_org] = u
                KEGG2Seq[ID_org] = s
            
            # If UniProt ID but no sequence --> Retrieve sequence directly
            else:
                Uni_id = KEGG2UniProt[ID]
                url = f"https://rest.uniprot.org/uniprotkb/{Uni_id}.fasta"
                response = requests.get(url)

                if not response.ok:
                    print(f'No UniProt fasta was found for protein {ID}')
                    Not_found.append(f'No UniProt fasta was found for protein {ID}')
                    continue

                fasta = response.text
                s = "".join(fasta.splitlines()[1:])  # Retrieves sequence only
                KEGG2Seq[ID] = s

                
    # NORMAL GENES
    if not ID.startswith('pWW0'):

        # Check if already processed (UniProt ID + Sequence)
        if ID in KEGG2Seq and ID in KEGG2UniProt and KEGG2Seq[ID] and KEGG2UniProt[ID]:
                continue
        
        if ID not in KEGG2UniProt or not KEGG2UniProt[ID]:  
            url = f'https://rest.kegg.jp/get/ppu:{ID}'
            response = requests.get(url)

            if not response.ok:  # Check access
                print(f'There was an error searching for protein {ID} in KEGG database')
                Not_found.append(f'There was an error searching for protein {ID} in KEGG database')
                continue
                
            info = response.text

            # Retrieve UniProt ID 
            r = re.findall(r'UniProt:\s+(\S+)', info)
            
            # Check if UniProt ID is found
            if not r:
                print(f'No UniProt ID was found for protein {ID}')
                Not_found.append(f'No UniProt ID was found for protein {ID}')
                continue
            
            Uni_id = r[0]  
            Uni_ids.append(Uni_id)

            KEGG2UniProt[ID] = Uni_id
            used_m.append(ID)

            # Retrieve sequence
            subseqs = re.findall(r'\s{8}[A-Z]+$', info, re.M) 

            if not subseqs:
                print(f'No sequence was found for protein {ID} in UniProt')
                Not_found.append(f'No sequence was found for protein {ID} in UniProt')
                continue

            s = ''.join(subseqs).replace(' ', '')  # Avoid white spaces
            KEGG2Seq[ID] = s

        else:
            # Retrieve sequence directly
            Uni_id = KEGG2UniProt[ID]
            url = f"https://rest.uniprot.org/uniprotkb/{Uni_id}.fasta"
            response = requests.get(url)

            if not response.ok:
                print(f'No UniProt fasta was found for protein {ID}')
                Not_found.append(f'No UniProt fasta was found for protein {ID}')
                continue

            fasta = response.text
            s = "".join(fasta.splitlines()[1:])  # Retrieves sequence only
            KEGG2Seq[ID] = s

# Remove all complexes which were already processed
complexes = new_complexes

print(complexes)
print(unique_complexes)
print('Completed\n')


# ===========================================================================================
# 3. Generate .fasta files for complexes
print('\nGenerating complex fasta files ...')

for c in complexes:
    write = True
    full_c = '-'.join(c)
    ruta = os.path.join(path, f'{full_c}.fasta')
    
    for ID in c:
        if ID not in KEGG2UniProt or ID not in KEGG2Seq:
            print(f'{full_c}.fasta file could not be written')
            print(f'No information for protein {ID} found')
            Not_found.append(f'No information for protein {ID} found')
            write = False
            break
    
    if write:
        with open(ruta, 'w') as file:
            for ID in c:            	
                file.write(f'>{KEGG2UniProt[ID]}\n')
                file.write(KEGG2Seq[ID])
                file.write('\n')

print('Completed\n')


# ===========================================================================================
# 4. Generate monomer PDB files
print('\nSearching for AlphaFold models ...')

for mon in monomers:

    filename = f"{mon}.pdb"
    filepath = os.path.join(path2, mon, filename)

    if os.path.exists(filepath):
        print(f"File {filename} already exists!")
        continue
    
    if mon not in KEGG2UniProt:
        continue
    
    Uni_id = KEGG2UniProt[mon]
    url = f'https://alphafold.ebi.ac.uk/files/AF-{Uni_id}-F1-model_v4.pdb'
    response = requests.get(url)

    if not response.ok:
        print(f'There was an error searching for AlphaFold model {Uni_id} ({mon})')
        Not_found.append(f'There was an error searching for AlphaFold model ({mon})')

        # Save .fasta file instead
        if mon in KEGG2Seq:
            ruta = os.path.join(path3, f'{mon}.fasta')
            with open(ruta, 'w') as file:
                file.write(f'>{Uni_id}\n')
                file.write(KEGG2Seq[mon])
                file.write('\n')
        else:
            print(f'No sequence found for protein {mon}')
            Not_found.append(f'No sequence found for protein {mon}')
        continue

    # Save PDB file
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w') as file:
        file.write(response.text)
    
print('Completed\n')


# ===========================================================================================
# 5. Extra files
print('Generating extra files ...')

# Error file
err_filename = 'err_MetabProc.txt'
with open(err_filename, 'w') as efile:
    for error in Not_found:
        efile.write(error + '\n')

# Writing dictionary
with open(id_rel_filename, 'w') as file:
    json.dump(KEGG2UniProt, file, indent=4)

print('Program finished successfully\n')
input('Press any key to close')
