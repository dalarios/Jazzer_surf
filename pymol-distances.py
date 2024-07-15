import pymol2
import csv
import os
import glob
import pandas as pd
import numpy as np
import cupy as cp
from concurrent.futures import ThreadPoolExecutor

# Function to convert three-letter amino acid code to single-letter code
def three_to_one(three_letter_code):
    aa_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_dict.get(three_letter_code.upper(), '?')

# Function to get complete atom names for common organic molecules and periodic elements from H to Kr
def get_full_atom_name(atom_symbol):
    atom_dict = {
        'H': 'Hydrogen', 'HE': 'Helium', 'LI': 'Lithium', 'BE': 'Beryllium', 'B': 'Boron',
        'C': 'Carbon', 'N': 'Nitrogen', 'O': 'Oxygen', 'F': 'Fluorine', 'NE': 'Neon',
        'NA': 'Sodium', 'MG': 'Magnesium', 'AL': 'Aluminium', 'SI': 'Silicon', 'P': 'Phosphorus',
        'S': 'Sulfur', 'CL': 'Chlorine', 'K': 'Potassium', 'CA': 'Calcium', 'MN': 'Manganese',
        'FE': 'Iron', 'CO': 'Cobalt', 'NI': 'Nickel', 'CU': 'Copper', 'ZN': 'Zinc',
        # Specific atoms in amino acids
        'CA': 'Carbon alpha', 'CB': 'Carbon beta', 'CG': 'Carbon gamma', 'CD': 'Carbon delta', 
        'CE': 'Carbon epsilon', 'CZ': 'Carbon zeta', 'ND1': 'Nitrogen delta 1', 'NE2': 'Nitrogen epsilon 2', 
        'CD2': 'Carbon delta 2', 'CE1': 'Carbon epsilon 1', 'CG2': 'Carbon gamma 2',
        'OG': 'Oxygen gamma', 'OG1': 'Oxygen gamma 1', 'OD1': 'Oxygen delta 1', 'OD2': 'Oxygen delta 2', 
        'OE1': 'Oxygen epsilon 1', 'OE2': 'Oxygen epsilon 2', 'ND2': 'Nitrogen delta 2',
        'NE': 'Nitrogen epsilon', 'NZ': 'Nitrogen zeta', 'SG': 'Sulfur gamma', 'CE2': 'Carbon epsilon 2', 'CD1': 'Carbon delta 1',
        'NH1': 'Nitrogen eta 1', 'NH2': 'Nitrogen eta 2',
        # Specific atoms in ATP
        'PA': 'Phosphorus alpha', 'PB': 'Phosphorus beta', 'PG': 'Phosphorus gamma',
        'O1A': 'Oxygen 1 alpha', 'O2A': 'Oxygen 2 alpha', 'O3A': 'Oxygen 3 alpha',
        'O1B': 'Oxygen 1 beta', 'O2B': 'Oxygen 2 beta', 'O3B': 'Oxygen 3 beta',
        'O1G': 'Oxygen 1 gamma', 'O2G': 'Oxygen 2 gamma', 'O3G': 'Oxygen 3 gamma',
        'N1': 'Nitrogen 1', 'N3': 'Nitrogen 3', 'N6': 'Nitrogen 6', 'N7': 'Nitrogen 7', 'N9': 'Nitrogen 9',
        'C1': 'Carbon 1', 'C2': 'Carbon 2', 'C3': 'Carbon 3', 'C4': 'Carbon 4', 'C5': 'Carbon 5',
        'C6': 'Carbon 6', 'C7': 'Carbon 7', 'C8': 'Carbon 8', 'C9': 'Carbon 9',
        # Prime atoms (commonly in nucleotides and carbohydrates)
        "C1'": "Carbon 1 prime", "C2'": "Carbon 2 prime", "C3'": "Carbon 3 prime",
        "C4'": "Carbon 4 prime", "C5'": "Carbon 5 prime", "O2'": "Oxygen 2 prime",
        "O3'": "Oxygen 3 prime", "O4'": "Oxygen 4 prime", "O5'": "Oxygen 5 prime"
    }
    return atom_dict.get(atom_symbol.upper(), 'Unknown')

def analyze_interchain_interactions(protein_file, distance_threshold=10):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(protein_file)
        pymol.cmd.select("all_atoms", "all")
        
        all_atoms = pymol.cmd.get_model("all_atoms").atom

        atom_coords = np.array([[atom.coord[0], atom.coord[1], atom.coord[2]] for atom in all_atoms])
        atom_coords_gpu = cp.array(atom_coords)
        interactions = []

        for i, atom1 in enumerate(all_atoms):
            if i % 100 == 0:  # Reduce output frequency
                print(f"Processing atom {i}/{len(all_atoms)}")
            if atom1.chain != '':  # Ensure valid chain identifier
                distances = cp.linalg.norm(atom_coords_gpu - atom_coords_gpu[i], axis=1)
                indices = cp.where((distances > 0) & (distances <= distance_threshold))[0]

                for j in indices.get():
                    if j < len(all_atoms):
                        atom2 = all_atoms[j]
                        if atom1.chain != atom2.chain:
                            distance = distances[j].get()
                            interactions.append((atom1.chain, atom1.resn, atom1.resi, atom1.name,
                                                 atom2.chain, atom2.resn, atom2.resi, atom2.name, distance))
        return interactions

def process_files_in_folder(folder_path, distance_threshold, output_csv=None):
    all_interactions = []

    def process_file(filename):
        print(f"Processing file: {filename}")
        interactions = analyze_interchain_interactions(filename, distance_threshold=distance_threshold)
        return [(os.path.splitext(os.path.basename(filename))[0],) + interaction for interaction in interactions]

    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor() as executor:
        results = executor.map(process_file, glob.glob(os.path.join(folder_path, "*.pdb")) + glob.glob(os.path.join(folder_path, "*.cif")))

    for result in results:
        all_interactions.extend(result)

    # Create a DataFrame from the interactions
    df = pd.DataFrame(all_interactions, columns=["Protein File", "Chain 1", "Residue Name 1", "Residue ID 1", "Atom Name 1", "Chain 2", "Residue Name 2", "Residue ID 2", "Atom Name 2", "Distance (Å)"])

    # Add columns for amino acid single-letter code and complete atom names
    df['Residue 1 Single Letter'] = df['Residue Name 1'].apply(three_to_one)
    df['Residue 2 Single Letter'] = df['Residue Name 2'].apply(three_to_one)
    df['Atom 1 Complete Name'] = df['Atom Name 1'].apply(get_full_atom_name)
    df['Atom 2 Complete Name'] = df['Atom Name 2'].apply(get_full_atom_name)

    if output_csv is None:
        output_csv = os.path.join(folder_path, "interactions_all_chains.csv")

    # Save the DataFrame to a CSV file
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    df.to_csv(output_csv, index=False)

    print(f"DataFrame has been saved to {output_csv}")
    return df

# Example usage
folder_path = r"C:\Users\Admin\git\Jazzer_surf\3d_predictions\chimeras\test"  # Ensure this path is correct
output_csv = r"C:\Users\Admin\git\Jazzer_surf\3d_predictions\chimeras\interactions_all_chains.csv"  # Use a directory where you have write permissions
df = process_files_in_folder(folder_path, distance_threshold=5.0, output_csv=output_csv)
