import pymol2
import pandas as pd
import os

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
        'S': 'Sulfur', 'CL': 'Chlorine', 'AR': 'Argon', 'K': 'Potassium', 'CA': 'Calcium',
        'SC': 'Scandium', 'TI': 'Titanium', 'V': 'Vanadium', 'CR': 'Chromium', 'MN': 'Manganese',
        'FE': 'Iron', 'CO': 'Cobalt', 'NI': 'Nickel', 'CU': 'Copper', 'ZN': 'Zinc',
        'GA': 'Gallium', 'GE': 'Germanium', 'AS': 'Arsenic', 'SE': 'Selenium', 'BR': 'Bromine',
        'KR': 'Krypton',
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

def analyze_all_interactions(protein_file, distance_threshold=5.0):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(protein_file)
        pymol.cmd.select("all_atoms", "all")
        
        all_atoms = pymol.cmd.get_model("all_atoms").atom

        interactions = []

        for i, atom1 in enumerate(all_atoms):
            atom1_selection = f"chain {atom1.chain} and resi {atom1.resi} and name {atom1.name}"
            for atom2 in all_atoms[i+1:]:
                atom2_selection = f"chain {atom2.chain} and resi {atom2.resi} and name {atom2.name}"
                
                try:
                    distance = pymol.cmd.get_distance(atom1_selection, atom2_selection)
                    if distance <= distance_threshold:
                        interactions.append((atom1.chain, atom1.resn, atom1.resi, atom1.name,
                                             atom2.chain, atom2.resn, atom2.resi, atom2.name, distance))
                except pymol.CmdException as e:
                    print(f"Error measuring distance between {atom1_selection} and {atom2_selection}: {e}")

        return interactions

def process_file(protein_file, distance_threshold=5.0, output_csv=None):
    interactions = analyze_all_interactions(protein_file, distance_threshold)

    # Create a DataFrame from the interactions
    df = pd.DataFrame(interactions, columns=["Chain1", "Residue Name1", "Residue ID1", "Atom Name1", "Chain2", "Residue Name2", "Residue ID2", "Atom Name2", "Distance (Å)"])

    # Add columns for amino acid single-letter code and complete atom names
    df['Residue Single Letter1'] = df['Residue Name1'].apply(three_to_one)
    df['Complete Atom Name1'] = df['Atom Name1'].apply(get_full_atom_name)
    df['Residue Single Letter2'] = df['Residue Name2'].apply(three_to_one)
    df['Complete Atom Name2'] = df['Atom Name2'].apply(get_full_atom_name)

    if output_csv is None:
        output_csv = os.path.splitext(protein_file)[0] + "_interactions.csv"

    # Save the DataFrame to a CSV file
    df.to_csv(output_csv, index=False)

    print(f"DataFrame has been saved to {output_csv}")
    return df

# Example usage
folder_path = "/Users/dalarios/git/Jazzer_surf/3d_predictions/chimeras/test/a_atp_tub.cif"
df = process_file(folder_path)

