from pymol import cmd, stored
import pandas as pd

# Dictionaries for single-letter amino acid codes and full atom names
aa_dict = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

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

def three_to_one(three_letter_code):
    return aa_dict.get(three_letter_code.strip().upper(), '?')

def get_full_atom_name(atom_symbol):
    return atom_dict.get(atom_symbol.strip().upper(), 'Unknown')

def log_interactions():
    # Select interacting atoms in chain A or B within 5 Å of any atom in chain C or D or ATP
    cmd.select("interacting_atoms", "(chain A within 5 of chain C) or (chain A within 5 of chain D) or (chain B within 5 of chain C) or (chain B within 5 of chain D) or (chain A within 5 of resn ATP) or (chain B within 5 of resn ATP) or (chain A within 5 of chain B)")
    
    interactions = []

    # Iterate over selected atoms and log interactions
    stored.list = []
    cmd.iterate("interacting_atoms", "stored.list.append((chain, resi, resn, name, index))")

    for chain, resi, resn, name, index in stored.list:
        # Interactions with chains C and D
        cmd.select("interaction_near", f"(chain C within 5 of index {index}) or (chain D within 5 of index {index})")
        stored.interaction_list = []
        cmd.iterate("interaction_near", "stored.interaction_list.append((chain, resi, resn, name, index))")
        for inter_chain, inter_resi, inter_resn, inter_name, inter_index in stored.interaction_list:
            distance = cmd.get_distance(f"index {index}", f"index {inter_index}")
            if distance <= 5.0:  # Only log interactions within 5 Å or less
                interactions.append((chain, resi, resn, name, inter_name, inter_resn, inter_chain, distance))
        
        # Interactions with ATP
        cmd.select("atp_near", f"(resn ATP within 5 of index {index})")
        stored.atp_list = []
        cmd.iterate("atp_near", "stored.atp_list.append((chain, resi, resn, name, index))")
        for atp_chain, atp_resi, atp_resn, atp_name, atp_index in stored.atp_list:
            distance = cmd.get_distance(f"index {index}", f"index {atp_index}")
            if distance <= 5.0:  # Only log interactions within 5 Å or less
                interactions.append((chain, resi, resn, name, atp_name, atp_resn, atp_chain, distance))

        # Interactions between chain A and chain B
        cmd.select("ab_near", f"(chain B within 5 of index {index})" if chain == 'A' else f"(chain A within 5 of index {index})")
        stored.ab_list = []
        cmd.iterate("ab_near", "stored.ab_list.append((chain, resi, resn, name, index))")
        for ab_chain, ab_resi, ab_resn, ab_name, ab_index in stored.ab_list:
            distance = cmd.get_distance(f"index {index}", f"index {ab_index}")
            if distance <= 5.0:  # Only log interactions within 5 Å or less
                interactions.append((chain, resi, resn, name, ab_name, ab_resn, ab_chain, distance))

    # Create a DataFrame from the interactions
    df = pd.DataFrame(interactions, columns=["chain", "resi", "resn", "atom_name", "interacting_atom", "interacting_resn", "interacting_chain", "distance (angstroms)"])

    # Add columns for amino acid single-letter code and complete atom names
    df['residue_one_letter'] = df['resn'].apply(three_to_one)
    df['full_atom_name'] = df['atom_name'].apply(get_full_atom_name)
    df['interacting_full_atom_name'] = df['interacting_atom'].apply(get_full_atom_name)

    # Save the DataFrame to a CSV file
    df.to_csv("interactions.csv", index=False)
    print(f"DataFrame has been saved to interactions.csv")

# Run the logging function
log_interactions()
