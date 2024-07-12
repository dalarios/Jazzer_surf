import pymol2
import csv
import os

# User-specific variables
user = "jazzeruncal"
base_path = f'/Users/{user}/git/Jazzer_surf/data/interactionsMT'
base_prot_path = f"/Users/{user}/git/Jazzer_surf/3d_predictions/MTanalysis/"

# Chain identifiers for kinesin, alpha-tubulin, and beta-tubulin
kinesin_chain = "A"  # Set the correct chain identifier for kinesin
alpha_chain = "C"    # Set the correct chain identifier for alpha-tubulin
beta_chain = "D"     # Set the correct chain identifier for beta-tubulin

def find_contacts(protein_file, chain1, chain2, distance_threshold=4.0, output_csv=None):
    if output_csv is None:
        output_csv = f"{base_path}/contacts.csv"

    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(protein_file)

        # Select residues around the chains of interest within the distance threshold
        pymol.cmd.select("chain1_near_chain2", f"byres (chain {chain2} around {distance_threshold}) and chain {chain1}")
        pymol.cmd.select("chain2_near_chain1", f"byres (chain {chain1} around {distance_threshold}) and chain {chain2}")

        # Get atoms for the selected residues
        chain1_atoms = pymol.cmd.get_model("chain1_near_chain2").atom
        chain2_atoms = pymol.cmd.get_model("chain2_near_chain1").atom

        interactions = []

        # Measure distances
        for atom1 in chain1_atoms:
            atom1_selection = f"id {atom1.index}"
            for atom2 in chain2_atoms:
                atom2_selection = f"id {atom2.index}"
                
                try:
                    distance = pymol.cmd.get_distance(atom1_selection, atom2_selection)
                    if distance <= distance_threshold:
                        interactions.append((atom1.resn, atom1.resi, atom1.name, atom1.chain,
                                             atom2.resn, atom2.resi, atom2.name, atom2.chain, distance))
                except pymol.CmdException as e:
                    print(f"Error measuring distance between {atom1_selection} and {atom2_selection}: {e}")

        # Write interactions to CSV
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        with open(output_csv, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Residue Name 1", "Residue ID 1", "Atom Name 1", "Chain 1",
                             "Residue Name 2", "Residue ID 2", "Atom Name 2", "Chain 2", "Distance (Å)"])
            for interaction in interactions:
                writer.writerow(interaction)

        print(f"Interactions have been saved to {output_csv}")

# Process each CIF file in the directory
for filename in os.listdir(base_prot_path):
    if filename.endswith(".cif"):
        protein_file = os.path.join(base_prot_path, filename)
        kinesin_name = filename.replace(".cif", "")

        # Find contacts between kinesin and alpha-tubulin
        alpha_output_csv = f"{base_path}/{kinesin_name}_kinesin_alpha.csv"
        find_contacts(protein_file, kinesin_chain, alpha_chain, output_csv=alpha_output_csv)

        # Find contacts between kinesin and beta-tubulin
        beta_output_csv = f"{base_path}/{kinesin_name}_kinesin_beta.csv"
        find_contacts(protein_file, kinesin_chain, beta_chain, output_csv=beta_output_csv)
