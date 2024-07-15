import pymol2
import csv
import os

def analyze_atp_binding(protein_file, ligand_name="ATP", distance_threshold=5.0, output_csv=None):
    if output_csv is None:
        output_csv = "/Users/dalarios/git/Jazzer_surf/data/interactions.csv"
    
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(protein_file)
        pymol.cmd.select("ligand", f"resn {ligand_name}")
        pymol.cmd.select("binding_residues", f"byres (ligand around {distance_threshold})")
        
        atp_atoms = pymol.cmd.get_model("ligand").atom
        binding_residues = pymol.cmd.get_model("binding_residues").atom

        interactions = []

        for atp_atom in atp_atoms:
            atp_selection = f"chain {atp_atom.chain} and resi {atp_atom.resi} and name {atp_atom.name}"
            for residue_atom in binding_residues:
                residue_selection = f"chain {residue_atom.chain} and resi {residue_atom.resi} and name {residue_atom.name}"
                
                try:
                    distance = pymol.cmd.get_distance(atp_selection, residue_selection)
                    if distance <= distance_threshold:
                        interactions.append((atp_atom.resn, atp_atom.resi, atp_atom.name,
                                             residue_atom.resn, residue_atom.resi, residue_atom.name, distance))
                except pymol.CmdException as e:
                    print(f"Error measuring distance between {atp_selection} and {residue_selection}: {e}")

        # Write interactions to CSV
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        with open(output_csv, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["ATP Residue Name", "ATP Residue ID", "ATP Atom Name", "Binding Residue Name", "Binding Residue ID", "Binding Atom Name", "Distance (Å)"])
            for interaction in interactions:
                writer.writerow(interaction)

        print(f"Interactions have been saved to {output_csv}")

# Example usage
protein_file = "/Users/dalarios/git/Jazzer_surf/3d_predictions/chimeras/fold_2c_motor_only_atp_model_0.pdb"
analyze_atp_binding(protein_file)