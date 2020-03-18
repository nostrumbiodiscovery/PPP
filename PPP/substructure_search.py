import numpy as np

CONSTRAINT = '{{ "type": "constrainAtomToPosition", "springConstant": 10, "equilibriumDistance": 0.0, "constrainThisAtom": "{}:{}:{}" }},'

def search_smart_pattern(ligand, smiles):
    from rdkit import Chem
    # Get smiles pattern
    patt = Chem.MolFromSmarts(smiles)
    Chem.Kekulize(patt)
    # Get ligand
    m = Chem.MolFromPDBFile(ligand)
    Chem.Kekulize(m)
    # Substructure search
    atom_idxs = np.array(m.GetSubstructMatch(patt))
    try:
        atoms_selected = np.array(m.GetAtoms())[atom_idxs]
    except IndexError:
        raise ValueError("No substructure was recognised with {}  as the smile pattern and {} as the molecule".format(
        smiles, Chem.MolToSmiles(m)))
    # Retrieve atom names and residue names of the common atoms
    atom_names = [atom.GetMonomerInfo().GetName().replace(" ", "_") for atom in atoms_selected]
    chain_names = [atom.GetPDBResidueInfo().GetChainId() for atom in atoms_selected]
    residue_numbers = [atom.GetPDBResidueInfo().GetResidueNumber() for atom in atoms_selected]
    # Get coords
    conformer = m.GetConformer()
    coords = np.array(conformer.GetPositions())[atom_idxs]
    return atom_names, chain_names, coords, residue_numbers

def extract_constraint_from_smiles(ligand, smiles):
    from rdkit import Chem
    m = Chem.MolFromPDBFile(ligand)
    atom_names, chain_names, coords, residue_numbers = search_smart_pattern(ligand, smiles)
    constants = [CONSTRAINT.format(chain_name, resname, atom_name) for atom_name, chain_name, coord, resname in zip(atom_names, chain_names, coords, residue_numbers)]
    return constants



if __name__ == "__main__":
    ligand = "ligand.pdb"
    smiles = "CCC([O-])O"
    constraints = extract_constraint_from_smiles(ligand, smiles)
