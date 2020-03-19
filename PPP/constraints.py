#!/bin/python
import sys
import os
import argparse
from collections import defaultdict
import PPP.substructure_search as ss

AMINOACIDS = ["VAL", "ASN", "GLY", "LEU", "ILE",
              "SER", "ASP", "LYS", "MET", "GLN",
              "TRP", "ARG", "ALA", "THR", "PRO",
              "PHE", "GLU", "HIS", "HIP", "TYR",
              "CYS", "HID", "HIE", "GLN"]

TER_CONSTR = 5

BACK_CONSTR = 0.5

CONSTR_ATOM = '''{{ "type": "constrainAtomToPosition", "springConstant": {0}, "equilibriumDistance": 0.0, "constrainThisAtom": "{1}:{2}:{3}" }},'''

CONSTR_DIST = '''{{ "type": "constrainAtomsDistance", "springConstant": {}, "equilibriumDistance": {}, "constrainThisAtom": "{}:{}:{}", "toThisOtherAtom": "{}:{}:{}" }},'''

CONSTR_CALPHA = '''{{ "type": "constrainAtomToPosition", "springConstant": {2}, "equilibriumDistance": 0.0, "constrainThisAtom": "{0}:{1}:_CA_" }},'''

class ConstraintBuilder(object):

    def __init__(self, pdb, ligand, gaps, metals, smiles):
        self.pdb = pdb
        self.ligand = ligand
        self.gaps = gaps
        self.metals = metals
        self.smiles = smiles
        self.residues = defaultdict(list)

    def _add_atom_id_to_dict(self, chain, atom_id, interval):
        # Only adds a new atom_id if it fullfills the interval condition
        for already_present_id in self.residues[chain]:
            if (abs(already_present_id - atom_id) < interval):
                break
        else:
            self.residues[chain].append(atom_id)

    def parse_atoms(self, interval=10):
        self.residues = defaultdict(list)
        initial_residue_found = False
        start = True
        with open(self.pdb, "r") as pdb:
            last_CA = {"atom_id": None, "chain": None}
            for line in pdb:
                atom_name = line[16:21].strip()
                atom_type = line[11:16].strip()
                # If you find something is not an atom pass line
                try:
                    atom_id = int(line[22:26].strip())
                except ValueError:
                    continue 
                chain = line[20:22].strip()
                if ((line.startswith("ATOM")) and (atom_name in AMINOACIDS) and (atom_type == "CA")):
                    try:
                        last_CA["chain"] = chain
                        last_CA["atom_id"] = atom_id
                        self._add_atom_id_to_dict(chain, atom_id, interval)
                    except TypeError:
                        pass
        terminal_id = last_CA["atom_id"]
        terminal_chain = last_CA["chain"]
        self.residues[terminal_chain].append(terminal_id)

    def build_constraint(self, BACK_CONSTR=BACK_CONSTR, TER_CONSTR=TER_CONSTR):

        init_constr = ['''"constraints":[''', ]

        terminal_constr = []
        back_constr = []
        # Backbone constraints
        for chain, atom_ids in self.residues.items():
            # Sort by ids to identify initial and final atoms
            atom_ids.sort()

            # Add constraints to terminal CA
            terminal_constr.append(CONSTR_CALPHA.format(chain, atom_ids[0], TER_CONSTR))
            terminal_constr.append(CONSTR_CALPHA.format(chain, atom_ids[-1], TER_CONSTR))

            # Add constraints to mid CA
            for atom_id in atom_ids[1:-1]:
                back_constr.append(CONSTR_CALPHA.format(chain, atom_id, BACK_CONSTR))

        # Gaps constraints
        gaps_constr = self.gaps_constraints()

        # Metal constraints
        metal_constr = self.metal_constraints()

        smiles_constr = self.constrain_smiles()

        final_constr = ["],"]

        terminal_constr[-1] = terminal_constr[-1].strip(',')

        constraints = init_constr + back_constr + gaps_constr + metal_constr + smiles_constr + terminal_constr + final_constr

        return constraints

    def gaps_constraints(self):
        #self.gaps = {}
        gaps_constr = []
        for chain, residues in self.gaps.items():
            gaps_constr = [CONSTR_ATOM.format(TER_CONSTR, chain, terminal, "_CA_") for terminals in residues for terminal in terminals]
        return gaps_constr

    def metal_constraints(self):

        metal_constr = []
        for metal, ligands in self.metals.items():
            metal_name, chain, metnum = metal.split(" ")
            for ligand in ligands:
                ligand_info, bond_lenght = ligand
                resname, resnum, chain, ligname = ligand_info.split(" ")
                metal_constr.append(CONSTR_DIST.format(TER_CONSTR, bond_lenght, chain, resnum, ligname, chain, metnum, metal_name))
        return metal_constr

    def constrain_smiles(self):
        if self.smiles:
            return ss.extract_constraint_from_smiles(self.ligand, self.smiles)
        else:
            return []
        


def retrieve_constraints(complex_pdb, ligand_pdb, gaps, metal, constrain_smiles=False, 
    back_constr=BACK_CONSTR, ter_constr=TER_CONSTR, interval=10):
    constr = ConstraintBuilder(complex_pdb, ligand_pdb, gaps, metal, constrain_smiles)
    constr.parse_atoms(interval=interval)
    constraints = constr.build_constraint(back_constr, ter_constr)
    return constraints

def parseargs():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('pdb', type=str, help='pdb to create the contraints on')
    parser.add_argument('conf', help='Control file to fill in. It need to templetazide with $CONSTRAINTS')
    parser.add_argument('--interval', type=int, help="Every how many CA to constraint")
    parser.add_argument('--ca', type=float, help="Constraint value to use on backbone CA", default=BACK_CONSTR)
    parser.add_argument('--terminal', type=float, help="Constraint value to use on terminal CA", default=TER_CONSTR)
    args = parser.parse_args()
    return os.path.abspath(args.pdb), os.path.abspath(args.conf), args.interval, args.conf, args.ca, args.terminal

if __name__ == "__main__":
    pdb, conf, interval, conf, back_constr, ter_constr = parseargs()
    constraints = retrieve_constraints(pdb, {}, {}, back_constr, ter_constr, interval)
