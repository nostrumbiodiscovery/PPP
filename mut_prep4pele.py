#!/bin/python
import sys
from prody import *

from adjustments_module import WritingAtomNames, FixStructureResnames, FixAtomNames, SolveClashes
from checks_module import CheckMutation, CheckClashes
from checks_module import CheckStructure, CheckforGaps
from global_processes import ParseArguments, FindInitialAndFinalResidues, PDBwriter, RenumberStructure
from hydrogens_addition import FixStructure
from mutational_module import Mutate

__author__ = 'jelisa'

"""
This block adds the hid, hip, and hie residues to prody. Otherwise this would consider
this aminoacids as heteroatoms.
"""
addNonstdAminoacid('HID', 'aromatic', 'basic', 'cyclic', 'large', 'polar', 'surface')
addNonstdAminoacid('HIE', 'aromatic', 'basic', 'cyclic', 'large', 'polar', 'surface')
addNonstdAminoacid('HIP', 'aromatic', 'basic', 'cyclic', 'large', 'polar', 'surface')
addNonstdAminoacid('CYT', 'neutral', 'acyclic', 'medium', 'polar', 'buried')
addNonstdAminoacid('LYN', 'neutral', 'acyclic', 'large', 'polar', 'buried')

def main(args):
    if args is None:
        sys.exit()

    try:
        initial_structure = parsePDB(args.input_pdb)
    except IOError:
        print "The file '{}' isn't a valid file\nCheck that it does exist and try again.".format(args.input_pdb)
        sys.exit()
    initial_residue, final_residue = FindInitialAndFinalResidues(initial_structure)
    # ff_parameters = ReadForceFieldParameters(args.force_field)

    print "* Checking for gaps."
    gaps, not_gaps = CheckforGaps(initial_structure)
    if gaps is None and not_gaps is None:
        print "WARNING: Problems when checking for gaps, so don't trust the existence of gaps."
        gaps, not_gaps = {}, {}
    print "* Checking for insertion codes."
    insertion_codes = [icode for icode in initial_structure.getIcodes() if icode]
    if insertion_codes:
        print " *The structure will be renumbered starting from 1 for each chain."
        structure2use = RenumberStructure(initial_structure, gaps, not_gaps)
    else:
        structure2use = initial_structure
    print "* Checking and Fixing the Residues Names:"
    structure2use = FixStructureResnames(structure2use, args.make_unique)
    print "* Checking and fixing the Atoms Names:"
    structure2use = FixAtomNames(structure2use, gaps, not_gaps)
    print "* Checking the structure for missing atoms:"
    residues2fix, residues2remove = CheckStructure(structure2use, gaps, not_gaps, args.remove_terminal_missing)
    if residues2fix:
        print '* Placing the missing atoms:'
        structure2use = FixStructure(structure2use, residues2fix)
    print args.mutation

    if not args.mutation:
        print 'Writing the structure to {}'.format(args.output_pdb[0])
        if args.make_unique:
            ligand_chain = structure2use.select("chain {}".format(args.make_unique))
            if ligand_chain:
                not_proteic_ligand = structure2use.select("chain {}".format(args.make_unique)).hetero
            else:
                not_proteic_ligand = None
            PDBwriter(args.output_pdb[0], WritingAtomNames(structure2use), args.make_unique, residues2remove,
                      args.no_gaps_ter, not_proteic_ligand, gaps, not_gaps)
        else:
            not_proteic_ligand = None
            PDBwriter(args.output_pdb[0], WritingAtomNames(structure2use), args.make_unique, residues2remove,
                      args.no_gaps_ter, not_proteic_ligand, gaps, not_gaps)
        sys.exit()
    else:
        clashes = []
        mutated_structure = None
        for mutation, output_file in zip(args.mutation, args.output_pdb):
            print '* Checking the mutation:'
            print " Mutation: {0[ini_resname]} {0[chain]} {0[resnum]} {0[fin_resname]}".format(mutation)
            correct_mutation = CheckMutation(structure2use, mutation)
            if not correct_mutation:
                exit_message = "The mutation was incorrect, check your parameters.\n" \
                               "The checked structure will be written to {}".format(output_file)
                PDBwriter(output_file, WritingAtomNames(structure2use), args.make_unique, gaps, args.no_gaps_ter, not_gaps)
                continue
            else:
                print "Output_file name: {0}".format(output_file)
                mutated_structure, zmatrix = Mutate(structure2use, mutation)
                if not args.mutant_multiple:
                    if args.mutation[0]['fin_resname'] in ["ALA", "GLY"]:
                        print "The ALA and the GLY don't have any rotamer to try."
                    else:
                        print "Checking Clashes:"
                        try:
                            clashes = CheckClashes(mutated_structure, mutation, zmatrix,
                                                   initial_residue, final_residue)
                        except ValueError:
                            pass
                        else:
                            if not clashes:
                                print "Structure without clashes."
                            else:
                                mutated_structure = SolveClashes(mutated_structure, clashes,
                                                                 mutation, zmatrix,
                                                                 initial_residue, final_residue)
                    mutated_structure.setTitle("mutated structure")
                    PDBwriter(output_file, WritingAtomNames(mutated_structure), args.make_unique, gaps, args.no_gaps_ter,
                              not_gaps)
                else:
                    print "Multiple mutations at the same time are still under development."
                    structure2use = mutated_structure
        if args.mutant_multiple and mutated_structure is not None:
            PDBwriter(args.output_pdb, WritingAtomNames(mutated_structure), gaps, not_gaps, args.no_gaps_ter)


if __name__ == '__main__':
    arguments = ParseArguments()
    main(arguments)
