# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

#create parser (takes input data and builds readable data structure)

parser = PDBParser(PERMISSIVE=1)

#load PDB file

def PDB_read(PDB_structure_ID, filename):
    structure_id = PDB_structure_ID
    filename = filename
    structure = parser.get_structure(structure_id, filename)
    return structure

#select only CA atoms
class Select(object):
    """Select everything for PDB output (for use as a base class).

    Default selection (everything) during writing - can be used as base class
    to implement selective output. This selects which entities will be written out.
    """
    
    def accept_model(self, model):
        """Overload this to reject models for output."""
        return 1

    def accept_chain(self, chain):
        """Overload this to reject chains for output."""
        return 1

    def accept_residue(self, residue):
        """Overload this to reject residues for output."""
        return 1

    def accept_atom(self, atom):
        if atom.get_id()=="CA" or atom.get_id()=="N" or atom.get_id()=="C":
            return 1
        else:
            return 0

def PDB_write(structure, structure_id, file_name):
    """generates output file, only CA, N and C atoms are selected. In order to change the selected residues, atoms etc.,
    the Select class has to be changed"""
    io = PDBIO()
    io.set_structure(structure)
    io.save("../output_data/" + structure_id+file_name+".pdb", Select())
