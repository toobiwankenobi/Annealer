#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 22 15:51:59 2018

@author: tobiashoch
"""

from Bio.PDB.PDBIO import PDBIO
import random

#change coordinates of atoms randomly
def uniform_randomizer(structure, start_point, end_point):
    for model in structure.get_list():
            for chain in model.get_list():
                for residue in chain.get_list():
                    for atom in residue.get_list():
                        coords = atom.get_coord()
                        random_value = random.uniform(start_point, end_point)
                        coords[0] = coords[0]*random_value
                        coords[1] = coords[1]*random_value
                        coords[2] = coords[2]*random_value
    return structure
                    
def random_write(structure, structure_id):
    io = PDBIO()
    io.set_structure(structure)
    io.save(structure_id+"_backbone_random.pdb")


