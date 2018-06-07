#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 22:17:27 2018

@author: tobiashoch
"""
from __future__ import division
import numpy as np
import math

def length(structure):
    atom_coords = []
    residues = []
    atom_name = []
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                resinfo = []
                resinfo.append(str(residue.get_resname()))
                for atom in residue.get_list():
                    resinfo.append(str(atom.get_name()))
                    resinfo.append(list(atom.get_coord()))
                    atom_coords.append(atom.get_coord())
                    atom_name.append(atom.get_name())
                    residues.append(resinfo)

    bond_distance = []
    for atoms in range(1,len(atom_coords)):
            dist = math.sqrt(((atom_coords[atoms][0]-atom_coords[atoms-1][0])**2) + ((atom_coords[atoms][1]-atom_coords[atoms-1][1])**2) + ((atom_coords[atoms][2]-atom_coords[atoms-1][2])**2))
            bond_distance.append(dist)

    return atom_coords, bond_distance, residues

def angles(coords):
    degrees = []
    radians = []
    for atoms in range(0,len(coords)-2):
        a = coords[atoms+1] - coords[atoms]
        b = coords[atoms+2] - coords[atoms+1]
        A = math.sqrt((a[0]*a[0]) + (a[1]*a[1]) + (a[2]*a[2]))
        B = math.sqrt((b[0]**2) + (b[1]**2) + (b[2]**2))
        dot = np.dot(b,a)
        gamma_rad = math.acos(-dot/(A*B))
        gamma_deg = math.degrees(gamma_rad)
        radians.append(gamma_rad)
        degrees.append(gamma_deg)

    return radians, degrees

def dihedral(coordinates, residue):
    torsion_radians = [0]
    torsion_degrees = [0]
    for atom in range(0, len(coordinates)-3):
        b1 = coordinates[atom+1]-coordinates[atom]
        b2 = coordinates[atom+2]-coordinates[atom+1]
        b3 = coordinates[atom+3]-coordinates[atom+2]
        cross1 = np.cross(b1,b2)
        cross2 = np.cross(b2,b3)
        n1 = cross1/(math.sqrt(cross1[0]**2 + cross1[1]**2 + cross1[2]**2))
        n2 = cross2/(math.sqrt(cross2[0]**2 + cross2[1]**2 + cross2[2]**2))
        u1 = n2
        u3 = b2/(math.sqrt(b2[0]**2 + b2[1]**2 + b2[2]**2))
        u2 = np.cross(u3,u1)
        theta_radians = - math.atan2(np.dot(n1,u2), np.dot(n1,u1))
        theta_degrees = math.degrees(theta_radians)
        torsion_radians.append(theta_radians)
        torsion_degrees.append(theta_degrees)

    extend = [0,0]
    torsion_radians.extend(extend)
    torsion_degrees.extend(extend)

    return torsion_radians, torsion_degrees
