#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 12 00:28:26 2018

@author: tobiashoch
"""
from __future__ import division

import math

##class variables
#defining bond parameters r_eq, k
N_CA = [1.455, 300]
CA_C = [1.510, 300]
C_N = [1.325, 300]

#defining bond angle parameters r_theta, K_theta,
N_CA_C = [111, 300]
CA_C_N = [116, 300]
C_N_CA = [122, 300]

#defining dihedral angle parameters k, n, psi0
#psi N-CA-C-N"""
psi = [-0.3, 1, 0]

#phi, C-N-CA-C
phi = [-0.3, 1, 0]

#omega, CA-C-N-CA
omega = [67.0, 1, 180]

#reduced potential energy function
"""E(pot) = E(bond) + E(angle) + E(torsion)"""

def bond_energy(length):
    energy = 0
    for i in range(0, len(length)-2, 3):
        energy += (N_CA[1]/2)*((length[i]-N_CA[0])**2)
        energy += (CA_C[1]/2)*((length[i+1]-CA_C[0])**2)
        energy += (C_N[1]/2)*((length[i+2]-C_N[0])**2)
    return float(energy)

def angle_energy(angles):
    energy = 0
    for i in range(0, len(angles)-2, 3):
        energy += (N_CA_C[1]/2)*((angles[i]-N_CA_C[0])**2)
        energy += (CA_C_N[1]/2)*((angles[i+1]-CA_C_N[0])**2)
        energy += (C_N_CA[1]/2)*((angles[i+2]-C_N_CA[0])**2)
    return float(energy)

def torsion_energy(dihedral_deg):
    energy = 0
    for i in range(0, len(dihedral_deg), 3):
        energy += psi[0]*(1-math.cos((psi[1]*dihedral_deg[i])-psi[2]))

    for i in range(1, len(dihedral_deg), 3):
        energy += phi[0]*(1-math.cos((phi[1]*dihedral_deg[i])-phi[2]))

    for i in range(2, len(dihedral_deg), 3):
        energy += omega[0]*(1-math.cos((omega[1]*dihedral_deg[i])-omega[2]))

    return energy

def pot_energy(distance, angles_deg, dihedral_deg):
    bond = bond_energy(distance)
    angle = angle_energy(angles_deg)
    torsion = torsion_energy(dihedral_deg)
    E_pot = bond + angle + torsion

    return E_pot
