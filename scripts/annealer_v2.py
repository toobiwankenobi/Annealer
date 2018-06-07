#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue May 15 10:25:06 2018

@author: tobiashoch
"""

from __future__ import division
import random
import randomizer
import calculations as calc
import PDB_reader as pdb
import energy
import math
import timeit
import copy
import matplotlib.pyplot as plt
import numpy as np

#probability of acceptance function
def acceptance_probability(energy_best, energy_new, str_new, str_actual, temperature):
    k = 1 #this paramter has a great impact on the probability of accepting a higher-energy structre, be careful!
    prob_acceptance = math.exp(float(-(energy_new-energy_best)) / float(k*temperature))
    x = random.uniform(0,1)
    if prob_acceptance > x:
        return str_new
    else:
        return str_actual

#annealing function
def annealer(initial_structure, initial_energy):
    str_best = copy.deepcopy(initial_structure)
    energy_best = initial_energy
    print("initial energy = %f " % energy_best)
    str_actual = copy.deepcopy(initial_structure)
    list_all = []
    list_best = [[0, initial_energy]]


    #all simulation paramters
    starting_point = 0.999 #range of 3D motion applied to each atom
    ending_point = 1.001
    T_max = 25000.0 #starting temperature for simulated annealing
    T_min = 0.02
    T_new = 25000.0
    steps = 50 #number or repeating steps for each temperature
    temp_updates = 10 #number of temperature updates used during the simulated annealing
    print("for loop starts")
    counter = 0
    counter_list = []
    for i in range(int(T_max), int(T_min), -int(T_max/temp_updates)):
        for i in range(steps):
            #do random alterations on each residue
            str_new = randomizer.uniform_randomizer(copy.deepcopy(str_actual), starting_point, ending_point)

            #re-calculate structure paramters
            coordinates_new, distance_new, residue_new = calc.length(str_new)
            angles_rad_new, angles_deg_new = calc.angles(coordinates_new)
            dihedral_rad_new, dihedral_deg_new = calc.dihedral(coordinates_new, residue_new)

            #calculate new energy
            energy_new = energy.pot_energy(distance_new, angles_deg_new, dihedral_deg_new)
            list_all.append([counter, energy_new])

            #set new structure if energy is lower or set new structure with higher energy with
            #a given probability
            if energy_new < energy_best:
                print("new energy = %f, best energy = %f, accepted" % (energy_new, energy_best))
                str_best, energy_best = copy.deepcopy(str_new), copy.deepcopy(energy_new)
                list_best.append([counter, energy_best])
                str_actual = copy.deepcopy(str_new)

            str_actual = acceptance_probability(energy_best, energy_new, str_new, str_actual, T_new)
            counter+=1
            counter_list.append(counter)
        #set new temperature
        T_new += -(T_max/temp_updates)
        print("The new temperature has been set to %s" % T_new + " Kelvin")

    calc_coordinates, calc_distance, calc_residue = calc.length(str_best)
    calc_angles_rad, calc_angles_deg = calc.angles(calc_coordinates)
    calc_dihedral_rad, calc_dihedral_deg = calc.dihedral(calc_coordinates, calc_residue)

    energy_best_calc = energy.pot_energy(calc_distance, calc_angles_deg, calc_dihedral_deg)

    print("calculated energy best = ", energy_best_calc)
    print("returned energy best = ", energy_best)

#    pdb.PDB_write(str_actual, "2ID8", "_actual_structure")

    if energy_best == initial_energy:
        print("Structure has already lowest-energy conformation")

    elif energy_best < initial_energy:
        pdb.PDB_write(str_best, "2ID8", "_new_structure") #write new pdb file with lower energy
        print("Structure with lower energy has been created")

    return str_best, energy_best, str_actual, steps, temp_updates, starting_point, ending_point, list_best, list_all, counter_list

# =============================================================================
# Structure loading, structure manipulation, simulated annealing commands
# =============================================================================

#load original pdb file and select C_alpa, C' and N atoms only and create new pdb file
original_structure = pdb.PDB_read("2ID8", "../input_data/2id8.pdb")
pdb.PDB_write(original_structure, "2ID8", "_backbone")

#laod backbone
backbone = pdb.PDB_read("2ID8", "../output_data/2ID8_backbone.pdb")

#change position of single atoms and create new pdb file
random_structure = randomizer.uniform_randomizer(backbone, 0.95, 1.05) #exampe: plus/minus 5 percent motion to each atom
pdb.PDB_write(random_structure, "2ID8", "_backbone_random")

#load and prepare inital state
initial_structure = pdb.PDB_read("2ID8", "../output_data/2ID8_backbone_random.pdb") #can be changed to the original instead of the randomized structure
initial_coordinates, initial_distance, initial_residue = calc.length(initial_structure)
initial_angles_rad, initial_angles_deg = calc.angles(initial_coordinates)

initial_dihedral_rad, initial_dihedral_deg = calc.dihedral(initial_coordinates, initial_residue)
initial_energy = energy.pot_energy(initial_distance, initial_angles_deg, initial_dihedral_deg)

#start simulated annealing
start = timeit.default_timer()
str_best, energy_best, str_actual, steps, updates, start_point, end_point, list_best, list_all, counter_list = annealer(initial_structure, initial_energy)
stop = timeit.default_timer()
time = stop-start

#create results output
with open('../output_data/simulation_results.txt', 'w') as f:
    print("Input structure = backbone_random.pdb", file=f) #has to be changed manually
    print("Number of steps for each temperature = %f" % steps, file=f)
    print("Number of temperature updates = %f" % updates, file=f)
    print("Changes in atom positions applied in percent = min: %f max: %f" % ((start_point-1)*100, (end_point-1)*100), file=f)
    print("Running time of the algorithm = %s" % time + " seconds", file=f)
    print("Input energy = %f" % initial_energy, file=f)
    print("New energy = %f" % energy_best, file=f)

#creating plot of energies
x_all = []
y_all = []
for i in range(0, len(list_all)):
    x_all.append(list_all[i][0])
    y_all.append(list_all[i][1])

x_best = []
y_best = []
for i in range(0, len(list_best)):
    x_best.append(list_best[i][0]+1)
    y_best.append(list_best[i][1])

plt.plot(counter_list, y_all, linestyle='-', mec='b', marker='o', color='grey', label="all energies")
plt.scatter(x_best, y_best, color='red', label="lowest energies", s=100)
plt.xlabel("Steps", fontsize=15)
plt.ylabel("Energy", fontsize=15)
plt.title("Conformational energy vs. simulation steps, motion range = ±" + str(round(((end_point-1)*100), 3)) + "%", fontsize=12, y=1.08)
plt.grid(True)
plt.legend()
plt.savefig("../output_data/energy_vs_steps.png")
plt.show()
