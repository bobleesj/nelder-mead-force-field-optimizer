import xyz_py as xyzp
import numpy as np
from scipy.optimize import minimize 
from numpy import sqrt
from numpy import around
import random as rand

# Do not use scientific notation
np.set_printoptions(suppress=True)

# Parse a N-atom argon file.
N = 3
# eps = 0.997 # kJ/mol
# alpha = 3.4 # Angstroms
eps = 1
alpha = 1
filename = "./xyz/argon_5_opt.xyz"

atom_labels, cooridnates = xyzp.load_xyz(filename)
print("Parsed Atom labels:\n", atom_labels, "\n") # List
print("Parsed atom coorindates:\n", cooridnates, "\n") # np.array

# Loop and find the distance between all interaction
# Atom interaction counts:
# 3 atoms: 1-2, 1-3, 2-3 (3)
# 4 atoms: 1-2, 1-3, 1-4, 2-3, 2-4, 3-4 (6)
# 5 atoms: 1-2, 1-3, 1-4, 1-5, 2-3, 2-4, 2,5, 3-4, 3-5, 4-5 (10)

i_list = list(range(1,N))
j_list = list(range(1,N))

interaction_count = 0
LJ_potential_sum = 0
for i in i_list:
    for j in j_list:
        if i <= j:
            interaction_count += 1
            ith = i
            jth = j+1 
            print("(i,j):", ith, jth)
            atom_i = cooridnates[ith - 1]
            atom_j = cooridnates[jth - 1]
            
            # Calculate the distance (r)
            r = np.sqrt(np.sum((atom_i-atom_j)**2, axis=0))

            # Calculate LJ potential
            LJ = ((alpha/r)**12 - (alpha/r)**6)

            # Print r and LJ
            print("r (Angstroms)", around(r, 3))
            print("LJ potential (kJ/mol):", round(LJ, 5), "\n")
            LJ_potential_sum += LJ
            
print("---OUTOUT---")
print("LJ potential sum (kJ/mol):",  4 * eps * LJ_potential_sum)

print("SUMMARY:", str(N), "atoms have", str(interaction_count), 
      "interactions total\n")
