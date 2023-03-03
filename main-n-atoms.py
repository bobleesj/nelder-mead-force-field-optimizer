import xyz_py as xyzp
import numpy as np
from scipy.optimize import minimize 
from numpy import sqrt
from numpy import around
import random as rand

# Global Variables
N = 2
ITERATION = 20
FILENAME = "./xyz/argon_" + str(N) + ".xyz"
# EPS = 0.997 # kJ/mol
# ALPHA = 3.4 # Angstroms
EPS = 1
ALPHA = 1
# Do not use scientific notation
np.set_printoptions(suppress=True)

atom_labels, cooridnates = xyzp.load_xyz(FILENAME)
print("Parsed Atom labels:\n", atom_labels, "\n") # List
print("Parsed atom coorindates:\n", cooridnates, "\n") # np.array

# Loop and find the distance between all interaction
i_list = list(range(1,N)) # (1 to N - 1)
j_list = list(range(1,N))

# LJ Potential Implementation
def LJPotential(params):
    LJ_sum = 0
    cooridnates = params.reshape(-1,3)
    for i in i_list:
        for j in j_list:
            if i <= j:
                ith = i
                jth = j+1 
                atom_i = cooridnates[ith - 1]
                atom_j = cooridnates[jth - 1]
                r = np.sqrt(np.sum((atom_i-atom_j)**2, axis=0))
                LJ = 4 * EPS * ((ALPHA/r)**12 - (ALPHA/r)**6)
                LJ_sum += LJ
    return LJ_sum
            
def minimizeNelderMead(initial_points):
    result = minimize(LJPotential, initial_points, method="nelder-mead", options={'maxiter': 100000})

    if result.success:
        cooridnate_list = result.x
        min_energy = result.fun
        # print(cooridnate_list)
        # print("\n")
        return (min_energy, cooridnate_list)
    else:
        raise ValueError(result.message)

for i in range(ITERATION):
    random_initial_points = np.random.random_sample(size = N * 3) * 6 - 5 #  
    min_energy, coordinate_list  = minimizeNelderMead(random_initial_points)
    print("Min energy:", min_energy)





# Atom interaction counts:
# 3 atoms: 1-2, 1-3, 2-3 (3)
# 4 atoms: 1-2, 1-3, 1-4, 2-3, 2-4, 3-4 (6)
# 5 atoms: 1-2, 1-3, 1-4, 1-5, 2-3, 2-4, 2,5, 3-4, 3-5, 4-5 (10)

