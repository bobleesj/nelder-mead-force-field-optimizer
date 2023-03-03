import xyz_py as xyzp
import numpy as np
from scipy.optimize import minimize 
from numpy import sqrt
from numpy import around
import random as rand

# Global Variables
N = 5
ITERATION = 100
FILENAME = "./xyz/argon_" + str(N)
# EPS = 0.997 # kJ/mol
# ALPHA = 3.4 # Angstroms
EPS = 1
ALPHA = 1

# Do not use scientific notation
np.set_printoptions(suppress=True)

# LJ Potential Implementation
def LJPotential(params):
    i_list = list(range(1,N)) # (1 to N - 1)
    j_list = list(range(1,N))
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
        return (min_energy, cooridnate_list)
    else:
        raise ValueError(result.message)


# Import 
atom_labels, cooridnates = xyzp.load_xyz(FILENAME + ".xyz")
intial_coordinates = np.array(cooridnates).flatten()
min_energy = 0
min_coordinate = []

for i in range(ITERATION):
    random_initial_points = np.random.random_sample(size = N * 3) * 6  #  
    energy, coordinate_list  = minimizeNelderMead(intial_coordinates + random_initial_points)
    if energy < min_energy:
        min_energy = energy
        min_coordinate = coordinate_list
        min_coordinate_reshape = min_coordinate.reshape(-1,3)
        print("\n*New min energy (kJ/mol):", around(min_energy, 3))
        print("*New min cooridnates:\n", around(min_coordinate_reshape, 3), "\n")
    
    print("Current energy (kJ/mol):", around(energy, 3))
    
print("\nMin energy for this run (kJ/mol):", around(min_energy, 3))
print("Min coordinates for this:\n", around(min_coordinate_reshape, 3), "\n")

xyzp.save_xyz(FILENAME + "_opt.xyz", atom_labels, min_coordinate_reshape)
f = open(FILENAME + "_opt.txt", "w")
f.writelines(str(min_energy))
f.close()