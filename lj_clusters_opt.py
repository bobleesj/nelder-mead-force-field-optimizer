
# Sangjoon (Bob) Lee
# CH393
# Spring 2023

import numpy as np
import random as rand
from algorithm.nelder_mead import minimize
import xyz_py as xyzp
import os

# Disable scientific notation for printing numpy arrays
np.set_printoptions(suppress=True)

'''
Summary - This script finds the global minimum of a Lennard-Jones (LJ) clusters
using the Nelder-Mead optimization algorithm.

The xyz coordinate values of Ar were acquired from the following source
http://doye.chem.ox.ac.uk/jon/structures/LJ.html
'''

# Global Variables
ITERATIONS = 30  # Number of iterations for the optimization algorithm
EPSILON = 1  # LJ potential parameter
ALPHA = 1  # LJ potential parameter

# USER INPUT:  Get file name from user input
while True:
    filename = input("Q1/2. Enter the file name Ex) ./xyz/argon/argon_2.xyz: ").strip()
    if not os.path.exists(filename):
        print(f"Error: File {filename} does not exist.")
    elif not filename.endswith('.xyz'):
        print(f"Error: File {filename} is not a valid xyz file.")
    else:
        break

# USER INPUT: Determine whether to add randomness to the initial posiiton of the molecules
isPerturbated = ""
while isPerturbated not in ["y", "n"]:
    isPerturbated = input("Q2/2. Would you like to add perturbation to initial positions? (y/n) ").strip().lower()

# Load initial coordinates from file
atom_labels, coordinates = xyzp.load_xyz(filename)
initial_coordinates = np.array(coordinates).flatten()

# Count the number of atoms
N = len(coordinates)
print(f"\n---LJ Potential for {N} Argon Atoms---")

def lj_potential(params):
    coordinates = params.reshape(-1, 3)  # Reshape the flattened array into a 2D array of coordinates
    lj_sum = 0  # Initialize the LJ potential energy to zero
    
    for i in range(N):
        for j in range(i+1, N):
            atom_i = coordinates[i]
            atom_j = coordinates[j]
            r = np.linalg.norm(atom_i - atom_j)
            lj = ((ALPHA / r)**12 - (ALPHA / r)**6)
            lj_sum += lj
    return 4 * EPSILON * lj_sum

print("Initial LJ potential", lj_potential(initial_coordinates), "\n")

min_energy = 0  # Initialize the minimum LJ potential energy to zero
min_coordinates = []  # Initialize the coordinates corresponding to the minimum LJ potential energy to an empty list

# Perform Nelder-Mead optimization for multiple random initial points
for i in range(ITERATIONS):
    if isPerturbated == "y":
        random_initial_points = np.random.random_sample(size=N * 3)
        result = minimize(lj_potential, initial_coordinates + random_initial_points, 0.01, 1000)
    else:
        result = minimize(lj_potential, initial_coordinates, 0.01, 1000)
    current_coordinates = result[0].reshape(-1, 3)
    current_energy = result[1]

    print(f"Iteration {i+1}: Min energy: {current_energy}")
    
    if current_energy < min_energy:
        min_energy = current_energy
        min_coordinates = current_coordinates

# Save optimized coordinates and energy to files
print("\nLowest energy found", min_energy)
print(min_coordinates, "\n")
xyzp.save_xyz(filename + "_opt.xyz", atom_labels, min_coordinates)

with open(filename + "_opt.txt", "w") as f:
    f.write(str(min_energy))