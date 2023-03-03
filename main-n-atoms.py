import xyz_py as xyzp
import numpy as np
from scipy.optimize import minimize 
from numpy import sqrt
from numpy import around
import random as rand

# Do not use scientific notation
np.set_printoptions(suppress=True)

# Parse a N-atom argon file.

N = 4
# 3 atoms: 1-2, 1-3, 2-3 (3)
# 4 atoms: 1-2, 1-3, 1-4, 2-3, 2-4, 3-4 (6)
# 5 atoms: 1-2, 1-3, 1-4, 1-5, 2-3, 2-4, 2,5, 3-4, 3-5, 4-5 (10)

filename = "./xyz/argon_" + str(N) + ".xyz"
atom_labels, cooridnates = xyzp.load_xyz(filename)
print("Parsed Atom labels:\n", atom_labels, "\n") # List
print("Parsed atom coorindates:\n", cooridnates, "\n") # np.array

# Loop and find the distance between all interaction

# Atom interaction counts:

i_list = list(range(1,N)) # 1 to 3 wher N = 4
j_list = list(range(1,N))
potential_count = 0

for i in i_list:
    for j in j_list:
        if i <= j:
            potential_count += 1
            ith = i
            jth = j+1 
            print("(i,j):", ith, jth)
            atom_i = cooridnates[ith - 1]
            atom_j = cooridnates[jth - 1]
            # Determine LJ Potential

            x_i = atom_i[0]
            y_i = atom_i[1]
            z_i = atom_i[2]
            x_j = atom_j[0]
            y_j = atom_j[1]
            z_j = atom_j[2]

            print(x_i, y_i, z_i)
            print(x_j, y_j, z_j, "\n")

            # Calculate the distance (r)
            r = np.sqrt(np.sum((atom_i-atom_j)**2, axis=0))
            print("Distance (r)", around(r, 3))

print("SUMMARY:", str(N), "atoms have", str(potential_count), 
      "interactions total\n")



# Save the global minimum coordinate
# Save the file









eps = 1 # kJ/mol
alpha = 1 # Angstroms

# # LJ Potential Implementation
# def LJPotential(params):
#     x1, y1, z1, x2, y2, z2, x3, y3, z3 = params
#     # Between 1 and 2
#     x_1 = x1 - x2
#     y_1 = y1 - y2
#     z_1 = z1 - z2
#     # Between 1 and 3
#     x_2 = x1 - x3
#     y_2 = y1 - y3
#     z_2 = z1 - z3
#     # Between 2 and 3
#     x_3 = x2 - x3
#     y_3 = y2 - y3
#     z_3 = z2 - z3

#     LJ_1 = ((alpha/sqrt((x_1)**2 + (y_1)**2 + (z_1)**2))**12
#                       - (alpha/sqrt((x_1)**2 + (y_1)**2 + (z_1)**2))**6)
#     LJ_2 = ((alpha/sqrt((x_2)**2 + (y_2)**2 + (z_2)**2))**12
#                       - (alpha/sqrt((x_2)**2 + (y_2)**2 + (z_2)**2))**6)
#     LJ_3 = ((alpha/sqrt((x_3)**2 + (y_3)**2 + (z_3)**2))**12
#                       - (alpha/sqrt((x_3)**2 + (y_3)**2 + (z_3)**2))**6)
#     return 4 * eps * (LJ_1 + LJ_2 + LJ_3)

# def minimizeNelderMead(initial_points):
#     result = minimize(LJPotential, initial_points, method="nelder-mead", options={'maxiter': 10000})

#     if result.success:
#         fitted_params = result.x
#         min_energy = result.fun
#         x1 = fitted_params[0]
#         y1 = fitted_params[1]
#         z1 = fitted_params[2]
#         x2 = fitted_params[3]
#         y2 = fitted_params[4]
#         z2 = fitted_params[5]
#         x3 = fitted_params[6]
#         y3 = fitted_params[7]
#         z3 = fitted_params[8]
#         print("\n")
#         return (min_energy)
#     else:
#         raise ValueError(result.message)

# # Generate a random set of coordinates for 3 atoms
# # [x1, y1, z1, x2, y2, z2, x3, y3, z3]

# for i in range(20):
#     random_initial_points = np.random.random_sample(size = 9) * 10 - 5 #  
#     min_energy = minimizeNelderMead(random_initial_points)
#     print("Min energy:", min_energy)
#     print("random init position", random_initial_points)

# 
