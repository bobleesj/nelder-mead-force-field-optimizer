import math
import numpy as np
import copy
import xyz_py as xyzp
from numpy import sqrt
from numpy import around
from scipy.optimize import minimize 
import random as rand
from numpy.linalg import norm
import numpy as np
import vg
from opt_algo import nelder_mead
 
# This code attempts to implement the flexible water model
# Potential is currently wrong. Need to fix again.


np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

# Global Variables
FILENAME = "./xyz/water_1"
k_l = 547.5 # kcal/(mol Å^-2)
l = 1.0 # Å
k_cub = -1.65 # Å^-1
k_theta = 49.9 # kcal/(mol rad^-1)
theta = 109.5 # deg

# Import pre-optimized structures
atom_labels, cooridnates = xyzp.load_xyz(FILENAME + ".xyz")
intial_coordinates = np.array(cooridnates).flatten()

# Import a water moelcule 

# U consists of the bonds adn the angles terms.

# Part I. Bonds
print(atom_labels)
print(cooridnates)

# Step 0. Check the math of the intramolecular at optimal oposition
# OH_1 = O - H1
# OH_2 = O - H2
# angle = vg.angle(OH_1, OH_2)
# angle_rad = angle * (math.pi/180)

# first_term = (r1-l)**2 + k_cub*(r1-l)**3
# second_term = (r2-l)**2 + k_cub*(r2-l)**3
# bonds_term = k_l*(first_term + second_term)

# angles_term = k_theta * (angle_rad - theta)**2

# U = bonds_term + angles_term



# Step 1. Determine the angle between the atoms (2 angles

# LJ Potential Implementation
def Ferguson(params):
    cooridnates = params.reshape(-1,3)
    O = cooridnates[0]
    H1 = cooridnates[1]
    H2 = cooridnates[2]
    r1 = np.sqrt(np.sum((O-H1)**2, axis=0))
    r2 = np.sqrt(np.sum((O-H2)**2, axis=0))

    OH_1 = O - H1
    OH_2 = O - H2
    angle = vg.angle(OH_1, OH_2)
    angle_rad = angle * (math.pi/180)
  
    first_term = (r1-l)**2 + k_cub*(r1-l)**3
    second_term = (r2-l)**2 + k_cub*(r2-l)**3

    bonds_term = k_l*(first_term + second_term)
    angles_term = k_theta * (angle - theta)**2

    U = bonds_term + angles_term
    
    print("r1:", np.round(r1, 4), "Å")
    print("r2:", np.round(r2, 4), "Å")
    print("angle:", np.round(angle, 3), "°")
    print("U:", np.round(U, 3), "kcal/mol")
    print("\n")

    return U
            
def minimizeNelderMead(initial_points):
    result = minimize(Ferguson, initial_points, method="nelder-mead", options={'maxiter': 1000})
    if result.success:
        cooridnate_list = result.x
        min_energy = result.fun
        return (min_energy, cooridnate_list)
    else:
        raise ValueError(result.message)

#energy, coordinate_list  = minimizeNelderMead(intial_coordinates)
#print(energy, coordinate_list)

# The moment it starts to increase the value, I should stop
res = nelder_mead(Ferguson, intial_coordinates, 0.00001, 10000)
print(res)

# Save the file
min_coordinate_reshape = res[0].reshape(-1,3)
xyzp.save_xyz(FILENAME + "_opt.xyz", atom_labels, min_coordinate_reshape)
# Import 
   

# min_coordinate_reshape = res[0].reshape(-1,3)
# xyzp.save_xyz(FILENAME + "_opt.xyz", atom_labels, min_coordinate_reshape)
# for i in range(5):
#     intial_coordinates = np.array(cooridnates).flatten()
#     energy, coordinate_list  = minimizeNelderMead(intial_coordinates)
#     if energy < min_energy:
#         min_energy = energy
#         min_coordinate = coordinate_list
#         min_coordinate_reshape = min_coordinate.reshape(-1,3)
#         print("\n*New min energy (kcal/mol):", around(min_energy, 3))
#         print("*New min cooridnates:\n", around(min_coordinate_reshape, 3), "\n")
    
#     print("Current energy (kJ/mol):", around(energy, 3))
# Check whether l refers to the lenght b/w O and H
# I calcualted twice since there are 2 bonds

# TRY TO OPTIMZIE THE STRUCTURE AND FIND THE MINIMUM.

# Now, let's use the neader mead algo to find the optmized structures.
# U: 578546.2738203785