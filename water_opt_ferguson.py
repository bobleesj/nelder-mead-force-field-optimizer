import math
import numpy as np
import xyz_py as xyzp
from numpy import sqrt
from numpy import around
from algorithm.nelder_mead import minimize
from util.water_info import calculate_molecule_properties
import random as rand
import vg
import os
# from nelder_mead import minimize
 
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

'''
Summary: This script performs a geometry optimization for water dimers using the
Ferguson Flexible Water Model. It uses the Nelder-Mead optimization algorithm
to find the lowest energy geometry based on the potential energy function
implemented in the Ferguson() function. The input is an XYZ file with the
initial geometry of the water dimer, and the user can also choose to add
perturbation to the initial coordinates for each iteration. The output is
the optimized geometry in an XYZ file and the lowest energy found in a text file.
'''

'''
The parameters are acquired from the following paper:
Parameterization and evaluation of a flexible water model
https://doi.org/10.1002/jcc.540160413
'''

# Global variables (Ferguson Flexible Water Model)
K_L = 547.5 # kcal/(mol Å^-2) - spring constant for bond lengths
L = 1.0 # Å / equlibrium bond length between O and H
K_CUB = -1.65 # Å^-1 - cubic force constant for bond lengths
K_THETA = 49.9 # kcal/(mol rad^-1) - force constant for bond angles
THETA = 109.5 # deg /equlibrium angle between O and H
A = 650000 # kcal mol^-1 Å^12 - Lennard-Jones potential constant
B = 625.47 # kcal mol^-1 Å^6 - Lennard-Jones potential constant
Q_H = 0.413 # kcal^0.5 Å^0.5 mol^-0.5 - partial charge on hydrogen atoms
Q_O = -0.826 # kcal^0.5 Å^0.5 mol^-0.5 - partial charge on oxygen atoms
K_VALUE = 332.0637133 # Coulomb's constant in kcal/mol-Å-e^2 units

# Neader Mead Parameters
ITERATION = 10
STEP_SIZE = 0.01
MAX_ITERATION = 10000

# USER INPUT (1/2):  Get file name from user input
while True:
    filename = input("Q1/2. Enter the file name Ex) ./xyz/water/tip4p_wales_n_2.xyz: ").strip()
    if not os.path.exists(filename):
        print(f"Error: File {filename} does not exist.")
    elif not filename.endswith('.xyz'):
        print(f"Error: File {filename} is not a valid xyz file.")
    else:
        break

# USER INPUT (2/2): Determine whether to add randomness to the initial posiiton of the molecules
isPerturbated = ""
while isPerturbated not in ["y", "n"]:
    isPerturbated = input("Q2/2. Would you like to add perturbation to the initial coordinates for every iteration? (y/n) ").strip().lower()

# Ferguson potential function
def Ferguson(params):
   
    coordinates = params.reshape(-1,3)
    O1 = coordinates[0]
    H11 = coordinates[1]
    H21 = coordinates[2]
    O2 = coordinates[3]
    H12 = coordinates[4]
    H22 = coordinates[5]

    # Hij where
    # i refers to the ith number of hydrogens
    # j refers to the jth number of oxygen
    
    # Calculate bond lengths
    r11 = np.sqrt(np.sum((O1-H11)**2, axis=0))
    r21 = np.sqrt(np.sum((O1-H21)**2, axis=0))
    r12 = np.sqrt(np.sum((O2-H12)**2, axis=0))
    r22 = np.sqrt(np.sum((O2-H22)**2, axis=0))

    # Calculate angles
    OH11 = O1 - H11
    OH21 = O1 - H21
    OH12 = O2 - H12
    OH22 = O2 - H22
    angle1 = vg.angle(OH11, OH21)
    angle2 = vg.angle(OH12, OH22)

    angle_rad1 = angle1 * (math.pi/180)
    angle_rad2 = angle2 * (math.pi/180)
  
    # Calculate bonds term
    first_term = (r11-L)**2 + K_CUB*(r11-L)**3
    second_term = (r21-L)**2 + K_CUB*(r21-L)**3
    third_term = (r12-L)**2 + K_CUB*(r12-L)**3
    fourth_term = (r22-L)**2 + K_CUB*(r22-L)**3
    bonds_term = K_L*(first_term + second_term + third_term + fourth_term)

    # Calculate angles term
    angles_term = K_THETA * ((angle_rad1 - THETA * math.pi/180)**2 + (angle_rad2 - THETA * math.pi/180)**2)

    # Calculate oxygen repulsion term
    rO = np.sqrt(np.sum((O1-O2)**2, axis=0))
    LJ_term = A/(rO**12) - B/(rO**6)
    
    # Calculate Coulombic term - OH interactions
    rO1H12 = np.sqrt(np.sum((O1-H12)**2, axis=0)) 
    rO1H22 = np.sqrt(np.sum((O1-H22)**2, axis=0))
    rO2H11 = np.sqrt(np.sum((O2-H11)**2, axis=0)) 
    rO2H21 = np.sqrt(np.sum((O2-H21)**2, axis=0))

    OH_interations = (Q_H*Q_O) * (1/rO1H12 + 1/rO1H22 + 1/rO2H11 + 1/rO2H21)

    # Calculate Coulombic term - HH interactionsb
    rH11H12 = np.sqrt(np.sum((H11-H12)**2, axis=0)) 
    rH21H12 = np.sqrt(np.sum((H21-H12)**2, axis=0)) 
    rH11H22 = np.sqrt(np.sum((H11-H22)**2, axis=0)) 
    rH21H22 = np.sqrt(np.sum((H21-H22)**2, axis=0))

    HH_interations = (Q_H*Q_H) * (1/rH11H12 + 1/rH21H12 + 1/rH11H22 + 1/rH21H22)
    OO_interaction = (Q_O*Q_O) * (1/rO)

    U = bonds_term + angles_term + LJ_term + K_VALUE * (OH_interations + HH_interations + OO_interaction)

    return U

# Import coordinates from XYZ file 
atom_labels, cooridnates = xyzp.load_xyz(filename)
initial_coordinates = np.array(cooridnates).flatten()
min_energy = 0
min_coordinate = []

# Perform Nelder-Mead optimization
for i in range(ITERATION):
    if isPerturbated == "y":
        random_initial_points = np.random.random_sample(size = 18) * 0.1 - 0.05
        result = minimize(Ferguson, initial_coordinates + random_initial_points, STEP_SIZE, MAX_ITERATION)
    else:
        result = minimize(Ferguson, initial_coordinates, STEP_SIZE, MAX_ITERATION)
    current_coordinates = result[0].reshape(-1, 3)
    current_energy = result[1]

    print(f"\nIteration {i+1}: energy converged: {current_energy} kcal/mol")
    calculate_molecule_properties(current_coordinates)
    
    if current_energy < min_energy:
        min_energy = current_energy
        min_coordinates = current_coordinates

print("\nThe lowest energy found was", min_energy, "kcal/mol")
print(min_coordinates, "\n")

# Save optimized coordinates and energy
xyzp.save_xyz(os.path.splitext(filename)[0]+'_opt.xyz', atom_labels, min_coordinates)

with open(os.path.splitext(filename)[0]+'_opt.txt', "w") as f:
    f.write(str(min_energy))

f.close()