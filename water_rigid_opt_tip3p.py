import math
import numpy as np
import copy
import xyz_py as xyzp
from numpy import sqrt
from numpy import around
from util.water_info import calculate_molecule_properties
from algorithm.nelder_mead import minimize
import random as rand
from numpy.linalg import norm
import vg
import os

'''
Use the TIP3P parameters from TransRot
'''
# from nelder_mead import minimize
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

# Global variables
K_L = 547.5 # kcal/(mol Å^-2) - spring constant for bond lengths
L = 0.9572 # Å - equlibrium bond length between O and H
K_CUB = -1.65 # Å^-1 - cubic force constant for bond lengths
K_THETA = 49.9 # kcal/(mol rad^-1) - force constant for bond angles
THETA = 104.52 # deg - equlibrium angle between O and H
A = 582002.6617 # kcal mol^-1 Å^12 - Lennard-Jones potential constant
B = 595.0550 # kcal mol^-1 Å^6 - Lennard-Jones potential constant
Q_H = 0.417	# kcal^0.5 Å^0.5 mol^-0.5 - partial charge on hydrogen atoms
Q_O = -0.834 # kcal^0.5 Å^0.5 mol^-0.5 - partial charge on oxygen atoms
K_VALUE = 332.0637133 # Coulomb's constant in kcal/mol-Å-e^2 units

# Neader Mead Parameters
ITERATION = 20
STEP_SIZE = 0.001
MAX_ITERATION = 10000
CONVERGENCE_CRITERIA = 0.00005

# Get file name from user input
while True:
    filename = input("Q1. Enter the file name Ex) ./xyz/water/tip4p_wales_n_2.xyz: ").strip()
    if not os.path.exists(filename):
        print(f"Error: File {filename} does not exist.")
    elif not filename.endswith('.xyz'):
        print(f"Error: File {filename} is not a valid xyz file.")
    else:
        break

# TIP3P force model implementation
def TIP3P_Potential(params, K_L_custom, K_THETA_custom):
    
    # Reshape coordinates
    coordinates = params.reshape(-1,3)
    # Define atoms
    O1 = coordinates[0]
    H11 = coordinates[1]
    H21 = coordinates[2]
    O2 = coordinates[3]
    H12 = coordinates[4]
    H22 = coordinates[5]
    
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
    bonds_term = K_L_custom * (first_term + second_term + third_term + fourth_term)
    angles_term = K_THETA_custom * ((angle_rad1 - THETA * math.pi/180)**2 + (angle_rad2 - THETA * math.pi/180)**2)

    # Calculate Oxygen Repulsion term
    rO = np.sqrt(np.sum((O1-O2)**2, axis=0))
    LJ_term = A/(rO**12) - B/(rO**6)
    
    # Calculate Coulombic term - OH interactions
    rO1H12 = np.sqrt(np.sum((O1-H12)**2, axis=0)) 
    rO1H22 = np.sqrt(np.sum((O1-H22)**2, axis=0))
    rO2H11 = np.sqrt(np.sum((O2-H11)**2, axis=0)) 
    rO2H21 = np.sqrt(np.sum((O2-H21)**2, axis=0))

    OH_interations = (Q_H*Q_O) * (1/rO1H12 + 1/rO1H22 + 1/rO2H11 + 1/rO2H21)

    # Calculate Coulombic term - HH interactions
    rH11H12 = np.sqrt(np.sum((H11-H12)**2, axis=0)) 
    rH21H12 = np.sqrt(np.sum((H21-H12)**2, axis=0)) 
    rH11H22 = np.sqrt(np.sum((H11-H22)**2, axis=0)) 
    rH21H22 = np.sqrt(np.sum((H21-H22)**2, axis=0))

    HH_interations = (Q_H*Q_H) * (1/rH11H12 + 1/rH21H12 + 1/rH11H22 + 1/rH21H22)
    OO_interaction = (Q_O*Q_O) * (1/rO)

    U = bonds_term + angles_term + LJ_term + K_VALUE * (OH_interations + HH_interations + OO_interaction)

    return U

# Main loop
atom_labels, cooridnates = xyzp.load_xyz(filename)
coordinates = np.array(cooridnates).flatten()
prev_energy = 0

for i in range(ITERATION):
    # Double K_L and K_THETA values per iteration
    K_L_custom = K_L * (2**i)
    K_THETA_custom = K_THETA * (2**i)

    # Minimize the TIP3P potential energy with the custom K_L and K_THETA values
    result = minimize(lambda params: TIP3P_Potential(params, K_L_custom, K_THETA_custom), coordinates, STEP_SIZE, MAX_ITERATION)
    coordinates = result[0]
    coordinates_3D = result[0].reshape(-1, 3)
    
    current_energy = result[1]

    # Check if the energy difference between consecutive iterations is below the convergence criteria
    if abs(prev_energy - current_energy) <= CONVERGENCE_CRITERIA:
        print("*Energy has converged\n")
        break

    prev_energy = current_energy
    
    # Print the results for the current iteration
    print(f"\nIteration {i+1}: K_L={K_L_custom} and K_THETA={K_THETA_custom}: energy converged: {prev_energy}, kcal/mol / {prev_energy* 4.184}, kJ/mol")
    calculate_molecule_properties(coordinates_3D)

# Save the optimized coordinates and energy to files
xyzp.save_xyz(os.path.splitext(filename)[0]+ '_tip3p_rigid_opt.xyz', atom_labels, coordinates_3D)
with open(os.path.splitext(filename)[0]+ '_tip3p_rigid_opt.txt', "w") as f:
    f.write(str(prev_energy))
print("\n")

f.close()