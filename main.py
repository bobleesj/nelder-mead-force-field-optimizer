import xyz_py as xyzp
import numpy as np
from scipy.optimize import minimize 
from numpy import sqrt
from numpy import around
import random as rand

# Do not use scientific notation
np.set_printoptions(suppress=True)

# Import XYZ file
def importFile(filename):
    atom_labels, cooridnates = xyzp.load_xyz(filename)
    print("Atom labels:\n", atom_labels, "\n") # List
    print("Atom coorindates:\n", cooridnates, "\n") # np.array
    return (atom_labels, cooridnates)

# Determine distance between the 2 Ar atoms given the coordinate
def calculateDistanceForTwoAtoms(atom_labels, coordinates):
    atom_count = xyzp.count_elements(atom_labels)
    if atom_count["Ar"] == 2:
        p1 = cooridnates[0]
        p2 = cooridnates[1]
        r = np.sqrt(np.sum((p1-p2)**2, axis=0))
        r_rounded = round(r, 2)
        print(r_rounded, "Å is the rounded distance between the unoptimized atoms. \n")
    return r

# Determine LJ given distance r, ε, α
def calculateLJPotential(r, esp, alpha):
    potential = 4 * eps * ((alpha/r)**12 - (alpha/r)**6)
    return potential

# Define ε, α
eps = 0.997 # kJ/mol
alpha = 3.4 # Angstroms
# ε, α values for Argon atoms are provided by LibreText Chemistry

# Calculate
atom_labels, cooridnates = importFile("argons.xyz")
r = calculateDistanceForTwoAtoms(atom_labels, cooridnates)
potential = calculateLJPotential(r, eps, alpha)
print(round(potential, 3), "(kJ/mol) is the LJ potential at d =", round(r, 3), "Å \n")


# LJ Potential Implementation
def LJPotential(params):
    x1, x2, y1, y2, z1, z2 = params
    x = x1 - x2
    y = y1 - y2
    z = z1 - z2
    return 4 * eps * ((alpha/sqrt((x)**2 + (y)**2 + (z)**2))**12
                      - (alpha/sqrt((x)**2 + (y)**2 + (z)**2))**6)

def minimizeNelderMead(initial_points):
    result = minimize(LJPotential, initial_points, method="nelder-mead")

    if result.success:
        fitted_params = result.x
        min_energy = result.fun
        x1 = fitted_params[0]
        x2 = fitted_params[1]
        y1 = fitted_params[2]
        y2 = fitted_params[3]
        z1 = fitted_params[4]
        z2 = fitted_params[5]
        r = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

        # Save the xyz file as a file
        coordinates_opt = [[x1, y1, z1], [x2, y2, z2]]
        xyzp.save_xyz("argons_opt.xyz", atom_labels, coordinates_opt)
        coordinates_opt_rounded = around(coordinates_opt, 3)
        return (r, coordinates_opt_rounded, min_energy)
    else:
        raise ValueError(result.message)

# Run minimization
# random_initial_position = [0, 10, 0, 10, 0, 10] # 
random_initial_points = np.random.random_sample(size = 6) * 10 # [x1, x2, y1, y2, z1, z2]
print(random_initial_points)

r, coordinates_opt_rounded, min_energy = minimizeNelderMead(random_initial_points)
print("r:", r)
print("optimized coordinates:", coordinates_opt_rounded)
print("min energy (kJ/mol):", min_energy)