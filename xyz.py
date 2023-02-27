import xyz_py as xyzp
import numpy as np
from scipy.optimize import minimize 

# Global variables
dist = 0

# LJ constants between 2 Argon (LibreText Chemistry)
eps = 0.997 # kJ/mol
alpha = 3.4 # Angstroms

def LJPotential(x):
  y = 4 * eps * ((alpha/x)**12 - (alpha/x)**6)
  return y

# Import xyz file
atom_labels, cooridnates = xyzp.load_xyz("argons.xyz")
atom_count = xyzp.count_elements(atom_labels)
print("Atom labels:\n", atom_labels, "\n") # List
print("Atom coorindates:\n", cooridnates, "\n") # np.array
print("Atom counts:\n", atom_count, "\n") # {"Ar": 2}

# Determine distance between the 2 Ar atoms
if atom_count["Ar"] == 2:
    p1 = cooridnates[0]
    p2 = cooridnates[1]
    squared_dist = np.sum((p1-p2)**2, axis=0)
    dist = np.sqrt(squared_dist)
    dist_rounded = round(dist, 2)
    print(dist_rounded, "Å is the rounded distance between the atoms. \n")
'''
Double-checked using https://www.calculatorsoup.com
(X1, Y1, Z1) = (0.26765, 2.64288, 0)
(X2, Y2, Z2) = (0.02059, -0.30971, 0)
d = 2.962908
'''

# Determine LJ for 2 Ar atoms
potential = LJPotential(dist)
potential_rounded = round(potential, 2)
print(potential_rounded, "(kJ/mol) is the LJ potential at d =", dist_rounded, "Å \n")

# How to 





