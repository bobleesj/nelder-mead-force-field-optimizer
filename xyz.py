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

print("Manually determining r")
x1 = 0.26765
y1 = 2.64288
z1 = 0
x2 = 0.02059
y2 = -0.30971
z2 = 0
r = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
print(r)

# Determine LJ for 2 Ar atoms
potential = LJPotential(dist)
potential_rounded = round(potential, 2)
print(potential_rounded, "(kJ/mol) is the LJ potential at d =", dist_rounded, "Å \n")

# Run more than 1 variable
def func(x):
  
  y = 4 * eps * ((alpha/x)**12 - (alpha/x)**6)
  return y

min = minimize(func, 2.5, method="nelder-mead")
# print(min)

def f(params):
    # print(params)  # <-- you'll see that params is a NumPy array
    x1, x2, y1, y2, z1, z2 = params # <-- for readability you may wish to assign names to the component variables
    # r = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    # r = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    return 4 * eps * ((alpha/np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2))**12 - (alpha/np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2))**6)

initial_guess = [0, 10, 0, 10, 0, 10]
result = minimize(f, initial_guess, method="nelder-mead")
if result.success:
    fitted_params = result.x
    x1 = fitted_params[0]
    x2 = fitted_params[1]
    y1 = fitted_params[2]
    y2 = fitted_params[3]
    z1 = fitted_params[4]
    z2 = fitted_params[5]
    print("x1:", x1)
    print("x2:", x2)
    print("y1:", y1)
    print("y2:", y2)
    print("z1:", z1)
    print("z2:", z2)

    r = np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    print(r)

    # Save the xyz file as a file. Use the library

    argon_coordinates_opt = [[x1, y1, z1], [x2, y2, z2]]

    xyzp.save_xyz("argons_opt.xyz", atom_labels, argon_coordinates_opt)
else:
    raise ValueError(result.message)

'''
def save_xyz
(
f_name: str, labels: list, coords: numpy.ndarray,
with_numbers: bool = False, verbose: bool = True, mask: list = [], atomic_numbers: bool = False)
'''

