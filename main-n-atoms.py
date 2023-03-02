import xyz_py as xyzp
import numpy as np
from scipy.optimize import minimize 
from numpy import sqrt
from numpy import around
import random as rand

# Generate a random set of coordinates for 3 atoms
# [x1, y1, z1, x2, y2, z2, x3, y3, z3]
random_initial_points = np.random.random_sample(size = 9) * 10 - 5 # 
x1 = random_initial_points[0]
y1 = random_initial_points[1]
z1 = random_initial_points[2]
x2 = random_initial_points[3]
y2 = random_initial_points[4]
z2 = random_initial_points[5]
x3 = random_initial_points[6]
y3 = random_initial_points[7]
z3 = random_initial_points[8]

# Do not use scientific notation
np.set_printoptions(suppress=True)

# Atom interaction counts:
# 3 atoms: 1-2, 1-3, 2-3 (3)
# 4 atoms: 1-2, 1-3, 1-4, 2-3, 2-4, 3-4 (6)
# 5 atoms: 1-2, 1-3, 1-4, 1-5, 2-3, 2-4, 2,5, 3-4, 3-5, 4-5 (10)
N = 5
i_list = list(range(1,N)) # 1 to 3 wher N = 4
j_list = list(range(1,N))
potential_count = 0

for i in i_list:
    for j in j_list:
        if i <= j:
            potential_count += 1
            i_atom = i
            j_atom = j+1 
            print("(i,j):", i_atom, j_atom)
            print("Let's calcjualte LJ \n")

            # Determine LJ Potential

print("SUMMARY:", str(N), "atoms have", str(potential_count), "interactions total! \n")