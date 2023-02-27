import numpy as np
import matplotlib.pyplot as plt

# Use Pyhton 3.9.7 for Anaconda to execute

from numpy.lib.npyio import NpzFile
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize 

# LJ constants between 2 Argon
# Ref: LibreText Chemistry
eps = 0.997 # kJ/mol
alpha = 3.4 # Angstroms

# LJ function
def LJPotential(x):
  y = 4 * eps * ((alpha/x)**12 - (alpha/x)**6)
  return y

def func(x):
  y = 4 * eps * ((alpha/x)**12 - (alpha/x)**6)
  return y

# x1 = 6, x2 = 4
# Compare x1 and x2, if x1 is lower, then reflect x2, vice versa

# Reflect and no extension
def reflect(array):
  sorted_array = np.sort(array)
  x1 = sorted_array[0]
  x2 = sorted_array[1]
  if func(x1) > func(x2):
    # print("local minima on the right - reflect across x2")
    # print("y_x1 - should be higher:", func(x1))
    # print("y_x2 - should be lower:", func(x2))
    r = x2 + abs(x2 - x1)
    if func(r) < func(x2):
      # print("y_r is lower than y_x2")
      return [x2, r]
    else:
      # print("y_r_ is higher than y_x2. Do not use r")
      return [x1, x2]

  elif func(x1) < func(x2): 
    # print("local minima on the left")
    # print("local minima on the right - reflect across x1")
    # print("y_x1 - should be lower:", func(x1))
    # print("y_x2 - should be higher:", func(x2))
    r = x1 - abs(x2 - x1)
    if func(r) < func(x1):
      # print("y_r is lower than y_x1")
      return [r, x1]
    else:
      # print("y_r_ is higher than y_x2. Do not use r")
      return [x1, x2]
    
# Right local minima
# print(reflect(4.2, 4.4)) # should return (4.0, 4.2)
# print(reflect(4.4, 4.2)) # should return (4.0, 4.2)
# print(reflect(3.9, 3.8)) # should return (3.8, 3.9)
# print(reflect(3.8, 3.9)) # should return (3.8, 3.9)

# # Left local minima test
# print(reflect([3.4, 3.5])) # should return (3.5, 3.6)
# print(reflect([3.5, 3.4])) # should return (3.5, 3.6)
# print(reflect([3.7, 3.8])) # should return (3.7, 3.8) x change
# print(reflect([3.8, 3.7])) # should return (3.7, 3.8) x change
print(reflect([3.70, 3.71])) # should return (3.70, 3.71) x change

# Next: Text whether it is required to extend further
def simplex(array):
  print("Input:", array)
  sorted_array = np.sort(array)
  reflected = reflect(sorted_array)
  extended = reflect(reflected)
  x1 = round(sorted_array[0], 6)
  x2 = round(sorted_array[1], 6)
  r = round(reflected[1], 6)
  e = round(extended[1], 6)
  print("x1, x2, r, e:", x1, x2, r, e)
  # Assume local minima is on the right side

  # Reflection and extension both applied
  if x2 != r and r != e:
    # print("Reflection and extension both were applied") 
    # print("Output (r, e):", r, e)
    return [r, e]

  if x2 != r and r == e:
    # print("Only Reflection was applied ")
    # print("Output (x2, r):", x2, r)  
    return [x2, r]

  else:
    # print("No reflection was applied - close to the minimum")  
    # print("Output (x1, x2):", x1, x2)
    return [x1, x2]
  
x = 3.700
increment = 0.001
x_1 = x + increment

for i in range(100):
  values = simplex([x, x_1])
  x = values[0]
  x_1 = values[1]

# Determine minimum using the built-in Nelder-mead method
def minimizeNelderMead(x0, func):
  min = minimize(func, x0, method="nelder-mead")
  min_x = str(round(min.x[0], 3))
  min_y = str(round(min.fun, 3))
  itr = str(min.nit)
  print("A minimum is found at (" + min_x + ", " + min_y + ")"
         + " after " + itr + " iterations.")

# Plot
def plotLJPotential():
  x = np.arange(2, 10, 0.1)
  y = LJPotential(x)
  y_upperbound = 10; x_upperbound = 10
  plt.xlabel("r (Å)")
  plt.ylabel("E (kJ/mol)")
  plt.axhline(y = 0, color = 'b', linestyle = 'dashed')
  plt.ylim(-2, y_upperbound)
  plt.xlim(2.5, x_upperbound)
  plt.title("LJ Potential as a function of distance between 2 Ar atoms")
  plt.plot(x, y)
  plt.grid()
  plt.show() 

# Execute functions
plotLJPotential()
minimizeNelderMead(1, LJPotential)