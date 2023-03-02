
import copy
import numpy as np
# x_1, .... x_n+1

# Global variables
x_1 = [4.424, 1.234]
x_2 = [2.243, 3.231]
x_3 = [4.424, 0.534]
x = [x_1, x_2, x_3]

def f(x_list):
    y_list = []
    for x, y in x_list:
        score = 0.26 * (x**2 + y**2) - 0.48 * x * y
        y_list.append(score)
    y_list = list(filter(None, y_list))
    print("x_list:\n", x_list, "\n")
    print("y_list:\n", np.around(y_list, 3), "\n")

    result = y_list + x_list
    print(result)
    # Now you have x_list and y_list
    # Combine them and sort
'''
# STEP 1. ORDER
# f(x_1) <= f(x_2) <= ... <= f(x_n+1)
'''

# Combine [[np.array([]), 3], [np.array([]), 3], [np.array([]), 3]]
# Sort 
res_list = f(x)
# print(np.around(res_list, 3))

# Step 2. 

