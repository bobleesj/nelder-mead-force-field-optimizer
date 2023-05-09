# Author: Sangjoon (Bob) Lee
# Spring 2023 CH393
# Research Problem III

import numpy as np
import copy

# Set the print format for numpy arrays
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

'''
I acknowledge that the original code was modified from the Nelder-Mead implementation found at this GitHub repository:
https://github.com/fchollet/nelder-mead/blob/master/nelder_mead.py

Additionally, the understanding of the Nelder-Mead method and its explanation is based on the information available on the Wikipedia page:
https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
'''

# Hyper-parameters
NO_IMPROVE_THR = 10e-8  # Threshold for considering an improvement in function value
NO_IMPROV_BREAK = 100  # Number of iterations without improvement before the algorithm terminates
ALPHA = 1.0  # Reflection coefficient
GAMMA = 2.0  # Expansion coefficient
RHO = -0.5  # Contraction coefficient
SIGMA = 0.2  # Reduction coefficient

def minimize(f, x_start, STEP=0.1, MAX_ITER=10000):
    """
    Nelder-Mead optimization algorithm.

    Args:
        f (function): Objective function to optimize.
        x_start (list): Initial point for optimization.
        step (float): Initial step size for creating the simplex.
        max_iter (int): Maximum number of iterations.

    Returns:
        list: Coordinates and function value of the best point found.
    """
    def reflection(x0, worst, alpha):
        return x0 + alpha * (x0 - worst)

    def expansion(x0, worst, gamma):
        return x0 + gamma * (x0 - worst)

    def contraction(x0, worst, rho):
        return x0 + rho * (x0 - worst)

    def reduction(best, simplex, sigma):
        return [best + sigma * (point - best) for point in simplex]

    no_improv = 0  # Counter for iterations with no improvement
    # Initialization
    dim = len(x_start)  # Determine the dimension of the problem
    prev_best = f(x_start)  # Calculate the initial function value
    res = [[x_start, prev_best]]  # Store the initial point and its function value in the results list

    # Create initial simplex
    for i in range(dim):
        x = copy.copy(x_start)  # Create a new point from the starting point
        x[i] = x[i] + STEP  # Modify the ith coordinate with the step value
        score = f(x)  # Calculate the function value for the new point
        res.append([x, score])  # Add the new point and its function value to the results list

    iters = 0  # Initialize the iteration counter

    # Main loop
    while True:
        # Sort the simplex points by function value (ascending order)
        res.sort(key=lambda x: x[1])
        best = res[0][1]  # Store the best function value

        # Break conditions
        if MAX_ITER and iters >= MAX_ITER:
            return res[0]  # Return the best point and its function value
        iters += 1

        # Check for improvement in the best function value
        if best < prev_best - NO_IMPROVE_THR:
            no_improv = 0  # Reset the no improvement counter
            prev_best = best  # Update the previous best function value
        else:
            no_improv += 1  # Increment the no improvement counter

        # Terminate the algorithm if no improvement for a specified number of iterations
        if no_improv > NO_IMPROV_BREAK:
            return res[0]

        # Calculate centroid of the current simplex (excluding the worst point)
        x0 = [0.0] * dim
        for tup in res[:-1]:
            for i, c in enumerate(tup[0]):
                x0[i] += c / (len(res) - 1)

        # REFLECTION
        xr = reflection(x0, res[-1][0], ALPHA) # Calculate the reflection point
        rscore = f(xr)  # Calculate the function value at the reflection point

        # If the reflected point's score is better than the worst point's score,
        # but worse than the best point's score, replace the worst point with the reflected point
        if res[0][1] <= rscore < res[-2][1]:
            del res[-1]  # Remove the worst point
            res.append([xr, rscore])  # Add the reflected point and its score to the results list
            continue

        # EXPANSION
        # If the reflected point's score is better than the best point, perform expansion
        if rscore < res[0][1]:
            xe = expansion(x0, res[-1][0], GAMMA)  # Calculate the expansion point
            escore = f(xe)  # Calculate the function value at the expansion point

            # If the expansion point's score is better than the reflected point's score,
            # replace the worst point with the expansion point
            if escore < rscore:
                del res[-1]  # Remove the worst point
                res.append([xe, escore])  # Add the expansion point and its score to the results list
                continue
            else:
                del res[-1]  # Remove the worst point
                res.append([xr, rscore])  # Add the reflected point and its score to the results list
                continue

        # CONTRACTION
        # If the reflected point's score is worse than the second worst point, perform contraction
        xc = contraction(x0, res[-1][0], RHO)  # Calculate the contraction point
        cscore = f(xc)  # Calculate the function value at the contraction point

        # If the contraction point's score is better than the worst point's score,
        # replace the worst point with the contraction point
        if cscore < res[-1][1]:
            del res[-1]  # Remove the worst point
            res.append([xc, cscore])  # Add the contraction point and its score to the results list
            continue

        # REDUCTION
        # If the contraction point's score is still worse than the worst point's score, perform reduction
        x1 = res[0][0]  # Get the best point
        red_points = reduction(x1, [point for point, score in res[1:]], SIGMA)  # Calculate the reduction points
        red_scores = [f(point) for point in red_points]  # Calculate the function values at the reduction points
        res = [[x1, best]] + list(zip(red_points, red_scores))  # Update the results list with the new simplex points2