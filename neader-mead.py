

import math
import numpy as np
import copy
import xyz_py as xyzp
from numpy import sqrt
from numpy import around
import random as rand

print("\nHello world")
np.set_printoptions(suppress=True)

# Global Variables
N = 2
ITERATION = 100
FILENAME = "./xyz/argon_" + str(N)
# EPS = 0.997 # kJ/mol
# ALPHA = 3.4 # Angstroms
EPS = 1
ALPHA = 1

# Import 
atom_labels, cooridnates = xyzp.load_xyz(FILENAME + ".xyz")
intial_coordinates = np.array(cooridnates).flatten()
print(intial_coordinates)
min_energy = 0
min_coordinate = []

def nelder_mead(f, x_start):
    # Step 2. Define dimension and initialize variables
    dim = len(x_start)
    prev_best = f(x_start)
    res = [[x_start, prev_best]]
    '''
    print(res)
    [[[1.0, 1.0, 1.0], 0.2273243]
    '''
    # Hyper-parameters
    step=0.1
    max_iter = 10000

    no_improv = 0
    no_improve_thr = 10e-8
    no_improv_break = 100
    alpha=1.
    gamma=2.
    rho=-0.5
    sigma=0.5

    for i in range(dim):
        x = copy.copy(x_start)
        x[i] = x[i] + step
        score = f(x)
        res.append([x, score]) 


    '''
    print(res)
    [[[1.0, 1.0, 1.0], 0.22732435670642046],
    [[1.1, 1.0, 1.0], 0.24076069582392554],
    [[1.0, 1.1, 1.0], 0.19084398750051138],
    [[1.0, 1.0, 1.1], 0.21649938733944804]]
    '''

    iters = 0
    
    while 1:
        # Step 3.1 Order the result
        res.sort(key=lambda x: x[1])
        '''print(res)
        [[[1.0, 1.1, 1.0], 0.19084398750051138],
         [[1.0, 1.0, 1.1], 0.21649938733944804],
         [[1.0, 1.0, 1.0], 0.22732435670642046],
         [[1.1, 1.0, 1.0], 0.24076069582392554]]
         '''
        best = res[0][1]
        '''print(best)
        0.19084398750051138
        '''
        
        # Step 3.2 Break after max_iter
        if max_iter and iters >= max_iter:
            return res[0]
        iters += 1

        # Step 3.3 Break after no improvement
        if best < prev_best - no_improve_thr:
            no_improv = 0
            '''
            print(prev_best)
            0.22732435670642046
            '''
            prev_best = best
            '''
            print(prev_best)
            0.19084398750051138
            '''
        else:
            no_improv += 1

        if no_improv > no_improv_break:
            return res[0]
        
        # Step 4. Determine centroid
        
        x0 = [0.] * dim
        '''
        print(x0)
        [0.0, 0.0, 0.0]
        '''
        '''
        print(len(res))
        4
        '''
        
        # loop through the top N best points
        for tup in res[:-1]:
            # print("Iterating through", tup)
            '''
            print(tup)
            [[1.0, 1.1, 1.0], 0.19084398750051138]
            '''
            for i, c in enumerate(tup[0]):
                '''
                print(i, c)
                0 1.0
                1 1.0
                2 1.0
                '''
                x0[i] += c / (len(res)-1)
                # For each iteration, you add divided by N-dim points
                '''
                # #Iterating through [[1.0, 1.1], 0.38168797500102275]
                # [0.5 0. ] 
                # [0.5  0.55] 

                # Iterating through [[1.0, 1.0], 0.4546487134128409]
                # [1.   0.55] 
                # [1.   1.05] 
                '''
        
        # Step 5. Reflection
        # print("Worst point:", res[-1][0]) # Worst point
        # print("Centroid:", x0)
        xr = x0 + alpha*(x0 - res[-1][0])
        # print("Relfect point:", xr)
        # print(res)
        # print(res[-1]) # Worst
        # print(res[-2]) # 2nd worst
        rscore = f(xr)
        if res[0][1] <= rscore < res[-2][1]:
            del res[-1] # delete the worst point
            res.append([xr, rscore])
            continue

        # Step 6. Expansion
        if rscore < res[0][1]: # if rscore is the best point
            xe = x0 + gamma*(x0 - res[-1][0])
            escore = f(xe)
            if escore < rscore:
                del res[-1] # delete the worst poin
                res.append([xe, escore])
                continue
            else:
                del res[-1]
                res.append([xr, score])
                continue
        
        # Step 7. Contraction (when rscore is bad)
        xc = x0 + rho * (x0 - res[-1][0])
        cscore = f(xc)
        if cscore < res[-1][-1]:
            del res[-1]
            res.append([xc, cscore])
            continue

        # Step 8. Reduction
        x1 = res[0][0] # best point
        nres = []
        for tup in res:
            redx = x1 + sigma * (tup[0] - x1)
            score = f(redx)
            nres.append([redx, score])
        res = nres
        
        
# # Step 1. Define a function of 3 variables
def f(x):
    #return math.sin(x[0]) * math.cos(x[1]) * (1. / (abs(x[2]) + 1))
    # return math.sin(x[0]) * math.cos(x[1])
    i_list = list(range(1,N)) # (1 to N - 1)
    j_list = list(range(1,N))
    LJ_sum = 0
    cooridnates = x.reshape(-1,3)
    for i in i_list:
        for j in j_list:
            if i <= j:
                ith = i
                jth = j+1 
                atom_i = cooridnates[ith - 1]
                atom_j = cooridnates[jth - 1]
                r = np.sqrt(np.sum((atom_i-atom_j)**2, axis=0))
                LJ = 4 * EPS * ((ALPHA/r)**12 - (ALPHA/r)**6)
                LJ_sum += LJ
    return LJ_sum

x_start = np.array(intial_coordinates)            
value = nelder_mead(f, x_start)
print(value)

LJ = f(np.array([0.26085988, 1.18855901, 0.1962963 , 0.5668863 , 0.08103074, 0.20925926]))
print(LJ)
#x_start = [1.0, 1.0, 1.0]
# x_start = np.array([1.0, 1.0])
# nelder_mead(f, x_start)



