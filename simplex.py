# # https://github.com/fchollet/nelder-mead/blob/master/nelder_mead.py


# import copy

# '''
#     Pure Python/Numpy implementation of the Nelder-Mead algorithm.
#     Reference: https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
# '''


# def nelder_mead(f, x_start,
#                 step=0.1, no_improve_thr=10e-6,
#                 no_improv_break=10, max_iter=3
#                 ,
#                 alpha=1., gamma=2., rho=-0.5, sigma=0.5):
#     '''
#         @param f (function): function to optimize, must return a scalar score
#             and operate over a numpy array of the same dimensions as x_start
#         @param x_start (numpy array): initial position
#         @param step (float): look-around radius in initial step
#         @no_improv_thr,  no_improv_break (float, int): break after no_improv_break iterations with
#             an improvement lower than no_improv_thr
#         @max_iter (int): always break after this number of iterations.
#             Set it to 0 to loop indefinitely.
#         @alpha, gamma, rho, sigma (floats): parameters of the algorithm
#             (see Wikipedia page for reference)
#         return: tuple (best parameter array, best score)
#     '''

#     # init
#     dim = len(x_start)
#     prev_best = f(x_start)
#     no_improv = 0
#     res = [[x_start, prev_best]]

#     print("dim:",dim)
#     print("prev_best:", prev_best)
#     print("no_improv:", no_improv)
#     print("res", res) # Prints res [[array([0., 0., 0.]), 0.0]]
#     # Notice there is one more point

#     for i in range(dim):
#         x = copy.copy(x_start)
#         # Step is 0.1 by default. It had 0.1 to all the initial values.
#         # x_start is 3-D

#         x[i] = x[i] + step
#         score = f(x) # It finds the score as 0.1 is added to each element iteratively
#         print("Score:", score)
#         res.append([x, score]) # res is an array of 3 points and 1 point
        
#     # Result
    
#         print("Result:", res)
#     # simplex iter
#     iters = 0
#     while 1:
#         # order
#         print("Result before sort", res, "\n")
#         res.sort(key=lambda x: x[1])
#         print("Result after sort",  res.sort(key=lambda x: x[1]), "\n")
#         best = res[0][1] # Accessing the 4th point
#         print("best:,", best, "\n")

#         # break after max_iter
#         if max_iter and iters >= max_iter:
#             return res[0]
#         iters += 1

#         # break after no_improv_break iterations with no improvement
#         print ('...best so far:', best)

#         # if best is lower, then make 
#         if best < prev_best - no_improve_thr:
#             no_improv = 0
#             prev_best = best
#             # prev_best is the new best 
#             # prev_best is what we just calcualted latest
#         else:
#             no_improv += 1

#         if no_improv >= no_improv_break:
#             return res[0]

#         # centroid
#         x0 = [0.] * dim # initial points? 
#         # res[:-1] does not include the last element
#         print("res[:-1]:", res[:-1], "\n") # get all x and 4th points
#         print("res:", res)
#         for tup in res[:-1]:
#             print("Print tup:", tup)
#             for i, c in enumerate(tup[0]): #tup[0] refers to the 3pts
#                 #Q. i -> index and c -> value 
#                 print("i:", i)
#                 print("c:", c)
#                 x0[i] += c / (len(res)-1) # We add the values then divided by 2:

#         # # reflection
#         # xr = x0 + alpha*(x0 - res[-1][0])
#         # rscore = f(xr)
#         # if res[0][1] <= rscore < res[-2][1]:
#         #     del res[-1]
#         #     res.append([xr, rscore])
#         #     continue

#         # # expansion
#         # if rscore < res[0][1]:
#         #     xe = x0 + gamma*(x0 - res[-1][0])
#         #     escore = f(xe)
#         #     if escore < rscore:
#         #         del res[-1]
#         #         res.append([xe, escore])
#         #         continue
#         #     else:
#         #         del res[-1]
#         #         res.append([xr, rscore])
#         #         continue

#         # # contraction
#         # xc = x0 + rho*(x0 - res[-1][0])
#         # cscore = f(xc)
#         # if cscore < res[-1][1]:
#         #     del res[-1]
#         #     res.append([xc, cscore])
#         #     continue

#         # # reduction
#         # x1 = res[0][0]
#         # nres = []
#         # for tup in res:
#         #     redx = x1 + sigma*(tup[0] - x1)
#         #     score = f(redx)
#         #     nres.append([redx, score])
#         # res = nres


# if __name__ == "__main__":
#     # test
#     import math
#     import numpy as np

#     def f(x):
#         return math.sin(x[0]) * math.cos(x[1]) * (1. / (abs(x[2]) + 1))

#     print (nelder_mead(f, np.array([1., 3., 1.])))