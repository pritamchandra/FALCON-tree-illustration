# # __________________________________________ IMPORTS

# Library imports
from sage.all import *
from scipy.linalg import solve, circulant, ldl
from numpy.linalg import inv
from numpy.fft import ifft
from copy import deepcopy

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# file imports
from importlib import reload
import tree_visualizer; reload(tree_visualizer)
import nested_list; reload(nested_list)
from ring_operations import *
from helper import *

# global options
pd.set_option('display.max_columns', None)

precision = 3
np.set_printoptions(legacy = '1.25', 
                    threshold = sys.maxsize,
                    precision = precision,
                    suppress = True)

# # _____________________________________________

# # __________________________________________ Parameters

n = 16 # must be a power of 2

# # ___________________________________________
# W.<y> = RR[] # field of real polynomials in the variable y
# the above syntax doesn't run in python
W = PolynomialRing(RR, 'y')

# # base vector 
a = list(np.random.randint(1, 10, n))
# a = list(range(0, n))
# a = [(-1)**i for i in range(n)]

p = W(a) # convert coefficients "a" to polynomial
F = Matrixify_polynomial(p)

tree_visualizer.Visualize.counter = 0
GSO_print_tree(F, None, [])
print()