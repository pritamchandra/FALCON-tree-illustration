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

# global options
pd.set_option('display.max_columns', None)

precision = 3
np.set_printoptions(legacy = '1.25', 
                    threshold = sys.maxsize,
                    precision = precision,
                    suppress = True)

# # _____________________________________________


def Print(A, round = precision, suppress_printing = False):
    ''' Parameters
            A : numpy matrix of polynomials, OR 
                a single polynomial
            suppress_printing : the function prints the string it returns
            precision : each coefficient of the polynomials show "precision"
                        number of decimal places

        Output
            a string corresponding to A using DataFrames in a well aligned format
    '''
    R = PolynomialRing(QQ, 'x')
    x = R.gens()[0]
    
    ''' rational polynomials are printed with fraction coefficients
        real polynomials are printed with coefficients with too many decimals
        (the rounding of the polynomial coefficients here happens only for 
        display)
    '''
    try:
        # "A" is a matrix of polynomials
        m, n = A.shape
        B = np.array([[Print(A[i, j], suppress_printing = True) for j in range(n)] \
                                                               for i in range(m)]) 
        C = [B[ : , 0]]
        for i in range(1, n):
            C.append([' ']*m); C.append(B[ : , i])
        C = np.array(C).T

        display_str = pd.DataFrame(C).to_string(index = False, \
                                                       header = False)
        if not suppress_printing: print(display_str)
        return display_str


    except AttributeError:
        # "A" is a single polynomial
        coeffs = np.round(A.list(), round)
        display_str = str(sum([coeffs[k]*(x**k) for k in range(coeffs.shape[0])]))
        
        if not suppress_printing: print(display_str)
        return display_str
    

def GSO_tree(A):
    ''' Parameters
            A : 2 x l numpy matrix with entries from R[x]/(x^d - 1)

        Output:
            T : binary tree, which is the falcon tree corresponding to A
                however, here it is computed in an unoptimized manner
                for the sake of illustration, by conducting GSO operations
                on polynomials (instead of optimized operations on fourier 
                coefficients).
    '''
    d = len(A[0, 0].list())
    L, B = GSO_2D(A)

    if d == 1: # base case
        return [L, []]

    L0 = GSO_tree(Matrixify_vector(B[0]))
    L1 = GSO_tree(Matrixify_vector(B[1]))

    return [L, [L0, L1]]

def GSO_print_tree(A, parent, pos):
    ''' Parameters
            A : 2 x l numpy matrix with entries from R[x]/(x^d - 1)
            parent : partial falcon tree that tracks the recursion call
            pos : nested list representing the position of the tree currently 
                    getting accessed. 

        Output:
            T : the computations are similar to that of T. But along with the 
                computation the tree-making process is also printed 
                for the sake of illustration.
    '''
    d = len(A[0, 0].list())
    L, B = GSO_2D(A)
    
    print("For the current matrix\n", Print(A), "\n")
    print("After orthogonalization, L = ", Print(L), ",")
    print("and the orthogonal matrix B is\n", Print(B), "\n")

    if parent is None: # the first call, the root of the tree is not created yet
        parent = nested_list.nestedList([np.array(L.list()), [None, None]])
    else: # every other call
        parent[pos] = nested_list.nestedList([np.array(L.list()), [None, None]])

    print("Falcon tree current form:")
    tree_visualizer.Visualize(parent, filename = "tree" + str(tree_visualizer.Visualize.counter))
    tree_visualizer.Visualize.counter += 1

    if d == 1: # base case
        print("This branch has reached the real numbers.")
        print("_"*70, "\n")
        return [np.array(L.list()), []]

    B0 = Matrixify_vector(B[0])
    B1 = Matrixify_vector(B[1])
    print("The rows of B are orthogonalized as B0 and B1 as follows.")
    print(Print(B0), "\n\n", Print(B1), "\n")
    print("B0 and B1 are input to the next run of the algorithm.")
    print("_"*70, "\n")

    # the position of the children as computed in this format 
    # [parent, [left child, right child]]
    L0 = GSO_print_tree(B0, parent, pos + [int(1), int(0)]) 
    L1 = GSO_print_tree(B1, parent, pos + [int(1), int(1)])

    return [L, [L0, L1]]

def Tree_formatting(T, invert_fourier = True, 
                       real_part = True, 
                       precision = precision):
    ''' Parameters
            T : a falcon binary tree in the format [root, [left, right]] where
                each node is a numpy array
            invert_fourier : if True, inverts the nodes of a tree which
                             are fourier transformed
            real_part : if True, keeps only the real part of the complex nodes 
                        of the tree
            precision : keeps only "precision" many decimal places in the components
                        of every node
            
        Output
            T : a falcon binary tree with the required formatting.
    '''
    L = T[0] 

    # apply the required formatting to the current root
    if invert_fourier: L = ifft(L)
    if real_part: L = L.real
    if precision: L = np.round(L, precision)

    if L.shape[0] == 1: # base case
        return [L]

    return [L, [Tree_formatting(T[1][0], invert_fourier, real_part, precision), 
                Tree_formatting(T[1][1], invert_fourier, real_part, precision)]]


def Factors(T):
    ''' Parameters
            T : falcon binary tree
                dimension of the root is n. 
                height of the tree is k.

        Output 
            factors : is a list of k lower triangular matrices of dimension 
                      n x n. Let A be the matrix whose falcon tree is T.
                      Then the L in the ldl-decomposition of A (reindexed) should
                      be factorable L = prod(factor[i]).

                      This is how T is a compact ldl-decomposition of A. 
    '''

    L = T[0]
    d = L.shape[0]

    if d == 1: # base case, dealing with numbers
        return [np.array([[   1, 0],
                          [L[0], 1]])]

    I = np.eye(d) # identity matrix
    Z = np.zeros((d, d)) # zero matrix

    # initialize factors with the first factor
    factors = [np.block([[             I, Z],
                         [circulant(L).T, I]])]

    # incorporate the factors arising from the subtrees T[1][0] and T[1][1]
    # into suitable matrices for T
    for factor1, factor2 in zip(Factors(T[1][0]), Factors(T[1][1])):
        factors.append(np.block([[factor1,       Z],
                                 [      Z, factor2]]))

    return factors

