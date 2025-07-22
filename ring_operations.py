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


# global options
pd.set_option('display.max_columns', None)

precision = 3
np.set_printoptions(legacy = '1.25', 
                    threshold = sys.maxsize,
                    precision = precision,
                    suppress = True)

# # _____________________________________________

def Conjugate(f):
    ''' Parameters
            f : polynomial in R_d

        Output
            f' : polynomial in R_d. f' is the conjugate of f, that is, 
                 for any root g of x^d - 1, it holds that f(g) is the 
                 complex conjugate of f'(g)
                 
                 if f = [f_0, f1, ..., f_(d - 1)] it turns out that
                    f' = [f_0, f_(d - 1), ..., f_1]
    '''

    if f == 0:
        return f # conjugate of 0

    if f.lift().degree() <= 1:
        return f # conjugate of a constant

    W = f.parent()
    F = f.list()
    F_conj = [F[0]] + list(reversed(F[1: ]))
    return W(F_conj) 

    
def Inner_product(F, G):
    ''' Parameters
            F, G : numpy vectors of same length with entries from R_d

        Output
            <F, G> : the ring inner product that is the sum(f, g').
    '''
    G_conj = np.array([Conjugate(g) for g in G])
    return np.dot(F, G_conj)


def Vectorize(f):
    ''' Parameters
            f : a polynomial in R_d
        
        Output
            f0, f1 : polynomials in R_(d/2) which are obtained by splitting f in a 
                     fourier fashion.
                     f = f0(x^2) + x f1(x^2)
    ''' 

    F = np.array(f.list())
    d = len(F)

    W = PolynomialRing(RR, 'x')
    x = W.gens()[0]

    R = W.quotient(x**(d//2) - 1)

    # R.<x> = QuotientRing(W, W.ideal(y^(d//2) - 1))
    
    f0 = R(F[[2*j for j in range(d//2)]]) # even terms
    f1 = R(F[[2*j + 1 for j in range(d//2)]]) # odd terms

    return f0, f1

def Matrixify_polynomial(f):
    ''' Parameters
            f : polynomial in R_d

        Output
            Mf : a numpy matrix in R_(d/2)^(2 x 2) such that it is the multiplication 
                 operator for f over R_(d/2). 
    '''

    f0, f1 = Vectorize(f)
    x = f0.parent().gens()[0]

    # one can calculate the Mf has the following form
    return np.array([[  f0, f1],
                     [x*f1, f0]])

def Matrixify_vector(F):
    ''' Parameters 
            F : a numpy vector of polynomials. F is in R_d^m

        Output
            MF : a numpy matrix in R_d^(2 x m) where the matrixify polynomial is 
                 applied to every coordinate.
    ''' 

    MF = []
    for f in F:
        MF.append(Matrixify_polynomial(f))

    return np.hstack(MF)

def GSO_2D(A):
    ''' Parameters
            A : 2 x 2 numpy array over a field of polynomials

        Output
            L : polynomials from the relevant field.
            B : 2 x 2 orthogonal matrix over the relevant field.
                       A = [1 0] @ B 
                           [L 1] 
                is the GSO of A.

        NOTE: for a row matrix [ - f - ], its GSO is given by
                               [ - g - ]
    
                [ - f - ] = [          1    0] [   f               ]
                [ - g - ]   [<g, f>/<f,f>    1] [g - (<g, f>/<f,f>)*f]

        NOTE: the inner product is conjugate linear, so <g, f> != <f, g>
    '''
    L = Inner_product(A[1], A[0]) / Inner_product(A[0], A[0])
    B = np.array([A[0], A[1] - L * A[0]])
    return L, B