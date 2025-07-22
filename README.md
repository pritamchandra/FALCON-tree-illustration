# FALCON-tree-illustration
FALCON is the NIST standard for post quantum cryptographic digital signature scheme. At the heart of it is the novel and delicate concept of Falcon Tree. These codes are not implementations (though they may easily be made so) but step by step illustration of what a FALCON tree is and how it is formed. 

> [!NOTE]
> On a high level a FALCON tree instantiates a very fast (quasi-linear compared to cubic) Gram-Schmidt Orthogonalization of a circulant matrix. 


> [!IMPORTANT]
> Most of the methods in this repository require SAGE.

Refer to [the exposition](An-Exposition-on-FFO.pdf) in this respository or Ducas and Prest's original work [Fast Fourier Orthogonalization](https://eprint.iacr.org/2015/1014.pdf).

Hop on to the [main](fti_main.py) and set your parameters.

```n``` : the dimension of the base ring, must be a power of 2. Then the underlying ring is $\frac{\mathbb R[x]}{(x^n - 1)}$.

```precision``` : an optional global parameters that decides how many decimal places to show on the outputs. This will need to be changed in the other files as well. 


## Procedure

Genere ```a``` a $n$ length vector of real numbers. Then
```python
W = PolynomialRing(RR, 'y') # initialize the ring R[x]
p = W(a) # convert coefficients "a" to polynomial
F = Matrixify_polynomial(p) # matrixify p
```

Once you have a matrix ```F``` as above, or a any other matrix ```F``` manually created, you can compute the FALCON tree of ```F```. There are two methods.

```python
T = GSO_tree(F)
T = GS_print_tree(F)
```
Both will return the FALCON tree associated to ```F```, except that ```GSO_print_tree``` will also illustrate the tree formation iteratively by printing the relevant intermediate steps. It will also generate images of binary trees that capture the iterative construction. Refer to [this](sample-output/) folder for an example output. 

Once the tree ```T``` is created the methods ```Tree_formatting(T)``` from the [helper](helper.py) file can be used with self-explanatory parameters to format the entries of the tree, which are vectors. This may include Fourier Transforming, or changing display precision. 

The factors corresponding to ```T``` can be generated using ```Factors(T)``` again in the [helper](helper.py) file. This will generate a list of unit lower triangular matrices use product is the $L$ part of the GSO.
