# Sparse rational function interpolation and fault-tolerant rational function reconstruction for parametric linear systems.

## Overview:

This repository solves systems `A(y) x = b(y)` where the entries of `A` and `b`
are polynomials in the parameters `y = (y1, ..., yn)` and returns the solution
`x_i(y) = N_i(y) / D_i(y)` as exact rational functions.

The algorithm is a black box variant, that is, it never builds the full symbolic matrix. It
evaluates the system over a prime field at a sequence of points, runs a
univariate rational function reconstruction along a fixed affine line, recovers
multivariate monomials with a fast Zippel style transposed Vandermonde solve and
then lifts the result from `Z_p` back to `Q` with integer rational
reconstruction.

The fault-tolerant version replaces the univariate reconstruction with **HFTRFR**
followed by **DFTRFR** which tolerates up to `E`
corrupt evaluations per reconstruction.

This is the working research repository. For the standalone four file
distribution of the base solver refer to
[cecm.sfu.ca/~mmonagan/code/KYsolve](https://www.cecm.sfu.ca/~mmonagan/code/KYsolve/). Please, note
that the above does not include the fault-tolerant versions. 

For a more detailed overview you can check out our paper at 
[https://www.cecm.sfu.ca/~mmonagan/papers/ICMS26.pdf]
