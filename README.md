This is a 2D HODLR2D Code to perform Matrix-Vector products of the form Ax.

HODLR2D is similar to FMM2D with the difference being: in HODLR2D the corner-neighbors
are considered to be part of interaction list and the interaction between them is low rank approximated.

It is a completely algebraic version. Compressions are made via ACA.

It is applicable to all HODLR2D-able and symmetric kernels.

The kernel that defines the matrix A is to be defined in "kernel.hpp" file.

The vector x is to be defined in  "kernel.hpp" file.

"ACA.hpp" file contains the ACA module.

"HODLR2DTree.hpp" file contains the algorithm.

Before running make sure Eigen and openmp library paths are specified in Makefile.

It takes these inputs at run time: number of Levels in tree, Number of particles in leaf along 1D, Half side length of the square domain, tol of ACA in powers of 10

To run it input in terminal:

make -f Makefile2D.mk clean

make -f Makefile2D.mk

./testHODLR2D 6 5 1 8

The output looks like:

Number of particles is: 102400

Time taken to create the tree is: 0.004739

Time taken to assemble is: 11.3064

Time taken to do Mat-Vec product is: 0.212237

Performing error calculation in box: 864

err: 2.15408e-05
