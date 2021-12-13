This is a 2D HODLR2D Code to perform Matrix-Vector products of the form Ax.

HODLR2D is similar to FMM2D with the difference being: in HODLR2D the corner-neighbors
are considered to be part of interaction list and the interaction between them is low rank approximated.

It is a completely algebraic version. Compressions are made via ACA.

It is applicable to all HODLR2D-able and symmetric kernels.

The kernel that defines the matrix A is to be defined in function "getMatrixEntry" of "userkernel" class in "kernel.hpp" file.

The vector x is to be defined in function "chargesFunction" of userkernel class in "kernel.hpp" file.

"ACA.hpp" file contains the ACA module.

"HODLR2DTree.hpp" file contains the algorithm.

Before running make sure Eigen and openmp library paths are specified in Makefile.

The algorithm is built on a KD Tree. It generates uniform tree - each leaf is at level nLevels and the number of particles in boxes at a given level differ by atmost 1.

It takes these inputs at run time: Number of particles in the domain, minimum number of particles in each leaf, tolerance of ACA in negative powers of 10

To run it input in terminal:

make -f Makefile2D.mk clean

make -f Makefile2D.mk

The output looks like:

./testHODLR2D 102400 16 8

nLevels: 6

Number of particles is: 102400

Time taken to create the tree is: 0.002087

Time taken to assemble is: 15.7373

Time taken to do Mat-Vec product is: 0.230789

Performing error calculation in box: 2732
err: 5.37262e-06
