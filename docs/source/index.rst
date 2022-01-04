.. role:: underline
    :class: underline

Welcome to HODLR2Dlib's Documentation!
**************************************

About :math:`\texttt{HODLR2Dlib}`:
==================================

HODLR2D is a new class of hierarchical matrices used to represent dense matrices arising in the discretization of elliptic PDEs in 2D. It is an extension of the HODLR idea in 2D. When the domain is 1D, the HODLR structure is such that, the sub-domain corresponding to the sub-matrices that are compressed share only a vertex. If we try to extend the HODLR in 2D, where we again compress the off-diagonal blocks, we see that (apart from compressing the interaction between the sub-domains that share a vertex), sub-domains that share an edge are also compressed. It is observed that the rank of the interaction corresponding to edge sharing boxes grows with N, the problem size.

HODLR2D is a new class of matrix representation, where all the sub-matrices except for the following are compressed

- the sub-matrices corresponding to the interaction between the edge sharing blocks
- the sub-matrices corresponding to the self interaction (diagonal blocks)

:math:`\texttt{HODLR2Dlib}` is a library consisting of fast matrix operations for HODLR2D matrices. In the current version, the operation available is matrix multiplication.

Low-rank approximation of the appropriate blocks is obtained using ACA. The domain is subdivided based on a KDTree. The algorithm has been parallelized using OpenMP.

The code is written in C++ and features an easy-to-use interface, where the user provides the following inputs:

- a ``kernel`` object which abstracts data of the matrix through a member function ``getMatrixEntry(int i, int j)`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix.

- locations of nodes in the domain through an Eigen matrix ``loc``

- the vector ``b`` to be multiplied to the matrix

The current release has the following capabilities:

- MatVecs: Obtains :math:`A x` at a cost of :math:`\mathcal{O}\left(N\right)`

Doc Contents
============
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   benchmarks

Other Links
===========

Learn more about :math:`\texttt{HODLR2Dlib}` by visiting the

* Code Repository: http://github.com/sivaramambikasaran/HODLR2
* Documentation: https://hodlr2.rtfd.io
