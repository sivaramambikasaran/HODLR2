.. role:: underline
    :class: underline

Welcome to HODLR2Dlib's Documentation!
**************************************

About :math:`\texttt{HODLR2Dlib}`:
==================================

HODLR2D is an extension of the HODLR idea in 2D. In 1D, the HODLR works fine since the sub-domain corresponding to the sub-matrices that are compressed share only a vertex. If we try to extend the HODLR in 2D, where we again compress the off-diagonal blocks, we see that (apart from compressing the interaction between the sub-domains that share a vertex), sub-domains that share an edge are also compressed.


:math:`\texttt{HODLR2Dlib}` is a library consisting of fast matrix operations for matrices based on HODLR2D structure. In the current version, the operation available is matrix multiplication.

Low-rank approximation of the appropriate blocks is obtained using ACA. The domain is subdivided based on a KDTree. The algorithm has been parallelized using OpenMP.

The code is written in C++ and features an easy-to-use interface, where the user provides input through a ``kernel`` object which abstracts data of the matrix through a member function ``getMatrixEntry(int i, int j)`` which returns the entry at the :math:`i^{\mathrm{th}}` row and :math:`j^{\mathrm{th}}` column of the matrix.

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

Learn more about :math:`\texttt{HODLRlib}` by visiting the

* Code Repository: http://github.com/sivaramambikasaran/HODLR
* Documentation: http://hodlrlib.rtfd.io
