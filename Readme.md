HODLR2lib is a library consisting of fast matrix operations for HODLR2 matrices. In the current version, the operation available is matrix multiplication.

HODLR2 is a new class of hierarchical matrices used to represent dense matrices arising in the discretization of elliptic PDEs in 2D. It is an extension of the HODLR idea in 2D. When the domain is 1D, the HODLR structure is such that, the sub-domain corresponding to the sub-matrices that are compressed share only a vertex. If we try to extend the HODLR in 2D, where we again compress the off-diagonal blocks, we see that (apart from compressing the interaction between the sub-domains that share a vertex), sub-domains that share an edge are also compressed. It is observed that the rank of the interaction corresponding to edge sharing boxes grows with N, the problem size.

HODLR2 is a new class of matrix representation, where all the sub-matrices except for the following are compressed

- the sub-matrices corresponding to the interaction between the edge sharing blocks
- the sub-matrices corresponding to the self interaction (diagonal blocks)

Low-rank approximation of the appropriate blocks is obtained using ACA. The domain is subdivided based on a KDTree. The algorithm has been parallelized using OpenMP. The code is written in C++ and features an easy-to-use interface.

For more details on the usage of the library, visit the [documentation](https://hodlr2.rtfd.io) page.
