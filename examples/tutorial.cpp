#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <map>

#include "HODLR2D.hpp"

double userkernel::getMatrixEntry(const unsigned i, const unsigned j) {
	pts2D r1 = particles_X[i];
	pts2D r2 = particles_X[j];
	double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
	double R	=	sqrt(R2);
	if (R < a) {
		return R/a;
	}
	else {
		return a/R;
	}
}

// double userkernel::getMatrixEntry(const unsigned i, const unsigned j) {
//   pts2D r1 = particles_X[i];
//   pts2D r2 = particles_X[j];
//   double R2	=	(r1.x-r2.x)*(r1.x-r2.x) + (r1.y-r2.y)*(r1.y-r2.y);
//   if (R2 < 1e-10) {
//     return 0.0;
//   }
//   else if (R2 < a*a) {
//     return 0.5*R2*log(R2)/a/a;
//   }
//   else {
//     return 0.5*log(R2);
//   }
// }

int main(int argc, char* argv[]) {
  unsigned N;
  unsigned MinParticlesInLeaf;
  int TOL_POW;
  if(argc < 4)
  {
      std::cout << "All arguments weren't passed to executable!" << std::endl;
      std::cout << "Using Default Arguments:" << std::endl;
      // Size of the Matrix in consideration:
      N          = 6400;
      // Size of Matrices at leaf level:
      MinParticlesInLeaf          = 200;
      // Tolerance of problem
      TOL_POW  = 12;
  }

  else
  {
	N    =       atoi(argv[1]);  //      Number of particles.
	MinParticlesInLeaf   =       atoi(argv[2]); // minimum particles in each leaf of KD Tree
	TOL_POW = atoi(argv[3]);
}

////////////////////////////
Eigen::VectorXd b(N); //vector definition
for (size_t i = 0; i < N; i++) {
  b(i) = 2*double(rand())/double(RAND_MAX)-1;
}
double L = 1.0;
unsigned Dimension = 2;
Eigen::MatrixXd loc(N,Dimension);
for (size_t j = 0; j < N; j++) {
  for (size_t k = 0; k < Dimension; k++) {
    loc(j,k) = 2.0*double(rand())/double(RAND_MAX)-1.0;
  }
}
double start, end;
///////////////////////// HODLR2 /////////////////////////////
start	=	omp_get_wtime();
HODLR2* hodlr2 = new HODLR2(N, MinParticlesInLeaf, TOL_POW, loc, L);
end	=	omp_get_wtime();
double timeCreateTreeAFMM = end-start;
std::cout << std::endl << "HODLR2 tree creation time: " << timeCreateTreeAFMM << std::endl;

std::cout << "========================= Problem Parameters =========================" << std::endl;
std::cout << "Matrix Size                        :" << N << std::endl;
std::cout << "Leaf Size                          :" << MinParticlesInLeaf << std::endl;
std::cout << "Tolerance                          :" << pow(10,-TOL_POW) << std::endl << std::endl;
unsigned nLevels = log(N/MinParticlesInLeaf)/log(4);
// std::cout << "nLevels: " << nLevels << std::endl;

start	=	omp_get_wtime();
hodlr2->assemble();
end	=	omp_get_wtime();
double timeAssemble = end-start;
std::cout << "========================= Assembly Time =========================" << std::endl;
std::cout << "Time for assemble in HODLR2 form    :" << timeAssemble << std::endl;


// What we are doing here is explicitly generating
// the matrix from its entries
start = omp_get_wtime();
Eigen::MatrixXd B = hodlr2->hodlr2dtree->K->getMatrix(0, 0, N, N);
end   = omp_get_wtime();
double exact_time = (end - start);
std::cout << "Time for direct matrix generation  :" << exact_time << std::endl;
std::cout << "Magnitude of Speed-Up              :" << (exact_time / timeAssemble) << std::endl << std::endl;

Eigen::VectorXd outputVec;
start = omp_get_wtime();
outputVec = hodlr2->computeMatVecProduct(b);
end   = omp_get_wtime();
double timeMatVecProduct = (end - start);
// hodlr2->evaluateError();

std::cout << "========================= Matrix-Vector Multiplication =========================" << std::endl;
std::cout << "Time for MatVec in HODLR form      :" << timeMatVecProduct << std::endl;

start = omp_get_wtime();
Eigen::VectorXd bSorted(N);
for (size_t i = 0; i < N; i++) {
  bSorted(i) = b(int(hodlr2->sorted_Properties[i]));
}
Eigen::MatrixXd r_exact_Sorted = B * bSorted;
Eigen::VectorXd r_exact(N);
for (size_t i = 0; i < N; i++) {
  r_exact(int(hodlr2->sorted_Properties[i])) = r_exact_Sorted(i);
}

end   = omp_get_wtime();
exact_time = (end - start);
std::cout << "Time for direct MatVec             :" << exact_time << std::endl;
std::cout << "Magnitude of Speed-Up              :" << (exact_time / timeMatVecProduct) << std::endl;
// Computing the relative error in the solution obtained:
std::cout << "Error in the solution is           :" << (outputVec-r_exact).norm() / (r_exact.norm()) << std::endl << std::endl;
}
