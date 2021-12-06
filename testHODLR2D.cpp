#include "kernel.hpp"
#include "ACA.hpp"
#include "HODLR2DTree.hpp"

int main(int argc, char* argv[]) {
	int nLevels		=	atoi(argv[1]);
	int nParticlesInLeafAlong1D	=	atoi(argv[2]); // assuming the particles are located at tensor product chebyNodes
	int L			=	atoi(argv[3]);
	int TOL_POW = atoi(argv[4]);
	double start, end;

	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	std::vector<pts2D> particles_X, particles_Y;
	userkernel* mykernel		=	new userkernel(particles_X, particles_Y);
	HODLR2DTree<userkernel>* A	=	new HODLR2DTree<userkernel>(mykernel, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);

	A->set_Uniform_Nodes();
	// A->set_Standard_Cheb_Nodes();
	A->createTree();
	A->assign_Tree_Interactions();

	end		=	omp_get_wtime();
	double timeCreateTree	=	(end-start);

	std::cout << std::endl << "Number of particles is: " << A->N << std::endl;
	std::cout << std::endl << "Time taken to create the tree is: " << timeCreateTree << std::endl;
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();

	A->assign_Center_Location();
	A->assignLeafChargeLocations();
	A->getNodes();
	A->assemble_M2L();
	int N = A->N;
	Eigen::VectorXd b(N);
	for (size_t i = 0; i < N; i++) {
		b(i) = A->K->chargesFunction(A->gridPoints[i]);
	}
	A->assignLeafCharges(b);

	end		=	omp_get_wtime();
	double timeAssignCharges=	(end-start);
	std::cout << std::endl << "Time taken to assemble is: " << timeAssignCharges << std::endl;
	/////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////
	start	=	omp_get_wtime();
	A->evaluate_M2M();
	A->evaluate_M2L();
	A->evaluate_L2L();
	A->evaluate_NearField();
	end		=	omp_get_wtime();
	double timeMatVecProduct=	(end-start);
	std::cout << std::endl << "Time taken to do Mat-Vec product is: " << timeMatVecProduct << std::endl;
	/////////////////////////////////////////////////////////////////////////
	srand(time(NULL));
	int nBox	=	rand()%A->nBoxesPerLevel[nLevels-1];
	std::cout << std::endl << "Performing error calculation in box: " << nBox << std::endl;
	std::cout << "err: " << A->compute_error(nBox) << std::endl;
}
