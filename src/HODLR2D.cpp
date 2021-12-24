#include "HODLR2D.hpp"

HODLR2::HODLR2(int N, int MinParticlesInLeaf, int TOL_POW, Eigen::MatrixXd& loc, double L) {
  this->N = N;
  this->L = L;
  this->MinParticlesInLeaf = MinParticlesInLeaf;
  this->TOL_POW = TOL_POW;
  nLevels      = log(N/MinParticlesInLeaf)/log(4);
  n_Dimension  =       2;  //      Dimension.
  n_Properties =       1;  //      Number of properties.
  locations    =       new double [N*n_Dimension];   //      Stores all the locations.
  properties   =       new double [N*n_Properties];  //      Stores all the properties.
  // Generates random locations and random values of property.
  // change the locations as per your need
  unsigned count_Location =       0;
  unsigned count_Property =       0;
  for (unsigned j=0; j<N; ++j) {
          for (unsigned k=0; k<n_Dimension; ++k) {
                  locations[count_Location] = loc(j,k);//       =       2*double(rand())/double(RAND_MAX)-1;//2*(int(rand())%2)-1;//double(rand())/double(RAND_MAX);
                  ++count_Location;
          }
          for (unsigned k=0; k<n_Properties; ++k) {
                  properties[count_Property]      =       j; //indices 0,1,2,3; after sorting this will be the indices of sorted array
                  ++count_Property;
          }
  }
  //////////////////////////////// KD Tree //////////////////////////////////////////
  std::vector<std::vector<int> > boxNumbers; // boxNumbers[nLevels] contains box numbers in N ordering. ex: [0 3 1 2]
  sorted_Locations        =       new double[N*n_Dimension];
  sorted_Properties       =       new double[N*n_Properties];
  std::vector<int> NumberOfParticlesInLeaves;// contains number of particels in each leaf in N ordering
  sort_KDTree(N, n_Dimension, locations, n_Properties, properties, MinParticlesInLeaf, nLevels, sorted_Locations, sorted_Properties, boxNumbers, NumberOfParticlesInLeaves); // Creates a KDTree given the locations. This KD Tree class generates a uniform tree - all leaves are at level nLevels. Number of particles in boxes at a given level differ by atmost 1.
  std::vector<pts2D> particles_X, particles_Y;
  mykernel		=	new userkernel(particles_X, particles_Y);
  hodlr2dtree	=	new HODLR2DTree(mykernel, int(N), int(nLevels), L, TOL_POW, sorted_Locations, boxNumbers[nLevels], NumberOfParticlesInLeaves);
}

void HODLR2::assemble() {
  hodlr2dtree->createTree();
  hodlr2dtree->assign_Tree_Interactions();
  hodlr2dtree->assign_Center_Location();
  hodlr2dtree->assignLeafChargeLocations();
  hodlr2dtree->getNodes();
  hodlr2dtree->assemble_M2L();
}

Eigen::VectorXd HODLR2::computeMatVecProduct(Eigen::VectorXd inputVecUnsorted) {
  Eigen::VectorXd outputVec;
  Eigen::VectorXd inputVec(N); //converting array to Eigen vector
  for (size_t i = 0; i < N; i++) {
    inputVec(i) = inputVecUnsorted(int(sorted_Properties[i]));
  }
  hodlr2dtree->assignLeafCharges(inputVec);
  hodlr2dtree->evaluate_M2M();
  hodlr2dtree->evaluate_M2L();
  hodlr2dtree->evaluate_L2L();
  hodlr2dtree->evaluate_NearField();
  outputVec = Eigen::VectorXd(N);
  Eigen::VectorXd outputVecSorted;
  outputVecSorted = hodlr2dtree->getMatVecProductOutput();//potential sorted
  outputVec = Eigen::VectorXd(N);
  for (size_t i = 0; i < N; i++) {
    outputVec(int(sorted_Properties[i])) = outputVecSorted(i);
  }
  return outputVec;
}

void HODLR2::evaluateError() {
  srand(time(NULL));
  std::cout << std::endl << "Performing error calculation in box: " << std::endl;
  for (size_t i = 0; i < 20; i++) {
    int nBox	=	rand()%hodlr2dtree->nBoxesPerLevel[nLevels];
    std::cout << "nBox: " << nBox << "	err: " << hodlr2dtree->compute_error(nBox) << std::endl;
  }
}

HODLR2::~HODLR2() {
  delete locations;
  delete properties;
  delete sorted_Locations;
  delete sorted_Properties;
  delete hodlr2dtree;
}
