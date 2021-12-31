#ifndef _HODLR2D_HPP__
#define _HODLR2D_HPP__

#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <map>

#include "kernel.hpp"
#include "ACA.hpp"
#include "HODLR2DBox.hpp"
#include "HODLR2DTree.hpp"
#include "KDTree.hpp"

class HODLR2 {
public:
	int N;
	int MinParticlesInLeaf;
	int TOL_POW;
	unsigned nLevels;
	unsigned n_Dimension;  //      Dimension.
	unsigned n_Properties;  //      Number of properties.
	double* locations;   //      Stores all the locations.
	double* properties;  //      Stores all the properties.
	double* sorted_Locations;
	double* sorted_Properties;
	userkernel* mykernel;
	HODLR2DTree* hodlr2dtree;
	HODLR2(int N, int MinParticlesInLeaf, int TOL_POW, Eigen::MatrixXd& loc);
  void assemble();
  Eigen::VectorXd computeMatVecProduct(Eigen::VectorXd inputVecUnsorted);
  void evaluateError();
  ~HODLR2();
};

#endif
