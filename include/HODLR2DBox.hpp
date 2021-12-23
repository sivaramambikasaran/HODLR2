#ifndef _HODLR2DBox_HPP__
#define _HODLR2DBox_HPP__

#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

// struct pts2D {
// 	double x,y;
// };
#include "definitions.hpp"

class HODLR2DBox {
public:
	int boxNumber;
	int parentNumber;
	int childrenNumbers[4];
	int neighborNumbers[4];
	int cornerNumbers[4];
	int innerNumbers[12];
	int outerNumbers[12];
	HODLR2DBox ();

	pts2D center;
	Eigen::VectorXd charges;

	std::map<int, Eigen::MatrixXd> M2L;
	Eigen::MatrixXd L2P;					//	Transfer from multipoles of 4 children to multipoles of parent.
	std::vector<pts2D> chebNodes;
	std::vector<int> chargeLocations;

	Eigen::VectorXd outgoing_charges;//equivalent densities {f_{k}^{B,o}}
	Eigen::VectorXd incoming_potential;//check potentials {u_{k}^{B,i}}
	std::vector<int> incoming_chargePoints;//equivalent points {y_{k}^{B,i}}
	std::vector<int> incoming_checkPoints;//check points {x_{k}^{B,i}}
	Eigen::VectorXd potential;//used only at leaf level
};

#endif
