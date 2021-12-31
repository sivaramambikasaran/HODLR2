#ifndef _HODLR2DTree_HPP__
#define _HODLR2DTree_HPP__

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "definitions.hpp"
#include "kernel.hpp"
#include "ACA.hpp"
#include "HODLR2DBox.hpp"

class HODLR2DTree {
public:
	userkernel* K;
	int nLevels;			//	Number of levels in the tree.
	int N;					//	Number of particles.
	double smallestBoxSize;	//	This is L/2.0^(nLevels).

	std::vector<int> nBoxesPerLevel;			//	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;				//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<std::vector<HODLR2DBox> > tree;	//	The tree storing all the information.

	int nParticlesInLeafAlong1D;
  int nParticlesInLeaf;
  std::vector<double> Nodes1D;
	std::vector<pts2D> Nodes;
  std::vector<pts2D> gridPoints; //all particles in domain
	int TOL_POW;
	double* locations;
	std::vector<int> boxNumbers;
	std::vector<int> NumberOfParticlesInLeaves;

// public:
	HODLR2DTree(userkernel* K, int N, int nLevels, int TOL_POW, double* locations, std::vector<int>& boxNumbers, std::vector<int>& NumberOfParticlesInLeaves);
  void shift_Nodes(double radius, pts2D center, std::vector<pts2D> &particle_loc);
	void createTree();
	void assign_Child0_Interaction(int j, int k);//	Assigns the interactions for child0 of a box
	void assign_Child1_Interaction(int j, int k);//	Assigns the interactions for child1 of a box
	void assign_Child2_Interaction(int j, int k);//	Assigns the interactions for child2 of a box
	void assign_Child3_Interaction(int j, int k);//	Assigns the interactions for child3 of a box
	void assign_Box_Interactions(int j, int k);//	Assigns the interactions for the children of a box
	void assign_Level_Interactions(int j);//	Assigns the interactions for the children all boxes at a given level
	void assign_Tree_Interactions();//	Assigns the interactions for the children all boxes in the tree
	void assign_Center_Location();
	void assignLeafChargeLocations();
	void getNodes();
	void getParticlesFromChildren_incoming(int j, int k, std::vector<int>& searchNodes);
	void getNodes_incoming_box(int j, int k, int& n_rows, int& n_cols, int& ComputedRank);
	void getNodes_incoming_level(int j);
	void assemble_M2L();
	void assignLeafCharges(Eigen::VectorXd &charges);
	void evaluate_M2M();
	void evaluate_M2L();
	void evaluate_L2L();
	void evaluate_NearField(); // evaluating at chargeLocations
	double compute_error(int nBox); // evaluating at chargeLocations
	Eigen::VectorXd getMatVecProductOutput();
};

#endif
