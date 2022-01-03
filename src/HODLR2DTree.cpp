#include "HODLR2DTree.hpp"
// #include "../include/HODLR2DTree.hpp"

HODLR2DTree::HODLR2DTree(userkernel* K, int N, int nLevels, int TOL_POW, double* locations, std::vector<int>& boxNumbers, std::vector<int>& NumberOfParticlesInLeaves) {
  this->K					=	K;
  this->nLevels		=	nLevels;
  this->locations = locations;
  this->boxNumbers= boxNumbers;
  this->NumberOfParticlesInLeaves = NumberOfParticlesInLeaves;
  // this->nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
  // this->nParticlesInLeaf = nParticlesInLeafAlong1D*nParticlesInLeafAlong1D;
  this->TOL_POW = TOL_POW;
  nBoxesPerLevel.push_back(1);
  boxRadius.push_back(1.0);
  for (int k=1; k<=nLevels; ++k) {
    nBoxesPerLevel.push_back(4*nBoxesPerLevel[k-1]);
    boxRadius.push_back(0.5*boxRadius[k-1]);
  }
  this->smallestBoxSize	=	boxRadius[nLevels];
  K->a					=	smallestBoxSize;
  this->N					=	N;
}

void HODLR2DTree::shift_Nodes(double radius, pts2D center, std::vector<pts2D> &particle_loc) {
  for (size_t i = 0; i < Nodes.size(); i++) {
    pts2D temp;
    temp.x = Nodes[i].x*radius + center.x;
    temp.y = Nodes[i].y*radius + center.y;
    particle_loc.push_back(temp);
  }
}

void HODLR2DTree::createTree() {
  //	First create root and add to tree
  HODLR2DBox root;
  root.boxNumber		=	0;
  root.parentNumber	=	-1;
  #pragma omp parallel for
  for (int l=0; l<4; ++l) {
    root.childrenNumbers[l]	=	l;
  }
  #pragma omp parallel for
  for (int l=0; l<4; ++l) {
    root.neighborNumbers[l]	=	-1;
  }
  #pragma omp parallel for
  for (int l=0; l<12; ++l) {
    root.innerNumbers[l]	=	-1;
  }
  #pragma omp parallel for
  for (int l=0; l<12; ++l) {
    root.outerNumbers[l]	=	-1;
  }
  #pragma omp parallel for
  for (int l=0; l<4; ++l) {
    root.cornerNumbers[l]	=	-1;
  }
  std::vector<HODLR2DBox> rootLevel;
  rootLevel.push_back(root);
  tree.push_back(rootLevel);

  for (int j=1; j<=nLevels; ++j) {
    std::vector<HODLR2DBox> level;
    for (int k=0; k<nBoxesPerLevel[j]; ++k) {
      HODLR2DBox box;
      box.boxNumber		=	k;
      box.parentNumber	=	k/4;
      for (int l=0; l<4; ++l) {
        box.childrenNumbers[l]	=	4*k+l;
      }
      level.push_back(box);
    }
    tree.push_back(level);
  }
}

//	Assigns the interactions for child0 of a box
void HODLR2DTree::assign_Child0_Interaction(int j, int k) {
  int nL	=	j+1;
  int nC	=	4*k;
  int nN;

  //	Assign siblings
  {

    tree[nL][nC].neighborNumbers[1]	=	nC+1;
    tree[nL][nC].cornerNumbers[2]	=	nC+2;
    tree[nL][nC].neighborNumbers[2]	=	nC+3;
  }

  //	Assign children of parent's zeroth neighbor
  {

    nN	=	tree[j][k].neighborNumbers[0];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].cornerNumbers[1]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's first neighbor
  {
    nN	=	tree[j][k].neighborNumbers[1];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[5]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's second neighbor
  {

    nN	=	tree[j][k].neighborNumbers[2];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[6]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's third neighbor
  {

    nN	=	tree[j][k].neighborNumbers[3];
    if (nN!=-1) {
      tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].cornerNumbers[3]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[3];
    }
  }
}


//	Assigns the interactions for child1 of a box
void HODLR2DTree::assign_Child1_Interaction(int j, int k) {
  int nL	=	j+1;
  int nC	=	4*k+1;
  int nN;

  //	Assign siblings
  {

    tree[nL][nC].neighborNumbers[3]	=	nC-1;
    tree[nL][nC].neighborNumbers[2]	=	nC+1;
    tree[nL][nC].cornerNumbers[3]	=	nC+2;
  }

  //	Assign children of parent's zeroth neighbor
  {

    nN	=	tree[j][k].neighborNumbers[0];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].neighborNumbers[0]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].cornerNumbers[0]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's first neighbor
  {

    nN	=	tree[j][k].neighborNumbers[1];
    if (nN != -1) {
      tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[5]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].cornerNumbers[2]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's second neighbor
  {

    nN	=	tree[j][k].neighborNumbers[2];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[8]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's third neighbor
  {

    nN	=	tree[j][k].neighborNumbers[3];
    if (nN != -1){
      tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[9]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[9]	=	tree[j][nN].childrenNumbers[3];
    }
  }
}

//	Assigns the interactions for child2 of a box
void HODLR2DTree::assign_Child2_Interaction(int j, int k) {
  int nL	=	j+1;
  int nC	=	4*k+2;
  int nN;

  //	Assign siblings
  {

    tree[nL][nC].cornerNumbers[0]	=	nC-2;
    tree[nL][nC].neighborNumbers[0]	=	nC-1;
    tree[nL][nC].neighborNumbers[3]	=	nC+1;
  }

  //	Assign children of parent's zeroth neighbor
  {

    nN	=	tree[j][k].neighborNumbers[0];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[0]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[0]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's first neighbor
  {

    nN	=	tree[j][k].neighborNumbers[1];
    if (nN != -1) {
      tree[nL][nC].cornerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].neighborNumbers[1]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's second neighbor
  {

    nN	=	tree[j][k].neighborNumbers[2];
    if (nN != -1) {
      tree[nL][nC].cornerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[8]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's third neighbor
  {

    nN	=	tree[j][k].neighborNumbers[3];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[11]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].outerNumbers[10]	=	tree[j][nN].childrenNumbers[3];
    }
  }
}

//	Assigns the interactions for child3 of a box
void HODLR2DTree::assign_Child3_Interaction(int j, int k) {
  int nL	=	j+1;
  int nC	=	4*k+3;
  int nN;

  //	Assign siblings
  {

    tree[nL][nC].neighborNumbers[0]	=	nC-3;
    tree[nL][nC].cornerNumbers[1]	=	nC-2;
    tree[nL][nC].neighborNumbers[1]	=	nC-1;
  }

  //	Assign children of parent's zeroth neighbor
  {

    nN	=	tree[j][k].neighborNumbers[0];
    if (nN != -1) {
      tree[nL][nC].outerNumbers[1]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[2]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[1]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's first neighbor
  {

    nN	=	tree[j][k].neighborNumbers[1];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[3]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].outerNumbers[3]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].outerNumbers[4]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[4]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's second neighbor
  {

    nN	=	tree[j][k].neighborNumbers[2];
    if (nN != -1) {
      tree[nL][nC].neighborNumbers[2]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].cornerNumbers[2]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].innerNumbers[6]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[7]	=	tree[j][nN].childrenNumbers[3];
    }
  }

  //	Assign children of parent's third neighbor
  {

    nN	=	tree[j][k].neighborNumbers[3];
    if (nN != -1) {
      tree[nL][nC].innerNumbers[11]	=	tree[j][nN].childrenNumbers[0];
      tree[nL][nC].cornerNumbers[0]	=	tree[j][nN].childrenNumbers[1];
      tree[nL][nC].neighborNumbers[3]	=	tree[j][nN].childrenNumbers[2];
      tree[nL][nC].innerNumbers[10]	=	tree[j][nN].childrenNumbers[3];
    }
  }
}

//	Assigns the interactions for the children of a box
void HODLR2DTree::assign_Box_Interactions(int j, int k) {
  this->assign_Child0_Interaction(j,k);
  this->assign_Child1_Interaction(j,k);
  this->assign_Child2_Interaction(j,k);
  this->assign_Child3_Interaction(j,k);
}

//	Assigns the interactions for the children all boxes at a given level
void HODLR2DTree::assign_Level_Interactions(int j) {
  #pragma omp parallel for
  for (int k=0; k<nBoxesPerLevel[j]; ++k) {
    this->assign_Box_Interactions(j,k);
  }
}

//	Assigns the interactions for the children all boxes in the tree
void HODLR2DTree::assign_Tree_Interactions() {
  for (int j=0; j<nLevels; ++j) {
    this->assign_Level_Interactions(j);
  }
}

void HODLR2DTree::assign_Center_Location() {
  int J;
  tree[0][0].center.x	=	0.0;
  tree[0][0].center.y	=	0.0;
  for (int j=0; j<nLevels; ++j) {
    J	=	j+1;
    double shift	=	0.5*boxRadius[j];
    #pragma omp parallel for
    for (int k=0; k<nBoxesPerLevel[j]; ++k) {
      tree[J][4*k].center.x		=	tree[j][k].center.x-shift;
      tree[J][4*k+1].center.x	=	tree[j][k].center.x+shift;
      tree[J][4*k+2].center.x	=	tree[j][k].center.x+shift;
      tree[J][4*k+3].center.x	=	tree[j][k].center.x-shift;

      tree[J][4*k].center.y		=	tree[j][k].center.y-shift;
      tree[J][4*k+1].center.y	=	tree[j][k].center.y-shift;
      tree[J][4*k+2].center.y	=	tree[j][k].center.y+shift;
      tree[J][4*k+3].center.y	=	tree[j][k].center.y+shift;
    }
  }
}

void HODLR2DTree::assignLeafChargeLocations() {
  for (size_t i = 0; i < N*2; i+=2) {
    pts2D temp;
    temp.x = locations[i];
    temp.y = locations[i+1];
    gridPoints.push_back(temp);
  }
  int startIndex = 0;
  for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    int boxNum = boxNumbers[k];
    int NumParticles = NumberOfParticlesInLeaves[k];
    for (size_t i = 0; i < NumParticles; i++) {
      tree[nLevels][boxNum].chargeLocations.push_back(startIndex+i);
    }
    startIndex += NumParticles;
  }
  K->particles_X = gridPoints;//object of base class FMM_Matrix
  K->particles_Y = gridPoints;
}

void HODLR2DTree::assignBoxChargeLocations() {
  for (size_t j = nLevels-1; j >= 1; j--) {
    for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
      int J = j+1;
      for (size_t c = 0; c < 4; c++) {
        int K = 4*k+c;
        tree[j][k].chargeLocations.insert(tree[j][k].chargeLocations.end(), tree[J][K].chargeLocations.begin(), tree[J][K].chargeLocations.end());
      }
    }
  }
}

void HODLR2DTree::assignBoxCharges() {
  for (size_t j = nLevels-1; j >= 1; j--) {
    for (size_t k = 0; k < nBoxesPerLevel[j]; k++) {
      int J = j+1;
      int vecLength;
      for (size_t c = 0; c < 4; c++) {
        int K = 4*k+c;
        vecLength += tree[J][K].charges.size();
      }
      tree[j][k].charges = Eigen::VectorXd(vecLength);
      int offset = 0;
      for (size_t c = 0; c < 4; c++) {
        int K = 4*k+c;
        tree[j][k].charges.segment(offset, tree[J][K].charges.size()) = tree[J][K].charges;
        offset += tree[J][K].charges.size();
      }
    }
  }
}

void HODLR2DTree::compress() {
  for (int j=1; j<=nLevels; j++) {
    this->compressPerLevel(j);
  }
}

void HODLR2DTree::compressPerLevel(int j) { //LFR; box interactions
  for (int k=0; k<nBoxesPerLevel[j]; ++k) {
    this->compressPerBox(j, k);
  }
}

void HODLR2DTree::compressPerBox(int j, int k) {
  std::vector<int> boxA_Nodes;
  boxA_Nodes.insert(boxA_Nodes.end(), tree[j][k].chargeLocations.begin(), tree[j][k].chargeLocations.end());
  std::vector<int> IL_Nodes;//indices
  for(int cn=0; cn<4; ++cn) {
    if(tree[j][k].cornerNumbers[cn] != -1) {
      int kIL = tree[j][k].cornerNumbers[cn];
      IL_Nodes.insert(IL_Nodes.end(), tree[j][kIL].chargeLocations.begin(), tree[j][kIL].chargeLocations.end());
    }
  }
  for(int in=0; in<12; ++in) {
    if(tree[j][k].innerNumbers[in] != -1) {
      int kIL = tree[j][k].innerNumbers[in];
      IL_Nodes.insert(IL_Nodes.end(), tree[j][kIL].chargeLocations.begin(), tree[j][kIL].chargeLocations.end());
    }
  }
  for(int on=0; on<12; ++on) {
    if(tree[j][k].outerNumbers[on] != -1) {
      int kIL = tree[j][k].outerNumbers[on];
      IL_Nodes.insert(IL_Nodes.end(), tree[j][kIL].chargeLocations.begin(), tree[j][kIL].chargeLocations.end());
    }
  }
  std::vector<int> row_indices, col_indices;
  row_indices = boxA_Nodes;
  col_indices = IL_Nodes;
  tree[j][k].numberOfChargesInIL = IL_Nodes.size();
  int ComputedRank;
  std::vector<int> row_bases, col_bases;
  LowRank* LR		=	new LowRank(K, TOL_POW, row_indices, col_indices);
  LR->ACA_only_nodes(row_bases, col_bases, ComputedRank);
  if(ComputedRank > 0) {
    for (int r = 0; r < row_bases.size(); r++) {
      tree[j][k].rows_basis.push_back(boxA_Nodes[row_bases[r]]);
    }
    for (int c = 0; c < col_bases.size(); c++) {
      tree[j][k].col_basis.push_back(IL_Nodes[col_bases[c]]);
    }
  }
}

void HODLR2DTree::segmentParentPotentialToChildren(int j, int k) {
  int J = j+1;
  int offset = 0;
  for (size_t c = 0; c < 4; c++) {
    int K = 4*k+c;
    tree[J][K].potential += tree[j][k].potential.segment(offset, tree[J][K].charges.size());
    offset += tree[J][K].charges.size();
  }
}

void HODLR2DTree::fastMultiply() {
  for (int j=1; j<=nLevels; j++) {
    for (int k=0; k<nBoxesPerLevel[j]; ++k) {
      tree[j][k].potential = Eigen::VectorXd::Zero(tree[j][k].charges.size());
    }
  }
  for (int j=1; j<=nLevels; j++) {
    this->fastMultiplyPerLevel(j);
  }
}

void HODLR2DTree::fastMultiplyPerLevel(int j) { //LFR; box interactions
  for (int k=0; k<nBoxesPerLevel[j]; ++k) {
    this->fastMultiplyPerBox(j, k);
    if(j <= nLevels-1) {
      segmentParentPotentialToChildren(j, k); //(j,k)- parent identity
    }
  }
}

void HODLR2DTree::fastMultiplyPerBox(int j, int k) {
  Eigen::VectorXd chargeVector(tree[j][k].numberOfChargesInIL);
  int offset = 0;
  for(int cn=0; cn<4; ++cn) {
    if(tree[j][k].cornerNumbers[cn] != -1) {
      int kIL = tree[j][k].cornerNumbers[cn];
      chargeVector.segment(offset, tree[j][kIL].charges.size()) = tree[j][kIL].charges;
      offset += tree[j][kIL].charges.size();
    }
  }
  for(int in=0; in<12; ++in) {
    if(tree[j][k].innerNumbers[in] != -1) {
      int kIL = tree[j][k].innerNumbers[in];
      chargeVector.segment(offset, tree[j][kIL].charges.size()) = tree[j][kIL].charges;
      offset += tree[j][kIL].charges.size();
    }
  }
  for(int on=0; on<12; ++on) {
    if(tree[j][k].outerNumbers[on] != -1) {
      int kIL = tree[j][k].outerNumbers[on];
      chargeVector.segment(offset, tree[j][kIL].charges.size()) = tree[j][kIL].charges;
      offset += tree[j][kIL].charges.size();
    }
  }
  if(tree[j][k].rows_basis.size() > 0) {
    std::vector<int> IL_Nodes;//indices
    for(int cn=0; cn<4; ++cn) {
      if(tree[j][k].cornerNumbers[cn] != -1) {
        int kIL = tree[j][k].cornerNumbers[cn];
        IL_Nodes.insert(IL_Nodes.end(), tree[j][kIL].chargeLocations.begin(), tree[j][kIL].chargeLocations.end());
      }
    }
    for(int in=0; in<12; ++in) {
      if(tree[j][k].innerNumbers[in] != -1) {
        int kIL = tree[j][k].innerNumbers[in];
        IL_Nodes.insert(IL_Nodes.end(), tree[j][kIL].chargeLocations.begin(), tree[j][kIL].chargeLocations.end());
      }
    }
    for(int on=0; on<12; ++on) {
      if(tree[j][k].outerNumbers[on] != -1) {
        int kIL = tree[j][k].outerNumbers[on];
        IL_Nodes.insert(IL_Nodes.end(), tree[j][kIL].chargeLocations.begin(), tree[j][kIL].chargeLocations.end());
      }
    }
    Eigen::VectorXd temp(tree[j][k].rows_basis.size());
    for (size_t i = 0; i < tree[j][k].rows_basis.size(); i++) {
      Eigen::MatrixXd eachRowOfAr = K->getMatrix(tree[j][k].rows_basis[i], IL_Nodes);
      Eigen::MatrixXd dumb = eachRowOfAr*chargeVector;
      temp(i) = dumb(0,0);
    }
    Eigen::MatrixXd D = K->getMatrix(tree[j][k].rows_basis, tree[j][k].col_basis);
    Eigen::PartialPivLU<Eigen::MatrixXd> D_LU = Eigen::PartialPivLU<Eigen::MatrixXd>(D);
    Eigen::VectorXd temp2 = D_LU.solve(temp);
    for (size_t i = 0; i < tree[j][k].charges.size(); i++) {
      Eigen::MatrixXd  eachRowOfAc= K->getMatrix(tree[j][k].chargeLocations[i], tree[j][k].col_basis);
      Eigen::MatrixXd dumb = eachRowOfAc*temp2;
      tree[j][k].potential(i) += dumb(0,0);
    }
  }
}

void HODLR2DTree::assignLeafCharges(Eigen::VectorXd &charges) {
  int start = 0;
  for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    int boxNum = boxNumbers[k];
    tree[nLevels][boxNum].charges	=	charges.segment(start, tree[nLevels][boxNum].chargeLocations.size());
    start += tree[nLevels][boxNum].chargeLocations.size();
  }
}

void HODLR2DTree::evaluate_NearField() { // evaluating at chargeLocations
  #pragma omp parallel for
  for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    #pragma omp parallel for
    for (size_t n = 0; n < 4; n++) {
      int nn = tree[nLevels][k].neighborNumbers[n];
      if(nn != -1) {
        Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][nn].chargeLocations);
        tree[nLevels][k].potential += R*tree[nLevels][nn].charges;
      }
    }
    Eigen::MatrixXd R = K->getMatrix(tree[nLevels][k].chargeLocations, tree[nLevels][k].chargeLocations);
    tree[nLevels][k].potential += R*tree[nLevels][k].charges; //self Interaction
  }
}

Eigen::VectorXd HODLR2DTree::getMatVecProductOutput() { // evaluating at chargeLocations
  Eigen::VectorXd r(N);
  int offset = 0;
  for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    int boxNum = boxNumbers[k];
    r.segment(offset,tree[nLevels][boxNum].potential.size()) = tree[nLevels][boxNum].potential;
    offset += tree[nLevels][boxNum].potential.size();
  }
  return r;
}

double HODLR2DTree::compute_error(int nBox) { // evaluating at chargeLocations
  Eigen::VectorXd true_potential = Eigen::VectorXd::Zero(tree[nLevels][nBox].chargeLocations.size());
  for (size_t k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    Eigen::MatrixXd R = K->getMatrix(tree[nLevels][nBox].chargeLocations, tree[nLevels][k].chargeLocations);
    true_potential += R*tree[nLevels][k].charges;
  }
  Eigen::VectorXd errVec = true_potential - tree[nLevels][nBox].potential;
  double error = errVec.norm()/true_potential.norm();
  return error;
}
