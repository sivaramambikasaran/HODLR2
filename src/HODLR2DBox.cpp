// #include "../include/HODLR2DBox.hpp"
#include "HODLR2DBox.hpp"

HODLR2DBox::HODLR2DBox () {
	boxNumber		=	-1;
	parentNumber	=	-1;
	for (int l=0; l<4; ++l) {
		childrenNumbers[l]	=	-1;
	}
	for (int l=0; l<4; ++l) {
		neighborNumbers[l]	=	-1;
	}
	for (int l=0; l<4; ++l) {
		cornerNumbers[l]	=	-1;
	}
	for (int l=0; l<12; ++l) {
		innerNumbers[l]		=	-1;
	}
	for (int l=0; l<12; ++l) {
		outerNumbers[l]		=	-1;
	}
}
