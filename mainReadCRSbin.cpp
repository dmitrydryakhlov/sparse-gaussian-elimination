#pragma once
#define __DFP754_H_INCLUDED
#include "mkl_types.h"
#include <mkl.h>
#include <mkl_pardiso.h>
#include "mkl_spblas.h"
#include <stdio.h>
#include "utils.h"
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <random>

int main(int argc, char* argv[]) {
	int error = 0;
	float CRSTStart, NodeTStart, BlockTStart, MKLTStart, CalcSNodesTStart,
		BlockTStartPrl, BlockFullTStart, BlockFullPrlTStart;
	float CRSTFinish, NodeTFinish, BlockTFinish, MKLTFinish, CalcSNodesTFinish,
		BlockTFinishPrl, BlockFullTFinish, BlockFullPrlTFinish;
	long long *I, *IU, *IL, *J, *JU, *JL, *colU, *colL, *indxL, *indxU;
	double *valCOO, *valUCOO, *valLCOO, *valUCRS, *valLCRS, *valLCCS, *valUCCS;
	long long *UpRowCCS, *LowRowCCS, *UpIndxCCS, *LowIndxCCS;
	long long *SNodesLow, *SNodesUp, NodesNLow, NodesNUp;

	long long *CRSErrors, *BlockErrors, *BlockPrlErrors, *NodeErrors, *MKLErrors, *BlockFullErrors, *BlockFullPrlErrors;
	double  *CRSAbsErrors, *BlockAbsErrors, *BlockPrlAbsErrors, *BlockFullAbsErrors, *NodeAbsErrors, *MKLAbsErrors, *BlockFullPrlAbsErrors;


	double **e, **y, **xCRS, **yCRS, **xBlock, **yBlock, **xMKL, **yMKL, **xNode, **yNode;
	double **b, **bMKL, **bCRS, **bNode, **bBlock, **useless, **xBlockPrl, **yBlockPrl, **bBlockPrl,
		**xBlockFull, **yBlockFull, **bBlockFull, **xBlockFullPrl, **yBlockFullPrl, **bBlockFullPrl;
	long long N, M, nz, nzU, nzL;

	int intNzU, intN;

	if (ReadMatrixFromBinaryFile(argv[1], intNzU, intN, &indxU, &colU, &valUCRS) != 0)
		exit(1);

	printf("\nreaded!\n");
	
	printf("N = %d --- nz = %d ---  \n", intN, intNzU);

	for (int i = 0; i < intN; i++) {
		printf("%d -- ", indxU[i]);
		printf("%d -- ", colU[i]);
		printf("%lf -- \n", valUCRS[i]);
	}

	

	return 0;
}
