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
	float CRSTStart, NodeTStart, BlockTStart, MKLTStart, CalcSNodesTStart;
	float CRSTFinish, NodeTFinish, BlockTFinish, MKLTFinish, CalcSNodesTFinish;

	long long *IU, *IL, *JU, *JL, *colU, *colL, *indxL, *indxU;
	double *valUCOO, *valLCOO, *valUCRS, *valLCRS, *valLCCS, *valUCCS;
	long long *UpRowCCS, *LowRowCCS, *UpIndxCCS, *LowIndxCCS;
	long long *SNodesLow, *SNodesUp, NodesNLow, NodesNUp;
	double *xCRSUp, *xBlockUp, *xMKLUp, *xNodeUp, *xCRSLow, *xBlockLow, *xMKLLow, *xNodeLow;
	double *bMKLLow, *bCRSLow, *bNodeLow, *bBlockLow, *bMKLUp, *bCRSUp, *bNodeUp, *bBlockUp;
	long long *fullRowL, *fullRowU;

	long long N = 60000;
	long long blockSizeU = 100;
	long long blockSizeL = 100;
	long long fullness = 2;

	long long nzL = calcNzL(N, blockSizeL, &fullRowL, fullness);
	long long nzU = calcNzU(N, blockSizeU, &fullRowU, fullness);

	mallocMatrixCOO(&IL, &JL, &valLCOO, nzL);
	mallocMatrixCOO(&IU, &JU, &valUCOO, nzU);
	generateBigBlockMatrixL(IL, JL, valLCOO, nzL, N, blockSizeL, fullRowL);
	generateBigBlockMatrixU(IU, JU, valUCOO, nzU, N, blockSizeU, fullRowU);

	//mm_write_mtx_crd("mtx/M100_10_2L.mtx", N, N, nzL, IL, JL, valLCOO);
	//mm_write_mtx_crd("mtx/M100_10_2U.mtx", N, N, nzU, IU, JU, valUCOO);

	COOtoCRS(N, nzU, IU, JU, valUCOO, &indxU, &colU, &valUCRS);
	COOtoCRS(N, nzL, IL, JL, valLCOO, &indxL, &colL, &valLCRS);

	saveBinCRS("bin/M60000_100_2L.bin", N, indxL, colL, valLCRS);
	saveBinCRS("bin/M60000_100_2U.bin", N, indxU, colU, valUCRS);

	printf("nzL = %d  |  nzU = %d  |  N = %d  | blockSize = %d \n", nzL, nzU, N, blockSizeL);
	printf("done!\n");

	free(IL);
	free(JL);
	free(IU);
	free(JU);
	free(valLCOO);
	free(valUCOO);

	free(indxU);
	free(indxL);
	free(colU);
	free(colL);
	free(valUCRS);
	free(valLCRS);
	return 0;
}