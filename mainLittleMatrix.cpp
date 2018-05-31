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
	long long N = 12;
	long long nzU = 36;
	long long nzL = 38;

	mallocVectors(&xCRSLow, &xCRSUp, &xBlockLow, &xBlockUp, &xMKLLow, &xMKLUp, &xNodeLow, &xNodeUp,
		&bCRSLow, &bCRSUp, &bBlockLow, &bBlockUp, &bMKLLow, &bMKLUp, &bNodeLow, &bNodeUp, N);

	mallocMatrixCOO(&IL, &JL, &valLCOO, nzL);
	mallocMatrixCOO(&IU, &JU, &valUCOO, nzU);

	makeBlockMatrix12x12LowCOORandom(IL, JL, valLCOO, nzL);
	makeBlockMatrix12x12UpCOORandom(IU, JU, valUCOO, nzU);

	COOtoCRS(N, nzU, IU, JU, valUCOO, &indxU, &colU, &valUCRS);
	COOtoCRS(N, nzL, IL, JL, valLCOO, &indxL, &colL, &valLCRS);
	COOtoCCS(N, nzU, IU, JU, valUCOO, &UpIndxCCS, &UpRowCCS, &valUCCS);
	COOtoCCS(N, nzL, IL, JL, valLCOO, &LowIndxCCS, &LowRowCCS, &valLCCS);

	for (long long i = 0; i < N; i++) {
		bNodeLow[i] = 3 * i + 1;
		bNodeUp[i] = 3 * i + 1;
		bBlockLow[i] = 3 * i + 1;
		bBlockUp[i] = 3 * i + 1;
		bMKLLow[i] = 3 * i + 1;
		bMKLUp[i] = 3 * i + 1;
		bCRSLow[i] = 3 * i + 1;
		bCRSUp[i] = 3 * i + 1;
	}

	/// CRS находим y из L*y=bx /// находим х из U*x=y
	printf("gaussBackLow started\n");
	CRSTStart = clock();
	gaussBackLow(0, N, xCRSLow, bCRSLow, valLCRS, colL, indxL);
	gaussBackUp(N, 0, xCRSUp, bCRSUp, valUCRS, colU, indxU);
	CRSTFinish = clock();
	printf("gaussBackUp finished\n");

	/// CalcSuperNodesLow & CalcSuperNodesUp
	printf("calcSuperNodes started\n");
	CalcSNodesTStart = clock();
	CalcSuperNodesLowCCS(valLCCS, LowRowCCS, LowIndxCCS, &SNodesLow, NodesNLow, nzL, N);
	CalcSuperNodesUpCCS(valUCCS, UpRowCCS, UpIndxCCS, &SNodesUp, NodesNUp, nzU, N);
	CalcSNodesTFinish = clock();
	printf("calcSuperNodes finished\n");

	/// PrintSuperNodes
	printf("\n NodesNlow / N : %d/%d \n", NodesNLow - 1, N);
	printf("\n NodesNUp / N : %d/%d \n", NodesNUp - 1, N);

	for (int i = 0; i < NodesNLow; i++) {
		printf("%d ", SNodesLow[i]);
	}printf("\n");

	for (int i = 0; i < NodesNUp; i++) {
		printf("%d ", SNodesUp[i]);
	}printf("\n");

	/// NodeSolverLow & NodeSolverUp
	printf("NodeSolverLow started\n");
	NodeTStart = clock();
	nodeSolverLow(valLCRS, colL, indxL, SNodesLow, NodesNLow, xNodeLow, bNodeLow, N);
	nodeSolverUp(valUCRS, colU, indxU, SNodesUp, NodesNUp, xNodeUp, bNodeUp, N);
	NodeTFinish = clock();
	printf("NodeSolverUp finished\n");

	/// BlockSolverLow & BlockSolverUp
	printf("BlockSolverLow started\n");
	BlockTStart = clock();
	blockSolverLowCCS(valLCRS, colL, indxL, valLCCS, LowRowCCS, LowIndxCCS, SNodesLow, NodesNLow, xBlockLow, bBlockLow, N);
	blockSolverUpCCS(valUCRS, colU, indxU, valUCCS, UpRowCCS, UpIndxCCS, SNodesUp, NodesNUp, xBlockUp, bBlockUp, N);
	BlockTFinish = clock();
	printf("BlockSolverUp finished\n");


	printf("\nNodesLow/N : %d/%d\n", NodesNLow - 1, N);
	for (int i = 0; i < NodesNLow; i++) {
		printf("%d  ", SNodesLow[i]);
	}printf("\n\n");
	printf("\nNodesUp/N : %d/%d\n", NodesNUp - 1, N);
	for (int i = 0; i < NodesNUp; i++) {
		printf("%d  ", SNodesUp[i]);
	}
	/// /////////////////////////////////////////////////////////////////////////
	/// MKL /////////////////////////////////////////////////////////////////////
	const long long  MKLn = N;
	long long	MKLerror = 0;
	const long long MKLn_short = MKLn;
	long long* colU_short;
	long long* colL_short;
	long long* indxU_short;
	long long* indxL_short;
	MKLPrepare(&colU_short, &colL_short, &indxU_short, &indxL_short,
		MKLn_short, indxU, indxL, bMKLLow, colU, colL, bMKLUp, N);
	printf("/n%d\n", N);
	MKLTStart = clock();
	mkl_dcsrtrsv("L", "N", "N", &MKLn_short, valLCRS, indxL_short, colL_short, bMKLUp, xMKLLow);
	mkl_dcsrtrsv("U", "N", "N", &MKLn_short, valUCRS, indxU_short, colU_short, bMKLUp, xMKLUp);
	MKLTFinish = clock();
	printf("/n%d\n", N);
	printf("MKL error: %lld\n", MKLerror);
	printf("\n Size: %d x %d , nzL = %d , nzU = %d\n", N, N, nzL, nzU);

	printf("\n\n [CRS]\n");
	printVectorF(xCRSLow, N);
	printVectorF(xCRSUp, N);
	printf("[Node]\n");
	printVectorF(xNodeLow, N);
	printVectorF(xNodeUp, N);
	printf("[Block]\n");
	printVectorF(xBlockLow, N);
	printVectorF(xBlockUp, N);
	printf("[MKL]\n");
	printVectorF(xMKLLow, N);
	printVectorF(xMKLUp, N);

	mallocVectors(&xCRSLow, &xCRSUp, &xBlockLow, &xBlockUp, &xMKLLow, &xMKLUp, &xNodeLow, &xNodeUp,
		&bCRSLow, &bCRSUp, &bBlockLow, &bBlockUp, &bMKLLow, &bMKLUp, &bNodeLow, &bNodeUp, N);

	free(xMKLLow); free(xMKLUp); free(bMKLLow);
	free(xCRSLow);	free(xCRSUp); free(bCRSLow); free(xNodeLow); free(xNodeUp); free(bNodeLow);
	free(xBlockLow); free(xBlockUp); free(bBlockUp); free(bCRSUp); free(bBlockLow); free(bMKLUp); free(bNodeUp);
	free(IU); free(JU); free(colU); free(colL);
	free(valUCOO); free(valUCRS);	free(valLCRS); free(indxU);	free(indxL); free(valLCCS); free(valUCCS);
	free(UpRowCCS); free(LowRowCCS); free(UpIndxCCS); free(LowIndxCCS);
	return 0;
}