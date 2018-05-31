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
	long long *I, *IU, *IL, *J, *JU, *JL, *colU, *colL, *indxL, *indxU;
	double *valCOO, *valUCOO, *valLCOO, *valUCRS, *valLCRS, *valLCCS, *valUCCS;
	long long *UpRowCCS, *LowRowCCS, *UpIndxCCS, *LowIndxCCS;
	long long *SNodesLow, *SNodesUp, NodesNLow, NodesNUp;
	double *e, *y, *xCRS, *yCRS, *xBlock, *yBlock, *xMKL, *yMKL, *xNode, *yNode;
	double *b, *bMKL, *bCRS, *bNode, *bBlock, *useless;
	long long N, M, nz, nzU, nzL;

	if (readMTX(argv[1], &I, &J, &valCOO, &M, &N, &nz) != 0)
		exit(0);
	checkAndFillDiag(&I, &J, nz, N, &valCOO);
	mallocVectors(&e, &y, &b, &xCRS, &yCRS, &bCRS, &xMKL, &yMKL,
		&bMKL, &xNode, &yNode, &bNode, &xBlock, &yBlock, &bBlock, &useless, N);
	randVector(e, N);
	cutUpperTriangleCOO(nz, I, J, valCOO, &nzL, &IL, &JL, &valLCOO);
	printf("Finished cutting coordinate matrix\n");
	transposeCOO(nzL, IL, JL, valLCOO, nzU, &IU, &JU, &valUCOO);


	printf("nzU = %d\n nzL = %d\n N = %d\n\n", nzU, nzL, N);

	/// Convert to CRS/CCS
	COOtoCRS(N, nzU, IU, JU, valUCOO, &indxU, &colU, &valUCRS);
	COOtoCRS(N, nzL, IL, JL, valLCOO, &indxL, &colL, &valLCRS);
	COOtoCCS(N, nzU, IU, JU, valUCOO, &UpIndxCCS, &UpRowCCS, &valUCCS);
	COOtoCCS(N, nzL, IL, JL, valLCOO, &LowIndxCCS, &LowRowCCS, &valLCCS);
	//saveBinCRS(argv[2], N, row, col, valCRS);
	//printf("Finished saving CRS matrix\n");

	/// Mult U*e = y /// L*y = b
	matrixMultVector(N, e, y, valUCRS, colU, indxU);
	matrixMultVector(N, y, b, valLCRS, colL, indxL);

	for (int i = 0; i < N; i++)
		bCRS[i] = bNode[i] = bBlock[i] = bMKL[i] = b[i];

	/// CRS находим y из L*y=bx /// находим х из U*x=y
	printf("gaussBackLow started\n");
	CRSTStart = clock();
	gaussBackLow(0, N, yCRS, bCRS, valLCRS, colL, indxL);
	gaussBackUp(N, 0, xCRS, yCRS, valUCRS, colU, indxU);
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

	/// NodeSolverLow & NodeSolverUp
	printf("NodeSolverLow started\n");
	NodeTStart = clock();
	nodeSolverLow(valLCRS, colL, indxL, SNodesLow, NodesNLow, yNode, bNode, N);
	nodeSolverUp(valUCRS, colU, indxU, SNodesUp, NodesNUp, xNode, yNode, N);
	NodeTFinish = clock();
	printf("NodeSolverUp finished\n");

	/// BlockSolverLow & BlockSolverUp
	printf("BlockSolverLow started\n");
	BlockTStart = clock();
	blockSolverLowCCS(valLCRS, colL, indxL, valLCCS, LowRowCCS, LowIndxCCS, SNodesLow, NodesNLow, yBlock, bBlock, N);
	blockSolverUpCCS(valUCRS, colU, indxU, valUCCS, UpRowCCS, UpIndxCCS, SNodesUp, NodesNUp, xBlock, yBlock, N);
	BlockTFinish = clock();
	printf("BlockSolverUp finished\n");

	///MKL Prepare
	const long long  MKLn = N;
	long long	MKLerror = 0;
	const long long MKLn_short = MKLn;
	long long* colU_short;
	long long* colL_short;
	long long* indxU_short;
	long long* indxL_short;
	MKLPrepare(&colU_short, &colL_short, &indxU_short, &indxL_short,
		MKLn_short, indxU, indxL, bMKL, colU, colL, b, N);

	/// MKLSolverUp &MKLSolverLow
	MKLTStart = clock();
	mkl_dcsrtrsv("L", "N", "N", &MKLn_short, valLCRS, indxL_short, colL_short, b, yMKL);
	mkl_dcsrtrsv("U", "N", "N", &MKLn_short, valUCRS, indxU_short, colU_short, yMKL, xMKL);
	MKLTFinish = clock();

	/// CheckResult
	printf("\n Size: %d x %d , nz = %d , nzL = %d , nzU = %d\n", N, N, nz, nzL, nzU);
	printf("\n\n [CRS]\n");
	long long CRSerrors = CheckSolv(N, e, xCRS);
	printf("[MKL]\n");
	long long MKLErrors = CheckSolv(N, e, xMKL);
	printf("[Node]\n");
	long long NodeErrors = CheckSolv(N, e, xNode);
	printf("[Block]\n");
	long long BlockErrors = CheckSolv(N, e, xBlock);

	double CRSAbsErrors = absError(N, e, xCRS);
	double MKLAbsErrors = absError(N, e, xMKL);
	double NodeAbsErrors = absError(N, e, xNode);
	double BlockAbsErrors = absError(N, e, xBlock);

	/// PrintResult
	printf("\n\n CalcNodesTime: %f\n\n", (CalcSNodesTFinish - CalcSNodesTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	if (CRSerrors == 0) {
		printf("\n[CRS] All is Ok! absErrors = %.10f \n Time is %f\n",
			CRSAbsErrors, (CRSTFinish - CRSTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n[CRS] Errors: %d/%d \n Time is %f\n",
			CRSerrors, N, (CRSTFinish - CRSTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	if (MKLErrors == 0) {
		printf("\n[MKL] All is Ok! absErrors = %.10f \n Time is %f\n",
			MKLAbsErrors, (MKLTFinish - MKLTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n[MKL] Errors: %d/%d \n Time is %f\n",
			MKLErrors, N, (MKLTFinish - MKLTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	if (NodeErrors == 0) {
		printf("\n[Node] All is Ok! absErrors = %.10f \n Time is %f\n",
			NodeAbsErrors, (NodeTFinish - NodeTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n[Node] Errors: %d/%d \n Time is %f\n",
			NodeErrors, N, (NodeTFinish - NodeTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	if (BlockErrors == 0) {
		printf("\n[Block] All is Ok! absErrors = %.10f \n Time is %f\n",
			BlockAbsErrors, (BlockTFinish - BlockTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n[Block] Errors: %d/%d \n Time is %f\n",
			BlockErrors, N, (BlockTFinish - BlockTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	/// FreeMem
	free(e); free(y); free(b); free(xMKL); free(yMKL); free(bMKL);
	free(xCRS);	free(yCRS);	free(bCRS);	free(xNode); free(yNode); free(bNode);
	free(xBlock); free(yBlock);	free(bBlock); free(I); free(IU); free(J);
	free(JU); free(colU); free(colL); free(valCOO); free(valUCOO); free(valUCRS);
	free(valLCRS); free(indxU);	free(indxL); free(valLCCS); free(valUCCS);
	free(UpRowCCS); free(LowRowCCS); free(UpIndxCCS); free(LowIndxCCS);

	return 0;
}