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

	//if (readMTX(argv[1], &I, &J, &valCOO, &M, &N, &nz) != 0)
		//exit(0);

	if (readMTX(argv[1], &IL, &JL, &valLCOO, &M, &N, &nzL) != 0)
		exit(1);

	if (readMTX(argv[2], &IU, &JU, &valUCOO, &M, &N, &nzU) != 0)
		exit(1);

	M = 2000;
	
	//checkAndFillDiag(&I, &J, nz, N, &valCOO);
	//mallocVectors(&e, &y, &b, &xCRS, &yCRS, &bCRS, &xMKL, &yMKL,
		//&bMKL, &xNode, &yNode, &bNode, &xBlock, &yBlock, &bBlock, &useless, N);

	
	CRSErrors = (long long*)malloc((N) * sizeof(long long));
	BlockErrors = (long long*)malloc((N) * sizeof(long long));
	BlockFullErrors = (long long*)malloc((N) * sizeof(long long));
	BlockFullPrlErrors = (long long*)malloc((N) * sizeof(long long));
	BlockPrlErrors = (long long*)malloc((N) * sizeof(long long));
	NodeErrors = (long long*)malloc((N) * sizeof(long long));
	MKLErrors = (long long*)malloc((N) * sizeof(long long));

	CRSAbsErrors = (double*)malloc((N) * sizeof(double));
	BlockAbsErrors = (double*)malloc((N) * sizeof(double));
	BlockFullAbsErrors = (double*)malloc((N) * sizeof(double));
	BlockFullPrlAbsErrors = (double*)malloc((N) * sizeof(double));
	BlockPrlAbsErrors = (double*)malloc((N) * sizeof(double));
	NodeAbsErrors = (double*)malloc((N) * sizeof(double));
	MKLAbsErrors = (double*)malloc((N) * sizeof(double));

	mallocMatrixNxM(&e, N, M);
	mallocMatrixNxM(&y, N, M);
	mallocMatrixNxM(&b, N, M);
	mallocMatrixNxM(&xCRS, N, M);
	mallocMatrixNxM(&yCRS, N, M);
	mallocMatrixNxM(&bCRS, N, M);
	mallocMatrixNxM(&xMKL, N, M);
	mallocMatrixNxM(&yMKL, N, M);
	mallocMatrixNxM(&bMKL, N, M);
	mallocMatrixNxM(&xNode, N, M);
	mallocMatrixNxM(&yNode, N, M);
	mallocMatrixNxM(&bNode, N, M);
	mallocMatrixNxM(&xBlock, N, M);
	mallocMatrixNxM(&yBlock, N, M);
	mallocMatrixNxM(&bBlock, N, M);
	mallocMatrixNxM(&xBlockPrl, N, M);
	mallocMatrixNxM(&yBlockPrl, N, M);
	mallocMatrixNxM(&bBlockPrl, N, M);
	mallocMatrixNxM(&xBlockFull, N, M);
	mallocMatrixNxM(&yBlockFull, N, M);
	mallocMatrixNxM(&bBlockFull, N, M);
	mallocMatrixNxM(&xBlockFullPrl, N, M);
	mallocMatrixNxM(&yBlockFullPrl, N, M);
	mallocMatrixNxM(&bBlockFullPrl, N, M);


	for (long long i = 0; i < M; i++) {
		randVector(e[i], N);
		//randVector(bCRS[i], N);
		//randVector(bMKL[i], N);
		//randVector(bNode[i], N);
		//randVector(bBlock[i], N);
	}

	//randVector(e, N);
	//cutUpperTriangleCOO(nz, I, J, valCOO, &nzL, &IL, &JL, &valLCOO);
	//printf("Finished cutting coordinate matrix\n");
	//transposeCOO(nzL, IL, JL, valLCOO, nzU, &IU, &JU, &valUCOO);


	printf("nzU = %d\n nzL = %d\n N = %d M = %d\n\n", nzU, nzL, N, M);

	/// Convert to CRS/CCS
	COOtoCRS(N, nzU, IU, JU, valUCOO, &indxU, &colU, &valUCRS);
	COOtoCRS(N, nzL, IL, JL, valLCOO, &indxL, &colL, &valLCRS);
	COOtoCCS(N, nzU, IU, JU, valUCOO, &UpIndxCCS, &UpRowCCS, &valUCCS);
	COOtoCCS(N, nzL, IL, JL, valLCOO, &LowIndxCCS, &LowRowCCS, &valLCCS);
	//saveBinCRS(argv[2], N, row, col, valCRS);
	//printf("Finished saving CRS matrix\n");

	/// Mult U*e = y /// L*y = b
	for (long long i = 0; i < M; i++) {
		matrixMultVector(N, e[i], y[i], valUCRS, colU, indxU);
		matrixMultVector(N, y[i], b[i], valLCRS, colL, indxL);
	}
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			bCRS[i][j] = bNode[i][j] = bBlockFull[i][j] = bBlockFullPrl[i][j] = bBlockPrl[i][j] = bBlock[i][j] = bMKL[i][j] = b[i][j];

	/// CRS находим y из L*y=bx /// находим х из U*x=y
	printf("gaussBackLow started\n");
	CRSTStart = clock();
	for (int i = 0; i < M; i++) {
		gaussBackLow(0, N, yCRS[i], bCRS[i], valLCRS, colL, indxL);
		gaussBackUp(N, 0, xCRS[i], yCRS[i], valUCRS, colU, indxU);
	}
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
	//printf("NodeSolverLow started\n");
	//NodeTStart = clock();
	//for (int i = 0; i < M; i++) {
	//	nodeSolverLow(valLCRS, colL, indxL, SNodesLow, NodesNLow, yNode[i], bNode[i], N);
	//	nodeSolverUp(valUCRS, colU, indxU, SNodesUp, NodesNUp, xNode[i], yNode[i], N);
	//}
	//NodeTFinish = clock();
	//printf("NodeSolverUp finished\n");

	/// BlockSolverLow & BlockSolverUp
	printf("BlockSolverLow started\n");
	BlockTStart = clock();
	for (int i = 0; i < M; i++) {
		blockSolverLowCCS(valLCRS, colL, indxL, valLCCS, LowRowCCS, LowIndxCCS, SNodesLow, NodesNLow, yBlock[i], bBlock[i], N);
		blockSolverUpCCS(valUCRS, colU, indxU, valUCCS, UpRowCCS, UpIndxCCS, SNodesUp, NodesNUp, xBlock[i], yBlock[i], N);
	}
	BlockTFinish = clock();
	printf("BlockSolverUp finished\n");

	/// BlockSolverLowFull & BlockSolverUpFull
	printf("BlockSolverLowFull started\n");
	BlockFullTStart = clock();
	blockSolverLowFull(valLCRS, colL, indxL, valLCCS, LowRowCCS, LowIndxCCS, SNodesLow, NodesNLow, yBlockFull, bBlockFull, N, M);
	blockSolverUpFull(valUCRS, colU, indxU, valUCCS, UpRowCCS, UpIndxCCS, SNodesUp, NodesNUp, xBlockFull, yBlockFull, N, M);
	BlockFullTFinish = clock();
	printf("BlockSolverUpFull finished\n");

	/// BlockSolverLowPrl & BlockSolverUpPrl
	printf("BlockSolverLowPrl started\n");
	BlockTStartPrl = clock();
#pragma omp parallel for 
	for (int i = 0; i < M; i++) {
		blockSolverLowCCS(valLCRS, colL, indxL, valLCCS, LowRowCCS, LowIndxCCS, SNodesLow, NodesNLow, yBlockPrl[i], bBlockPrl[i], N);
		blockSolverUpCCS(valUCRS, colU, indxU, valUCCS, UpRowCCS, UpIndxCCS, SNodesUp, NodesNUp, xBlockPrl[i], yBlockPrl[i], N);
	}
	BlockTFinishPrl = clock();
	printf("BlockSolverUpPrl finished\n");

	/// BlockSolverLowFullPrl & BlockSolverUpFullPrl
	printf("BlockSolverLowFullPrl started\n");
	BlockFullPrlTStart = clock();
	blockSolverLowFullPrl(valLCRS, colL, indxL, valLCCS, LowRowCCS, LowIndxCCS, SNodesLow, NodesNLow, yBlockFullPrl, bBlockFullPrl, N, M);
	blockSolverUpFullPrl(valUCRS, colU, indxU, valUCCS, UpRowCCS, UpIndxCCS, SNodesUp, NodesNUp, xBlockFullPrl, yBlockFullPrl, N, M);
	BlockFullPrlTFinish = clock();
	printf("BlockSolverUpFullPrl finished\n");


	///MKL Prepare
	/*const long long  MKLn = N;
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
	for (int i = 0; i < N; i++) {
		mkl_dcsrtrsv("L", "N", "N", &MKLn_short, valLCRS, indxL_short, colL_short, b, yMKL);
		mkl_dcsrtrsv("U", "N", "N", &MKLn_short, valUCRS, indxU_short, colU_short, yMKL, xMKL);
	}
	MKLTFinish = clock();*/

	/// CheckResult
	double CRSAbsSumError = 0.0f, MKLAbsSumError = 0.0f, NodeAbsSumError = 0.0f, BlockAbsSumError = 0.0f,
		BlockPrlAbsSumError = 0.0f, BlockFullAbsSumError = 0.0f, BlockFullPrlAbsSumError = 0.0f;
	long long CRSSumErrors = 0, MKLSumErrors = 0, NodeSumErrors = 0, BlockSumErrors = 0,
		BlockPrlSumErrors = 0, BlockFullSumErrors = 0, BlockFullPrlSumErrors = 0;
	printf("\n Size: %d x %d , nzL = %d , nzU = %d\n", N, N, nzL, nzU);
	printf("\n\n [CRS]\n");
	for (int i = 0; i < M; i++) {
		CRSErrors[i] = CheckSolv(N, e[i], xCRS[i]);
		CRSSumErrors += CRSErrors[i];
	}
	//printf("[MKL]\n");
	//for (int i = 0; i < N; i++) {
	//	MKLErrors[i] = CheckSolv(N, e[i], xMKL[i]);
	// MKLSumErrors += MKLErrors[i];
	//}
	//printf("[Node]\n");
	//for (int i = 0; i < N; i++) {
	//	NodeErrors[i] = CheckSolv(N, e[i], xNode[i]);
	//	NodeSumErrors += NodeErrors[i];
	//}
	printf("[Block]\n");
	for (int i = 0; i < M; i++) {
		BlockErrors[i] = CheckSolv(N, e[i], xBlock[i]);
		BlockSumErrors += BlockErrors[i];
	}
	printf("[BlockPrl]\n");
	for (int i = 0; i < M; i++) {
		BlockPrlErrors[i] = CheckSolv(N, e[i], xBlockPrl[i]);
		BlockPrlSumErrors += BlockErrors[i];
	}
	printf("[BlockFull]\n");
	for (int i = 0; i < M; i++) {
		BlockFullErrors[i] = CheckSolv(N, e[i], xBlockFull[i]);
		BlockFullSumErrors += BlockFullErrors[i];
	}
	printf("[BlockFullPrl]\n");
	for (int i = 0; i < M; i++) {
		BlockFullPrlErrors[i] = CheckSolv(N, e[i], xBlockFullPrl[i]);
		BlockFullPrlSumErrors += BlockFullPrlErrors[i];
	}

	for (int i = 0; i < M; i++) {
		CRSAbsErrors[i] = absError(N, e[i], xCRS[i]);
		CRSAbsSumError += CRSAbsErrors[i];
		//MKLAbsErrors[i] = absError(N, e[i], xMKL[i]);
		//MKLAbsSumError += MKLAbsErrors[i];
		//NodeAbsErrors[i] = absError(N, e[i], xNode[i]);
		//NodeAbsSumError += NodeAbsErrors[i];
		BlockAbsErrors[i] = absError(N, e[i], xBlock[i]);
		BlockAbsSumError += BlockAbsErrors[i];
		BlockPrlAbsErrors[i] = absError(N, e[i], xBlockPrl[i]);
		BlockPrlAbsSumError += BlockPrlAbsErrors[i];
		BlockFullAbsErrors[i] = absError(N, e[i], xBlockFull[i]);
		BlockFullAbsSumError += BlockFullAbsErrors[i];
		BlockFullPrlAbsErrors[i] = absError(N, e[i], xBlockFullPrl[i]);
		BlockFullPrlAbsSumError += BlockFullPrlAbsErrors[i];
	}

	/// PrintResult
	printf("\n\n CalcNodesTime: %f\n\n", (CalcSNodesTFinish - CalcSNodesTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	if (CRSSumErrors == 0) {
		printf("\n[CRS] All is Ok! absErrors = %.10f \n Time is %f\n",
			CRSAbsSumError, (CRSTFinish - CRSTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n[CRS] Errors: %d/%d \n Time is %f\n",
			CRSSumErrors, N, (CRSTFinish - CRSTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	//if (MKLSumErrors == 0) {
	//	printf("\n[MKL] All is Ok! absErrors = %.10f \n Time is %f\n",
	//		MKLAbsSumError, (MKLTFinish - MKLTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	//}
	//else {
	//	printf("\n[MKL] Errors: %d/%d \n Time is %f\n",
	//		MKLSumErrors, N, (MKLTFinish - MKLTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	//}

	//if (NodeSumErrors == 0) {
	//	printf("\n[Node] All is Ok! absErrors = %.10f \n Time is %f\n",
	//		NodeAbsSumError, (NodeTFinish - NodeTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	//}
	//else {
	//	printf("\n[Node] Errors: %d/%d \n Time is %f\n",
	//		NodeSumErrors, N, (NodeTFinish - NodeTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	//}

	if (BlockSumErrors == 0) {
		printf("\n[Block] All is Ok! absErrors = %.10f \n Time is %f\n",
			BlockAbsSumError, (BlockTFinish - BlockTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n[Block] Errors: %d/%d \n Time is %f\n",
			BlockSumErrors, N, (BlockTFinish - BlockTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	if (BlockPrlSumErrors == 0) {
		printf("\n[BlockPrl] All is Ok! absErrors = %.10f \n Time is %f\n",
			BlockPrlAbsSumError, (BlockTFinishPrl - BlockTStartPrl) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n[BlockPrl] Errors: %d/%d \n Time is %f\n",
			BlockPrlSumErrors, N, (BlockTFinishPrl - BlockTStartPrl) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	if (BlockFullSumErrors == 0) {
		printf("\n[BlockFull] All is Ok! absErrors = %.10f \n Time is %f\n",
			BlockFullAbsSumError, (BlockFullTFinish - BlockFullTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n[BlockFull] Errors: %d/%d \n Time is %f\n",
			BlockFullSumErrors, N, (BlockFullTFinish - BlockFullTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	if (BlockFullPrlSumErrors == 0) {
		printf("\n[BlockFullPrl] All is Ok! absErrors = %.10f \n Time is %f\n",
			BlockFullPrlAbsSumError, (BlockFullPrlTFinish - BlockFullPrlTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n[BlockFullPrl] Errors: %d/%d \n Time is %f\n",
			BlockFullPrlSumErrors, N, (BlockFullPrlTFinish - BlockFullPrlTStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	/// FreeMem
	free(CRSErrors);
	free(BlockErrors);
	free(BlockFullErrors);
	free(BlockFullPrlErrors);
	free(NodeErrors);
	free(MKLErrors);

	free(CRSAbsErrors);
	free(BlockAbsErrors);
	free(BlockFullAbsErrors);
	free(BlockFullPrlAbsErrors);
	free(NodeAbsErrors);
	free(MKLAbsErrors);

	free(e); free(y); free(b);
	free(xMKL); free(yMKL); free(bMKL);
	free(xCRS);	free(yCRS);	free(bCRS);
	free(xNode); free(yNode); free(bNode);
	free(xBlock); free(yBlock);	free(bBlock);
	free(xBlockFull); free(yBlockFull);	free(bBlockFull);
	free(xBlockFullPrl); free(yBlockFullPrl); free(bBlockFullPrl);

	free(IU); free(JU); free(colU);
	free(IL); free(JL); free(colL);
	free(valUCOO); free(valLCOO); free(valUCRS); free(valUCCS);
	free(valLCRS); free(valLCCS); free(indxU);	free(indxL);
	free(UpRowCCS); free(LowRowCCS); free(UpIndxCCS); free(LowIndxCCS);

	return 0;
}
