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
	double **xCRSUp, **xBlockUp, **xMKLUp, **xNodeUp, **xCRSLow, **xBlockLow, **xMKLLow, **xNodeLow;
	double **bMKLLow, **bCRSLow, **bNodeLow, **bBlockLow, **bMKLUp, **bCRSUp, **bNodeUp, **bBlockUp;
	double **bBlockFullLow, **bBlockFullUp, **xBlockFullUp, **xBlockFullLow;
	long long *fullRowL, *fullRowU;

	long long N = 5000, M;
	long long blockSizeU = 30;
	long long blockSizeL = 30;

	/*mallocVectors(&xCRSLow, &xCRSUp, &xBlockLow, &xBlockUp, &xMKLLow, &xMKLUp, &xNodeLow, &xNodeUp,
		&bCRSLow, &bCRSUp, &bBlockLow, &bBlockUp, &bMKLLow, &bMKLUp, &bNodeLow, &bNodeUp, N);*/
	mallocMatrixNxN(&xCRSLow, N);
	mallocMatrixNxN(&xCRSUp, N);
	mallocMatrixNxN(&xBlockLow, N);
	mallocMatrixNxN(&xBlockUp, N);
	mallocMatrixNxN(&xMKLLow, N);
	mallocMatrixNxN(&xMKLUp, N);
	mallocMatrixNxN(&xNodeLow, N);
	mallocMatrixNxN(&xNodeUp, N);
	mallocMatrixNxN(&xBlockFullLow, N);
	mallocMatrixNxN(&xBlockFullUp, N);

	mallocMatrixNxN(&bCRSLow, N);
	mallocMatrixNxN(&bCRSUp, N);
	mallocMatrixNxN(&bBlockLow, N);
	mallocMatrixNxN(&bBlockUp, N);
	mallocMatrixNxN(&bBlockFullLow, N);
	mallocMatrixNxN(&bBlockFullUp, N);
	mallocMatrixNxN(&bMKLLow, N);
	mallocMatrixNxN(&bMKLUp, N);
	mallocMatrixNxN(&bNodeLow, N);
	mallocMatrixNxN(&bNodeUp, N);

	long long nzL = calcNzL(N, blockSizeL, &fullRowL);
	long long nzU = calcNzU(N, blockSizeU, &fullRowU);

	mallocMatrixCOO(&IL, &JL, &valLCOO, nzL);
	mallocMatrixCOO(&IU, &JU, &valUCOO, nzU);
	generateBigBlockMatrixL(IL, JL, valLCOO, nzL, N, blockSizeL, fullRowL);
	generateBigBlockMatrixU(IU, JU, valUCOO, nzU, N, blockSizeU, fullRowU);

	//makeBlockMatrix12x12LowCOORandom(IL, JL, valLCOO, nzL);
	//makeBlockMatrix12x12UpCOORandom(IU, JU, valUCOO, nzU);

	mm_write_mtx_crd("myMatrix5_5L.mtx", N, N, nzL, IL, JL, valLCOO);
	mm_write_mtx_crd("myMatrix5_5U.mtx", N, N, nzU, IU, JU, valUCOO);

	free(IL);
	free(JL);
	free(IU);
	free(JU);
	free(valLCOO);
	free(valUCOO);

	if (readMTX("myMatrix5_5L.mtx", &IL, &JL, &valLCOO, &N, &M, &nzL) != 0)
		exit(0);
	if (readMTX("myMatrix5_5U.mtx", &IU, &JU, &valUCOO, &N, &M, &nzU) != 0)
		exit(0);

	COOtoCRS(N, nzU, IU, JU, valUCOO, &indxU, &colU, &valUCRS);
	COOtoCRS(N, nzL, IL, JL, valLCOO, &indxL, &colL, &valLCRS);
	COOtoCCS(N, nzU, IU, JU, valUCOO, &UpIndxCCS, &UpRowCCS, &valUCCS);
	COOtoCCS(N, nzL, IL, JL, valLCOO, &LowIndxCCS, &LowRowCCS, &valLCCS);

	for (long long i = 0; i < N; i++) {
		randVector(bNodeLow[i], N);
		/*randVector(bNodeUp[i], N);
		randVector(bBlockLow[i], N);
		randVector(bBlockUp[i], N);
		randVector(bBlockFullLow[i], N);
		randVector(bBlockFullUp[i], N);
		randVector(bMKLLow[i], N);
		randVector(bMKLUp[i], N);
		randVector(bCRSLow[i], N);
		randVector(bCRSUp[i], N);*/
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			bNodeUp[i][j] = bBlockLow[i][j] = bBlockUp[i][j] = bBlockFullLow[i][j] = bBlockFullUp[i][j] = bCRSLow[i][j] = bCRSUp[i][j] = bNodeLow[i][j];
		}
	}

	/// CRS находим y из L*y=bx /// находим х из U*x=y
	printf("gaussBackLow started\n");
	CRSTStart = clock();
	for (int i = 0; i < N; i++) {
		gaussBackLow(0, N, xCRSLow[i], bCRSLow[i], valLCRS, colL, indxL);
		gaussBackUp(N, 0, xCRSUp[i], bCRSUp[i], valUCRS, colU, indxU);
	}
	CRSTFinish = clock();
	printf("gaussBackUp finished\n");

	/// CalcSuperNodesLow & CalcSuperNodesUp
	printf("calcSuperNodes started\n");
	CalcSNodesTStart = clock();
	for (long long i = 0; i < N; i++) {
		CalcSuperNodesLowCCS(valLCCS, LowRowCCS, LowIndxCCS, &SNodesLow, NodesNLow, nzL, N);
		CalcSuperNodesUpCCS(valUCCS, UpRowCCS, UpIndxCCS, &SNodesUp, NodesNUp, nzU, N);
		CalcSNodesTFinish = clock();
	}
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
	for (long long i = 0; i < N; i++) {
		nodeSolverLow(valLCRS, colL, indxL, SNodesLow, NodesNLow, xNodeLow[i], bNodeLow[i], N);
		nodeSolverUp(valUCRS, colU, indxU, SNodesUp, NodesNUp, xNodeUp[i], bNodeUp[i], N);
	}
	NodeTFinish = clock();
	printf("NodeSolverUp finished\n");

	/// BlockSolverLow & BlockSolverUp
	printf("BlockSolverLow started\n");
	BlockTStart = clock();
	for (long long i = 0; i < N; i++) {
		blockSolverLowCCS(valLCRS, colL, indxL, valLCCS, LowRowCCS, LowIndxCCS, SNodesLow, NodesNLow, xBlockLow[i], bBlockLow[i], N);
		blockSolverUpCCS(valUCRS, colU, indxU, valUCCS, UpRowCCS, UpIndxCCS, SNodesUp, NodesNUp, xBlockUp[i], bBlockUp[i], N);
	}
	BlockTFinish = clock();
	printf("BlockSolverUp finished\n");

	/// BlockSolverLowFull & BlockSolverUpFull
	printf("BlockSolverLowFull started\n");
	BlockTStart = clock();
	blockSolverLowFull(valLCRS, colL, indxL, valLCCS, LowRowCCS, LowIndxCCS, SNodesLow, NodesNLow, xBlockFullLow, bBlockFullLow, N);
	blockSolverUpFull(valUCRS, colU, indxU, valUCCS, UpRowCCS, UpIndxCCS, SNodesUp, NodesNUp, xBlockFullUp, bBlockFullUp, N);
	BlockTFinish = clock();
	printf("BlockSolverUpFull finished\n");

	/*printf("\nNodesLow/N : %d/%d\n", NodesNLow - 1, N);
	for (int i = 0; i < NodesNLow; i++) {
		printf("%d  ", SNodesLow[i]);
	}printf("\n\n");
	printf("\nNodesUp/N : %d/%d\n", NodesNUp - 1, N);
	for (int i = 0; i < NodesNUp; i++) {
		printf("%d  ", SNodesUp[i]);
	}*/
	/// /////////////////////////////////////////////////////////////////////////
	/// MKL /////////////////////////////////////////////////////////////////////
	/*const long long  MKLn = N;
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
	*/
	printf("\n\n [CRS]\n");
	for (long long i = 0; i < N; i++) {
		printVectorF(xCRSLow[i], N);
		printVectorF(xCRSUp[i], N);
	}
	printf("[Node]\n");
	for (long long i = 0; i < N; i++) {
		printVectorF(xNodeLow[i], N);
		printVectorF(xNodeUp[i], N);
	}
	printf("[Block]\n");
	for (long long i = 0; i < N; i++) {
		printVectorF(xBlockLow[i], N);
		printVectorF(xBlockUp[i], N);
	}
	printf("[BlockFull]\n");
	for (long long i = 0; i < N; i++) {
		printVectorF(xBlockFullLow[i], N);
		printVectorF(xBlockFullUp[i], N);
	}
	//printf("[MKL]\n");
	for (long long i = 0; i < N; i++) {
		//printVectorF(xMKLLow[i], N);
		//printVectorF(xMKLUp[i], N);
	}
	long long errorsU = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if ((fabs(xNodeUp[i][j] - xBlockUp[i][j]) < 0.01) && (fabs(xBlockFullUp[i][j] - xCRSUp[i][j]) < 0.01)) {
				printf("[U] %lf --- %lf --- %lf --- %lf\n", xNodeUp[i][j], xBlockUp[i][j], xCRSUp[i][j], xBlockFullUp[i][j]);
			}
			else {
				errorsU++;
				printf("[U] %lf --- %lf --- %lf --- %lf\n", xNodeUp[i][j], xBlockUp[i][j], xCRSUp[i][j], xBlockFullUp[i][j]);
			}
		}
	}

	long long errorsL = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if ((fabs(xNodeLow[i][j] - xBlockLow[i][j]) < 0.01) && (fabs(xBlockFullLow[i][j] - xCRSLow[i][j]) < 0.01)) {
				printf("[L] %lf --- %lf --- %lf --- %lf\n", xNodeLow[i][j], xBlockLow[i][j], xCRSLow[i][j], xBlockFullLow[i][j]);
			}
			else {
				errorsL++;
				printf("[L] %lf --- %lf --- %lf --- %lf\n", xNodeLow[i][j], xBlockLow[i][j], xCRSLow[i][j], xBlockFullLow[i][j]);
			}
		}
	}
	printf("\n errorU: %d / %d\n", errorsU, N*N);
	printf("\n errorL: %d / %d\n", errorsL, N*N);

	for (long long i = 0; i < N; i++) {
		free(xMKLLow[i]);
		free(xMKLUp[i]);
		free(bMKLLow[i]);
		free(xCRSLow[i]);
		free(xCRSUp[i]);
		free(bCRSLow[i]);
		free(xNodeLow[i]);
		free(xNodeUp[i]);
		free(bNodeLow[i]);
		free(xBlockLow[i]);
		free(xBlockUp[i]);
		free(bBlockUp[i]);
		free(bCRSUp[i]);
		free(bBlockLow[i]);
		free(bMKLUp[i]);
		free(bNodeUp[i]);
		free(bBlockFullLow[i]);
		free(bBlockFullUp[i]);
		free(xBlockFullLow[i]);
		free(xBlockFullUp[i]);
	}
	free(xMKLLow); free(xMKLUp); free(bMKLLow);
	free(xCRSLow);	free(xCRSUp); free(bCRSLow); free(xNodeLow); free(xNodeUp); free(bNodeLow);
	free(xBlockLow); free(xBlockUp); free(bBlockUp); free(bCRSUp); free(bBlockLow);
	free(xBlockFullLow); free(xBlockFullUp); free(bBlockFullUp); free(bBlockFullLow);
	free(bMKLUp);
	free(bNodeUp);
	free(IU); free(JU); free(colU); free(colL);
	free(valUCOO); free(valUCRS);	free(valLCRS); free(indxU);	free(indxL); free(valLCCS); free(valUCCS);
	free(UpRowCCS); free(LowRowCCS); free(UpIndxCCS); free(LowIndxCCS);
	return 0;
}