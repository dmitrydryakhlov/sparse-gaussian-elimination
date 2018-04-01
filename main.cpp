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
	float MyTimeStart;
	float MyTimeFinish;
	float MKLTimeStart;
	float MKLTimeFinish;
	long long *I, *IU, *J, *JU, *indxU, *colU, *IL, *JL, *colL, *indxL;
	double *val, *valU, *valCrsU, *valCrsL, *valCooU, *valL, valCooL;
	double *e, *y, *bx, *MKLbx, *by, *MKLby, *x, *MKLx, *MKLy;
	long long N, M, nz, nzU, nzL;
	long long i, j;
	if (readMTX(argv[1], &I, &J, &val, &M, &N, &nz) != 0)
		exit(0);

	checkAndFillDiag(&I, &J, nz, N, &val);
	mallocVectors(&y, &e, &bx, &MKLbx, &by, &MKLby, &x, &MKLx, &MKLy, N);
	randVector(e, N);
	cutUpperTriangleCOO(nz, I, J, val, &nzL, &IL, &JL, &valL);
	printf("Finished cutting coordinate matrix\n");

	transposeCOO(nzL, IL, JL, valL, nzU, &IU, &JU, &valU);
	COOtoCRS(N, nzU, IU, JU, valU, &indxU, &colU, &valCrsU);
	COOtoCRS(N, nzL, IL, JL, valL, &indxL, &colL, &valCrsL);

	//saveBinCRS(argv[2], N, row, col, valCRS);
	//printf("Finished saving CRS matrix\n");

	matrixMultVector(N, e, by, valCrsU, colU, indxU);// U*x = by
	matrixMultVector(N, by, bx, valCrsL, colL, indxL);//L*by = bx

	MyTimeStart = clock();
	//	gaussBackLow(N, y, bx, valCrsL, colL, indxL);//находим y из L*y=bx
	//	gaussBackUp(N, x, y, valCrsU, colU, indxU);//находим х из U*x=y
	MyTimeFinish = clock();

	const long long  MKLn = N;
	long long	MKLerror = 0;
	const long long MKLn_short = MKLn;
	long long* colU_short;
	long long* colL_short;
	long long* indxU_short;
	long long* indxL_short;
	MKLPrepare(&colU_short, &colL_short, &indxU_short, &indxL_short,
		MKLn_short, indxU, indxL, MKLbx, colU, colL, bx, N);

	MKLTimeStart = clock();
	mkl_dcsrtrsv("L", "N", "N", &MKLn_short, valCrsL, indxL_short, colL_short, bx, MKLy);
	mkl_dcsrtrsv("U", "N", "N", &MKLn_short, valCrsU, indxU_short, colU_short, MKLy, MKLx);
	MKLTimeFinish = clock();

	printf("MKL error: %lld\n", MKLerror);
	printf("\n Size: %d x %d , nz = %d , nzL = %d , nzU = %d\n", N, N, nz, nzL, nzU);
	int errors = CheckSolv(N, e, x);
	int MKLErrors = CheckSolv(N, e, MKLx);
	double absErrors = absError(N, e, x);
	double MKLAbsErrors = absError(N, e, MKLx);

	if (errors == 0) {
		printf("\n All is Ok! absErrors = %.10f \n Time is %f\n",
			absErrors, (MyTimeFinish - MyTimeStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n Errors: %d/%d  , Time is %f\n",
			errors, N, (MyTimeFinish - MyTimeStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	if (MKLErrors == 0) {
		printf("\n[MKL] All is Ok! absErrors = %.10f \n Time is %f\n",
			MKLAbsErrors, (MKLTimeFinish - MKLTimeStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}
	else {
		printf("\n Errors: %d/%d  , MKLTime is %f\n",
			MKLErrors, N, (MKLTimeFinish - MKLTimeStart) / (float)CLOCKS_PER_SEC * 1000.0f);
	}

	freeMem(&y, &e, &bx, &MKLbx, &by, &MKLby, &x, &MKLx, &MKLy, &I, &IU,
		&J, &JU, &colU, &colL, &val, &valU, &valCrsU, &valCrsL, &indxU, &indxL);
	return 0;
}