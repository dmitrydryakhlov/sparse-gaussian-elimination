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

	long long zeroDiagCount = countZeroDiag(I, J, nz, N);  //сколько диагональных элементов - нули
	long long *addDiag = (long long*)malloc((zeroDiagCount + 1) * sizeof(long long)); //
	long long *Inew = (long long*)malloc((nz + zeroDiagCount + 1) * sizeof(long long));
	long long *Jnew = (long long*)malloc((nz + zeroDiagCount + 1) * sizeof(long long));
	double *valNew = (double*)malloc((nz + zeroDiagCount + 1) * sizeof(double));
	getZerosDiagNumbers(I, J, nz, N, zeroDiagCount, addDiag); //получаем номера строк, где нули на главной диагонали
	fillDiag(I, J, val, Inew, Jnew, valNew, N, nz, nz + zeroDiagCount, addDiag);//избавляемя от нулей на гл диагонали
	free(I); I = Inew;
	free(J); J = Jnew;
	free(val); val = valNew;
	nz = nz + zeroDiagCount; //корректируем количество ненулевых элементов в итоговой матрице

	y = (double*)malloc((N) * sizeof(double));
	e = (double*)malloc((N) * sizeof(double));
	bx = (double*)malloc((N) * sizeof(double));
	MKLbx = (double*)malloc((N) * sizeof(double));
	by = (double*)malloc((N) * sizeof(double));
	MKLby = (double*)malloc((N) * sizeof(double));
	x = (double*)malloc((N) * sizeof(double));
	MKLx = (double*)malloc((N) * sizeof(double));
	MKLy = (double*)malloc((N) * sizeof(double));


	for (i = 0; i < N; i++)
		e[i] = (double)rand() / 3.0;


	//транспонирование
	//transposeCOO(nz, I, J);  
	//printf("Finished transposing coordinate matrix\n");
	//cutLowerTriangleCOO(nz, I, J, val, &nzU, &IU, &JU, &valU);

	cutUpperTriangleCOO(nz, I, J, val, &nzL, &IL, &JL, &valL);
	printf("Finished cutting coordinate matrix\n");

	nzU = nzL;
	IU = (long long*)malloc((nzU) * sizeof(long long));
	JU = (long long*)malloc((nzU) * sizeof(long long));
	valU = (double*)malloc((nzU) * sizeof(double));

	transposeCOO(nzL, IL, JL, valL, nzU, IU, JU, valU);

	COOtoCRS(N, nzU, IU, JU, valU, &indxU, &colU, &valCrsU);

	COOtoCRS(N, nzL, IL, JL, valL, &indxL, &colL, &valCrsL);

	//saveBinCRS(argv[2], N, row, col, valCRS);
	//printf("Finished saving CRS matrix\n");

	matrixMultVector(N, e, by, valCrsU, colU, indxU);// U*x = by
	matrixMultVector(N, by, bx, valCrsL, colL, indxL);//L*by = bx

	MyTimeStart = clock() / (float)CLOCKS_PER_SEC;
	gaussBackLow(N, y, bx, valCrsL, colL, indxL);//находим y из L*y=bx
	gaussBackUp(N, x, y, valCrsU, colU, indxU);//находим х из U*x=y
	MyTimeFinish = clock() / (float)CLOCKS_PER_SEC;

	//////////////MKL////////////////////////////////////////////////////////
	const long long  MKLn = N;
	long long	MKLerror = 0;

	const long long MKLn_short = MKLn;

	long long *colU_short = new long long[indxU[MKLn]];
	long long *colL_short = new long long[indxL[MKLn]];

	long long *indxU_short = new long long[MKLn + 1];
	long long *indxL_short = new long long[MKLn + 1];

	for (int l = 0; l < N; l++) {
		MKLbx[l + 1] = bx[l];
	}
	for (int l = 0; l < indxU[MKLn]; l++) {
		colU[l] += 1;
		colU_short[l] = colU[l];
	}
	for (int l = 0; l < indxL[MKLn]; l++) {
		colL[l] += 1;
		colL_short[l] = colL[l];
	}
	for (int l = 0; l < MKLn + 1; l++) {
		indxU[l] += 1;
		indxU_short[l] = indxU[l];
	}
	for (int l = 0; l < MKLn + 1; l++) {
		indxL[l] += 1;
		indxL_short[l] = indxL[l];
	}

	MKLTimeStart = clock() / (float)1000;
	mkl_dcsrtrsv("L", "N", "N", &MKLn_short, valCrsL, indxL_short, colL_short, bx, MKLy);
	mkl_dcsrtrsv("U", "N", "N", &MKLn_short, valCrsU, indxU_short, colU_short, MKLy, MKLx);
	MKLTimeFinish = clock() / (float)1000;

	/// ////////////////MKL//////////////////////////////////////////////////////////

	printf("MKL error: %lld\n", MKLerror);
	
	for (int l = 0; l < N; l++) {
		//printf(" MKLx[%d] = %lf,  x[%d] = %f, e[%d] = %f \n", l, MKLx[l], l, x[l], l, e[l]);
	}

	printf("\n Size: %d x %d , nz = %d , nzL = %d , nzU = %d\n", N, N, nz, nzL, nzU);
	int errors = CheckSolv(N, e, x);
	int MKLErrors = CheckSolv(N, e, MKLx);
	double absErrors = absError(N, e, x);
	double MKLAbsErrors = absError(N, e, MKLx);
	if (errors == 0) { printf("\n All is Ok! absErrors = %.10f \n Time is %f\n", absErrors, (MyTimeFinish - MyTimeStart) * 1000); }
	else { printf("\n Errors: %d/%d  , Time is %f\n", errors, N, (MyTimeFinish - MyTimeStart) * 1000); }


	if (MKLErrors == 0) { printf("\n[MKL] All is Ok! absErrors = %.10f \n Time is %f\n", MKLAbsErrors, (MKLTimeFinish - MKLTimeStart) * 1000); }
	else { printf("\n Errors: %d/%d  , MKLTime is %f\n", MKLErrors, N, (MKLTimeFinish - MKLTimeStart) * 1000); }

	free(e);
	free(y);
	free(bx);
	free(by);
	free(I);
	free(IU);
	free(J);
	free(JU);
	free(val);
	free(valU);
	free(valCrsU);
	free(valCrsL);
	free(colL);
	free(colU);
	free(MKLbx);
	free(MKLby);
	free(MKLx);
	free(MKLy);
	free(indxU);
	free(indxL);
	return 0;
}
