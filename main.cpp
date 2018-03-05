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
//#include "type.h"

int main(int argc, char* argv[]) {
	int error = 0;
	float MyTimeStart;
	float MyTimeFinish;
	float MKLTimeStart;
	float MKLTimeFinish;

	long long *I, *IU, *J, *JU, *indxU, *colU, *IL, *JL, *colL, *indxL;
	double *val, *valU, *valCrsU, *valCrsL, *valCooU, *valL, valCooL;
	double *e, *y, *bx, *by, *x, *MKLx;
	long long N, M, nz, nzU, nzL;
	long long i, j;
	if (readMTX(argv[1], &I, &J, &val, &M, &N, &nz) != 0)
		exit(0);

	long long zeroDiagCount = countZeroDiag(I, J, nz, N);  //сколько диагональных элементов - нули
	long long *addDiag = (long long*)malloc((zeroDiagCount+1) * sizeof(long long)); //
	long long *Inew = (long long*)malloc((nz + zeroDiagCount+1) * sizeof(long long));
	long long *Jnew = (long long*)malloc((nz + zeroDiagCount+1) * sizeof(long long));
	double *valNew = (double*)malloc((nz + zeroDiagCount+1) * sizeof(double));
	getZerosDiagNumbers(I, J, nz, N, zeroDiagCount, addDiag); //получаем номера строк, где нули на главной диагонали
	fillDiag(I, J, val, Inew, Jnew, valNew, N, nz, nz + zeroDiagCount, addDiag);//избавляемя от нулей на гл диагонали
	free(I); I = Inew;
	free(J); J = Jnew;
	free(val); val = valNew;
	nz = nz + zeroDiagCount; //корректируем количество ненулевых элементов в итоговой матрице

	y = (double*)malloc((N) * sizeof(double));
	e = (double*)malloc((N) * sizeof(double));
	bx = (double*)malloc((N) * sizeof(double));
	by = (double*)malloc((N) * sizeof(double));
	x = (double*)malloc((N) * sizeof(double));
	MKLx = (double*)malloc((N) * sizeof(double));
	

	for (i = 0; i < N; i++) {
		e[i] = (double)rand()/3.0;
		//bx[i] = (double)rand();
	}
	
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
	//printf("Finished cutting coordinate matrix\n");
	
	COOtoCRS(N, nzU, IU, JU, valU, &indxU, &colU, &valCrsU);
	printf("Finished converting coordinate matrix to CRS matrix\n");

	COOtoCRS(N, nzL, IL, JL, valL, &indxL, &colL, &valCrsL);
	printf("Finished converting coordinate matrix to CRS matrix\n");

	//saveBinCRS(argv[2], N, row, col, valCRS);
	//printf("Finished saving CRS matrix\n");

	
	matrixMultVector(N, e, by, valCrsU, colU, indxU);// U*x = by
	matrixMultVector(N, by, bx, valCrsL, colL, indxL);//L*by = bx

	MyTimeStart = clock() / (float)CLOCKS_PER_SEC;
	gaussBackLow(N, y, bx, valCrsL, colL, indxL);//находим y из L*y=bx
	gaussBackUp (N, x, y, valCrsU, colU, indxU);//находим х из U*x=y
	MyTimeFinish = clock() /(float)CLOCKS_PER_SEC;


	//////////////MKL////////////////
	void *pt[64];
	pt[0] = 0;
	const long long maxfct = 1;
	const long long MKLnz = nzU;
	const long long	mnum = 1;
	const long long	mtype = 11; // мб 1? матрица симметрична и положительно определена
	const long long	phase = 333; // backward solver
	long long	*perm = new long long[nzL+1];
	perm[60] = 1;// перестановка
	const long long	nrhs = 1; // число правых частей
	const long long	msglvl = 1; // вывод статистической информации
	long long	MKLerror = 0; // код ошибки
	long long iparm[64];
	for (i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	//long long p = 0;
	//iparm[60] = p;// используются параметры по умолчанию
	iparm[0] = 1;
	iparm[1] = 0;//2 
	iparm[3] = 0;
	iparm[4] = 0;
	iparm[5] = 0;
	iparm[6] = 0;
	iparm[7] = 2;
	iparm[8] = 0;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[11] = 0;
	iparm[12] = 1; 
	iparm[13] = 0;
	iparm[14] = 0;
	iparm[15] = 0;
	iparm[16] = 0;
	iparm[17] = -1;
	iparm[18] = -1; 
	iparm[19] = 0;

	iparm[35] = 1;
	iparm[20] = 1;
	iparm[24] = 1;
	//iparm[33] = 1;
	//iparm[34] = 1;
	/*
	iparm[21] = 0;
	iparm[22] = 0;
	iparm[23] = 0;
	iparm[25] = 0;
	iparm[26] = 0;
	iparm[27] = 0;
	iparm[28] = 0;
	iparm[29] = 0;
	iparm[30] = 0;
	iparm[31] = 0;
	iparm[32] = 0;
	iparm[35] = 0;
	iparm[36] = 0;
	iparm[37] = 0;
	iparm[38] = 0;
	iparm[39] = 0;
	iparm[40] = 0;
	iparm[41] = 0;
	iparm[42] = 0;
	iparm[43] = 0;
	iparm[44] = 0;
	iparm[45] = 0;
	iparm[46] = 0;
	iparm[47] = 0;
	iparm[48] = 0;
	iparm[49] = 0;
	iparm[50] = 0;
	iparm[51] = 0;
	iparm[52] = 0;
	iparm[53] = 0;
	iparm[54] = 0;
	iparm[55] = 0;
	iparm[56] = 0;
	iparm[57] = 0;
	iparm[58] = 0;
	iparm[59] = 2;
	iparm[60] = 0;
	iparm[61] = 0;
	iparm[62] = 0;
	iparm[63] = 0;
	*/
	for (i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}

	/*for (int i = nzU+1; i > 0; i--) {
		valCrsU[i] = valCrsU[i - 1];
		colU[i] = colU[i - 1];
	}
	for (int i = N + 1; i > 0; i--) {
		indxU[i] = indxU[i - 1];
	}*/ //index 1 or 0

	MKLx = NULL; // ????????????????????????????????????????
	MKLTimeStart = clock() / (float)1000;
	for (i = 0; i < 3; i++) {
		iparm[11] = i;
		PARDISO(pt,
			&maxfct,
			&mnum,
			&mtype,
			&phase,
			&MKLnz,
			valCrsU,
			indxU,
			colU,
			perm,
			&nrhs,
			iparm,
			&msglvl,
			bx,
			MKLx,
			&MKLerror);
	}
	MKLTimeFinish = clock() / (float)1000;
	///////////////////MKL//////////////////////

	printf("\n Size: %d x %d , nz = %d , nzL = %d , nzU = %d", N, N, nz, nzL, nzU);
	int errors = CheckSolv(N, e, x);
	double absErrors = absError(N, e, x);
	if (errors == 0) { printf("\n All is Ok! absErrors = %.10f \n Time is %f\n", absErrors, (MyTimeFinish-MyTimeStart)*1000); }
	else { printf("\n Errors: %d/%d  , Time is %f\n", errors, N, (MyTimeFinish - MyTimeStart) * 1000); }
	//printf("\n All is Ok! MyTime is %f (ms)", (MyTimeFinish - MyTimeStart)*1000);
	//printf("\n All is Ok! MKLTime is %f (ms)", (MKLTimeFinish - MKLTimeStart) * 1000);


	free(e);
	free(IU);
	free(JU);
	free(valU);
	free(y);
	free(by);
	free(bx);
	free(I);
	free(J);
	free(val);
	free(indxU);
	free(colU);
	free(valCrsU);
	free(indxL);
	free(colL);
	free(valCrsL);
	return 0;
}
