#pragma once
#include <stdio.h>
#include "utils.h"
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <random>
//#include <mkl.h>
//#include <mkl_pardiso.h>
//#include "type.h"

int main(int argc, char* argv[]) {
	int error = 0;
	float MyTimeStart;
	float MyTimeFinish;
	float MKLTimeStart;
	float MKLTimeFinish;

	int *I, *IU, *J, *JU, *indxU, *colU, *IL, *JL, *colL, *indxL;
	double *val, *valU, *valCrsU, *valCrsL, *valCooU, *valL, valCooL;
	double *e, *y, *bx, *by, *x, *MKLx;
	int N, M, nz, nzU, nzL;
	int i, j;
	if (readMTX(argv[1], &I, &J, &val, &M, &N, &nz) != 0)
		exit(0);

	int zeroDiagCount = countZeroDiag(I, J, nz, N);  //сколько диагональных элементов - нули
	int *addDiag = (int*)malloc(zeroDiagCount * sizeof(int)); //
	int *Inew = (int*)malloc((nz + zeroDiagCount) * sizeof(int));
	int *Jnew = (int*)malloc((nz + zeroDiagCount) * sizeof(int));
	double *valNew = (double*)malloc((nz + zeroDiagCount) * sizeof(double));
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
	IU = (int*)malloc((nzU) * sizeof(int));
	JU = (int*)malloc((nzU) * sizeof(int));
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
	//void *pt[64];
	//pt[0] = 0;
	//const long long maxfct = 1;
	//const long long MKLnz = nzU;
	//	mnum = 1;
	//	mtype = 2; // матрица симметрична и положительно определена
	//	phase = 11; // переупорядочивание и символическая факторизация
	//	*perm = new long long[nzL]; // перестановка
	//	nrhs = 1; // число правых частей
	//	msglvl = 1; // вывод статистической информации
	//	MKLerror = 0; // код ошибки
	//const long long iparm[64];
	//long long iparm[0] = 0; // используются параметры по умолчанию
				  // используются переупорядочивание методом вложенных сечений
	//double *MKLx = NULL;
	//MKLTimeStart = clock() / (float)1000;
	//PARDISO(pt,
	//	&maxfct,
	//	&mnum,
	//	&mtype,
	//	&phase,
	//	MKLnz,
	//	valCrsU,
	//	indxU,
	//	colU,
	//	perm,
	//	&nrhs,
	//	iparm,
	//	&msglvl,
	//	bx,
	//	MKLx,
	//	&MKLerror);
	//MKLTimeFinish = clock() / (float)1000;
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
