#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MM_MAX_LINE_LENGTH 1025
#define MM_PREMATURE_EOF 12


long long readMTX(const char* fileName, long long **I, long long **J, double **val, long long *M, long long *N, long long *nz);
long long mm_read_mtx_crd_size(FILE *f, long long *M, long long *N, long long *nz);
long long COOtoCRS(long long n, long long nz, long long *I, long long *J, double *valCOO, long long **indx, long long **col, double **valCrs);
long long saveBinCRS(const char* fileName, long long n, long long *row, long long *col, double *val);
long long cutLowerTriangleCOO(long long nz, long long *I, long long *J, double *val, long long *nzU, long long **IU, long long **JU, double **valU);
long long cutUpperTriangleCOO(long long nz, long long *I, long long *J, double *val, long long *nzL, long long **IL, long long **JL, double **valL);
long long transposeCOO(long long nz, long long *I, long long *J, double *val, long long &nzL, long long **IL, long long **JL, double **valL);
long long transposeCOO(long long nz, long long *I, long long *J);
void printSparceCOO(long long n, long long nz, long long *I, long long *J, double  *val);
long long countZeroDiag(long long *I, long long *J, long long nz, long long N);
void getZerosDiagNumbers(long long *I, long long *J, long long nz, long long N, long long count, long long* addDiag);
void fillDiag(long long *I, long long *J, double *val, long long *Inew, long long *Jnew, double *valNew,
	long long N, long long nz, long long nzNew, long long *addDiag);
int CheckSolv(long long n, double* x, double*xCheck);
double absError(long long n, double* x, double * xCheck);
void gaussBackLow(long long n, double* y, double* b, double *valCrsL, long long* colL, long long* indxL);
void gaussBackUp(long long n, double* x, double* y, double *valCrsU, long long* colU, long long* indxU);
void matrixMultVector(long long n, double* x, double* xCheck, double *valCrs, long long* col, long long* indx);
void checkAndFillDiag(long long** I, long long** J, long long &nz, long long N, double** val);
void mallocMatrixCCS(double** AVal, long long** ARow, long long** AIndx, long long nz, long long N);
void makeBlock6x6LowCCS(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx, long long NzL, long long N);
void makeBlock6x6UpCCS(double *MatrixUpVal, long long *MatrixUpRow, long long *MatrixUpIndx, long long NzU, long long N);
void mallocVectors(double** y, double** e, double** bx, double** xCCS, double** yCCS, double** bCCS,
	double** MKLbx, double** by, double** MKLby, double** x, double** MKLx, double** MKLy, long long N);
void randVector(double* b, long long N); 
void MKLPrepare(long long ** colU_short, long long ** colL_short, long long ** indxU_short,
	long long ** indxL_short, const long long MKLn, long long* indxU, long long* indxL,
	double* MKLbx, long long* colU, long long* colL, double* bx, long long N);
void freeMem(double** y, double** e, double** bx,
	double** MKLbx, double** by, double** MKLby, double** x, double** MKLx, double** MKLy,
	long long** I, long long** IU, long long** J, long long** JU, long long** colU, long long** colL,
	double** val, double** valU, double** valCrsU, double** valCrsL, long long** indxU, long long ** indxL,
	double** valLowCCS, double** valUpCCS, long long** UpRowCCS, long long** LowRowCCS, long long** UpIndxCCS, long long** LowIndxCCS);



void CalcSuperNodesLowCCS(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx,
	long long **SNodesLow, long long &NodesNLow, long long NzL, long long N);

void CalcSuperNodesUpCCS(double *MatrixUpVal, long long *MatrixUpRow, long long *MatrixUpIndx,
	long long **SNodesUp, long long &NodesNUp, long long NzU, long long N);

void changeBLowCCS(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx,
	long long indexX, long long N, double* x, double* b);
void changeBUpCCS(double *MatrixUpVal, long long *MatrixUpRow, long long *MatrixUpIndx,
	long long indexX, long long N, double* x, double* b);
void blockSolverLowCCS(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx,
	long long* SNodesLow, long long NodesNLow, double* xLow, double* bLow, long long N);
void blockSolverUpCCS(double *MatrixUpVal, long long *MatrixUpRow, long long *MatrixUpIndx,
	long long* SNodesUp, long long NodesNUp, double* xUp, double* bUp, long long N);

void mallocMatrixCOO(long long** I, long long** J, double** COOVal, long long NzL);

void COOtoCCS(long long N, long long nz, long long *I, long long *J, double *valCOO, long long **indx, long long **row, double **valCcs);

void makeMatrix6x6COO(long long* I, long long* J, double* COOVal, long long NzL);

void mallocMatrixCOO(long long ** I, long long ** J, double ** COOVal, long long NzL);
