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
long long transposeCOO(long long nz, long long *I, long long *J, double *val, long long nzL, long long *IL, long long *JL, double *valL);
long long transposeCOO(long long nz, long long *I, long long *J );
void printmatrixSparceCOO(long long n, long long nz, long long *I, long long *J, double  *val);
long long countZeroDiag(long long *I, long long *J, long long nz, long long N);
void getZerosDiagNumbers(long long *I, long long *J, long long nz, long long N, long long count, long long* addDiag);
void fillDiag(long long *I, long long *J, double *val, long long *Inew, long long *Jnew, double *valNew,
	long long N, long long nz, long long nzNew, long long *addDiag);
long long CheckSolv(long long n, double* x, double*xCheck);
double absError(long long n, double* x, double * xCheck);
void gaussBackLow(long long n, double* y, double* b, double *valCrsL, long long* colL, long long* indxL);
void gaussBackUp(long long n, double* x, double* y, double *valCrsU, long long* colU, long long* indxU);
void matrixMultVector(long long n, double* x, double* xCheck, double *valCrs, long long* col, long long* indx);