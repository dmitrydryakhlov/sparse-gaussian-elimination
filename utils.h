#pragma once
#include <stdio.h>
#include <stdlib.h>
#define MM_MAX_LINE_LENGTH 1025
#define MM_PREMATURE_EOF 12


int readMTX(const char* fileName, int **I, int **J, double **val, int *M, int *N, int *nz);
int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz);
int COOtoCRS(int n, int nz, int *I, int *J, double *valCOO, int **indx, int **col, double **valCrs);
int saveBinCRS(const char* fileName, int n, int *row, int *col, double *val);
int cutLowerTriangleCOO(int nz, int *I, int *J, double *val, int *nzU, int **IU, int **JU, double **valU);
int cutUpperTriangleCOO(int nz, int *I, int *J, double *val, int *nzL, int **IL, int **JL, double **valL);
int transposeCOO(int nz, int *I, int *J);
void printmatrixSparceCOO(int n, int nz, int *I, int *J, double  *val);
int countZeroDiag(int *I, int *J, int nz, int N);
void getZerosDiagNumbers(int *I, int *J, int nz, int N, int count, int* addDiag);
void fillDiag(int *I, int *J, double *val, int *Inew, int *Jnew, double *valNew,
	int N, int nz, int nzNew, int *addDiag);
int CheckSolv(int n, double* x, double*xCheck, double * solution);
void gaussBackLow(int n, double* y, double* b, double *valCrsL, int* colL, int* indxL);
void gaussBackUp(int n, double* x, double* y, double *valCrsU, int* colU, int* indxU);
void matrixMultVector(int n, double* x, double* xCheck, double *valCrs, int* col, int* indx);