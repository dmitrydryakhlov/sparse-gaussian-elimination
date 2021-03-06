#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MM_MAX_LINE_LENGTH 1025
#define MM_PREMATURE_EOF 12


long long readMTX(const char* fileName, long long **I, long long **J, double **val, long long *M, long long *N, long long *nz);
long long mm_write_mtx_crd(char fname[], long long M, long long N, long long nz, long long I[], long long J[], double val[]);
long long mm_read_mtx_crd_size(FILE *f, long long *M, long long *N, long long *nz);
long long mm_write_mtx_crd_size(FILE *f, int M, int N, int nz);
long long saveBinCRS(const char* fileName, long long n, long long *row, long long *col, double *val);
int ReadMatrixFromBinaryFile(char* matrixName, long long &nz, long long &N, long long **indx, long long **col, double** value);

long long COOtoCRS(long long n, long long nz, long long *I, long long *J, double *valCOO, long long **indx, long long **col, double **valCrs);
void CRStoCOO(long long N, long long nz, long long* I, long long *J, double* valCOO,
	long long* indxCRS, long long* colCRS, double* valCRS);
void COOtoCCS(long long N, long long nz, long long *I, long long *J,
	double *valCOO, long long **indx, long long **row, double **valCcs);

long long cutLowerTriangleCOO(long long nz, long long *I, long long *J, double *val, long long *nzU, long long **IU, long long **JU, double **valU);
long long cutUpperTriangleCOO(long long nz, long long *I, long long *J, double *val, long long *nzL, long long **IL, long long **JL, double **valL);

long long transposeCOO(long long nz, long long *I, long long *J, double *val, long long &nzL, long long **IL, long long **JL, double **valL);
long long transposeCOO(long long nz, long long *I, long long *J);

void printSparceCOO(long long n, long long nz, long long *I, long long *J, double  *val);
void printVectorF(double* b, long long N);
void printVectorI(long long* b, long long N);

long long countZeroDiag(long long *I, long long *J, long long nz, long long N);
void getZerosDiagNumbers(long long *I, long long *J, long long nz, long long N, long long count, long long* addDiag);
void checkAndFillDiag(long long** I, long long** J, long long &nz, long long N, double** val);
void fillDiag(long long *I, long long *J, double *val, long long *Inew, long long *Jnew, double *valNew,
	long long N, long long nz, long long nzNew, long long *addDiag);

long long CheckSolv(long long n, double* x, double*xCheck);
double absError(long long n, double* x, double * xCheck);

void gaussBackLow(long long startNode, long long finishNode, double* y, double* b, double *valCrsL, long long* colL, long long* indxL);
void gaussBackUp(long long startNode, long long finishNode, double* x, double* y, double *valCrsU, long long* colU, long long* indxU);

void nodeSolverLow(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx,
	long long* SNodesLow, long long NodesNLow, double* xLow, double* bLow, long long N);
void nodeSolverUp(double *MatrixUpVal, long long *MatrixUpRow, long long *MatrixUpIndx,
	long long* SNodesUp, long long NodesNUp, double* xUp, double* bUp, long long N);

void matrixMultVector(long long n, double* x, double* xCheck, double *valCrs, long long* col, long long* indx);
void mallocMatrixCCS(double** AVal, long long** ARow, long long** AIndx, long long nz, long long N);

void makeBlock6x6LowCCS(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx, long long NzL, long long N);
void makeBlock6x6UpCCS(double *MatrixUpVal, long long *MatrixUpRow, long long *MatrixUpIndx, long long NzU, long long N);

void mallocVectors(double** x, double** y, double** b,
	double** xCRS, double** yCRS, double** bCRS,
	double** xMKL, double** yMKL, double** bMKL,
	double** xBlock, double** yBlock, double** bBlock, long long N);
void mallocMatrixNxM(double*** A, long long N, long long M);

void randVector(double* b, long long N);

void MKLPrepare(long long ** colU_short, long long ** colL_short, long long ** indxU_short,
	long long ** indxL_short, const long long MKLn, long long* indxU, long long* indxL,
	double* MKLbx, long long* colU, long long* colL, double* bx, long long N);
void freeMem(double** xUp, double** xLow, double** xCCSUp, double** xCCSLow,
	double** bUp, double** bLow, double** bCCSUp, double** bCCSLow,
	double** bMKLLow, double** bMKLUp, double** xMKLLow, double** xMKLUp,
	long long ** IU, long long ** JU, long long ** colU, long long ** colL,
	double ** valU, double ** valCrsU, double ** valCrsL, long long ** indxU,
	long long ** indxL, double** valLowCCS, double** valUpCCS, long long** UpRowCCS,
	long long** LowRowCCS, long long** UpIndxCCS, long long** LowIndxCCS);



void CalcSuperNodesLowCCS(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx,
	long long **SNodesLow, long long &NodesNLow, long long NzL, long long N);
void CalcSuperNodesUpCCS(double *MatrixUpVal, long long *MatrixUpRow, long long *MatrixUpIndx,
	long long **SNodesUp, long long &NodesNUp, long long NzU, long long N);

void changeBLowCCS(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx,
	long long SNodesLowl, long long SNodesLowl_1, long long N, double* x, double* bLow);
void changeBUpCCS(double *MatrixUpVal, long long *MatrixUpRow, long long *MatrixUpIndx,
	long long SNodesUpl, long long SNodesUpl_1, long long N, double* x, double* b);

void blockSolverLowCCS(double *MatrixLowValCRS, long long *MatrixLowColCRS, long long *MatrixLowIndxCRS,
	double *MatrixLowValCCS, long long *MatrixLowRowCCS, long long *MatrixLowIndxCCS,
	long long* SNodesLow, long long NodesNLow, double* xLow, double* bLow, long long N);
void blockSolverUpCCS(double *MatrixUpValCRS, long long *MatrixUpColCRS, long long *MatrixUpIndxCRS,
	double *MatrixUpValCCS, long long *MatrixUpRowCCS, long long *MatrixUpIndxCCS,
	long long* SNodesUp, long long NodesNUp, double* xUp, double* bUp, long long N);

void blockSolverLowCRS(double *MatrixLowValCRS, long long *MatrixLowColCRS, long long *MatrixLowIndxCRS,
	double *MatrixLowValCCS, long long *MatrixRowColCCS, long long *MatrixLowIndxCCS,
	long long* SNodesLow, long long NodesNLow, double* xLow, double* bLow, long long N);

void blockSolverUpCRS(double *MatrixUpValCRS, long long *MatrixUpColCRS, long long *MatrixUpIndxCRS,
	double *MatrixUpValCCS, long long *MatrixUpColCCS, long long *MatrixUpIndxCCS,
	long long* SNodesUp, long long NodesNUp, double* xUp, double* bUp, long long N);

void blockSolverLowFullCCS(double *MatrixLowValCRS, long long *MatrixLowColCRS, long long *MatrixLowIndxCRS,
	double *MatrixLowValCCS, long long *MatrixLowRowCCS, long long *MatrixLowIndxCCS,
	long long* SNodesLow, long long NodesNLow, double** xLow, double** bLow, long long N, long long M);
void blockSolverUpFullCCS(double *MatrixUpValCRS, long long *MatrixUpColCRS, long long *MatrixUpIndxCRS,
	double *MatrixUpValCCS, long long *MatrixUpRowCCS, long long *MatrixUpIndxCCS,
	long long* SNodesUp, long long NodesNUp, double** xUp, double** bUp, long long N, long long M);

void blockSolverLowFullCCSPrl(double *MatrixLowValCRS, long long *MatrixLowColCRS, long long *MatrixLowIndxCRS,
	double *MatrixLowValCCS, long long *MatrixLowRowCCS, long long *MatrixLowIndxCCS,
	long long* SNodesLow, long long NodesNLow, double** xLow, double** bLow, long long N, long long M);
void blockSolverUpFullCCSPrl(double *MatrixUpValCRS, long long *MatrixUpColCRS, long long *MatrixUpIndxCRS,
	double *MatrixUpValCCS, long long *MatrixUpRowCCS, long long *MatrixUpIndxCCS,
	long long* SNodesUp, long long NodesNUp, double** xUp, double** bUp, long long N, long long M);

void blockSolverLowFullCRS(double *MatrixLowValCRS, long long *MatrixLowColCRS, long long *MatrixLowIndxCRS,
	double *MatrixLowValCCS, long long *MatrixLowRowCCS, long long *MatrixLowIndxCCS,
	long long* SNodesLow, long long NodesNLow, double** xLow, double** bLow, long long N, long long M);
void blockSolverUpFullCRS(double *MatrixUpValCRS, long long *MatrixUpColCRS, long long *MatrixUpIndxCRS,
	double *MatrixUpValCCS, long long *MatrixUpRowCCS, long long *MatrixUpIndxCCS,
	long long* SNodesUp, long long NodesNUp, double** xUp, double** bUp, long long N, long long M);

void blockSolverLowFullCRSPrl(double *MatrixLowValCRS, long long *MatrixLowColCRS, long long *MatrixLowIndxCRS,
	double *MatrixLowValCCS, long long *MatrixLowRowCCS, long long *MatrixLowIndxCCS,
	long long* SNodesLow, long long NodesNLow, double** xLow, double** bLow, long long N, long long M);
void blockSolverUpFullCRSPrl(double *MatrixUpValCRS, long long *MatrixUpColCRS, long long *MatrixUpIndxCRS,
	double *MatrixUpValCCS, long long *MatrixUpRowCCS, long long *MatrixUpIndxCCS,
	long long* SNodesUp, long long NodesNUp, double** xUp, double** bUp, long long N, long long M);

void mallocMatrixCOO(long long** I, long long** J, double** COOVal, long long NzL);

void makeMatrix6x6COO(long long* I, long long* J, double* COOVal, long long NzL);
void mallocMatrixCOO(long long ** I, long long ** J, double ** COOVal, long long NzL);

void makeBlockMatrix8x8LowCOO(long long * I, long long * J, double * COOVal, long long NzL);
//void makeBlockMatrix8x8UpCOO(long long * I, long long * J, float * COOVal, long long NzL);

void makeBlock8x8LowCCS(double *MatrixLowValCCS, long long *MatrixLowRowCCS, long long *MatrixLowIndxCCS, long long NzL, long long N);
void makeBlock8x8UpCCS(double *MatrixUpValCCS, long long *MatrixUpRowCCS, long long *MatrixUpIndxCCS, long long NzU, long long N);

void makeBlock8x8LowCRS(double *MatrixLowValCRS, long long *MatrixLowColCRS, long long *MatrixLowIndxCRS, long long NzL, long long N);
void makeBlock8x8UpCRS(double *MatrixUpValCRS, long long *MatrixUpColCRS, long long *MatrixUpIndxCRS, long long NzU, long long N);

void makeBlockMatrix12x12LowCOO(long long * I, long long * J, double * COOVal, long long NzL);
void makeBlockMatrix12x12UpCOO(long long * I, long long * J, double * COOVal, long long NzL);

void makeBlockMatrix12x12LowCOORandom(long long * I, long long * J, double * COOVal, long long NzL);
void makeBlockMatrix12x12UpCOORandom(long long * I, long long * J, double * COOVal, long long NzU);

void generateBigBlockMatrixL(long long * I, long long * J, double * COOVal, long long NzL, long long N, long long blockSize, long long* fullRowL);
void generateBigBlockMatrixU(long long * I, long long * J, double * COOVal, long long NzU, long long N, long long blockSize, long long* fullRowU);

long long calcNzL(long long N, long long blockSizeL, long long **fullRowL, long long fullness);
long long calcNzU(long long N, long long blockSizeU, long long **fullRowU, long long fullness);


void gaussBackBlockUp(long long startNode, long long finishNode, double * x,
	double * y, double * valCrsU, long long * colU, long long * indxU);

void gaussBackBlockLow(long long startNode, long long finishNode, double * x,
	double * y, double * valCrsL, long long * colL, long long * indxL);

void gaussBackBlockUpFull(long long startNode, long long finishNode, double ** x,
	double ** y, double * valUCRS, long long * colUCRS, long long * indxUCRS,
	double * valUCCS, long long * rowUCCS, long long * indxUCCS, long long N, long long M);

void gaussBackBlockLowFull(long long startNode, long long finishNode, double ** x,
	double ** y, double * valUCRS, long long * colLCRS, long long * indxLCRS,
	double * valLCCS, long long * rowLCCS, long long * indxUCCS, long long N, long long M);

void gaussBackBlockUpFullPrl(long long startNode, long long finishNode, double ** x,
	double ** y, double * valUCRS, long long * colUCRS, long long * indxUCRS,
	double * valUCCS, long long * rowUCCS, long long * indxUCCS, long long N, long long M);

void gaussBackBlockLowFullPrl(long long startNode, long long finishNode, double ** x,
	double ** y, double * valUCRS, long long * colLCRS, long long * indxLCRS,
	double * valLCCS, long long * rowLCCS, long long * indxUCCS, long long N, long long M);