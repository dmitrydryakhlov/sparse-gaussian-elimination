#include "utils.h"
#include <omp.h>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>

long long mm_read_mtx_crd_size(FILE *f, long long *M, long long *N, long long *nz)
{
	char line[MM_MAX_LINE_LENGTH];
	long long num_items_read;

	/* set return null parameter values, in case we exit with errors */
	*M = *N = *nz = 0;

	/* now continue scanning until you reach the end-of-comments */
	do
	{
		if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
			return MM_PREMATURE_EOF;
	} while (line[0] == '%');

	/* line[] is either blank or has M,N, nz */
	if (sscanf(line, "%d %d %d", M, N, nz) == 3)
		return 0;

	else
		do
		{
			num_items_read = fscanf(f, "%d %d %d", M, N, nz);
			if (num_items_read == EOF) return MM_PREMATURE_EOF;
		} while (num_items_read != 3);

		return 0;
}

long long mm_write_mtx_crd_size(FILE *f, int M, int N, int nz) {
	if (fprintf(f, "%d %d %d\n", M, N, nz) != 3)
		return 1;
	else
		return 0;
}

long long readMTX(const char * fileName, long long ** I, long long ** J, double ** val, long long * M, long long * N, long long * nz)
{
	FILE *f;
	long long i;
	long long ret_code;

	if ((f = fopen(fileName, "r")) == NULL) {
		printf("Can't open file %s\n", fileName);
		return 1;
	}
	else
		printf("Matrix: %s\n", fileName);

	if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) != 0)
		return 1;

	(*I) = (long long *)malloc((*nz) * sizeof(long long));
	(*J) = (long long *)malloc((*nz) * sizeof(long long));
	(*val) = (double *)malloc((*nz) * sizeof(double));

	for (i = 0; i < (*nz); i++) {
		int x, y;
		double value;
		fscanf(f, "%d %d %lf\n", &x, &y, &value);
		(*I)[i] = (long long)(x - 1);
		(*J)[i] = (long long)(y - 1);
		(*val)[i] = (long long)value;
	}
	if (f != stdin) fclose(f);
	return 0;
}

long long mm_write_mtx_crd(char fname[], long long M, long long N, long long nz,
	long long I[], long long J[], double val[]) {
	FILE *f;
	long long i;

	if ((f = fopen(fname, "w")) == NULL)
		return 1;

	fprintf(f, "%d %d %d\n", M, N, nz);

	for (i = 0; i < nz; i++)
		fprintf(f, "%d %d %20.16g\n", I[i] + 1, J[i] + 1, val[i]);

	if (f != stdout) fclose(f);
	return 0;
}

long long COOtoCRS(long long n, long long nz, long long *I, long long *J,
	double *valCOO, long long **indx, long long **col, double **valCrs)
{
	long long i;
	long long *places;

	(*indx) = (long long*)malloc((n + 1) * sizeof(long long));
	(*col) = (long long*)malloc((nz) * sizeof(long long));
	(*valCrs) = (double*)malloc((nz) * sizeof(double));
	places = (long long*)malloc((n + 1) * sizeof(long long));

	for (i = 0; i < n + 1; i++) {
		(*indx)[i] = 0;
	}
	for (i = 0; i < nz; i++) {
		(*indx)[I[i] + 1] += 1;
	}
	for (i = 0; i < n; i++) {
		(*indx)[i + 1] += (*indx)[i];
		places[i] = (*indx)[i];
	}
	for (i = 0; i < nz; i++) {
		(*col)[places[I[i]]] = J[i];
		(*valCrs)[places[I[i]]] = valCOO[i];
		places[I[i]]++;
	}
	//sort each column
	for (i = 0; i < n; i++) {
		long long start, finish, j, k;
		start = (*indx)[i];
		finish = (*indx)[i + 1];
		//need to write bubble sort
		for (j = start; j < finish; j++)
			for (k = start; k < finish - j - 1; k++)
				if ((*col)[k] > (*col)[k + 1]) {
					long long tempCol = (*col)[k];
					(*col)[k] = (*col)[k + 1];
					(*col)[k + 1] = tempCol;

					double tempVal = (*valCrs)[k];
					(*valCrs)[k] = (*valCrs)[k + 1];
					(*valCrs)[k + 1] = tempVal;
				}
	}
	free(places);
	return 0;
}

void printVectorF(double* b, long long N) {
	printf("Vector double:\n");
	for (long long i = 0; i < N; i++) {
		printf("%.3lf   ", b[i]);
	}
	printf("\n\n");
}

void printVectorI(long long* b, long long N) {
	printf("Vector int:\n");
	for (long long i = 0; i < N; i++) {
		printf("%d   ", b[i]);
	}
	printf("\n\n");
}

long long saveBinCRS(const char * fileName, long long n, long long * row, long long * col, double * val)
{
	FILE *f;
	f = fopen(fileName, "wb");
	fwrite(&n, sizeof(long long), 1, f);
	fwrite(&(row[n]), sizeof(long long), 1, f);
	fwrite(col, sizeof(long long), row[n], f);
	fwrite(row, sizeof(long long), n + 1, f);
	fwrite(val, sizeof(double), row[n], f);
	fclose(f);
	return 0;
}

long long cutLowerTriangleCOO(long long nz, long long * I, long long * J, double * val, long long * nzU, long long ** IU, long long ** JU, double ** valU)
{
	long long i, j;
	(*nzU) = 0;
	for (i = 0; i < nz; i++) {
		if (J[i] >= I[i])
			(*nzU)++;
	}
	(*IU) = (long long*)malloc((*nzU) * sizeof(long long));
	(*JU) = (long long*)malloc((*nzU) * sizeof(long long));
	(*valU) = (double*)malloc((*nzU) * sizeof(double));

	j = 0;
	for (i = 0; i < nz; i++) {
		if (J[i] >= I[i]) {
			(*IU)[j] = I[i];
			(*JU)[j] = J[i];
			(*valU)[j] = val[i];
			j += 1;
		}
	}
	return 0;
}

long long cutUpperTriangleCOO(long long nz, long long * I, long long * J, double * val, long long * nzL, long long ** IL, long long ** JL, double ** valL)
{
	long long i, j;
	(*nzL) = 0;
	for (i = 0; i < nz; i++) {
		if (J[i] <= I[i])
			(*nzL)++;
	}
	(*IL) = (long long*)malloc((*nzL) * sizeof(long long));
	(*JL) = (long long*)malloc((*nzL) * sizeof(long long));
	(*valL) = (double*)malloc((*nzL) * sizeof(double));

	j = 0;
	for (i = 0; i < nz; i++) {
		if (J[i] <= I[i]) {
			(*IL)[j] = I[i];
			(*JL)[j] = J[i];
			(*valL)[j] = val[i];
			j += 1;
		}
	}
	return 0;
}

long long transposeCOO(long long nz, long long * I, long long * J, double * val, long long &nzT, long long ** IT, long long ** JT, double ** valT)
{
	nzT = nz;
	(*IT) = (long long*)malloc((nzT) * sizeof(long long));
	(*JT) = (long long*)malloc((nzT) * sizeof(long long));
	(*valT) = (double*)malloc((nzT) * sizeof(double));

	long long i;
	for (i = 0; i < nz; i++) {
		(*IT)[i] = J[i];
		(*JT)[i] = I[i];
		(*valT)[i] = val[i];
	}
	return 0;
}

long long transposeCOO(long long nz, long long * I, long long * J)
{
	long long i;
	for (i = 0; i < nz; i++) {
		long long temp;
		temp = I[i];
		I[i] = J[i];
		J[i] = temp;
	}
	return 0;
}

void mallocVectors(double** x, double** y, double** b,
	double** xCRS, double** yCRS, double** bCRS,
	double** xMKL, double** yMKL, double** bMKL,
	double** xNode, double** yNode, double** bNode,
	double** xBlock, double** yBlock, double** bBlock,
	double** bNodeUp, long long N) {

	(*x) = (double*)malloc((N) * sizeof(double));
	(*y) = (double*)malloc((N) * sizeof(double));
	(*b) = (double*)malloc((N) * sizeof(double));
	(*xCRS) = (double*)malloc((N) * sizeof(double));
	(*yCRS) = (double*)malloc((N) * sizeof(double));
	(*bCRS) = (double*)malloc((N) * sizeof(double));
	(*xMKL) = (double*)malloc((N + 1) * sizeof(double));
	(*yMKL) = (double*)malloc((N + 1) * sizeof(double));
	(*bMKL) = (double*)malloc((N + 1) * sizeof(double));
	(*xNode) = (double*)malloc((N) * sizeof(double));
	(*yNode) = (double*)malloc((N) * sizeof(double));
	(*bNode) = (double*)malloc((N) * sizeof(double));
	(*xBlock) = (double*)malloc((N) * sizeof(double));
	(*yBlock) = (double*)malloc((N) * sizeof(double));
	(*bBlock) = (double*)malloc((N) * sizeof(double));
	(*bNodeUp) = (double*)malloc((N) * sizeof(double));
}

void randVector(double* b, long long N) {
	for (long long i = 0; i < N; i++) {
		b[i] = rand() % 10 / 3.0f + 1;
		// b[i] = i + 1;
	}
}

void MKLPrepare(long long ** colU_short, long long ** colL_short, long long ** indxU_short,
	long long ** indxL_short, const long long MKLn, long long* indxU, long long* indxL,
	double* MKLbx, long long* colU, long long* colL, double* bx, long long N)
{
	*colU_short = new long long[indxU[MKLn]];
	*colL_short = new long long[indxL[MKLn]];
	*indxU_short = new long long[MKLn + 1];
	*indxL_short = new long long[MKLn + 1];

	for (long long l = 0; l < N; l++) {
		MKLbx[l + 1] = bx[l];
	}
	for (long long l = 0; l < indxU[MKLn]; l++) {
		colU[l] += 1;
		(*colU_short)[l] = colU[l];
	}
	for (long long l = 0; l < indxL[MKLn]; l++) {
		colL[l] += 1;
		(*colL_short)[l] = colL[l];
	}
	for (long long l = 0; l < MKLn + 1; l++) {
		indxU[l] += 1;
		(*indxU_short)[l] = indxU[l];
	}
	for (long long l = 0; l < MKLn + 1; l++) {
		indxL[l] += 1;
		(*indxL_short)[l] = indxL[l];
	}
}

void freeMem(double** xUp, double** xLow, double** xCCSUp, double** xCCSLow,
	double** bUp, double** bLow, double** bCCSUp, double** bCCSLow,
	double** bMKLLow, double** bMKLUp, double** xMKLLow, double** xMKLUp,
	long long ** IU, long long ** JU, long long ** colU, long long ** colL,
	double ** valU, double ** valCrsU, double ** valCrsL, long long ** indxU,
	long long ** indxL, double** valLowCCS, double** valUpCCS, long long** UpRowCCS,
	long long** LowRowCCS, long long** UpIndxCCS, long long** LowIndxCCS) {
	free(*(xUp));
	free((*xLow));

	free(*(xCCSUp));
	free((*xCCSLow));
	free((*bUp));
	free((*bLow));
	free((*bCCSLow));
	free((*bMKLLow));
	free((*bMKLUp));
	free((*xMKLLow));
	free((*xMKLUp));

	free((*IU));
	free((*JU));
	free((*valU));
	free((*valCrsU));
	free((*valCrsL));
	free((*colL));
	free((*colU));
	free((*indxU));
	free((*indxL));
	free((*valUpCCS));
	free((*UpRowCCS));
	free((*LowRowCCS));
	free((*UpIndxCCS));
	free((*LowIndxCCS));
	free((*valLowCCS));
}

void mallocMatrix(double*** a, long long N) {
	(*a) = new double *[N];
	for (long long i = 0; i < N; i++) {
		(*a)[i] = new double[N];
	}
}

void printSparceCOO(long long n, long long nz, long long * I, long long * J, double * val)
{
	long long i;
	for (i = 0; i < nz; i++)
		printf("%lf  %d  %d \n", I[i], J[i], val[i]);
}

long long countZeroDiag(long long * I, long long * J, long long nz, long long N)
{
	long long count = 0;
	for (long long i = 0; i < nz; i++)
		if (I[i] == J[i])
			count++;
	return N - count;
}

void getZerosDiagNumbers(long long * I, long long * J, long long nz, long long N, long long count, long long * addDiag)
{
	long long ind = 0;
	long long* indexes = (long long*)malloc(N * sizeof(long long));
	for (long long i = 0; i < N; i++)
	{
		indexes[i] = 0;
	}
	for (long long i = 0; i < nz; i++)
		if (I[i] == J[i])
			indexes[I[i]] = 1;
	for (long long i = 0; i < N; i++)
		if (indexes[i] == 0) {
			addDiag[ind] = i;
			ind += 1;
		}
}

void fillDiag(long long * I, long long * J, double * val, long long * Inew, long long * Jnew, double * valNew, long long N, long long nz, long long nzNew, long long * addDiag) {
	for (long long i = 0; i < nz; i++) {
		Inew[i] = I[i];
		Jnew[i] = J[i];
		valNew[i] = val[i];
	}
	for (long long i = nz; i < nzNew; i++) {
		Inew[i] = addDiag[i - nz];
		Jnew[i] = addDiag[i - nz];
		valNew[i] = 1.0;
	}
	for (long long i = nz; i < nzNew; i++) {
		long long j = i;
		while (j > 0 && Inew[j] <= Inew[j - 1]) {
			long long temp = Inew[j];
			Inew[j] = Inew[j - 1];
			Inew[j - 1] = temp;

			temp = Jnew[j];
			Jnew[j] = Jnew[j - 1];
			Jnew[j - 1] = temp;

			double tempD = valNew[j];
			valNew[j] = valNew[j - 1];
			valNew[j - 1] = tempD;

			j--;
		}
	}
}

long long CheckSolv(long long n, double * x, double * xCheck) {
	double eps = 0.1;
	long long i, errors = 0;
	for (i = 0; i < n; i++) {
		// if (i % 100000 == 0)printf(" i = %d , x[%d] = %lf, xCheck[%d] = %lf, \n\n", i, i, x[i], i, xCheck[i]);
		if (((x[i] - xCheck[i]) > eps) || ((x[i] - xCheck[i]) < -eps)) {
			errors++;
		}
	}
	return errors;
}

double absError(long long n, double * x, double * xCheck)
{
	double abserror = 0.0;
	long long i;
	for (i = 0; i < n; i++) {
		//prlong longf("%.10f --- %.10f\n", x[i], xCheck[i]);
		if (fabs(x[i] - xCheck[i]) > fabs(abserror)) {
			abserror = fabs(x[i] - xCheck[i]);
		}
	}
	return abserror;
}

void gaussBackLow(long long startNode, long long finishNode, double * y,
	double * b, double * valCrsL, long long * colL, long long * indxL) {
	long long j, k, index_top, index_low;
	double sum;
	y[startNode] = b[startNode] / valCrsL[indxL[startNode + 1] - 1];
	for (k = startNode; k < finishNode; k++) {///k = startNode + 1   / k = 1
		sum = 0;
		index_top = indxL[k + 1] - 1;
		index_low = indxL[k];
		for (j = index_top - 1; j >= index_low; j--) {
			sum += valCrsL[j] * y[colL[j]];
		}
		y[k] = (b[k] - sum) / valCrsL[indxL[k + 1] - 1];
	}
	return;
}

void gaussBackUp(long long startNode, long long finishNode, double * x,
	double * y, double * valCrsU, long long * colU, long long * indxU) {
	long long j, k, index_top, index_low;
	double sum;
	x[startNode - 1] = y[startNode - 1] / valCrsU[indxU[startNode] - 1];
	for (k = startNode - 1; k >= finishNode; k--) {
		sum = 0;
		index_top = indxU[k + 1];
		index_low = indxU[k];
		for (j = index_low + 1; j < index_top; j++) {
			sum += valCrsU[j] * x[colU[j]];
		}
		x[k] = (y[k] - sum) / valCrsU[index_low];
	}
	return;
}

void matrixMultVector(long long N, double * x, double * xCheck, double * valCrs, long long * col, long long * indx)
{
	long long i, j, index_top, index_low;
	double sum;
	for (i = 0; i < N; i++) {
		sum = 0;
		index_low = indx[i]; index_top = indx[i + 1];
		for (j = index_low; j < index_top; j++) {
			sum += valCrs[j] * x[col[j]];
		}
		xCheck[i] = sum;
	}
}

void checkAndFillDiag(long long** I, long long** J, long long &nz, long long N, double** val) {
	long long zeroDiagCount = countZeroDiag((*I), (*J), nz, N);
	long long *addDiag = (long long*)malloc((zeroDiagCount + 1) * sizeof(long long)); //
	long long *Inew = (long long*)malloc((nz + zeroDiagCount + 1) * sizeof(long long));
	long long *Jnew = (long long*)malloc((nz + zeroDiagCount + 1) * sizeof(long long));
	double *valNew = (double*)malloc((nz + zeroDiagCount + 1) * sizeof(double));
	getZerosDiagNumbers((*I), (*J), nz, N, zeroDiagCount, addDiag); //получаем номера строк, где нули на главной диагонали
	fillDiag((*I), (*J), (*val), Inew, Jnew, valNew, N, nz, nz + zeroDiagCount, addDiag);//избавляемя от нулей на гл диагонали
	free(*I); (*I) = Inew;
	free(*J); (*J) = Jnew;
	free(*val); (*val) = valNew;
	nz = nz + zeroDiagCount;
	printf("nzf = %d \n", nz);
}

void mallocMatrixCCS(double** AVal, long long** ACol, long long** AIndx, long long nz, long long N) {
	(*AVal) = new double[nz];
	(*ACol) = new long long[nz];
	(*AIndx) = new long long[N + 1];
}

void makeBlock6x6LowCCS(double * MatrixLowVal, long long * MatrixLowRow, long long * MatrixLowIndx, long long NzL, long long N)
{
	for (long long i = 0; i < NzL; i++)
		MatrixLowVal[i] = 1;

	MatrixLowRow[0] = 0; MatrixLowRow[6] = 5; MatrixLowRow[12] = 4;
	MatrixLowRow[1] = 1; MatrixLowRow[7] = 2; MatrixLowRow[13] = 5;
	MatrixLowRow[2] = 4; MatrixLowRow[8] = 3; MatrixLowRow[14] = 5;
	MatrixLowRow[3] = 5; MatrixLowRow[9] = 5;
	MatrixLowRow[4] = 1; MatrixLowRow[10] = 3;
	MatrixLowRow[5] = 4; MatrixLowRow[11] = 5;

	MatrixLowIndx[0] = 0;
	MatrixLowIndx[1] = 4;
	MatrixLowIndx[2] = 7;
	MatrixLowIndx[3] = 10;
	MatrixLowIndx[4] = 12;
	MatrixLowIndx[5] = 14;
	MatrixLowIndx[6] = 16;
}

void makeBlock6x6UpCCS(double * MatrixUpVal, long long * MatrixUpRow, long long * MatrixUpIndx, long long NzU, long long N)
{
	for (long long i = 0; i < NzU; i++)
		MatrixUpVal[i] = 1;

	MatrixUpRow[0] = 0; MatrixUpRow[6] = 0; MatrixUpRow[12] = 3;
	MatrixUpRow[1] = 0; MatrixUpRow[7] = 1; MatrixUpRow[13] = 4;
	MatrixUpRow[2] = 1; MatrixUpRow[8] = 4; MatrixUpRow[14] = 5;
	MatrixUpRow[3] = 2; MatrixUpRow[9] = 0;
	MatrixUpRow[4] = 2; MatrixUpRow[10] = 1;
	MatrixUpRow[5] = 3; MatrixUpRow[11] = 2;

	MatrixUpIndx[0] = 0;
	MatrixUpIndx[1] = 1;
	MatrixUpIndx[2] = 3;
	MatrixUpIndx[3] = 4;
	MatrixUpIndx[4] = 6;
	MatrixUpIndx[5] = 9;
	MatrixUpIndx[6] = 16;
}

void CalcSuperNodesLowCCS(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx,
	long long **SNodesLow, long long &NodesNLow, long long NzL, long long N) {
	long long *TmpNodes = new long long[N];
	TmpNodes[0] = 0;
	NodesNLow = 1;
	long long lowIndx = 0;
	long long upIndx = 0;
	for (long long j = 0; j < N - 1; j++) { //по всем столбцам/////////////////
		if ((MatrixLowIndx[j + 1] - MatrixLowIndx[j] - 1) == (MatrixLowIndx[j + 2] - MatrixLowIndx[j + 1])) {

			lowIndx = MatrixLowIndx[j];
			upIndx = MatrixLowIndx[j + 1];
			if ((upIndx - lowIndx != 1) && (MatrixLowRow[lowIndx + 1] == MatrixLowRow[lowIndx] + 1)) {
				for (long long i = lowIndx + 1; i < upIndx; i++) {
					if (MatrixLowRow[i] != MatrixLowRow[i + (upIndx - lowIndx - 1)]) {
						TmpNodes[NodesNLow] = j + 1;
						NodesNLow++;
						break;
					}
				}
			}
			else {
				TmpNodes[NodesNLow] = j + 1;
				NodesNLow++;
			}
		}
		else {
			TmpNodes[NodesNLow] = j + 1;
			NodesNLow++;
		}
	}
	NodesNLow++;//for last = arr.length
	(*SNodesLow) = new long long[NodesNLow];
	for (long long i = 0; i < NodesNLow - 1; i++)
		(*SNodesLow)[i] = TmpNodes[i];
	(*SNodesLow)[NodesNLow - 1] = N;
	delete[]TmpNodes;
}

void CalcSuperNodesUpCCS(double * MatrixUpVal, long long * MatrixUpRow, long long * MatrixUpIndx,
	long long ** SNodesUp, long long & NodesNUp, long long NzU, long long N) {

	long long *TmpNodes = new long long[N];
	TmpNodes[0] = N;
	NodesNUp = 1;
	long long lowIndx = 0;
	long long upIndx = 0;
	for (long long j = N - 1; j > 0; j--) { //по всем столбцам
		if ((MatrixUpIndx[j + 1] - MatrixUpIndx[j] - 1) == (MatrixUpIndx[j] - MatrixUpIndx[j - 1])) {
			upIndx = MatrixUpIndx[j + 1];
			lowIndx = MatrixUpIndx[j];
			for (long long i = upIndx - 1; i > lowIndx; i--) {
				if (MatrixUpRow[i - 1] != MatrixUpRow[i - (upIndx - lowIndx)]) {
					TmpNodes[NodesNUp] = j;
					NodesNUp++;
					break;
				}
			}
		}
		else {
			TmpNodes[NodesNUp] = j;
			NodesNUp++;
		}
	}
	NodesNUp++;//for last = arr.length
	(*SNodesUp) = new long long[NodesNUp];
	for (long long i = 0; i < NodesNUp - 1; i++)
		(*SNodesUp)[i] = TmpNodes[i];

	(*SNodesUp)[NodesNUp - 1] = 0;
	delete[]TmpNodes;
	return;
}

void gaussBackBlockLow(long long startNode, long long finishNode, double * y,
	double * b, double * valCrsL, long long * colL, long long * indxL) {
	long long j, k, index_top, index_low;
	double sum;
	long long m = 1;
	y[startNode] = b[startNode] / valCrsL[indxL[startNode + 1] - 1];
	for (k = startNode + 1; k < finishNode; k++) {// по строкам нода
		sum = 0;
		index_top = indxL[k + 1] - 1;
		index_low = indxL[k + 1] - 1 - m;
		for (j = index_top - 1; j >= index_low; j--) {
			sum += valCrsL[j] * y[colL[j]];
		}
		y[k] = (b[k] - sum) / valCrsL[indxL[k + 1] - 1];
		m++;
	}
	return;
}

void gaussBackBlockUp(long long startNode, long long finishNode, double * x,
	double * y, double * valCrsU, long long * colU, long long * indxU) {
	long long j, k, index_top, index_low;
	double sum;
	long long m = 2;
	x[startNode - 1] = y[startNode - 1] / valCrsU[indxU[startNode - 1]];
	for (k = startNode - 2; k >= finishNode; k--) {// по строкам нода
		sum = 0;
		index_top = indxU[k] + m;
		index_low = indxU[k] + 1;
		for (j = index_low; j < index_top; j++) {
			sum += valCrsU[j] * x[colU[j]];
		}
		if (valCrsU[indxU[k]] == 0) {
			printf("valCrsU[indxU[%d]] = %lf\n", k, valCrsU[indxU[k]]);
		}
		x[k] = (y[k] - sum) / valCrsU[indxU[k]];
		m++;
	}
	return;
}

void changeBLowCCS(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx,
	long long SNodesLowl, long long SNodesLowl_1, long long N, double* x, double* bLow) {
	for (int j = SNodesLowl; j < SNodesLowl_1; j++) {//по всем столбам
		for (int i = MatrixLowIndx[j]; i < MatrixLowIndx[j + 1]; i++) {//вдоль столбца
			bLow[MatrixLowRow[i]] -= MatrixLowVal[i] * x[j];
		}
	}
}

void changeBUpCCS(double *MatrixUpVal, long long *MatrixUpRow, long long *MatrixUpIndx,
	long long SNodesUpl, long long SNodesUpl_1, long long N, double* x, double* bUp) {
	for (int j = SNodesUpl_1; j < SNodesUpl; j++) {//по всем столбам
		for (int i = MatrixUpIndx[j]; i < MatrixUpIndx[j + 1]; i++) {//вдоль столбца
			bUp[MatrixUpRow[i]] -= MatrixUpVal[i] * x[j];
		}
	}
}

void blockSolverLowCCS(double *MatrixLowValCRS, long long *MatrixLowColCRS, long long *MatrixLowIndxCRS,
	double *MatrixLowValCCS, long long *MatrixRowColCCS, long long *MatrixLowIndxCCS,
	long long* SNodesLow, long long NodesNLow, double* xLow, double* bLow, long long N) {
	for (long long l = 0; l < NodesNLow - 1; l++) { //по всем нодам
		gaussBackBlockLow(SNodesLow[l], SNodesLow[l + 1], xLow, bLow, MatrixLowValCRS, MatrixLowColCRS, MatrixLowIndxCRS);
		changeBLowCCS(MatrixLowValCCS, MatrixRowColCCS, MatrixLowIndxCCS, SNodesLow[l], SNodesLow[l + 1], N, xLow, bLow);
	}
}

void blockSolverUpCCS(double *MatrixUpValCRS, long long *MatrixUpColCRS, long long *MatrixUpIndxCRS,
	double *MatrixUpValCCS, long long *MatrixUpColCCS, long long *MatrixUpIndxCCS,
	long long* SNodesUp, long long NodesNUp, double* xUp, double* bUp, long long N) {

	for (long long l = 0; l < NodesNUp - 1; l++) { //по всем нодам
		gaussBackBlockUp(SNodesUp[l], SNodesUp[l + 1], xUp, bUp, MatrixUpValCRS, MatrixUpColCRS, MatrixUpIndxCRS);
		changeBUpCCS(MatrixUpValCCS, MatrixUpColCCS, MatrixUpIndxCCS, SNodesUp[l], SNodesUp[l + 1], N, xUp, bUp);
	}
}

void nodeSolverLow(double *MatrixLowVal, long long *MatrixLowRow, long long *MatrixLowIndx,
	long long* SNodesLow, long long NodesNLow, double* xLow, double* bLow, long long N) {
	for (long long l = 0; l < NodesNLow - 1; l++) { //по всем нодам
		gaussBackLow(SNodesLow[l], SNodesLow[l + 1], xLow, bLow, MatrixLowVal, MatrixLowRow, MatrixLowIndx);
	}
	return;
}

void nodeSolverUp(double *MatrixUpVal, long long *MatrixUpRow, long long *MatrixUpIndx,
	long long* SNodesUp, long long NodesNUp, double* xUp, double* bUp, long long N) {
	for (long long l = 0; l < NodesNUp - 1; l++) { //по всем нодам
		gaussBackUp(SNodesUp[l], SNodesUp[l + 1], xUp, bUp, MatrixUpVal, MatrixUpRow, MatrixUpIndx);
	}
	return;
}

void COOtoCCS(long long N, long long nz, long long *I, long long *J, double *valCOO, long long **indx, long long **row, double **valCcs) {
	long long i;
	long long *places;

	(*indx) = (long long *)malloc((N + 1) * sizeof(long long));
	(*row) = (long long *)malloc((nz) * sizeof(long long));
	(*valCcs) = (double *)malloc((nz) * sizeof(double));
	places = (long long *)malloc((N + 1) * sizeof(long long));

	for (i = 0; i < N + 1; i++) {
		(*indx)[i] = 0;
	}
	for (i = 0; i < nz; i++) {
		(*indx)[J[i] + 1] += 1;
	}
	for (i = 0; i < N; i++) {
		(*indx)[i + 1] += (*indx)[i];
		places[i] = (*indx)[i];
	}
	for (i = 0; i < nz; i++) {
		(*row)[places[J[i]]] = I[i];
		(*valCcs)[places[J[i]]] = valCOO[i];
		places[J[i]]++;
	}
	//sort each row
	for (i = 0; i < N; i++) {
		long long start, finish, j, k;
		start = (*indx)[i];
		finish = (*indx)[i + 1];
		//need to write bubble sort
		for (j = start; j < finish; j++)
			for (k = start; k < finish - j - 1; k++)
				if ((*row)[k] > (*row)[k + 1]) {
					long long tempCol = (*row)[k];
					(*row)[k] = (*row)[k + 1];
					(*row)[k + 1] = tempCol;

					double tempVal = (*valCcs)[k];
					(*valCcs)[k] = (*valCcs)[k + 1];
					(*valCcs)[k + 1] = tempVal;
				}
	}
	free(places);
}

void makeMatrix6x6COO(long long * I, long long * J, double * COOVal, long long NzL) {
	I[0] = 0;	J[0] = 0;	COOVal[0] = 1;
	I[1] = 0;	J[1] = 1;	COOVal[1] = 1;
	I[2] = 0;	J[2] = 4;	COOVal[2] = 1;
	I[3] = 0;	J[3] = 5;	COOVal[3] = 1;
	I[4] = 1;	J[4] = 1;	COOVal[4] = 1;
	I[5] = 1;	J[5] = 4;	COOVal[5] = 1;
	I[6] = 1;	J[6] = 5;	COOVal[6] = 1;
	I[7] = 2;	J[7] = 2;	COOVal[7] = 1;
	I[8] = 2;	J[8] = 3;	COOVal[8] = 1;
	I[9] = 2;	J[9] = 5;	COOVal[9] = 1;
	I[10] = 3;	J[10] = 3;	COOVal[10] = 1;
	I[11] = 3;	J[11] = 5;	COOVal[11] = 1;
	I[12] = 4;	J[12] = 4;	COOVal[12] = 1;
	I[13] = 4;	J[13] = 5;	COOVal[13] = 1;
	I[14] = 5;	J[14] = 5;	COOVal[14] = 1;
}

void mallocMatrixCOO(long long ** I, long long ** J, double ** COOVal, long long NzL) {
	(*I) = new long long[NzL];
	(*J) = new long long[NzL];
	(*COOVal) = new double[NzL];
}

void makeBlockMatrix8x8LowCOO(long long * I, long long * J, double * COOVal, long long NzL) {
	I[0] = 0;	J[0] = 0;	COOVal[0] = 1;
	I[1] = 1;	J[1] = 0;	COOVal[1] = 1;
	I[2] = 1;	J[2] = 1;	COOVal[2] = 1;
	I[3] = 2;	J[3] = 0;	COOVal[3] = 1;
	I[4] = 2;	J[4] = 1;	COOVal[4] = 1;
	I[5] = 2;	J[5] = 2;	COOVal[5] = 1;
	I[6] = 3;	J[6] = 2;	COOVal[6] = 1;
	I[7] = 3;	J[7] = 3;	COOVal[7] = 1;
	I[8] = 4;	J[8] = 2;	COOVal[8] = 1;
	I[9] = 4;	J[9] = 3;	COOVal[9] = 1;
	I[10] = 4;	J[10] = 4;	COOVal[10] = 1;
	I[11] = 5;	J[11] = 0;	COOVal[11] = 1;
	I[12] = 5;	J[12] = 1;	COOVal[12] = 1;
	I[13] = 5;	J[13] = 5;	COOVal[13] = 1;
	I[14] = 6;	J[14] = 2;	COOVal[14] = 1;
	I[15] = 6;	J[15] = 3;	COOVal[15] = 1;
	I[16] = 6;	J[16] = 6;	COOVal[16] = 1;
	I[17] = 7;	J[17] = 4;	COOVal[17] = 1;
	I[18] = 7;	J[18] = 6;	COOVal[18] = 1;
	I[19] = 7;	J[19] = 7;	COOVal[19] = 1;
}

void makeBlockMatrix8x8UpCOO(long long * I, long long * J, double * COOVal, long long NzL) {
	I[0] = 0;	J[0] = 0;	COOVal[0] = 1;
	I[1] = 0;	J[1] = 1;	COOVal[1] = 1;
	I[2] = 0;	J[2] = 2;	COOVal[2] = 1;
	I[3] = 0;	J[3] = 5;	COOVal[3] = 1;
	I[4] = 1;	J[4] = 1;	COOVal[4] = 1;
	I[5] = 1;	J[5] = 2;	COOVal[5] = 1;
	I[6] = 1;	J[6] = 5;	COOVal[6] = 1;
	I[7] = 2;	J[7] = 2;	COOVal[7] = 1;
	I[8] = 2;	J[8] = 3;	COOVal[8] = 1;
	I[9] = 2;	J[9] = 4;	COOVal[9] = 1;
	I[10] = 2;	J[10] = 6;	COOVal[10] = 1;
	I[11] = 3;	J[11] = 3;	COOVal[11] = 1;
	I[12] = 3;	J[12] = 4;	COOVal[12] = 1;
	I[13] = 3;	J[13] = 6;	COOVal[13] = 1;
	I[14] = 4;	J[14] = 4;	COOVal[14] = 1;
	I[15] = 4;	J[15] = 7;	COOVal[15] = 1;
	I[16] = 5;	J[16] = 5;	COOVal[16] = 1;
	I[17] = 6;	J[17] = 6;	COOVal[17] = 1;
	I[18] = 6;	J[18] = 7;	COOVal[18] = 1;
	I[19] = 7;	J[19] = 7;	COOVal[19] = 1;
}

void makeBlock8x8LowCCS(double * MatrixLowValCCS, long long * MatrixLowRowCCS, long long * MatrixLowIndxCCS, long long NzL, long long N) {
	for (long long i = 0; i < NzL; i++)
		MatrixLowValCCS[i] = 1;

	MatrixLowRowCCS[0] = 0; MatrixLowRowCCS[10] = 6;
	MatrixLowRowCCS[1] = 1; MatrixLowRowCCS[11] = 3;
	MatrixLowRowCCS[2] = 2; MatrixLowRowCCS[12] = 4;
	MatrixLowRowCCS[3] = 5; MatrixLowRowCCS[13] = 6;
	MatrixLowRowCCS[4] = 1; MatrixLowRowCCS[14] = 4;
	MatrixLowRowCCS[5] = 2; MatrixLowRowCCS[15] = 7;
	MatrixLowRowCCS[6] = 5; MatrixLowRowCCS[16] = 5;
	MatrixLowRowCCS[7] = 2; MatrixLowRowCCS[17] = 6;
	MatrixLowRowCCS[8] = 3; MatrixLowRowCCS[18] = 7;
	MatrixLowRowCCS[9] = 4; MatrixLowRowCCS[19] = 7;

	MatrixLowIndxCCS[0] = 0;
	MatrixLowIndxCCS[1] = 4;
	MatrixLowIndxCCS[2] = 7;
	MatrixLowIndxCCS[3] = 11;
	MatrixLowIndxCCS[4] = 14;
	MatrixLowIndxCCS[5] = 16;
	MatrixLowIndxCCS[6] = 17;
	MatrixLowIndxCCS[7] = 19;
	MatrixLowIndxCCS[8] = 20;
}

void makeBlock8x8UpCCS(double * MatrixUpValCCS, long long * MatrixUpRowCCS, long long * MatrixUpIndxCCS, long long NzU, long long N) {
	for (long long i = 0; i < NzU; i++)
		MatrixUpValCCS[i] = 1;

	MatrixUpRowCCS[0] = 0; MatrixUpRowCCS[10] = 4;
	MatrixUpRowCCS[1] = 0; MatrixUpRowCCS[11] = 0;
	MatrixUpRowCCS[2] = 1; MatrixUpRowCCS[12] = 1;
	MatrixUpRowCCS[3] = 0; MatrixUpRowCCS[13] = 5;
	MatrixUpRowCCS[4] = 1; MatrixUpRowCCS[14] = 2;
	MatrixUpRowCCS[5] = 2; MatrixUpRowCCS[15] = 3;
	MatrixUpRowCCS[6] = 2; MatrixUpRowCCS[16] = 6;
	MatrixUpRowCCS[7] = 3; MatrixUpRowCCS[17] = 4;
	MatrixUpRowCCS[8] = 2; MatrixUpRowCCS[18] = 6;
	MatrixUpRowCCS[9] = 3; MatrixUpRowCCS[19] = 7;

	MatrixUpIndxCCS[0] = 0;
	MatrixUpIndxCCS[1] = 1;
	MatrixUpIndxCCS[2] = 3;
	MatrixUpIndxCCS[3] = 6;
	MatrixUpIndxCCS[4] = 8;
	MatrixUpIndxCCS[5] = 11;
	MatrixUpIndxCCS[6] = 14;
	MatrixUpIndxCCS[7] = 17;
	MatrixUpIndxCCS[8] = 20;
}

void makeBlock8x8LowCRS(double * MatrixLowValCRS, long long * MatrixLowColCRS, long long * MatrixLowIndxCRS, long long NzL, long long N) {
	for (long long i = 0; i < NzL; i++)
		MatrixLowValCRS[i] = 1;

	MatrixLowColCRS[0] = 0; MatrixLowColCRS[10] = 4;
	MatrixLowColCRS[1] = 0; MatrixLowColCRS[11] = 0;
	MatrixLowColCRS[2] = 1; MatrixLowColCRS[12] = 1;
	MatrixLowColCRS[3] = 0; MatrixLowColCRS[13] = 5;
	MatrixLowColCRS[4] = 1; MatrixLowColCRS[14] = 2;
	MatrixLowColCRS[5] = 2; MatrixLowColCRS[15] = 3;
	MatrixLowColCRS[6] = 2; MatrixLowColCRS[16] = 6;
	MatrixLowColCRS[7] = 3; MatrixLowColCRS[17] = 4;
	MatrixLowColCRS[8] = 2; MatrixLowColCRS[18] = 6;
	MatrixLowColCRS[9] = 3; MatrixLowColCRS[19] = 7;

	MatrixLowIndxCRS[0] = 0;
	MatrixLowIndxCRS[1] = 1;
	MatrixLowIndxCRS[2] = 3;
	MatrixLowIndxCRS[3] = 6;
	MatrixLowIndxCRS[4] = 8;
	MatrixLowIndxCRS[5] = 11;
	MatrixLowIndxCRS[6] = 14;
	MatrixLowIndxCRS[7] = 17;
	MatrixLowIndxCRS[8] = 20;
}

void makeBlock8x8UpCRS(double * MatrixUpValCRS, long long * MatrixUpColCRS, long long * MatrixUpIndxCRS, long long NzU, long long N) {
	for (long long i = 0; i < NzU; i++)
		MatrixUpValCRS[i] = 1;

	MatrixUpColCRS[0] = 0; MatrixUpColCRS[10] = 6;
	MatrixUpColCRS[1] = 1; MatrixUpColCRS[11] = 3;
	MatrixUpColCRS[2] = 2; MatrixUpColCRS[12] = 4;
	MatrixUpColCRS[3] = 5; MatrixUpColCRS[13] = 6;
	MatrixUpColCRS[4] = 1; MatrixUpColCRS[14] = 4;
	MatrixUpColCRS[5] = 2; MatrixUpColCRS[15] = 7;
	MatrixUpColCRS[6] = 5; MatrixUpColCRS[16] = 5;
	MatrixUpColCRS[7] = 2; MatrixUpColCRS[17] = 6;
	MatrixUpColCRS[8] = 3; MatrixUpColCRS[18] = 7;
	MatrixUpColCRS[9] = 4; MatrixUpColCRS[19] = 7;

	MatrixUpIndxCRS[0] = 0;
	MatrixUpIndxCRS[1] = 4;
	MatrixUpIndxCRS[2] = 7;
	MatrixUpIndxCRS[3] = 11;
	MatrixUpIndxCRS[4] = 14;
	MatrixUpIndxCRS[5] = 16;
	MatrixUpIndxCRS[6] = 17;
	MatrixUpIndxCRS[7] = 19;
	MatrixUpIndxCRS[8] = 20;
}

void makeBlockMatrix12x12LowCOO(long long * I, long long * J, double * COOVal, long long NzL) {
	I[0] = 0;	J[0] = 0;	COOVal[0] = 1;
	I[1] = 1;	J[1] = 0;	COOVal[1] = 1;
	I[2] = 1;	J[2] = 1;	COOVal[2] = 1;
	I[3] = 2;	J[3] = 2;	COOVal[3] = 1;
	I[4] = 3;	J[4] = 0;	COOVal[4] = 1;
	I[5] = 3;	J[5] = 1;	COOVal[5] = 1;
	I[6] = 3;	J[6] = 3;	COOVal[6] = 1;
	I[7] = 4;	J[7] = 0;	COOVal[7] = 1;
	I[8] = 4;	J[8] = 1;	COOVal[8] = 1;
	I[9] = 4;	J[9] = 3;	COOVal[9] = 1;
	I[10] = 4;	J[10] = 4;	COOVal[10] = 1;
	I[11] = 5;	J[11] = 3;	COOVal[11] = 1;
	I[12] = 5;	J[12] = 4;	COOVal[12] = 1;
	I[13] = 5;	J[13] = 5;	COOVal[13] = 1;
	I[14] = 6;	J[14] = 6;	COOVal[14] = 1;
	I[15] = 7;	J[15] = 0;	COOVal[15] = 1;
	I[16] = 7;	J[16] = 1;	COOVal[16] = 1;
	I[17] = 7;	J[17] = 3;	COOVal[17] = 1;
	I[18] = 7;	J[18] = 4;	COOVal[18] = 1;
	I[19] = 7;	J[19] = 5;	COOVal[19] = 1;
	I[20] = 7;	J[20] = 7;	COOVal[20] = 1;
	I[21] = 8;	J[21] = 8;	COOVal[21] = 1;
	I[22] = 9;	J[22] = 3;	COOVal[22] = 1;
	I[23] = 9;	J[23] = 4;	COOVal[23] = 1;
	I[24] = 9;	J[24] = 5;	COOVal[24] = 1;
	I[25] = 9;	J[25] = 6;	COOVal[25] = 1;
	I[26] = 9;	J[26] = 8;	COOVal[26] = 1;
	I[27] = 9;	J[27] = 9;	COOVal[27] = 1;
	I[28] = 10;	J[28] = 0;	COOVal[28] = 1;
	I[29] = 10;	J[29] = 1;	COOVal[29] = 1;
	I[30] = 10;	J[30] = 6;	COOVal[30] = 1;
	I[31] = 10;	J[31] = 8;	COOVal[31] = 1;
	I[32] = 10;	J[32] = 9;	COOVal[32] = 1;
	I[33] = 10;	J[33] = 10;	COOVal[33] = 1;
	I[34] = 11;	J[34] = 8;	COOVal[34] = 1;
	I[35] = 11;	J[35] = 9;	COOVal[35] = 1;
	I[36] = 11;	J[36] = 10;	COOVal[36] = 1;
	I[37] = 11;	J[37] = 11;	COOVal[37] = 1;
}

void makeBlockMatrix12x12UpCOO(long long * I, long long * J, double * COOVal, long long NzU) {
	I[0] = 0;	J[0] = 0;	COOVal[0] = 1;
	I[1] = 0;	J[1] = 1;	COOVal[1] = 1;
	I[2] = 0;	J[2] = 5;	COOVal[2] = 1;
	I[3] = 0;	J[3] = 6;	COOVal[3] = 1;
	I[4] = 0;	J[4] = 7;	COOVal[4] = 1;
	I[5] = 1;	J[5] = 1;	COOVal[5] = 1;
	I[6] = 1;	J[6] = 3;	COOVal[6] = 1;
	I[7] = 1;	J[7] = 4;	COOVal[7] = 1;
	I[8] = 2;	J[8] = 2;	COOVal[8] = 1;
	I[9] = 2;	J[9] = 8;	COOVal[9] = 1;
	I[10] = 2;	J[10] = 9;	COOVal[10] = 1;
	I[11] = 2;	J[11] = 10;	COOVal[11] = 1;
	I[12] = 2;	J[12] = 11;	COOVal[12] = 1;
	I[13] = 3;	J[13] = 3;	COOVal[13] = 1;
	I[14] = 3;	J[14] = 4;	COOVal[14] = 1;
	I[15] = 4;	J[15] = 4;	COOVal[15] = 1;
	I[16] = 5;	J[16] = 5;	COOVal[16] = 1;
	I[17] = 5;	J[17] = 6;	COOVal[17] = 1;
	I[18] = 5;	J[18] = 7;	COOVal[18] = 1;
	I[19] = 6;	J[19] = 6;	COOVal[19] = 1;
	I[20] = 6;	J[20] = 7;	COOVal[20] = 1;
	I[21] = 6;	J[21] = 8;	COOVal[21] = 1;
	I[22] = 6;	J[22] = 9;	COOVal[22] = 1;
	I[23] = 6;	J[23] = 10;	COOVal[23] = 1;
	I[24] = 6;	J[24] = 11;	COOVal[24] = 1;
	I[25] = 7;	J[25] = 7;	COOVal[25] = 1;
	I[26] = 8;	J[26] = 8;	COOVal[26] = 1;
	I[27] = 8;	J[27] = 9;	COOVal[27] = 1;
	I[28] = 8;	J[28] = 10;	COOVal[28] = 1;
	I[29] = 8;	J[29] = 11;	COOVal[29] = 1;
	I[30] = 9;	J[30] = 9;	COOVal[30] = 1;
	I[31] = 9;	J[31] = 10;	COOVal[31] = 1;
	I[32] = 9;	J[32] = 11;	COOVal[32] = 1;
	I[33] = 10;	J[33] = 10;	COOVal[33] = 1;
	I[34] = 10;	J[34] = 11;	COOVal[34] = 1;
	I[35] = 11;	J[35] = 11;	COOVal[35] = 1;
}

void makeBlockMatrix12x12LowCOORandom(long long * I, long long * J, double * COOVal, long long NzL) {
	I[0] = 0;	J[0] = 0;	COOVal[0] = 1;
	I[1] = 1;	J[1] = 0;	COOVal[1] = 2;
	I[2] = 1;	J[2] = 1;	COOVal[2] = 3;
	I[3] = 2;	J[3] = 2;	COOVal[3] = 4;
	I[4] = 3;	J[4] = 0;	COOVal[4] = 5;
	I[5] = 3;	J[5] = 1;	COOVal[5] = 6;
	I[6] = 3;	J[6] = 3;	COOVal[6] = 7;
	I[7] = 4;	J[7] = 0;	COOVal[7] = 8;
	I[8] = 4;	J[8] = 1;	COOVal[8] = 9;
	I[9] = 4;	J[9] = 3;	COOVal[9] = 10;
	I[10] = 4;	J[10] = 4;	COOVal[10] = 11;
	I[11] = 5;	J[11] = 3;	COOVal[11] = 12;
	I[12] = 5;	J[12] = 4;	COOVal[12] = 13;
	I[13] = 5;	J[13] = 5;	COOVal[13] = 14;
	I[14] = 6;	J[14] = 6;	COOVal[14] = 15;
	I[15] = 7;	J[15] = 0;	COOVal[15] = 16;
	I[16] = 7;	J[16] = 1;	COOVal[16] = 17;
	I[17] = 7;	J[17] = 3;	COOVal[17] = 18;
	I[18] = 7;	J[18] = 4;	COOVal[18] = 19;
	I[19] = 7;	J[19] = 5;	COOVal[19] = 20;
	I[20] = 7;	J[20] = 7;	COOVal[20] = 21;
	I[21] = 8;	J[21] = 8;	COOVal[21] = 22;
	I[22] = 9;	J[22] = 3;	COOVal[22] = 23;
	I[23] = 9;	J[23] = 4;	COOVal[23] = 24;
	I[24] = 9;	J[24] = 5;	COOVal[24] = 25;
	I[25] = 9;	J[25] = 6;	COOVal[25] = 26;
	I[26] = 9;	J[26] = 8;	COOVal[26] = 27;
	I[27] = 9;	J[27] = 9;	COOVal[27] = 28;
	I[28] = 10;	J[28] = 0;	COOVal[28] = 29;
	I[29] = 10;	J[29] = 1;	COOVal[29] = 30;
	I[30] = 10;	J[30] = 6;	COOVal[30] = 31;
	I[31] = 10;	J[31] = 8;	COOVal[31] = 32;
	I[32] = 10;	J[32] = 9;	COOVal[32] = 33;
	I[33] = 10;	J[33] = 10;	COOVal[33] = 34;
	I[34] = 11;	J[34] = 8;	COOVal[34] = 35;
	I[35] = 11;	J[35] = 9;	COOVal[35] = 36;
	I[36] = 11;	J[36] = 10;	COOVal[36] = 37;
	I[37] = 11;	J[37] = 11;	COOVal[37] = 38;
}

void makeBlockMatrix12x12UpCOORandom(long long * I, long long * J, double * COOVal, long long NzU) {
	I[0] = 0;	J[0] = 0;	COOVal[0] = 1;
	I[1] = 0;	J[1] = 1;	COOVal[1] = 2;
	I[2] = 0;	J[2] = 5;	COOVal[2] = 3;
	I[3] = 0;	J[3] = 6;	COOVal[3] = 4;
	I[4] = 0;	J[4] = 7;	COOVal[4] = 5;
	I[5] = 1;	J[5] = 1;	COOVal[5] = 6;
	I[6] = 1;	J[6] = 3;	COOVal[6] = 7;
	I[7] = 1;	J[7] = 4;	COOVal[7] = 8;
	I[8] = 2;	J[8] = 2;	COOVal[8] = 9;
	I[9] = 2;	J[9] = 8;	COOVal[9] = 10;
	I[10] = 2;	J[10] = 9;	COOVal[10] = 11;
	I[11] = 2;	J[11] = 10;	COOVal[11] = 12;
	I[12] = 2;	J[12] = 11;	COOVal[12] = 13;
	I[13] = 3;	J[13] = 3;	COOVal[13] = 14;
	I[14] = 3;	J[14] = 4;	COOVal[14] = 15;
	I[15] = 4;	J[15] = 4;	COOVal[15] = 16;
	I[16] = 5;	J[16] = 5;	COOVal[16] = 17;
	I[17] = 5;	J[17] = 6;	COOVal[17] = 18;
	I[18] = 5;	J[18] = 7;	COOVal[18] = 19;
	I[19] = 6;	J[19] = 6;	COOVal[19] = 20;
	I[20] = 6;	J[20] = 7;	COOVal[20] = 21;
	I[21] = 6;	J[21] = 8;	COOVal[21] = 22;
	I[22] = 6;	J[22] = 9;	COOVal[22] = 23;
	I[23] = 6;	J[23] = 10;	COOVal[23] = 24;
	I[24] = 6;	J[24] = 11;	COOVal[24] = 25;
	I[25] = 7;	J[25] = 7;	COOVal[25] = 26;
	I[26] = 8;	J[26] = 8;	COOVal[26] = 27;
	I[27] = 8;	J[27] = 9;	COOVal[27] = 28;
	I[28] = 8;	J[28] = 10;	COOVal[28] = 29;
	I[29] = 8;	J[29] = 11;	COOVal[29] = 30;
	I[30] = 9;	J[30] = 9;	COOVal[30] = 31;
	I[31] = 9;	J[31] = 10;	COOVal[31] = 32;
	I[32] = 9;	J[32] = 11;	COOVal[32] = 33;
	I[33] = 10;	J[33] = 10;	COOVal[33] = 34;
	I[34] = 10;	J[34] = 11;	COOVal[34] = 35;
	I[35] = 11;	J[35] = 11;	COOVal[35] = 36;
}

void generateBigBlockMatrixL(long long * I, long long * J, double * COOVal, long long NzL, long long N,
	long long blockSize, long long* fullRowL) {
	for (long long i = 0; i < NzL; i++) {
		COOVal[i] = 1;
	}

	lldiv_t blockCount = div(N, blockSize);
	// lldiv_t.quot - целая часть
	// lldiv_t.rem - остаток

	long long k = 0;
	for (long long blockI = 0; blockI < blockCount.quot; blockI++) {
		for (long long i = 0; i < blockSize; i++) {
			if (fullRowL[i + blockI*blockSize] == 1) {
				for (long long j = 0; j <= blockI * blockSize + i; j++) {
					I[k] = blockI * blockSize + i;
					J[k] = j;
					k++;
				}
			}
			else {
				for (long long j = 0; j <= i; j++) {
					I[k] = blockI * blockSize + i;
					J[k] = blockI * blockSize + j;
					k++;
				}
			}
		}
	}
	for (long long i = N - blockCount.rem; i < N; i++) {
		if (fullRowL[i] == 1) {
			for (long long j = 0; j <= i; j++) {
				I[k] = i;
				J[k] = j;
				k++;
			}
		}
		else {
			for (long long j = N - blockCount.rem; j <= i; j++) {
				I[k] = i;
				J[k] = j;
				k++;
			}
		}
	}

	for (int i = 0; i < k; i++) {
		printf("[%d]: %d   %d  \n", i, I[i], J[i]);
	}
	return;
}

void generateBigBlockMatrixU(long long * I, long long * J, double * COOVal,
	long long NzU, long long N, long long blockSize, long long* fullRowU) {

	for (long long i = 0; i < NzU; i++) {
		COOVal[i] = 1;
	}

	lldiv_t blockCount = div(N, blockSize);

	long long k = 0;
	for (long long blockI = 0; blockI < blockCount.quot; blockI++) {
		for (long long i = 0; i < blockSize; i++) {
			if (fullRowU[i + blockI*blockSize] == 1) {
				for (long long j = i + blockI*blockSize; j < N; j++) {
					I[k] = blockI * blockSize + i;
					J[k] = j;
					k++;
				}
			}
			else {
				for (long long j = i; j < blockSize; j++) {
					I[k] = blockI * blockSize + i;
					J[k] = blockI * blockSize + j;
					k++;
				}
			}
		}
	}

	for (long long i = N - blockCount.rem; i < N; i++) {
		for (long long j = i; j < N; j++) {
			I[k] = i;
			J[k] = j;
			k++;
		}
	}

	for (int i = 0; i < k; i++) {
		printf("[%d]: %d   %d  \n", i, I[i], J[i]);
	}
	return;
}

long long calcNzL(long long N, long long blockSizeL, long long **fullRowL)
{
	(*fullRowL) = (long long *)malloc(N * sizeof(long long));

	for (int i = 0; i < N; i++) {
		(*fullRowL)[i] = rand() % 2;
	}
	(*fullRowL)[6] = 1;
	for (int i = 0; i < N; i++) {
		printf("FullRowL[%d] = %d\n", i, (*fullRowL)[i]);
	}

	lldiv_t blockCount = div(N, blockSizeL);
	long long k = 0;
	for (long long blockI = 0; blockI < blockCount.quot; blockI++) {
		for (long long i = 0; i < blockSizeL; i++) {
			if ((*fullRowL)[i + blockI*blockSizeL] == 1) {
				for (long long j = 0; j <= blockI * blockSizeL + i; j++) {
					k++;
				}
			}
			else {
				for (long long j = 0; j <= i; j++) {
					k++;
				}
			}
		}
	}
	for (long long i = N - blockCount.rem; i < N; i++) {
		if ((*fullRowL)[i] == 1) {
			for (long long j = 0; j <= i; j++) {
				k++;
			}
		}
		else {
			for (long long j = N - blockCount.rem; j <= i; j++) {
				k++;
			}
		}
	}
	return k;
}

long long calcNzU(long long N, long long blockSizeU, long long **fullRowU) {
	lldiv_t blockCount = div(N, blockSizeU);

	(*fullRowU) = (long long *)malloc(N * sizeof(long long));

	for (int i = 0; i < N; i++) {
		(*fullRowU)[i] = rand() % 2;
	}
	(*fullRowU)[6] = 1;
	for (int i = 0; i < N; i++) {
		printf("FullRowU[%d] = %d\n", i, (*fullRowU)[i]);
	}

	long long k = 0;
	for (long long blockI = 0; blockI < blockCount.quot; blockI++) {
		for (long long i = 0; i < blockSizeU; i++) {
			if (fullRowU[i + blockI*blockSizeU]) {
				for (long long j = i; j < N; j++) {
					k++;
				}
			}
			else {
				for (long long j = i; j < blockSizeU; j++) {
					k++;
				}
			}
		}
	}

	for (long long i = N - blockCount.rem; i < N; i++) {
		for (long long j = i; j < N; j++) {
			k++;
		}
	}
	return k;
}
