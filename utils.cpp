#include "utils.h"

int readMTX(const char * fileName, int ** I, int ** J, double ** val, int * M, int * N, int * nz)
{
	FILE *f;
	int i;
	int ret_code;

	if ((f = fopen(fileName, "r")) == NULL) {
		printf("Can't open file %s\n", fileName);
		return 1;
	}
	else
		printf("Matrix: %s\n", fileName);

	if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz)) != 0)
		return 1;

	(*I) = (int *)malloc((*nz) * sizeof(int));
	(*J) = (int *)malloc((*nz) * sizeof(int));
	(*val) = (double *)malloc((*nz) * sizeof(double));

	for (i = 0; i < (*nz); i++)
	{
		fscanf(f, "%d %d %lg\n", &(*I)[i], &(*J)[i], &(*val)[i]);
		(*I)[i]--;  /* adjust from 1-based to 0-based */
		(*J)[i]--;
	}
	if (f != stdin) fclose(f);
	return 0;
}

int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz)
{
	char line[MM_MAX_LINE_LENGTH];
	int num_items_read;

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

int COOtoCRS(int n, int nz, int *I, int *J, double *valCOO, int **indx, int **col, double **valCrs)
{
	int i;
	int *places;

	(*indx) = (int*)malloc((n + 1) * sizeof(int));
	(*col) = (int*)malloc((nz) * sizeof(int));
	(*valCrs) = (double*)malloc((nz) * sizeof(double));
	places = (int*)malloc(n * sizeof(int));

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
		int start, finish, j, k;
		start = (*indx)[i];
		finish = (*indx)[i + 1];
		//need to write bubble sort
		for (j = start; j < finish; j++)
			for (k = start; k < finish - j - 1; k++)
				if ((*col)[k] > (*col)[k + 1]) {
					int tempCol = (*col)[k];
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

int saveBinCRS(const char * fileName, int n, int * row, int * col, double * val)
{
	FILE *f;
	f = fopen(fileName, "wb");
	fwrite(&n, sizeof(int), 1, f);
	fwrite(&(row[n]), sizeof(int), 1, f);
	fwrite(col, sizeof(int), row[n], f);
	fwrite(row, sizeof(int), n + 1, f);
	fwrite(val, sizeof(double), row[n], f);
	fclose(f);
	return 0;
}

int cutLowerTriangleCOO(int nz, int * I, int * J, double * val, int * nzU, int ** IU, int ** JU, double ** valU)
{
	int i, j;
	(*nzU) = 0;
	for (i = 0; i < nz; i++) {
		if (J[i] >= I[i])
			(*nzU)++;
	}
	(*IU) = (int*)malloc((*nzU) * sizeof(int));
	(*JU) = (int*)malloc((*nzU) * sizeof(int));
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

int cutUpperTriangleCOO(int nz, int * I, int * J, double * val, int * nzL, int ** IL, int ** JL, double ** valL)
{
	int i, j;
	(*nzL) = 0;
	for (i = 0; i < nz; i++) {
		if (J[i] <= I[i])
			(*nzL)++;
	}
	(*IL) = (int*)malloc((*nzL) * sizeof(int));
	(*JL) = (int*)malloc((*nzL) * sizeof(int));
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

int transposeCOO(int nz, int * I, int * J)
{
	int i;
	for (i = 0; i < nz; i++) {
		int temp;
		temp = I[i];
		I[i] = J[i];
		J[i] = temp;
	}
	return 0;
}

void printmatrixSparceCOO(int n, int nz, int * I, int * J, double * val)
{
	int i;
	for (i = 0; i < nz; i++)
		printf("%lf  %d  %d \n", I[i], J[i], val[i]);
}

int countZeroDiag(int * I, int * J, int nz, int N)
{
	int count = 0;
	for (int i = 0; i < nz; i++)
		if (I[i] == J[i])
			count++;
	return N - count;
}

void getZerosDiagNumbers(int * I, int * J, int nz, int N, int count, int * addDiag)
{
	int ind = 0;
	int* indexes = (int*)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++)
	{
		indexes[i] = 0;
	}
	for (int i = 0; i < nz; i++)
		if (I[i] == J[i])
			indexes[I[i]] = 1;
	for (int i = 0; i < N; i++)
		if (indexes[i] == 0) {
			addDiag[ind] = i;
			ind += 1;
		}
}

void fillDiag(int * I, int * J, double * val, int * Inew, int * Jnew, double * valNew, int N, int nz, int nzNew, int * addDiag) {
	for (int i = 0; i < nz; i++) {
		Inew[i] = I[i];
		Jnew[i] = J[i];
		valNew[i] = val[i];
	}
	for (int i = nz; i < nzNew; i++) {
		Inew[i] = addDiag[i - nz];
		Jnew[i] = addDiag[i - nz];
		valNew[i] = 1.0;
	}
	for (int i = nz; i < nzNew; i++) {
		int j = i;
		while (j > 0 && Inew[j] <= Inew[j - 1]) {
			int temp = Inew[j];
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

int CheckSolv(int n, double * x, double * xCheck, double * solution)
{
	float eps = 0.001f;
	int i, errors = 0;
	for (i = 0; i < n; i++) {
		if (((x[i] - xCheck[i]) > eps) || ((x[i] - xCheck[i]) < -eps))
			errors++;
		//printf("%f", x[i] - xCheck[i]);
		printf(" i = %d , x[%d] = %f, xCheck[%d] = %f, solutin[%d] = %f \n", i, i, x[i], i, xCheck[i], i, solution[i]);
	}
	return errors;
}

void gaussBackLow(int N, double * y, double * b, double * valCrsL, int * colL, int * indxL)
{
	int j, k, index_top, index_low;
	double sum;
	y[0] = b[0] / valCrsL[indxL[0]];
	for (k = 1; k < N; k++) {
		sum = 0;
		index_top = indxL[k + 1] - 1;
		index_low = indxL[k];
		for (j = index_top - 1; j >= index_low; j--) {
			sum += valCrsL[j] * y[colL[j]];
		}
		y[k] = (b[k] - sum) / valCrsL[indxL[k+1]-1];
	}
}

void gaussBackUp(int N, double * x, double * y, double * valCrsU, int * colU, int * indxU)
{
	int j, k, index_top, index_low;
	double sum;
	x[N - 1] = y[N - 1] / valCrsU[indxU[N] - 1];
	for (k = N - 2; k >= 0; k--) {
		sum = 0;
		index_top = indxU[k + 1];
		index_low = indxU[k];
		for (j = index_low + 1; j < index_top; j++) {
			sum += valCrsU[j] * x[colU[j]];
		}
		x[k] = (y[k] - sum) / valCrsU[index_low];
	}
}

void matrixMultVector(int N, double * x, double * xCheck, double * valCrs, int * col, int * indx)
{
	int i, j, index_top, index_low;
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
