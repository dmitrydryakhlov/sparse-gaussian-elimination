#include "utils.h"

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

	for (i = 0; i < (*nz); i++)
	{
		fscanf(f, "%d %d %lg\n", &(*I)[i], &(*J)[i], &(*val)[i]);
		(*I)[i]--;  /* adjust from 1-based to 0-based */
		(*J)[i]--;
	}
	if (f != stdin) fclose(f);
	return 0;
}

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

long long COOtoCRS(long long n, long long nz, long long *I, long long *J, double *valCOO, long long **indx, long long **col, double **valCrs)
{
	long long i;
	long long *places;

	(*indx) = (long long*)malloc((n + 1) * sizeof(long long));
	(*col) = (long long*)malloc((nz) * sizeof(long long));
	(*valCrs) = (double*)malloc((nz) * sizeof(double));
	places = (long long*)malloc((n+1) * sizeof(long long));

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

long long transposeCOO(long long nz, long long * I, long long * J, double * val, long long nzT, long long * IT, long long * JT, double * valT)
{
	long long i;
	for (i = 0; i < nz; i++) {
		IT[i] = J[i];
		JT[i] = I[i];
		valT[i] = val[i];
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

void printmatrixSparceCOO(long long n, long long nz, long long * I, long long * J, double * val)
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

long long CheckSolv(long long n, double * x, double * xCheck)
{
	double eps = 0.001;
	long long i, errors = 0;
	for (i = 0; i < n; i++) {
		if (((x[i] - xCheck[i]) > eps) || ((x[i] - xCheck[i]) < -eps))
			errors++;
		//prlong longf("%f", x[i] - xCheck[i]);
		//prlong longf(" i = %d , x[%d] = %f, xCheck[%d] = %f, \n", i, i, x[i], i, xCheck[i], i);
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

void gaussBackLow(long long N, double * y, double * b, double * valCrsL, long long * colL, long long * indxL)
{
	long long j, k, index_top, index_low;
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

void gaussBackUp(long long N, double * x, double * y, double * valCrsU, long long * colU, long long * indxU)
{
	long long j, k, index_top, index_low;
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
