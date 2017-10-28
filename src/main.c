#include <stdio.h>

#include "../include/checker.h"
#include "../include/reader.h"
#include "../include/solver.h"
#include "../include/utils.h"

int main(int argc, char* argv[]) {
	int error = 0;
	int *I, *IU, *J, *JU;
	double *val, *valU;
	int *row, *col;
	double *valCRS;
	double *y, *b, *b1, *x;

	int M, N, nz, nzU;
	if (argc < 3) {
		printHelp();
		exit(0);
	}
	if (readMTX(argv[1], &I, &J, &val, &M, &N, &nz) != 0)
		exit(0);

	transposeCOO(nz, I, J);
	printf("Finished transposing coordinate matrix\n");

	cutLowerTriangleCOO(nz, I, J, val, &nzU, &IU, &JU, &valU);
	printf("Finished cutting coordinate matrix\n");

	COOtoCRS(N, nzU, IU, JU, val, &row, &col, &valCRS);
	printf("Finished converting coordinate matrix to CRS matrix\n");

	saveBinCRS(argv[2], N, row, col, valCRS);
	printf("Finished saving CRS matrix\n");
	int i, j = 0;
	printf("\n Value \n");
	for (i = 0; i < nzU; i++) {
		printf("%f   ", valCRS[i]);
	}
	printf("\n col \n");
	for (i = 0; i < nzU; i++) {
		printf("%d   ", col[i]);
	}
	printf("\n row\n");
	for (i = 0; i < N + 1; i++) {
		printf("%d   ", row[i]);
	}
	printf("%d", N);

	b = (double*)malloc((N) * sizeof(double));
	y = (double*)malloc((N) * sizeof(double));

	printf("\n b \n");
	for (i = 0; i < N; i++) {
		b[i] = i;
		printf("%f   ", b[i]);
	}

	int index1, index2;
	for (i = N - 1; i >= 0; i--) {
		double sum = 0;
		index1 = row[i];
		index2 = row[i + 1];
		for (j = index1+1; j < index2; j++) {
			sum += val[j] * y[j];
		}
		y[i] = (b[i] - sum) / val[row[i]];
	}
	printf("\n y: \n");
	for (i = 0; i < N; i++) {
		printf("%f   ", y[i]);
	}
	free(y);
	free(b);
	free(I);
	free(J);
	free(val);
	free(row);
	free(col);
	free(valCRS);
	return 0;
}