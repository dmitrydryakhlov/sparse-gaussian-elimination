#include <stdio.h>

#include "../include/checker.h"
#include "../include/reader.h"
#include "../include/solver.h"
#include "../include/utils.h"

int main(int argc, char* argv[]) {
	int m = 100;
	int n = 100;
	int nzz = 1000;
	int Col[100];
	int Row[100];
	float Val[100];
	float B[100];
	float X[100];

	//ReadMatrix();
	SolverMatrix(n, m, nzz, Col, Row, Val, B, X);
	CheckSolv(n, m, Col, Row, Val, X, B);
	printf("All: success\n");

	int error = 0;
	int *I, *IU, *J, *JU;
	double *val, *valU;
	int *row, *col;
	double *valCRS;

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

	free(I);
	free(J);
	free(val);
	free(row);
	free(col);
	free(valCRS);

	return 0;
}