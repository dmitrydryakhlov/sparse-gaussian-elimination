#include "utils.h"

int main(int argc, char* argv[]) {
	int error = 0;
	int *I, *IU, *J, *JU, *indxU, *colU, *IL, *JL, *colL, *indxL;
	double *val, *valU, *valCrsU, *valCrsL, *valCooU, *valL, valCooL;
	double *y, *b, *bCheck, *x, *yCheck;
	int N, M, nz, nzU, nzL;
	int i, j;
	if (readMTX(argv[1], &I, &J, &val, &M, &N, &nz) != 0)
		exit(0);

	int zeroDiagCount = countZeroDiag(I, J, nz, N);
	int *addDiag = (int*)malloc(zeroDiagCount * sizeof(int));
	int *Inew = (int*)malloc((nz + zeroDiagCount) * sizeof(int));
	int *Jnew = (int*)malloc((nz + zeroDiagCount) * sizeof(int));
	double *valNew = (double*)malloc((nz + zeroDiagCount) * sizeof(double));
	getZerosDiagNumbers(I, J, nz, N, zeroDiagCount, addDiag);
	fillDiag(I, J, val, Inew, Jnew, valNew, N, nz, nz + zeroDiagCount, addDiag);
	free(I); I = Inew;
	free(J); J = Jnew;
	free(val); val = valNew;
	nz = nz + zeroDiagCount;

	y = (double*)malloc((N) * sizeof(double));
	b = (double*)malloc((N) * sizeof(double));
	bCheck = (double*)malloc((N) * sizeof(double));
	x = (double*)malloc((N) * sizeof(double));
	yCheck = (double*)malloc((N) * sizeof(double));
	
	for (i = 0; i < N; i++) {
		b[i] = (double)(i + 1);
		y[i] = (double)(i + 1);
		bCheck[i] = 0;
		x[i] = 0;
		yCheck[i] = 0;
	}

	transposeCOO(nz, I, J);
	printf("Finished transposing coordinate matrix\n");

	cutLowerTriangleCOO(nz, I, J, val, &nzU, &IU, &JU, &valU);
	printf("Finished cutting coordinate matrix\n");

	cutUpperTriangleCOO(nz, I, J, val, &nzL, &IL, &JL, &valL);
	printf("Finished cutting coordinate matrix\n");

	COOtoCRS(N, nzU, IU, JU, valU, &indxU, &colU, &valCrsU);
	printf("Finished converting coordinate matrix to CRS matrix\n");

	COOtoCRS(N, nzL, IL, JL, valL, &indxL, &colL, &valCrsL);
	printf("Finished converting coordinate matrix to CRS matrix\n");

	//saveBinCRS(argv[2], N, row, col, valCRS);
	//printf("Finished saving CRS matrix\n");

	gaussBackLow(N, y, b, valCrsL, colL, indxL);
	gaussBackUp (N, x, y, valCrsU, colU, indxU);
	matrixMultVector(N, x, yCheck, valCrsU, colU, indxU);
	matrixMultVector(N, y, bCheck, valCrsL, colL, indxL);

	int errors = CheckSolv(N, b, bCheck, y);
	if (errors == 0) { printf("\n All is Ok\n"); }
	else { printf("\n Errors: %d/%d \n", errors, N); }

	free(IU);
	free(JU);
	free(valU);
	free(y);
	free(b);
	free(bCheck);
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
