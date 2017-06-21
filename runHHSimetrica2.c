#include "matrix.h"

#define VERBOSE

int main(int argc, char *argv[]) {

	if(argc <= 1){
		printf("usage: %s inputFile\n", argv[0]);
		exit(1);
	}

	matrix *A, *HH, *Q, *R;

	loadMatrix(&A, argv[1]);
	HH = newMatrix(A->rows, A->cols);
	Q  = newMatrix(A->rows, A->rows);
	R  = newMatrix(A->rows, A->cols);

#ifdef VERBOSE
	printf("\nA =\n");
	printMatrix(A);
#endif	
	isSymmetric(A) ? printf("\nMatrix A is symmetric\n") : printf("\nMatrix A is Asymmetric\n");

	printf("\nHouseholder transformation: (results in a tridiagonal matrix)\n");
	householder(A, HH);

#ifdef VERBOSE
	printf("\nHH =\n");
	printMatrix(HH);
#endif

	printf("\nQR decomposition of HH matrix...\n");

	for (int i=1; i<=50; i++){
		qr(HH, Q, R);
		product(R, Q, HH);
	}

#ifdef VERBOSE
	printf("\nHH =\n");
	printMatrix(HH);

	printf("\nQ =\n");
	printMatrix(Q);

	printf("\nEigenvalues diagonal matrix:\n");
	printf("\nR =\n");
	printMatrix(R);

	matrix *v = newMatrix(R->rows, 1);
	diagonalToVector(R, v);
	printf("\nEigenvalues of matrix A:\n\n");
	printMatrix(v);
#endif

}
