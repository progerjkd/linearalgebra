#include "matrix.h"

int main() {

	matrix *A, *HH, *Q, *R;
	int n = 4;

	A  = newMatrix(n, n);
	HH = newMatrix(n, n);
	Q  = newMatrix(n, n);
	R  = newMatrix(n, n);

	setElement(A, 1, 1,   4);
	setElement(A, 1, 2,   2);
	setElement(A, 1, 3,   2);
	setElement(A, 1, 4,   1);
	setElement(A, 2, 1,   2);
	setElement(A, 2, 2,  -3);
	setElement(A, 2, 3,   1);
	setElement(A, 2, 4,   1);
	setElement(A, 3, 1,   2);
	setElement(A, 3, 2,   1);
	setElement(A, 3, 3,   3);
	setElement(A, 3, 4,   1);
	setElement(A, 4, 1,   1);
	setElement(A, 4, 2,   1);
	setElement(A, 4, 3,   1);
	setElement(A, 4, 4,   2);


	printf("\nA =\n");
	printMatrix(A);
	
	isSymmetric(A) ? printf("\nMatrix A is symmetric\n") : printf("\nMatrix A is Asymmetric\n");

	printf("\nHouseholder transformation: (results in a tridiagonal matrix)\n");
	householder(A, HH);

	printf("\nHH =\n");
	printMatrix(HH);

	printf("\nQR decomposition of HH matrix...\n");

	for (int i=1; i<=50; i++){
		qr(HH, Q, R);
		product(R, Q, HH);
	}

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

}
