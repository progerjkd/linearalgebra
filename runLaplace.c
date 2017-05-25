#include "matrix.h"

int main() {

	matrix *A, *x, *b;
	int n = 4;

	A = newMatrix(n, n);
	x = newMatrix(n, 1);
	b = newMatrix(n, 1);

	setElement(A, 1, 1,   2);
	setElement(A, 1, 2,   1);
	setElement(A, 1, 3,   6);
	setElement(A, 1, 4,   0);
	setElement(A, 2, 1,   5);
	setElement(A, 2, 2,   0);
	setElement(A, 2, 3,   1);
	setElement(A, 2, 4,   3);
	setElement(A, 3, 1,  -2);
	setElement(A, 3, 2,   5);
	setElement(A, 3, 3,   1);
	setElement(A, 3, 4,   8);
	setElement(A, 4, 1,  11);
	setElement(A, 4, 2,   4);
	setElement(A, 4, 3,  -2);
	setElement(A, 4, 4,  -7);

	setElement(b, 1, 1,  -4);
	setElement(b, 2, 1,   0);
	setElement(b, 3, 1,   7);
	setElement(b, 4, 1, -10);

	printMatrix(A);
	double det = laplace(A);
	printf("Determinante: %f\n", det);
}
