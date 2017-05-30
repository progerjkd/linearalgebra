#include "matrix.h"

int main() {

	matrix *A, *Q, *R;
	int m = 5, n = 3;

	A = newMatrix(m, n);
	Q = newMatrix(m, m);
	R = newMatrix(m, n);

	setElement(A, 1, 1,  12);
	setElement(A, 1, 2, -51);
	setElement(A, 1, 3,   4);
	setElement(A, 2, 1,   6);
	setElement(A, 2, 2, 167);
	setElement(A, 2, 3, -68);
	setElement(A, 3, 1,  -4);
	setElement(A, 3, 2,  24);
	setElement(A, 3, 3, -41);
	setElement(A, 4, 1,  -1);
	setElement(A, 4, 2,   1);
	setElement(A, 4, 3,   0);
	setElement(A, 5, 1,   2);
	setElement(A, 5, 2,   0);
	setElement(A, 5, 3,   3);

	householder(A, R, Q);


	matrix *prod = newMatrix(A->rows, A->cols);
	product(Q, R, prod);
	printf("\nQ * R =\n");
	printMatrix(prod);

}
