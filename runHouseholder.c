#include "matrix.h"

int main() {

	matrix *A, *HH;
	int n = 4;

	A = newMatrix(n, n);
	HH = newMatrix(n, n);

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


	householder(A, HH);
}
