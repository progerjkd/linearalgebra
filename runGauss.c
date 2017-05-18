#include "matrix.h"

int main() {

	matrix *A, *x, *b, *L, *U, *y;
	int n = 3;

	A = newMatrix(n, n);
	x = newMatrix(n, 1);
	b = newMatrix(n, 1);
	L = newMatrix(n, n);
	U = newMatrix(n, n);
	y = newMatrix(n, 1);
/*
	setElement(A, 1, 1,   5);
	setElement(A, 1, 2,   2);
	setElement(A, 1, 3,   3);
	setElement(A, 2, 1,  -2);
	setElement(A, 2, 2,   1);
	setElement(A, 2, 3,   2);
	setElement(A, 3, 1,   1);
	setElement(A, 3, 2,   4);
	setElement(A, 3, 3,  -1);

	setElement(b, 1, 1,  10);
	setElement(b, 2, 1,   4);
	setElement(b, 3, 1,   3);

	gauss(A, x, b, n);
*/

	setElement(A, 1, 1,   1);
	setElement(A, 1, 2,  -1);
	setElement(A, 1, 3,   2);
	setElement(A, 2, 1,   2);
	setElement(A, 2, 2,   1);
	setElement(A, 2, 3,  -1);
	setElement(A, 3, 1,  -2);
	setElement(A, 3, 2,  -5);
	setElement(A, 3, 3,   3);

	setElement(b, 1, 1,   2);
	setElement(b, 2, 1,   1);
	setElement(b, 3, 1,   3);

//	decomposicaoLU(A, n, L, U);


	gauss(A, x, b, n);
//	retroSubstituicao(U, x, y, n);

}
