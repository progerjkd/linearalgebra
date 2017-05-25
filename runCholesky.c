#include "matrix.h"

int main() {

	matrix *A, *x, *b, *S, *St, *y;
	int n = 3;

	A = newMatrix(n, n);
	x = newMatrix(n, 1);
	b = newMatrix(n, 1);
	S = newMatrix(n, n);
	St = newMatrix(n, n);
	y = newMatrix(n, 1);

	setElement(A, 1, 1,   4);
	setElement(A, 1, 2,  -2);
	setElement(A, 1, 3,   2);
	setElement(A, 2, 1,  -2);
	setElement(A, 2, 2,  10);
	setElement(A, 2, 3,  -7);
	setElement(A, 3, 1,   2);
	setElement(A, 3, 2,  -7);
	setElement(A, 3, 3,  30);

	setElement(b, 1, 1,   6);
	setElement(b, 2, 1,  -3);
	setElement(b, 3, 1,  78);

	decomposicaoCholesky(A, S, St, n);

/*
	printf("\ngauss(L, y, b, n):\n");
	gauss(L, y, b, n);
	printf("\nretroSubstituicao(U, x, y, n):\n");
	retroSubstituicao(U, x, y, n);
	
*/
	printf("\ngaussP(St, y, b, n):\n");
	gaussP(St, y, b, n);
	printf("\nretroSubstituicao(St, y, b, n):\n");
	retroSubstituicao(St, y, b, n);
	printf("\nretroSubstituicao(S, x, y, n):\n");
	retroSubstituicao(S, x, y, n);
}
