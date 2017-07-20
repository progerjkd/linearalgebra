#include "matrix.h"

int main() {

	matrix *A, *Q, *R;
	int m = 2, n = 2;

	A = newMatrix(m, n);
	Q = newMatrix(m, m);
	R = newMatrix(m, n);

	/*
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
	*/

	setElement(A, 1, 1,   57.3970);
	setElement(A, 1, 2,   -0.0047);
	setElement(A, 2, 1,   -0.0047);
	setElement(A, 2, 2,  0);

	_qr QR;
	QR = qr(A);


	printMatrix2(QR.Q, "QR.Q");
	printMatrix2(QR.R, "QR.R");
/*
	matrix *prod = newMatrix(A->rows, A->cols);
	product(Q, R, prod);
	printf("\nQ * R =\n");
	printMatrix(prod);
*/
	// falta fazer o restante do algoritmo para o cálculo de autovetores
	// aqui só está sendo computada a primeira decomposição QR
}
