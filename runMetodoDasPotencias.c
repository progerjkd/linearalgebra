#include "matrix.h"

int main() {

	matrix *A, *y;
	int n = 3;

	A = newMatrix(n, n);
	y = newMatrix(n, 1);


	setElement(A, 1, 1,   8);
	setElement(A, 1, 2,   1);
	setElement(A, 1, 3,   2);
	setElement(A, 2, 1,  -1);
	setElement(A, 2, 2,   5);
	setElement(A, 2, 3,   1);
	setElement(A, 3, 1,   0);
	setElement(A, 3, 2,   1);
	setElement(A, 3, 3,  90);

	setElement(y, 1, 1,   1);
	setElement(y, 2, 1,   1);
	setElement(y, 3, 1,   1);

	int tol=100;
	double lambda;
	matrix *z;

	metodoDasPotencias(A, y, (double)0.01, 6, &z, &lambda); 

	printf("\nMaior autovalor: %f\n", lambda);
	printf("Autovetor associado:\n");
	printMatrix(z);
}
