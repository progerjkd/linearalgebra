#include "matrix.h"

int main() {

	matrix *A, *b;
	int n = 3;

	A = newMatrix(n, n);
	b = newMatrix(n, 1);

	setElement(A, 1, 1,   2);
	setElement(A, 1, 2,   1);
	setElement(A, 1, 3,   3);
	setElement(A, 2, 1,   4);
	setElement(A, 2, 2,   2);
	setElement(A, 2, 3,  -5);
	setElement(A, 3, 1,   2);
	setElement(A, 3, 2,   1);
	setElement(A, 3, 3,   1);

	setElement(b, 1, 1,   4);
	setElement(b, 2, 1,  -5);
	setElement(b, 3, 1,   3);

	printf("\n-> Normas vetoriais\n");
	printf("b =\n");
	printMatrix(b);

	printf("\np-normas vetoriais:\n");
	for(int i=1; i<=10; i++){
		printf("||b||%d = %f\n", i, norma(b, i));
	}

	printf("\nNorma infinita: %f\n", normaInf(b));



	printf("\n-> Normas matriciais\n");
	printf("A =\n");
	printMatrix(A);

	printf("\np-normas matriciais:\n");
	for(int i=1; i<=10; i++){
		printf("||A||%d = %f\n", i, normaMatricial(A, i));
	}
	
	printf("\nNorma 1: %f\n", normaMatricial1(A));
	printf("Norma infinita: %f\n", normaMatricialInf(A));
}
