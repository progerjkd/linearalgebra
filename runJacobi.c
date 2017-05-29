#include "matrix.h"

int main() {

	matrix *A, *e, *V;
	int n = 3, sweeps;

	A = newMatrix(n, n);
	V = newMatrix(n, n);
	e = newMatrix(n, 1);


	setElement(A, 1, 1,  2.518);
	setElement(A, 1, 2,  0.459);
	setElement(A, 1, 3, -0.713);
	setElement(A, 2, 1,  0.459);
	setElement(A, 2, 2,  1.725);
	setElement(A, 2, 3,  0.282);
	setElement(A, 3, 1, -0.713);
	setElement(A, 3, 2,  0.282);
	setElement(A, 3, 3,  1.757);


	matrix *B;
	B = copyMatrix(A);
	sweeps = jacobi(A, e, V);


	printf("\n\n\nsweeps = %d\n", sweeps);

	printf("\nV'*A*V (should be diagonal):\n");
	matrix *Vt, *tmp;
	Vt = newMatrix(V->rows, V->cols);
	transpose(V, Vt);
	tmp = newMatrix(Vt->rows, B->cols);
	product(Vt, B, tmp);
	product(tmp, V, A);
	printMatrix(A);

	/*	
	printf("\nEigenvalues (should be the diagonal elements):\n");
	printMatrix(e);
	*/

	printf("\n");
	matrix *eigenvectors[V->rows+1];

	for(int i=1; i<=V->rows; i++){
		double eigenvalue;

		eigenvectors[i] = matrixColToVector(V, i);
		printf("Autovetor[%d]:\n", i);
		printMatrix(eigenvectors[i]);

		getElement(e, i, 1, &eigenvalue);
		printf("Autovalor = %f\n\n", eigenvalue);

	}

}

