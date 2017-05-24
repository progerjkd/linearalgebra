#include "matrix.h"
#include <math.h>

int decomposicaoCholesky(matrix *A, matrix *S,  matrix *St, int n) {

	int i, j, k;

	printf("\nEntrada - Decomposição de Cholesky\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nS =\n");
	printMatrix(S);

	printf("\nSt =\n");
	printMatrix(St);

	if(!isSymmetric(A)){
		printf("A matriz A não é simétrica. Abortando a decomposição de Cholesky...\n");
		return -1;
	}

	for (j=1; j<=n; j++){
		double A_jj;
		getElement(A, j, j, &A_jj);
		setElement(S, j, j, A_jj);

		for (i=1; i<=j-1; i++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			setElement(S, i, j, A_ij);

			for(k=1; k<=i-1; k++){
				double S_ij, S_ki, S_kj;
				getElement(S, i, j, &S_ij);
				getElement(S, k, i, &S_ki);
				getElement(S, k, j, &S_kj);

				setElement(S, i, j, S_ij - S_ki * S_kj);
			}
			double S_ij, S_ii, S_jj;
			getElement(S, i, j, &S_ij);
			getElement(S, i, i, &S_ii);
			setElement(S, i, j, S_ij / S_ii);

			getElement(S, j, j, &S_jj);
			getElement(S, i, j, &S_ij);
			setElement(S, j, j, S_jj - S_ij * S_ij);
		}

		double S_jj;
		getElement(S, j, j, &S_jj);
		setElement(S, j, j, sqrt(S_jj));
	}

	transpose(S, St);

	printf("\nSaída - Decomposição de Cholesky\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nS =\n");
	printMatrix(S);

	printf("\nSt =\n");
	printMatrix(St);
}