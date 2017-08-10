#include "matrix.h"

void decomposicaoLU(matrix *A, int n, matrix *L, matrix *U){

	int i, j, k;
	double soma;

	printf("\nEntrada - Decomposição LU\n");

	printf("\nA =\n");
	printMatrix(A);

	setElement(L, 1, 1, 1);
	setElement(L, 2, 2, 1);
	setElement(L, 3, 3, 1);

	printf("\nL =\n");
	printMatrix(L);

	printf("\nU =\n");
	printMatrix(U);

	for(j=1; j<=n; j++){
		for(i=1; i<=j; i++){
			soma = 0;
			for(k=1; k<=i-1; k++){
				double a_ik, a_kj;
				getElement(L, i, k, &a_ik);
				getElement(U, k, j, &a_kj);
				soma += a_ik * a_kj;
			}
			double a_ij;
			getElement(A, i, j, &a_ij);
			setElement(U, i, j, a_ij - soma);
		} // U

		for(i=j+1; i<=n; i++){
			soma = 0;
			for(k=1; k<=j-1; k++){
				double a_ik, a_kj;
				getElement(L, i, k, &a_ik);
				getElement(U, k, j, &a_kj);
				soma += a_ik * a_kj;
			}
			double a_ij, u_jj;
			getElement(A, i, j, &a_ij);
			getElement(U, j, j, &u_jj);
			setElement(L, i, j, (a_ij - soma)/u_jj );
		} // L
	}

	printf("\nSaída - Decomposição LU\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nL =\n");
	printMatrix(L);

	printf("\nU =\n");
	printMatrix(U);
}
