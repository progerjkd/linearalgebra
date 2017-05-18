#include "matrix.h"

int retroSubstituicao(matrix *A, matrix *x, matrix *b, int n) {

	// Implementar verificação de diagonal vazia
	// Caso ocorra algum elemento dadiagonal vazio o método torna-se inválido
	
	if(diagonalContainsZero(A)){
		printf("Diagonal da matriz contém valores nulos, retrosubstituição abortada...\n");
		printf("A =\n");
		printMatrix(A);
		return -1;
	}
		
	int i, j;

	printf("\nEntrada - Retrosubstuição\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nx =\n");
	printMatrix(x);

	printf("\nb =\n");
	printMatrix(b);



	double b_n, a_nn;
	getElement(b, n, 1, &b_n);
	getElement(A, n, n, &a_nn);
	setElement(x, n, 1, b_n / a_nn);
			
	for (i=n-1; i>=1; i--){
		double a_ii, b_i, soma = 0;
		for (j=i+1; j<=n; j++){
			double a_ij, x_j;
			getElement(A, i, j, &a_ij);
			getElement(x, j, 1, &x_j);
			soma += a_ij * x_j;
		}
		getElement(A, i, i, &a_ii);
		getElement(b, i, 1, &b_i);
		setElement(x, i, 1, (b_i - soma) / a_ii);
	}

	printf("\nSaída - Retrosubstuição\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nx =\n");
	printMatrix(x);

	printf("\nb =\n");
	printMatrix(b);
}
