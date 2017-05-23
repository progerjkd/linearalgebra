#include "matrix.h"
#include <math.h>

int gaussJordan(matrix *A, matrix *x, matrix *b, int n) {

	int i, j, k;
	double soma, alfa;

	printf("\nEntrada - Eliminação de Gauss-Jordan\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nx =\n");
	printMatrix(x);

	printf("\nb =\n");
	printMatrix(b);


	for(k=1; k<=n; k++){
		// escolha do pivô
		double max = 0;
		int kk = 0;
		for(i=k; i<=n; i++){
			//printf("%d,%d\n",i,k);
			double a_ik;
			getElement(A, i, k, &a_ik);
			if (fabs(a_ik) > max){
				max = fabs(a_ik);
				kk = i;
			}
		}
		//printf("max = %f\n", max);
		//printf("kk = %d\n", kk);

		if(kk == 0){
			printf("O sistema não tem solução única.\n");
		}
		matrixSwapLines(A, n, k, kk);	// swap lines k and k from matrix A of order n
		matrixSwapLines(b, 1, k, kk);	// swap lines k and k from vector b of n lines

		/*
		printf("\nPasso %d - inicio:\n \nA =\n",k);
		printMatrix(A);

		printf("\nx =\n");
		printMatrix(x);

		printf("\nb =\n");
		printMatrix(b);
		*/


		double a_kk;
		getElement(A, k, k, &a_kk);

		for(j=k; j<=n; j++){
			double a_kj;
			getElement(A, k, j, &a_kj);
			setElement(A, k, j, a_kj / a_kk);
		}

		double b_k;
		getElement(b, k, 1, &b_k);
		setElement(b, k, 1, b_k / a_kk);

		for(i=1; i<=n; i++){
			if(i != k){
				double a_ik;
				getElement(A, i, k, &a_ik);
				for(j=k; j<=n; j++){
					double a_ij, a_kj;
					getElement(A, i, j, &a_ij);
					getElement(A, k, j, &a_kj);
					setElement(A, i, j, a_ij - a_kj * a_ik);
				}
				double b_i, b_k;
				getElement(b, i, 1, &b_i);
				getElement(b, k, 1, &b_k);
				setElement(b, i, 1, b_i - b_k * a_ik);
			}
		}

	}


	printf("\nSaída - Eliminação de Gauss-Jordan\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nx =\n");
	printMatrix(x);

	printf("\nb =\n");
	printMatrix(b);
}
