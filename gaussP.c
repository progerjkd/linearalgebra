#include "matrix.h"
#include <math.h>

void gaussP(matrix *A, matrix *x, matrix *b, int n) {

	int i, j, k;
	double soma, alfa;

	printf("\nEntrada - Eliminação de Gauss com pivoteamento\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nx =\n");
	printMatrix(x);

	printf("\nb =\n");
	printMatrix(b);


	for(k=1; k<=n-1; k++){
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



		for(i=k+1; i<=n; i++){
			double m, a_ik, a_kk;
			getElement(A, i, k, &a_ik);
			getElement(A, k, k, &a_kk);
			m = a_ik/a_kk;
			//printf("m = a_%d%d/a_%d%d => %f = %f/%f\n", i, k, k, k, m, a_ik, a_kk);

			for(j=k+1; j<=n+1; j++){
				double a_ij, a_kj;
				getElement(A, i, j, &a_ij);
				getElement(A, k, j, &a_kj);

				setElement(A, i, j, a_ij - m * a_kj);
			}

			double b_i, b_k;
			getElement(b, i, 1, &b_i);
			getElement(b, k, 1, &b_k);

			setElement(b, i, 1, b_i - m * b_k);
		}

		/*
		printf("\nPasso %d - fim:\n \nA =\n",k);
		printMatrix(A);

		printf("\nx =\n");
		printMatrix(x);

		printf("\nb =\n");
		printMatrix(b);
		*/

	}


	printf("\nSaída - Eliminação de Gauss com pivoteamento\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nx =\n");
	printMatrix(x);

	printf("\nb =\n");
	printMatrix(b);
}
