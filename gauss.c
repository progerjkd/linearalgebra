#include "matrix.h"

void gauss(matrix *A, matrix *x, matrix *b, int n) {

	int i, j, k;
	double soma, alfa;

	printf("\nEntrada - Eliminação de Gauss\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nx =\n");
	printMatrix(x);

	printf("\nb =\n");
	printMatrix(b);


	for (j=1; j<=n-1; j++){
		for (i=j+1; i<=n; i++){
			double a_jj, a_ij;
			getElement(A, j, j, &a_jj);
			getElement(A, i, j, &a_ij);
			if (a_ij != 0 )
				alfa = (- a_jj / a_ij);
			else
				alfa = 0;

//			printf("a_jj: %f, a_ij: %f, alfa: %f\n", a_jj, a_ij, alfa);
			
			for (k=j+1; k<=n; k++){
				double a_ik, a_jk;
				getElement(A, i, k, &a_ik);
				getElement(A, j, k, &a_jk);
				setElement(A, i, k, alfa * a_ik + a_jk);
			}
			
			double b_i, b_j;
			getElement(b, i, 1, &b_i);
			getElement(b, j, 1, &b_j);

			setElement(b, i, 1, alfa * b_i + b_j);

		}

	}

	for (i=n; i>=1; i--){
		soma = 0;
		for (k=i+1; k<=n; k++){
			double a_ik, x_k;
			getElement(A, i, k, &a_ik);
			getElement(x, k, 1, &x_k);

			soma += a_ik * x_k;
		}

		double b_i, a_ii;
		getElement(b, i, 1, &b_i);
		getElement(A, i, i, &a_ii);

		setElement(x, i, 1, (b_i - soma) / a_ii);
	}

	printf("\nSaída - Eliminação de Gauss\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nx =\n");
	printMatrix(x);

	printf("\nb =\n");
	printMatrix(b);
}
