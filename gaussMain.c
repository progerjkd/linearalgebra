#include "matrix.h"

int main() {

	matrix *A, *x, *b;
	int i, j, k, n = 3;
	double soma, alfa;

	A = newMatrix(n, n);
	b = newMatrix(n, 1);
	x = newMatrix(n, 1);

	setElement(A, 1, 1,   1);
	setElement(A, 1, 2,   1);
	setElement(A, 1, 3,  -2);
	setElement(A, 2, 1,  -2);
	setElement(A, 2, 2,   2);
	setElement(A, 2, 3,  -3);
	setElement(A, 3, 1,   3);
	setElement(A, 3, 2,  -1);
	setElement(A, 3, 3,   2);

	setElement(b, 1, 1,   0);
	setElement(b, 2, 1,   2);
	setElement(b, 3, 1,  12);

	printf("A =\n");
	printMatrix(A);

	printf("x =\n");
	printMatrix(x);

	printf("b =\n");
	printMatrix(b);

	printf("\nCalculando...\n\n");

	for (j=1; j<=n-1; j++){
		for (i=j+1; i<=n; i++){
			double a_jj, a_ij;
			getElement(A, j, j, &a_jj);
			getElement(A, i, j, &a_ij);
			alfa = (- a_jj / a_ij);
			
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

	printf("A =\n");
	printMatrix(A);

	printf("x =\n");
	printMatrix(x);

	printf("b =\n");
	printMatrix(b);
/*
	for(i=1; i<=n; i++){
		for(j=1; j<=n; j++){
			double tmp;
			getElement(A, i, j, &tmp);
			printf("%d,%d: %f\n", i, j, tmp); 
		}
	}
	printf("\n");

	setElement(L, 1, 1, 1);
	setElement(L, 2, 2, 1);
	setElement(L, 3, 3, 1);

	printf("L=\n");
	printMatrix(L);

	printf("U=\n");
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
		}
	}

	printf("A=\n");
	printMatrix(A);

	printf("L=\n");
	printMatrix(L);

	printf("U=\n");
	printMatrix(U);
	*/
}

