#include "matrix.h"
#include <math.h>

int main() {

	matrix *A, *Hess, *Q, *R;
	int n = 4;

	A  = newMatrix(n, n);
	Hess = newMatrix(n, n);
	Q  = newMatrix(n, n);
	R  = newMatrix(n, n);

	setElement(A, 1, 1,   4);
	setElement(A, 1, 2,   6);
	setElement(A, 1, 3,   7);
	setElement(A, 1, 4,   2);
	setElement(A, 2, 1,   1);
	setElement(A, 2, 2,   2);
	setElement(A, 2, 3,   0);
	setElement(A, 2, 4,   1);
	setElement(A, 3, 1,  -2);
	setElement(A, 3, 2,   0);
	setElement(A, 3, 3,   3);
	setElement(A, 3, 4,  -2);
	setElement(A, 4, 1,   2);
	setElement(A, 4, 2,   1);
	setElement(A, 4, 3,  -2);
	setElement(A, 4, 4,  -1);


	printf("\nA =\n");
	printMatrix(A);

	isSymmetric(A) ? printf("\nMatrix A is symmetric\n") : printf("\nMatrix A is Asymmetric\n");

	printf("\nHouseholder transformation: (results in a Hessenberg Matrix)\n");

	householder(A, Hess);

	printf("\nHess =\n");
	printMatrix(Hess);

	printf("\nQR decomposition of Hess matrix...\n");

	for (int i=1; i<=50; i++){
		qr(Hess, Q, R);
		product(R, Q, Hess);
	}

//	printf("\nMatriz com dentes (ou não) abaixo da diagonal:\n");
	printf("\nHess =\n");
	printMatrix(Hess);

	double epsilon = 0.001;
	n = Hess->rows;
	int i = 1;

	printf("\nEigenvalues of matrix A:\n\n");
	while(i <= n){
		double Hesstmp;
		getElement(Hess, i+1, i, &Hesstmp);
		if((i<n) && (fabs(Hesstmp) > epsilon)){
			// seleciona a matriz 2x2 formada com os dentes e calcula os autovalores
			matrix *M2 = newMatrix(2, 2);
			subMatrix(Hess, i, i+1, i, i+1, M2);

			double M2_11, M2_22;
			getElement(M2, 1, 1, &M2_11);
			getElement(M2, 2, 2, &M2_22);


			_eqSegGrau x;
			x = eqSegGrau(1, -(M2_11 + M2_22), laplace(M2));
			printf("%10.4f %s%fi\n", creal(x.x1), (cimag(x.x1) >= 0) ? "+" : "", cimag(x.x1));
			printf("%10.4f %s%fi\n", creal(x.x2), (cimag(x.x2) >= 0) ? "+" : "", cimag(x.x2));
			i++;
		} else {
			double Hess_ii;
			getElement(Hess, i, i, &Hess_ii);
			printf("%10.4f\n", Hess_ii);
		}
		i++;
	}

/*	printf("\nR =\n");
	printMatrix(R);

	printf("\nMatrix diagonal de autovalores:\n");
	printf("Q =\n");
	printMatrix(Q);
*/
}
