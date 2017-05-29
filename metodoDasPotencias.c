#include "matrix.h"
#define LAMBDA "\u019b"

int metodoDasPotencias(matrix *A, matrix *z0, double tol, int maxit, matrix **_z, double *_lambda){

 	printf("\nEntrada - Método das potências\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nz[0] =\n");
	printMatrix(z0);

	printf("\ntol = %f\nmaxit = %d\n\n", tol, maxit);


	matrix *z[maxit];

	z[0] = z0;

	printf("Saída - Método das potências\n\n");
	// z[0]= z[0]/ ||z[0]||
//	printf("\n||z0|| = %f\n", norma(z[0], 2));
	matrix *tmp;
	tmp = newMatrix(z[0]->rows, z[0]->cols);
	divisionByScalar(z[0], norma(z[0], 2), tmp);
	z[0] = tmp;
	printf("k= 0\n");
	printf("z[0] =\n");
	printMatrix(z[0]);

	// lambda[0] = z[0]' * A * z[0]
	double lambda[maxit];
//	lambda = malloc(maxit * sizeof(double));
	tmp = newMatrix(A->rows, z[0]->cols);
	product(A, z[0], tmp);
	dotProduct(tmp, z[0], &lambda[0]);
	printf("%s[0]= %f\n", LAMBDA, lambda[0]);

	int k;	
	for(k=1; k<=maxit; k++){

		printf("\n\nk= %d\n", k);

		// q = A * z[k-1]
		matrix *q;
		q = newMatrix(A->rows, z[k-1]->cols);
		product(A, z[k-1], q);
//		printf("q =\n");
//		printMatrix(q);

		// z[k] = q / ||q||
//		printf("||q|| = %f\n", norma(q, 2));
		z[k] = newMatrix(q->rows, q->cols);
		divisionByScalar(q, norma(q, 2), z[k]);
		printf("z[%d] =\n", k);
		printMatrix(z[k]);

		// lambda[k] = z[k]' * A * z[k]
		tmp = newMatrix(A->rows, z[k]->cols);
		product(A, z[k], tmp);
		dotProduct(tmp, z[k], &lambda[k]);
		printf("%s[%d]= %f\n", LAMBDA, k, lambda[k]);

		if(lambda[k] - lambda[k-1] < tol) {
//			printf("\n\n%s[%d] - %s[%d] < %f: %f - %f < %f\n", LAMBDA, k, LAMBDA, k-1, tol, lambda[k], lambda[k-1], tol);
//			printf("%f < %f\nbreaking...\n", lambda[k] - lambda[k-1], tol);
			break;
		}
		// evita o autoincremento de k ao sair do loop
		if(k == maxit)
			break;
//		printf("lambda[%d] - lambda[%d] > %f = %f > %f\n", k, k-1, tol, lambda[k] - lambda[k-1], tol);
	}

	*_z = z[k];
	*_lambda = lambda[k];

}
