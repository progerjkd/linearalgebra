#include "matrix.h"
#include <math.h>

int metodoDaPotencia2(matrix *A, matrix *x, double tol, int maxit){


	int p, k = 1;

	printf("\nEntrada - Método das potências\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nx =\n");
	printMatrix(x);

	printf("\ntol = %f\nmaxit = %d\n\n", tol, maxit);

	printf("\nSaída - Método das potências\n");

	for(int i=1; i<=x->rows; i++){
		double x_i;
		getElement(x, i, 1, &x_i);
		x_i = fabs(x_i);

		if(x_i == normaInf(x)){
			p = i;
			break;
		}
	}

	matrix *tmp;
	tmp = newMatrix(x->rows, x->cols);
	double x_p;
	getElement(x, p, 1, &x_p);
	x_p = fabs(x_p);
	divisionByScalar(A, x_p, tmp);

	while(k <= maxit){
		matrix *y;
		y = newMatrix(A->rows, x->cols);
		product(A, x, y);

		double u;
		getElement(y, p, 1, &u);

		for(int i=1; i<=y->rows; i++){
			double y_i;
			getElement(y, i, 1, &y_i);
			y_i = fabs(y_i);

			if(y_i == normaInf(y)){
				p = i;
				break;
			}
		}


		double y_p;
		getElement(y, p, 1, &y_p);
		y_p = fabs(y_p);

		if(y_p == 0){
			printf("\n\nk = %d:\n", k);
			printf("y_p = 0\n");
			printf("autovalor = 0\n");
			printf("autovetor =\n");
			printMatrix(x);
			//return 1;
		}
		
		matrix *tmp2, *tmp3, *tmp4, *tmp5;
		tmp2 = newMatrix(y->rows, y->cols);
		tmp3 = newMatrix(y->rows, y->cols);
		tmp4 = newMatrix(y->rows, y->cols);
		tmp5 = newMatrix(y->rows, y->cols);
		divisionByScalar(y, y_p, tmp2);
/*printf("x =\n");
printMatrix(x);
printf("\ntmp2 =\n");
printMatrix(tmp2);
printf("\ntmp3 =\n");
printMatrix(tmp3);
*/
		subtraction(x, tmp2, tmp3);

		if (normaInf(tmp3) < tol){
			printf("\n\nk = %d:\n", k);
			printf("normaInf(tmp3) < tol: ");
			printf("%f < %f\n", normaInf(tmp3), tol);
			printf("autovalor = %f\n", u);
			printf("autovetor =\n");
			divisionByScalar(y, y_p, tmp4);
			printMatrix(tmp4);
			//return 1;

		}


		divisionByScalar(y, y_p, tmp5);
		x = copyMatrix(tmp5);
		k++;
	}
	printf("\nNúmero máximo de iterações atingido\n");


/*
	printf("\nA =\n");
	printMatrix(A);

	printf("\nx =\n");
	printMatrix(x);
*/
}
