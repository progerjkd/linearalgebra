#include "matrix.h"
#include <math.h>

// requerimento: matriz quadrada
// TODO: verificar se a matriz é nula

double laplace(matrix *A){

	double det = 0.0;

	if(!isSquare(A)){
		printf("Matriz A não é quadrada. Cancelando o cálculo do determinante\n");
		return -1;
	}

	if(A->rows == 1){	// 1x1
		getElement(A, 1, 1, &det);
	}
	else if(A->rows == 2){	// 2x2
		double A_11, A_12, A_21, A_22;

		getElement(A, 1, 1, &A_11);
		getElement(A, 1, 2, &A_12);
		getElement(A, 2, 1, &A_21);
		getElement(A, 2, 2, &A_22);

		det = A_11*A_22 - A_12*A_21;
	} else if(A->rows == 3){ // 3x3
		double A_11, A_12, A_13, A_21, A_22, A_23, A_31, A_32, A_33;
		getElement(A, 1, 1, &A_11);
		getElement(A, 1, 2, &A_12);
		getElement(A, 1, 3, &A_13);
		getElement(A, 2, 1, &A_21);
		getElement(A, 2, 2, &A_22);
		getElement(A, 2, 3, &A_23);
		getElement(A, 3, 1, &A_31);
		getElement(A, 3, 2, &A_32);
		getElement(A, 3, 3, &A_33);

		det = A_11 * A_22 * A_33 + A_12 * A_23 * A_31 + A_13 * A_21 * A_32 - A_13 * A_22 * A_31 - A_11 * A_23 * A_32 - A_12 * A_21 * A_33;
	} else{

		matrix *aux;
		int i_aux, j_aux, row, col, i;
		for(i=1; i<=A->rows; i++){
			double A_1i;
			getElement(A, 1, i, &A_1i);
			if(A_1i != 0){
				aux = newMatrix(A->rows - 1, A->cols - 1);
				i_aux = 1;
				j_aux = 1;

				for(row = 2; row<=A->rows; row++){
					for(col = 1; col<=A->cols; col++){
						if(col != i){
							double A_rc;
							getElement(A, row, col, &A_rc);
							setElement(aux, i_aux, j_aux, A_rc);
							j_aux++;
						}
					}

					i_aux++;
					j_aux = 1;
				}
				
				double A_0i;
				getElement(A, 0, i, &A_0i);

				printf("\n aux =\n");
				printMatrix(aux);

				det += pow(-1, i) * A_0i * laplace(aux);
			}
		}
	}

	return det;
}


/*


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

*/
