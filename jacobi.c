#include <stdio.h>
#include <stdbool.h>
#include "matrix.h"
#include <math.h>
#include <float.h>


int jacobi(matrix *A, matrix *e, matrix *V){

	if(!isSquare(A) || !isSquare(V)){
		printf("Matrix A or V are not square.\nQuitting.\n");
		return -1;
	}

	printf("\nEntrada - Método de Jacobi\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\ne =\n");
	printMatrix(e);

	printf("\nV =\n");
	printMatrix(V);

	bool changed;
	int sweeps=0, n=A->rows;

	for(int i=1; i<=n; i++){
		double A_ii;
		getElement(A, i, i, &A_ii);
		setElement(e, i, 1, A_ii);
	}
	identity(V);

	do{
		sweeps++;
		changed = false;
		int p, q;

		for(p=1; p<=n; p++)
			for(q=p+1; q<=n; q++){
				double a_pp, a_qq, a_pq, phi, c, s, a_pp1, a_qq1;

				getElement(e, p, 1, &a_pp);
				getElement(e, q, 1, &a_qq);
				getElement(A, p, q, &a_pq);
				phi = 0.5 * atan2(2 * a_pq, a_qq - a_pp);
				c = cos(phi);
				s = sin(phi);
				a_pp1 = c * c * a_pp-2 * s * c * a_pq+s * s * a_qq;
				a_qq1 = s * s * a_pp+2 * s * c * a_pq+c * c * a_qq;


				if(a_pp1 != a_pp || a_qq1 != a_qq){
					changed = true;
					setElement(e, p, 1, a_pp1);
					setElement(e, q, 1, a_qq1);
					setElement(A, p, q, 0);
					double a_ip, a_pi, a_iq, a_qi, v_ip, v_iq;

					for(int i=1; i<p; i++){
						getElement(A, i, p, &a_ip);
						getElement(A, i, q, &a_iq);
						setElement(A, i, p, c * a_ip - s * a_iq);
						setElement(A, i, q, c * a_iq + s * a_ip);
					}

					for(int i=p+1; i<q; i++){
						getElement(A, p, i, &a_pi);
						getElement(A, i, q, &a_iq);
						setElement(A, p, i, c * a_pi - s * a_iq);
						setElement(A, i, q, c * a_iq + s * a_pi);
					}

					for(int i=q+1; i<=n; i++){
						getElement(A, p, i, &a_pi);
						getElement(A, q, i, &a_qi);
						setElement(A, p, i, c * a_pi - s * a_qi);
						setElement(A, q, i, c * a_qi + s * a_pi);
					}

					for(int i=1; i<=n; i++){
						getElement(V, i, p, &v_ip);
						getElement(V, i, q, &v_iq);
						setElement(V, i, p, c * v_ip - s * v_iq);
						setElement(V, i, q, c * v_iq + s * v_ip);
					}
				}
			}
	}while(changed);

	printf("\nSaída - Método de Jacobi\n");

	printf("\nA = (should be left-triangular)\n");
	printMatrix(A);

	printf("\ne =\n");
	printMatrix(e);

	printf("\nV =\n");
	printMatrix(V);

//	printf("\nsweeps = %d\n", sweeps);

	return sweeps;


}
