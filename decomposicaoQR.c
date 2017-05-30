#include "matrix.h"
#include <math.h>

int decomposicaoQR(matrix *A, int n, matrix *L, matrix *U){

}

matrix* matrix_minor(matrix *x, int d){

	matrix *m = newMatrix(x->rows, x->cols);

	for(int i=1; i<=d; i++)
		setElement(m, i, i, 1);

	for(int i=d; i<=x->rows; i++)
		for(int j=d; j<=x->cols; j++){
			double x_ij;
			getElement(x, i, j, &x_ij);
			setElement(m, i, j, x_ij);
		}
	return m;
}

/* c = a + b * s */
matrix* vmadd(matrix *a, matrix *b, double s, matrix *c, int n){
	for (int i=1; i<=n; i++){
		double a_i, b_i;
		getElement(a, i, 1, &a_i);
		getElement(b, i, 1, &b_i);
		setElement(c, i, 1, a_i + s * b_i);
		setElement(c, i, 1, a_i + b_i * s);
	}

	return c;
}

/* y = x / d */
matrix* vdiv(matrix *x, double d, matrix *y, int n){
	for (int i=1; i<=n; i++){
		double x_i;
		getElement(x, i, 1, &x_i);
		setElement(y, i, 1, x_i / d);
	}

	return y;
}

/* m = I - v v^T */
matrix* vmul(matrix *v, int n){

	matrix *x = newMatrix(n, n);

	for(int i=1; i<=n; i++)
		for(int j=1; j<=n; j++){
			double v_i, v_j;
			getElement(v, i, 1, &v_i);
			getElement(v, j, 1, &v_j);
			setElement(x, i, j, -2 * v_i * v_j);
		}

	for (int i=1; i<=n; i++){
		double x_ii;
		getElement(x, i, i, &x_ii);
		setElement(x, i, i, x_ii + 1);
	}

	return x;
}

matrix* matrix_mul(matrix *x, matrix *y)
{       
	if (x->cols != y->rows)
		return 0;
	
	matrix *r = newMatrix(x->rows, y->cols);

	for(int i=1; i<=x->rows; i++)
		for (int j=1; j<=y->cols; j++)
			for (int k=1; k<=x->cols; k++){
				double x_ik, y_kj, r_ij;
				getElement(x, i, k, &x_ik);
				getElement(y, k, j, &y_kj);
				getElement(r, i, j, &r_ij);
				setElement(r, i, j, r_ij + x_ik * y_kj);
			}
	return r;
}


double vnorm(matrix *x, int n){

	double sum=0;
	for (int i=1; i<=n; i++){
		double x_i;
		getElement(x, i, 1, &x_i);
		sum += x_i * x_i;
	}
	return sqrt(sum);
}


void householder(matrix *M, matrix *R, matrix *Q){

	matrix *q[M->rows+1];
	matrix *z = M, *z1;

	matrix *Q_ADDR;
	matrix *R_ADDR;

	Q_ADDR = Q;
	R_ADDR = R;

	printf("\nEntrada - Decomposição QR\n");

	printf("A =\n");
	printMatrix(M);

	printf("\nQ =\n");
	printMatrix(Q);

	printf("\nR =\n");
	printMatrix(R);

	for(int k=1; k<=M->cols && k<=M->rows - 1; k++){
		matrix *e = newMatrix(M->rows, 1);
		matrix *x = newMatrix(M->rows, 1);
		double a;
		
		z1 = matrix_minor(z, k);

		if(!areEqual(z, M))
			deleteMatrix(z);
	
		z = z1;	

		x = matrixColToVector(z, k);
		a = vnorm(x, M->rows);

		double m_kk;
		getElement(M, k, k, &m_kk);
		if(m_kk > 0)
			a = -a;
		

		for(int i=1; i<=M->rows; i++)
			setElement(e, i, 1, (i == k) ? 1 : 0);

		vmadd(x, e, a, e, M->rows);
		vdiv(e, vnorm(e, M->rows), e, M->rows);
		q[k] = vmul(e, M->rows);
		z1 = matrix_mul(q[k], z);
		if(!areEqual(z, M))
			deleteMatrix(z);
		z = z1;
	}

	
	deleteMatrix(z);
	Q = q[1];	// *

	R = newMatrix(q[1]->rows, M->cols);
	R = matrix_mul(q[1], M);

	for(int i=2; i<=M->cols && i<=M->rows-1; i++){
		z1 = matrix_mul(q[i], Q);
		if(i>2)
			deleteMatrix(Q);
		
		Q = z1;
		deleteMatrix(q[i]);
	}
	deleteMatrix(q[1]);
	z = matrix_mul(Q, M);
	deleteMatrix(R);	//*
	R = z;

	matrix *Qt = newMatrix(Q->rows, Q->cols);
	transpose(Q, Qt);
	Q = copyMatrix(Qt);

	*Q_ADDR = *Q;
	*R_ADDR = *R;

	printf("\nSaída - Decomposição QR\n");

	printf("\nA =\n");
	printMatrix(M);

	printf("\nQ =\n");
        printMatrix(Q);

        printf("\nR =\n");
        printMatrix(R);

}
