#include "matrix.h"
#include <math.h>

#define VERBOSE

matrix* matrix_minor(matrix *x, int d){

	matrix *m = newMatrix(x->rows, x->cols);

	for(int i=1; i<=d; i++)
		setElement(m, i, i, 1);

	for(int i=d; i<=x->rows; i++)
		for(int j=d; j<=x->cols; j++)
			setElement(m, i, j, getElement2(x, i, j));
		
	return m;
}

/* c = a + b * s */
matrix* vmadd(matrix *a, matrix *b, double s, matrix *c, int n){
	for (int i=1; i<=n; i++)
		setElement(c, i, 1, getElement2(a, i, 1) + s * getElement2(b, i, 1));

	return c;
}

/* y = x / d */
matrix* vdiv(matrix *x, double d, matrix *y, int n){
	for (int i=1; i<=n; i++)
		setElement(y, i, 1, getElement2(x, i, 1) / d);

	return y;
}

/* m = I - v v^T */
matrix* vmul(matrix *v, int n){

	matrix *x = newMatrix(n, n);

	for(int i=1; i<=n; i++)
		for(int j=1; j<=n; j++)
			setElement(x, i, j, - 2 * getElement2(v, i, 1) * getElement2(v, j, 1));

	for (int i=1; i<=n; i++)
		setElement(x, i, i, getElement2(x, i, i) + 1);
	
	return x;
}

matrix* matrix_mul(matrix *x, matrix *y)
{       
	if (x->cols != y->rows)
		return 0;
	
	matrix *r = newMatrix(x->rows, y->cols);

	for(int i=1; i<=x->rows; i++)
		for (int j=1; j<=y->cols; j++)
			for (int k=1; k<=x->cols; k++)
				setElement(r, i, j, getElement2(r, i, j) + getElement2(x, i, k) * getElement2(y, k, j));

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


void qrold(matrix *M, matrix *Q, matrix *R){

	matrix *q[M->rows+1];
	matrix *z = M, *z1;

	matrix *Q_ADDR;
	matrix *R_ADDR;

	Q_ADDR = Q;
	R_ADDR = R;

#ifdef VERBOSE
	printf("\nEntrada - Decomposição QR\n");

	printf("A =\n");
	printMatrix(M);

	printf("\nQ =\n");
	printMatrix(Q);

	printf("\nR =\n");
	printMatrix(R);
#endif

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

#ifdef VERBOSE

	printf("\nSaída - Decomposição QR\n");

	printf("\nA =\n");
	printMatrix(M);

	printf("\nQ =\n");
        printMatrix(Q);

        printf("\nR =\n");
        printMatrix(R);
#endif

}
