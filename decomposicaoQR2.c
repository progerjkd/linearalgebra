#include "matrix.h"
#include <math.h>

//#define VERBOSE

matrix *qr2(matrix *A, double epsilon){

	int m = A->rows;
	matrix *D = newMatrix(A->rows, 1);
	matrix *B = copyMatrix(A);

	__DELETE_matrixAbs2__ = true;
	__DELETE_getMinElemVec__ = true;
	__DELETE_qr__ = true;
	__DELETE_productByScalar2__ = true;
	__DELETE_product2__ = true;
	__DELETE_sum2__ = true;

	while(m > 1){
		while(fabs(getElement2(B, m, m-1)) >= epsilon){
			// calculate shift
			// S=eig(B(m-1:m,m-1:m))
	
			matrix *s = newMatrix(2, 1);
			matrix *Beigenvectors = newMatrix(A->rows-1, A->cols-1);
			jacobi(subMatrix2(B, m-1, m, m-1, m), s, Beigenvectors); //TODO OTIMIZAR JACOBI
			deleteMatrix(Beigenvectors);


			// [j,k]=min([abs(B(m,m)*[1 1]'-S)])
			matrix *tmp_11 = newMatrix(2, 1);
			setElement(tmp_11, 1, 1, 1);
			setElement(tmp_11, 2, 1, 1);

			_minElemVec minElemVec = getMinElemVec(matrixAbs2(subtraction2(productByScalar2(tmp_11, getElement2(B, m, m)), s)), 1);
			double j = minElemVec.min;
			int k = minElemVec.pos;


			// QR factorization of B
			// [Q,U]=qr(B-S(k)*eye(m));


			_qr QR;

			__DELETE_subtraction2__ = true;
			QR = qr(subtraction2(B, productByScalar2(identity2(m), getElement2(s, k, 1))));
			__DELETE_subtraction2__ = false;

			// Calculate next B
			// B=U*Q+S(k)*eye(m)

			B = sum2(product2(QR.R, QR.Q), productByScalar2(identity2(m), getElement2(s, k, 1)));
			
		}

		// Place mth eigenvalue in A(m,m)
		for(int i=1; i<=m; i++)
			for(int j=1; j<=m; j++)
				setElement(A, i, j, getElement2(B, i, j));
			


		// Repeat process on the m-1 x m-1 submatrix of A
		// m=m-1;   
		// B=A(1:m,1:m);   

		m--;
		deleteMatrix(B);
		B = subMatrix2(A, 1, m, 1, m);

	}


	return diagonalToVector2(A);
}
