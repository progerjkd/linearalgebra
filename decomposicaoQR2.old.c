#include "matrix.h"
#include <math.h>

//#define VERBOSE

void qr2(matrix *A, double epsilon){

	int m = A->rows;
	matrix *D = newMatrix(A->rows, 1);
	matrix *B = copyMatrix(A);

	while(m > 1){
		printf("M = %d\n", m);
		double B_mm;
		getElement(B, m, m-1, &B_mm);
		while(fabs(B_mm) >= epsilon){
			// calculate shift
			// S=eig(B(m-1:m,m-1:m))

//			matrix Btmp = newMatrix(A->rows-1, A->cols-1);

			matrix *Beigenvectors = newMatrix(A->rows-1, A->cols-1);
			matrix *Bnew = newMatrix(2, 2);
			subMatrix(B, m-1, m, m-1, m, Bnew);

			printf("\nsubmatriz Bnew =\n");
			printMatrix(Bnew);
			matrix *s = newMatrix(2, 1);
			jacobi(Bnew, s, Beigenvectors);
			printf("\n S =\n");
			printMatrix(s);

			// [j,k]=min([abs(B(m,m)*[1 1]'-S)])
			matrix *tmp_11 = newMatrix(2, 1);
			matrix *tmp2 = newMatrix(2, 1);
			setElement(tmp_11, 1, 1, 1);
			setElement(tmp_11, 2, 1, 1);
			double B_tmp;
			getElement(B, m, m, &B_tmp);
			printf("B_mm = %f\n", B_tmp);

			productByScalar(tmp_11, B_tmp, tmp2);
			printf("\nB(m,m)  * [1 1] =\n");
			printMatrix(tmp2);
			matrix *tmp1 = newMatrix(2, 1);
			subtraction(tmp2, s, tmp1);
			printf("\nB(m,m)  * [1 1]' - S =\n");
			printMatrix(tmp1);

			matrixAbs(tmp1, tmp2);
			printf("abs = \n");
			printMatrix(tmp2);

			_minElemVec minElemVec = getMinElemVec(tmp2, 1);
			double j = minElemVec.min;
			int k = minElemVec.pos;

			printf("\nj = %f, k = %d\n", j, k);

			// QR factorization of B
			// [Q,U]=qr(B-S(k)*eye(m));

			double s_k;
			getElement(s, k, 1, &s_k);
			matrix *tmp3 = newMatrix(m, m);
			matrix *tmp4 = newMatrix(m, m);
			matrix *eye = newMatrix(m, m);
			identity(eye);
			productByScalar(eye, s_k, tmp3);
			subtraction(B, tmp3, tmp4);

			matrix *Q = newMatrix(m, m);
			matrix *U = newMatrix(m, m);
			qr(tmp4, Q, U);


			// Calculate next B
			// B=U*Q+S(k)*eye(m);
			
			matrix *UQ = newMatrix(m, m);
			product(U, Q, UQ);

			matrix *tmp5 = newMatrix(m, m);
			getElement(s, k, 1, &s_k);
			productByScalar(eye, s_k, tmp5);
			sum(UQ, tmp5, tmp4);

			B = copyMatrix(tmp4);
			printf("\n B =\n");
			printMatrix(B);

			getElement(B, m, m-1, &B_mm);
			printf("\nAQUI B(m, m-1) = %f\n", B_mm);
		}

		// Place mth eigenvalue in A(m,m)
		for(int i=1; i<=m; i++)
			for(int j=1; j<=m; j++){
				double B_ij;
				getElement(B, i, j, &B_ij);
				setElement(A, i, j, B_ij);
			}


		// Repeat process on the m-1 x m-1 submatrix of A
		// m=m-1;   
		// B=A(1:m,1:m);   

		m--;

		matrix *tmp6 = newMatrix(m, m);

		for(int i=1; i<=m; i++)
			for(int j=1; j<=m; j++){
				double B_ij;
				getElement(A, i, j, &B_ij);
				setElement(tmp6, i, j, B_ij);
			}

		B = tmp6;

	}

	matrix *d = newMatrix(A->rows, 1);
	diagonalToVector(A, d);
	printf("\n*** INSIDE QR2 ***\nA =\n");
	printMatrix(A);
	printf("\nd =\n");
	printMatrix(d);


//	return *d;
}
