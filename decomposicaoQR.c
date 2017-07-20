#include "matrix.h"
#include <math.h>

//#define VERBOSE

// Processo de Gram-Schmidt modificado

bool __DELETE_qr__ = false;

_qr qr(matrix *A){

	int m = A->rows;
	int n = A->cols;

	matrix *V = copyMatrix(A);

	_qr QR;
	QR.Q = newMatrix(m, n);
	QR.R = newMatrix(n, n);

#ifdef VERBOSE
	printf("\nEntrada - Decomposição QR\n");

	printMatrix2(A, "A");
	printMatrix2(QR.Q, "Q");
	printMatrix2(QR.R, "R");

#endif

	__DELETE_norma__ = true;
	__DELETE_dotProduct2__ = true;

	for(int i=1; i<=n; i++){
		// R(i,i)=norm(V(:,i))
		setElement(QR.R, i, i, norma(matrixColToVector(V, i), 2));
		
		// Q(:,i)=V(:,i)/R(i,i);
		for(int x=1; x<=m; x++)
			setElement(QR.Q, x, i, getElement2(V, x, i) / getElement2(QR.R, i, i));

		for(int j=i+1; j<=n; j++){
			// R(i,j)=Q(:,i)'*V(:,j);
			setElement(QR.R, i, j, dotProduct2(matrixColToVector(QR.Q, i), matrixColToVector(V, j)));

			// V(:,j)=V(:,j)-R(i,j)*Q(:,i);
			for(int y=1; y<=m; y++)
				setElement(V, y, j, getElement2(V, y, j) - getElement2(QR.R, i, j) * getElement2(QR.Q, y, i));
		}

	}
	
	deleteMatrix(V);


#ifdef VERBOSE

	printf("\nSaída - Decomposição QR\n");

	printMatrix2(A, "A");
	printMatrix2(QR.Q, "Q");
	printMatrix2(QR.R, "R");
#endif

	if(__DELETE_qr__)
		deleteMatrix(A);

	return QR;
	
}
