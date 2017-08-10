#include "matrix.h"


int main() {
  matrix * A, * Ac, * B, * c, * d, * M, * ct, * mdp, *tmp, *y, *x;
  double dp;
  int i, j, n = 3;

/*  double y[n+1];
  double x[n+1];

  x[1] = 2;
  x[2] = 3;
  x[3] = 4;
*/

  A = newMatrix(n, n);
  setElement(A, 1, 1, 9);
  setElement(A, 1, 2, 7);
  setElement(A, 1, 3, 5);
  setElement(A, 2, 1, 7);
  setElement(A, 2, 2, 3);
  setElement(A, 2, 3, 2);
  setElement(A, 3, 1, 5);
  setElement(A, 3, 2, 2);
  setElement(A, 3, 3, 0);

  printMatrix2(A, "A");

  tmp=vec(A);
  printf("\nElementos da triangular inferior:\n");
  printf("\nA.vec =\n");
  printMatrix(tmp);

  printf("\ny = A x + y\n");

  printf("\nEntrada - symmetricStorageGaxpy\n");

  y = newMatrix(n, 1);
  printMatrix2(y, "y");

  x = newMatrix(n, 1);
  setElement(x, 1, 1, 2.0);
  setElement(x, 2, 1, 3.0);
  setElement(x, 3, 1, 4.0);
  printMatrix2(x, "x");


/*  printf("\nvecGetLowerIndex(1,1,3): %d\n", vecGetLowerIndex(1,1,3));
  printf("\nvecGetLowerIndex(2,1,3): %d\n", vecGetLowerIndex(2,1,3));
  printf("\nvecGetLowerIndex(2,2,3): %d\n", vecGetLowerIndex(2,2,3));
  printf("\nvecGetLowerIndex(3,1,3): %d\n", vecGetLowerIndex(3,1,3));
  printf("\nvecGetLowerIndex(3,2,3): %d\n", vecGetLowerIndex(3,2,3));
  printf("\nvecGetLowerIndex(3,3,3): %d\n", vecGetLowerIndex(3,3,3));

  printf("\nvecGetLowerValue(tmp,1,1,3): %f\n", vecGetLowerValue(tmp,1,1,3));
  printf("\nvecGetLowerValue(tmp,2,1,3): %f\n", vecGetLowerValue(tmp,2,1,3));
  printf("\nvecGetLowerValue(tmp,2,2,3): %f\n", vecGetLowerValue(tmp,2,2,3));
  printf("\nvecGetLowerValue(tmp,3,1,3): %f\n", vecGetLowerValue(tmp,3,1,3));
  printf("\nvecGetLowerValue(tmp,3,2,3): %f\n", vecGetLowerValue(tmp,3,2,3));
  printf("\nvecGetLowerValue(tmp,3,3,3): %f\n", vecGetLowerValue(tmp,3,3,3));
*/
// double vecGetLowerValue(matrix * vec, int i, int j, int n)
// double vecGetUpperValue(matrix * vec, int i, int j, int n)
// int vecGetLowerIndex(int i, int j, int n)
// int vecGetUpperIndex(int i, int j, int n)


  for(j=1; j<=n; j++){
	  for(i=1; i<=j-1; i++){
		//y[i] = vecGetUpperValue(tmp, i, j, n) * x[j] + y[i];
		double x_j, y_i, gaxpy;
		getElement(x,j,1, &x_j);
		getElement(y,i,1, &y_i);
		gaxpy = vecGetUpperValue(tmp, i, j, n) * x_j + y_i;
//		printf("\n teste=%f x_j=%f y_i=%f gaxpy=%f i=%d\n", teste, x_j, y_i, gaxpy, i);
		setElement(y, i, 1, gaxpy);
	  }
	  for(i=j; i<=n; i++){
		//y[i] = vecGetLowerValue(tmp, i, j, n) * x[j] + y[i];
		double x_j, y_i, gaxpy;
                getElement(x,j,1, &x_j);
		getElement(y,i,1, &y_i);
		gaxpy = vecGetLowerValue(tmp, i, j, n) * x_j + y_i;
//		printf("\n teste=%f x_j=%f y_i=%f gaxpy=%f i=%d\n", teste, x_j, y_i, gaxpy, i);
		setElement(y, i, 1, gaxpy);

	  }
  }

  printf("\nSaída - symmetricStorageGaxpy\n");

  printMatrix2(y, "y");
  printMatrix2(x, "x");
  printMatrix2(A, "A");

  return 0;
}
