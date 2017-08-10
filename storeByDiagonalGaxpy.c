#include "matrix.h"

int diagPrint(matrix * mtx) {
  if (!isSquare(mtx)) return 0;
  int row, col;
  int vec_sequency = 1;
  int start = 0;
  int vec_size = 0;
  int x;
  matrix * diag;

  // calculate diag's size
  for (row = 1; row <= mtx->rows; row++){
	  int i = row;
	  int j = 1;  
	  while ((i <= mtx->rows) && (j <= mtx->cols)){
		  printf("%d,%d\n", i, j);
		  i++;
		  j++;
	  }
  }
  return 0;
}

matrix * diag(matrix * mtx) {
  if (!isSquare(mtx)) return 0;
  int row, col;
  int diag_sequency = 1;
  int diag_size = 0;
  matrix * diag;

  // calculate diag's size
  for (row = 1; row <= mtx->rows; row++){
	  int i = row;
	  int j = 1;  
	  while ((i <= mtx->rows) && (j <= mtx->cols)){
		  diag_size++;
		  i++;
		  j++;
	  }
  }

  diag = newMatrix(diag_size, 1);

  for (row = 1; row <= mtx->rows; row++){
	  int i = row;
	  int j = 1;  
	  while ((i <= mtx->rows) && (j <= mtx->cols)){
	    	  setElement(diag, diag_sequency, 1, ELEM(mtx, i, j));
		//printf("\n diag_sequency = %d, i = %d, j = %d, ELEM(mtx,i,j) = %f\n", diag_sequency, i, j, ELEM(mtx,i,j));
		  diag_sequency++;
		  i++;
		  j++;
	  }
  }

  return diag;
}


matrix * D(matrix * mtx, int index) {
  if (!isSquare(mtx)) return 0;
  matrix *D = newMatrix(mtx->rows, mtx->cols);

  // calculate diag's size
	  int i = 1;
	  int j = 1+index;  
	  while ((i <= mtx->rows) && (j <= mtx->cols)){
		  //printf("\ni = %d, j = %d, ELEM(mtx,i,j) = %f\n", i, j, ELEM(mtx,i,j));
  		  setElement(D, i, j, ELEM(mtx, i, j));
		  i++;
		  j++;
	  }

  return D;
}



int main() {
  matrix * A, * Ac, * B, * c, * d, * M, * ct, * mdp, *tmp, *y, *x, *tmp2, *A_diag;
  double dp;
  int i, j, k, n = 3;

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

  printf("\nDiagonais:\n");
  A_diag = diag(A);  
  printMatrix2(A_diag, "A.diag");

  printf("\ny = A x + y\n");

  printf("\nEntrada - storeByDiagonalGaxpy\n");

  y = newMatrix(n, 1);
  printMatrix2(y, "y");

  x = newMatrix(n, 1);
  setElement(x, 1, 1, 2.0);
  setElement(x, 2, 1, 3.0);
  setElement(x, 3, 1, 4.0);
  printMatrix2(x, "x");



  for(i=1; i<=n; i++){
	double x_i, y_i, A_diag_value, gaxpy;
	getElement(x, i, 1, &x_i);
	getElement(y, i, 1, &y_i);
	getElement(A_diag, i, 1, &A_diag_value);
	gaxpy = A_diag_value * x_i + y_i;
	setElement(y, i, 1, gaxpy);
	
  }

  for(k=1; k<=n-1; k++){
  	int t = n * k - k * (k-1) /2;
	for(i=1; i<=n-k; i++){

		double x_ik, y_i, A_diag_it_value, gaxpy;
		getElement(x, i+k, 1, &x_ik);
		getElement(y, i, 1, &y_i);
		getElement(A_diag, i+t, 1, &A_diag_it_value);
		gaxpy = A_diag_it_value * x_ik + y_i;
		setElement(y, i, 1, gaxpy);
	}
	for(i=1; i<=n-k; i++){

		double x_i, y_ik, A_diag_it_value, gaxpy;
		getElement(x, i, 1, &x_i);
		getElement(y, i+k, 1, &y_ik);
		getElement(A_diag, i+t, 1, &A_diag_it_value);
		gaxpy = A_diag_it_value * x_i + y_ik;
		setElement(y, i+k, 1, gaxpy);
	}
  }
  printf("\nSaída - storeByDiagonalGaxpy\n");

  printMatrix2(y, "y");
  printMatrix2(x, "x");
  printMatrix2(A, "A");

  return 0;
}
