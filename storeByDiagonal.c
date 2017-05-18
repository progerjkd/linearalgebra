#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int rows;
  int cols;
  double * data;
} matrix;

/* Creates a ``rows by cols'' matrix with all values 0.  
 * Returns NULL if rows <= 0 or cols <= 0 and otherwise a
 * pointer to the new matrix.
 */
matrix * newMatrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0) return NULL;

  // allocate a matrix structure
  matrix * m = (matrix *) malloc(sizeof(matrix));

  // set dimensions
  m->rows = rows;
  m->cols = cols;

  // allocate a double array of length rows * cols
  m->data = (double *) malloc(rows*cols*sizeof(double));
  // set all data to 0
  int i;
  for (i = 0; i < rows*cols; i++)
    m->data[i] = 0.0;

  return m;
}

/* Deletes a matrix.  Returns 0 if successful and -1 if mtx 
 * is NULL.
 */
int deleteMatrix(matrix * mtx) {
  if (!mtx) return -1;
  // free mtx's data
  assert (mtx->data);
  free(mtx->data);
  // free mtx itself
  free(mtx);
  return 0;
}

#define ELEM(mtx, row, col) \
  mtx->data[(col-1) * mtx->rows + (row-1)]

/* Copies a matrix.  Returns NULL if mtx is NULL.
 */
matrix * copyMatrix(matrix * mtx) {
  if (!mtx) return NULL;

  // create a new matrix to hold the copy
  matrix * cp = newMatrix(mtx->rows, mtx->cols);

  // copy mtx's data to cp's data
  memcpy(cp->data, mtx->data, 
         mtx->rows * mtx->cols * sizeof(double));

  return cp;
}

/* Sets the (row, col) element of mtx to val.  Returns 0 if
 * successful, -1 if mtx is NULL, and -2 if row or col are
 * outside of the dimensions of mtx.
 */
int setElement(matrix * mtx, int row, int col, double val) 
{
  if (!mtx) return -1;
  assert (mtx->data);
  if (row <= 0 || row > mtx->rows ||
      col <= 0 || col > mtx->cols)
    return -2;

  ELEM(mtx, row, col) = val;
  return 0;
}

/* Sets the reference val to the value of the (row, col) 
 * element of mtx.  Returns 0 if successful, -1 if either 
 * mtx or val is NULL, and -2 if row or col are outside of 
 * the dimensions of mtx.
 */
int getElement(matrix * mtx, int row, int col, 
               double * val) {
  if (!mtx || !val) return -1;
  assert (mtx->data);
  if (row <= 0 || row > mtx->rows ||
      col <= 0 || col > mtx->cols)
    return -2;

  *val = ELEM(mtx, row, col);
  return 0;
}

/* Sets the reference n to the number of rows of mtx.
 * Returns 0 if successful and -1 if mtx or n is NULL.
 */
int nRows(matrix * mtx, int * n) {
  if (!mtx || !n) return -1;
  *n = mtx->rows;
  return 0;
}

/* Sets the reference n to the number of columns of mtx.
 * Returns 0 if successful and -1 if mtx is NULL.
 */
int nCols(matrix * mtx, int * n) {
  if (!mtx || !n) return -1;
  *n = mtx->rows;
  return 0;
}

/* Prints the matrix to stdout.  Returns 0 if successful 
 * and -1 if mtx is NULL.
 */
int printMatrix(matrix * mtx) {
  if (!mtx) return -1;
  
  int row, col;
  for (row = 1; row <= mtx->rows; row++) {
    for (col = 1; col <= mtx->cols; col++) {
      // Print the floating-point element with
      //  - either a - if negative or a space if positive
      //  - at least 3 spaces before the .
      //  - precision to the hundredths place
      printf("% 6.2f ", ELEM(mtx, row, col));
    }
    // separate rows by newlines
    printf("\n");
  }
  return 0;
}

/* Writes the transpose of matrix in into matrix out.  
 * Returns 0 if successful, -1 if either in or out is NULL,
 * and -2 if the dimensions of in and out are incompatible.
 */
int transpose(matrix * in, matrix * out) {
  if (!in || !out) return -1;
  if (in->rows != out->cols || in->cols != out->rows)
    return -2;

  int row, col;
  for (row = 1; row <= in->rows; row++)
    for (col = 1; col <= in->cols; col++)
      ELEM(out, col, row) = ELEM(in, row, col);
  return 0;
}

/* Writes the sum of matrices mtx1 and mtx2 into matrix 
 * sum. Returns 0 if successful, -1 if any of the matrices 
 * are NULL, and -2 if the dimensions of the matrices are
 * incompatible.
 */
int sum(matrix * mtx1, matrix * mtx2, matrix * sum) {
  if (!mtx1 || !mtx2 || !sum) return -1;
  if (mtx1->rows != mtx2->rows ||
      mtx1->rows != sum->rows ||
      mtx1->cols != mtx2->cols ||
      mtx1->cols != sum->cols)
    return -2;

  int row, col;
  for (col = 1; col <= mtx1->cols; col++)
    for (row = 1; row <= mtx1->rows; row++)
      ELEM(sum, row, col) = 
        ELEM(mtx1, row, col) + ELEM(mtx2, row, col);
  return 0;
}

/* Writes the product of matrices mtx1 and mtx2 into matrix
 * prod.  Returns 0 if successful, -1 if any of the 
 * matrices are NULL, and -2 if the dimensions of the 
 * matrices are incompatible.
 */
int product(matrix * mtx1, matrix * mtx2, matrix * prod) {
  if (!mtx1 || !mtx2 || !prod) return -1;
  if (mtx1->cols != mtx2->rows ||
      mtx1->rows != prod->rows ||
      mtx2->cols != prod->cols)
    return -2;

  int row, col, k;
  for (col = 1; col <= mtx2->cols; col++)
    for (row = 1; row <= mtx1->rows; row++) {
      double val = 0.0;
      for (k = 1; k <= mtx1->cols; k++)
        val += ELEM(mtx1, row, k) * ELEM(mtx2, k, col);
      ELEM(prod, row, col) = val;
    }
  return 0;
}

int areEqual(matrix * mtx1, matrix * mtx2) {
  if (!mtx1 || !mtx2 ) return -1;
  if (mtx1->cols != mtx2->cols || mtx1->rows != mtx2->rows)  return -2;
  int row, col;
  // looks at positions below the diagonal
  for (row = 1; row <= mtx1->rows; row++) 
    for (col = 1; col <= mtx1->cols; col++)
      if (ELEM(mtx1, row, col) != ELEM(mtx2, row, col))
        return 0;
  return 1;
}

/* Writes the dot product of vectors v1 and v2 into 
 * reference prod.  Returns 0 if successful, -1 if any of
 * v1, v2, or prod are NULL, -2 if either matrix is not a 
 * vector, and -3 if the vectors are of incompatible 
 * dimensions.
 */
int dotProduct(matrix * v1, matrix * v2, double * prod) {
  if (!v1 || !v2 || !prod) return -1;
  if (v1->cols != 1 || v2->cols != 1) return -2;
  if (v1->rows != v2->rows) return -3;

  *prod = 0;
  int i;
  for (i = 1; i <= v1->rows; i++)
    *prod += ELEM(v1, i, 1) * ELEM(v2, i, 1);
  return 0;
}

int identity(matrix * m) {
  if (!m || m->rows != m->cols) return -1;
  int row, col;
  for (col = 1; col <= m->cols; col++)
    for (row = 1; row <= m->rows; row++)
      if (row == col) 
        ELEM(m, row, col) = 1.0;
      else 
        ELEM(m, row, col) = 0.0;
  return 0;
}

int isSquare(matrix * mtx) {
  return mtx && mtx->rows == mtx->cols;
}

int isDiagonal(matrix * mtx) {
  if (!isSquare(mtx)) return 0;
  int row, col;
  for (col = 1; col <= mtx->cols; col++)
    for (row = 1; row <= mtx->rows; row++)
      // if the element is not on the diagonal and not 0
      if (row != col && ELEM(mtx, row, col) != 0.0)
        // then the matrix is not diagonal
        return 0;
  return 1;
}

int isUpperTriangular(matrix * mtx) {
  if (!isSquare(mtx)) return 0;
  int row, col;
  // looks at positions below the diagonal
  for (col = 1; col <= mtx->cols; col++)
    for (row = col+1; row <= mtx->rows; row++) 
      if (ELEM(mtx, row, col) != 0.0)
        return 0;
  return 1;
}

int diagonal(matrix * v, matrix * mtx) {
  if (!v || !mtx ||
      v->cols > 1 || v->rows != mtx->rows ||
      mtx->cols != mtx->rows)
    return -1;
  int row, col;
  for (col = 1; col <= mtx->cols; col++)
    for (row = 1; row <= mtx->rows; row++)
      if (row == col) 
        ELEM(mtx, row, col) = ELEM(v, col, 1);
      else
        ELEM(mtx, row, col) = 0.0;
  return 0;
}

// INÍCIO DO MEU CÓDIGO

int isSymmetric(matrix * mtx) {
  if (!isSquare(mtx)) return 0;

  matrix * mtxt = newMatrix(mtx->rows, mtx->cols);
  transpose(mtx, mtxt);
  if (areEqual(mtx, mtxt))
	  return 1;
  else
	  return 0;
}

int printLowerTriangular(matrix * mtx) {
  if (!isSquare(mtx)) return 0;
  int row, col;
  // looks at positions below the diagonal
  for (col = 1; col <= mtx->cols; col++)
    for (row = 1; row <= col; row++){
    	printf("A[%d,%d] = %f\n", row, col, ELEM(mtx, row, col));
    }
}

int printLowerTriangularCol(matrix * mtx) {
  if (!isSquare(mtx)) return 0;
  int row, col, start;

  start = 0;
  for (col = 1; col <= mtx->cols; col++){
	  start++;
	  for (row = start; row <= mtx->rows; row++){
    	    printf("A[%d,%d] = %f\n", row, col, ELEM(mtx, row, col));
    	  }
  }
}

matrix * addElementMatrix(int rows, int cols, matrix * mtx) {
	mtx->rows = rows;
	mtx->cols = cols;
	mtx->data = (double *) malloc(rows*cols*sizeof(double));
	int i;
	for (i = 0; i < rows*cols; i++)
		mtx->data[i] = 0.0;
	return mtx;
}

matrix * vec(matrix * mtx) {
  if (!isSquare(mtx)) return 0;
  int row, col;
  int vec_sequency = 1;
  int start = 0;
  int vec_size = 0;
  matrix * vec;

  // calculate vec's size
  for (col = 1; col <= mtx->cols; col++){
	  start++;
	  for (row = start; row <= mtx->rows; row++){
		  vec_size++;
    	  }
  }

  vec = newMatrix(vec_size, 1);
  start = 0;

  for (col = 1; col <= mtx->cols; col++){
	  start++;
	  for (row = start; row <= mtx->rows; row++){
	    setElement(vec, vec_sequency, 1, ELEM(mtx, row, col));
	    vec_sequency++;
    	  }
  }

  return vec;
}

double vecGetLowerValue(matrix * vec, int i, int j, int n){

	int indice = (j - 1) * n - j * ( j - 1) / 2 + i;
	double value;

	getElement(vec, indice, 1, &value);

	return value;
}

double vecGetUpperValue(matrix * vec, int i, int j, int n){

	int indice = (i - 1) * n - i * ( i - 1) / 2 + j;
	double value;

	getElement(vec, indice, 1, &value);

	return value;
}


int vecGetLowerIndex(int i, int j, int n){

	int indice;
	indice = (j - 1) * n - j * ( j - 1) / 2 + i;

	return indice;
}

int vecGetUpperIndex(int i, int j, int n){

	int indice;
	indice = (i - 1) * n - i * ( i - 1) / 2 + j;

	return indice;
}

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
  printf("A =\n");
  printMatrix(A);

  y = newMatrix(n, 1);
  printf("\nVector y:\n");
  printMatrix(y);

  x = newMatrix(n, 1);
  setElement(x, 1, 1, 2.0);
  setElement(x, 2, 1, 3.0);
  setElement(x, 3, 1, 4.0);
  printf("\nVector x:\n");
  printMatrix(x);


 // tmp2=D(A,-2);
 // printf("\nD(A,-2)=\n");
 // printMatrix(tmp2);

  printf("\n COMECO ALGORITMO 2\n");


  A_diag = diag(A);  
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

  printf("\n Resultados:\n");
  printf("\nVector y:\n");
  printMatrix(y);

  printf("\nVector x:\n");
  printMatrix(x);

  printf("A =\n");
  printMatrix(A);

  printf("\n");

  return 0;
}
