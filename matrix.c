#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ELEM(mtx, row, col) \
  mtx->data[(col-1) * mtx->rows + (row-1)]

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
      //printf("% 6.2f ", ELEM(mtx, row, col));
      printf("% 10.4f ", ELEM(mtx, row, col));
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

// copy the vector *v into diagonal of *mtx
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

int diagonalContainsZero(matrix * mtx) {
  if (!isSquare(mtx)) return 0;
  int row, col;
  for (col = 1; col <= mtx->cols; col++)
    for (row = 1; row <= mtx->rows; row++)
      // if the element is on the diagonal and is 0
      if (row == col && ELEM(mtx, row, col) == 0.0)
        return 1;
  return 0;
}

//  My code begins here

int subtraction(matrix *A, matrix *B, matrix *sub){
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij;
			double B_ij;
			getElement(A, i, j, &A_ij);
			getElement(B, i, j, &B_ij);
			setElement(sub, i, j, A_ij - B_ij);
		}
	}
}

int division(matrix *A, matrix *B, matrix *div){
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij;
			double B_ij;
			getElement(A, i, j, &A_ij);
			getElement(B, i, j, &B_ij);
			setElement(div, i, j, A_ij / B_ij);
		}
	}
}

int productByScalar(matrix *A, double x, matrix *prod){
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			setElement(prod, i, j, A_ij * x);
		}
	}
}

int divisionByScalar(matrix *A, double x, matrix *div){
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			setElement(div, i, j, A_ij / x);
		}
	}
}

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

int matrixSwapLines(matrix *matrix, int n, int x, int y){

	double tmp, tmp2;
	for (int i=1; i<=n; i++){
		getElement(matrix, x, i, &tmp);
		getElement(matrix, y, i, &tmp2);
		setElement(matrix, x, i, tmp2);
		setElement(matrix, y, i, tmp);
	}
}

// copia o menor principal de ordem n da matriz A para a matriz B
int menorPrincipal(matrix *A, matrix *B, int n){
	
//	B = newMatrix(A->rows - 1, A->cols - 1);
	for(int i=1; i<=n; i++){
		for(int j=1; j<=n; j++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			setElement(B, i, j, A_ij);
		}
	}
}

double norma(matrix *A, int p){

	double sum = 0;
	for(int i=1; i<=A->rows; i++){
		double A_i;
		getElement(A, i, 1, &A_i);
		sum += pow(fabs(A_i), p);
	}

/*	printf("1/%d = %f\n", p, (double)1/p);
	printf("pow(1/%d, %f) = %f\n", p, sum, pow(sum, (double)1/p));
*/
	return pow(sum,(double)1/p);
}

double normaInf(matrix *A){

	int i;
	double maior;

	getElement(A, 1, 1, &maior);

	for(i=1; i<=A->rows; i++){
		double A_i;
		getElement(A, i, 1, &A_i);
		A_i = fabs(A_i);
		if(A_i > maior)
			maior = A_i;
	}

	return maior;
}

double normaMatricial(matrix *A, int p){

	double sum = 0;
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			sum += pow(fabs(A_ij), p);
		}
	}

/*	printf("1/%d = %f\n", p, (double)1/p);
	printf("pow(1/%d, %f) = %f\n", p, sum, pow(sum, (double)1/p));
*/
	return pow(sum,(double)1/p);
}

double normaMatricialInf(matrix *A){

	double maior = 0;

//	getElement(A, 1, 1, &maior);

	for(int i=1; i<=A->rows; i++){
		double sumRow = 0;

		for(int j=1; j<=A->cols; j++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			sumRow += fabs(A_ij);
		}
			if(sumRow > maior)
				maior = sumRow;
	}

	return maior;
}

double normaMatricial1(matrix *A){

	double maior = 0;

//	getElement(A, 1, 1, &maior);

	for(int j=1; j<=A->cols; j++){
		double sumCol = 0;

		for(int i=1; i<=A->rows; i++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			sumCol += fabs(A_ij);
		}
			if(sumCol > maior)
				maior = sumCol;
	}

	return maior;
}


// IMPLEMENTAR POSTO
// IMPLEMENTAR NORMA

// IMPLEMENTAR POSTO
// IMPLEMENTAR NORMA
