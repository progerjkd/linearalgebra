#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <sys/resource.h>
#include <errno.h>



// TO DO: DEFINIR TODOS OS HEADERS NO matrix.h e chamá-lo aqui

// AQUI OU NO HEADER?
#define ELEM(mtx, row, col) \
  mtx->data[(col-1) * mtx->rows + (row-1)]

// AQUI OU NO HEADER?
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

double getElement2(matrix * mtx, int row, int col) {
  if (!mtx) return -1;
  assert (mtx->data);
  if (row <= 0 || row > mtx->rows ||
      col <= 0 || col > mtx->cols)
    return -2;

  return ELEM(mtx, row, col);
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

int printMatrix2(matrix * mtx, char *name) {
  if (!mtx) return -1;
  
  printf("\n%s =\n\n", name);
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
  printf("\n");
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

matrix *transpose2(matrix * in) {
	if (!in){
		printf("transpose: input matrix is null\n"); 
		return -1;
	}

  matrix *out = newMatrix(in->cols, in->rows);

  int row, col;
  for (row = 1; row <= in->rows; row++)
    for (col = 1; col <= in->cols; col++)
      ELEM(out, col, row) = ELEM(in, row, col);
  return out;
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

bool __DELETE_sum2__ = false;
matrix *sum2(matrix * mtx1, matrix * mtx2) {
	if (!mtx1 || !mtx2) return -1;
	if (mtx1->rows != mtx2->rows || mtx1->cols != mtx2->cols )
		return -2;

	matrix *sum = newMatrix(mtx1->rows, mtx1->cols);
	int row, col;
	for (col = 1; col <= mtx1->cols; col++)
		for (row = 1; row <= mtx1->rows; row++)
			ELEM(sum, row, col) = ELEM(mtx1, row, col) + ELEM(mtx2, row, col);

	if(__DELETE_sum2__){
		deleteMatrix(mtx1);
		deleteMatrix(mtx2);
	}

	return sum;
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

bool __DELETE_product2__ = false;
matrix *product2(matrix * A, matrix * B) {
	if (!A || !B){
		printf("product: input matrix is null\n"); 
		return -1;
	}
	if (A->cols != B-> rows){
		printf("product: Incompatible matrix sizes\n"); 
		return -2;
	}

	matrix *prod = newMatrix(A->rows, B->cols);
	int row, col, k;
	for (col = 1; col <= B->cols; col++)
		for (row = 1; row <= A->rows; row++) {
			double val = 0.0;
			for (k = 1; k <= A->cols; k++)
				val += ELEM(A, row, k) * ELEM(B, k, col);
			ELEM(prod, row, col) = val;
		}

	if(__DELETE_product2__){
		deleteMatrix(A);
		deleteMatrix(B);
	}

	return prod;
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

bool __DELETE_dotProduct2__ = false;
double dotProduct2(matrix * v1, matrix * v2) {
	if (!v1 || !v2) return -1;
	if (v1->cols != 1 || v2->cols != 1) return -2;
	if (v1->rows != v2->rows) return -3;

	double dot = 0;
	int i;
	for (i = 1; i <= v1->rows; i++)
		dot += ELEM(v1, i, 1) * ELEM(v2, i, 1);

	if(__DELETE_dotProduct2__){
		deleteMatrix(v1);
		deleteMatrix(v2);
	}

	return dot;
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

matrix *identity2(int order) {
	if(order < 1){
		printf("error: Order of identity matrix must be >= 1\n");
		return -1;
	}
  matrix *eye = newMatrix(order, order);
  int row, col;
  for (col = 1; col <= eye->cols; col++)
    for (row = 1; row <= eye->rows; row++)
      if (row == col) 
        ELEM(eye, row, col) = 1.0;
      else 
        ELEM(eye, row, col) = 0.0;
  return eye;
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
int vectorToDiagonal(matrix * v, matrix * mtx) {
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

// subtracts matrices A - B
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

bool __DELETE_subtraction2__ = false;
matrix *subtraction2(matrix *A, matrix *B){
	if (!A || !B){
		printf("subtraction: input matrix is null\n"); 
		return -1;
	}
	if (A->cols != B-> cols || A->rows != B->rows){
		printf("subtraction: Incompatible matrix sizes\n"); 
		return -2;
	}

	matrix *sub = newMatrix(A->rows, A->cols);
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij = getElement2(A, i, j);
			double B_ij = getElement2(B, i, j);
			setElement(sub, i, j, A_ij - B_ij);
		}
	}

	if(__DELETE_subtraction2__){
		deleteMatrix(A);
		deleteMatrix(B);
	}

	return sub;
}

// divides A / B
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

matrix *division2(matrix *A, matrix *B){
	if (!A || !B){
		printf("division: input matrix is null\n"); 
		return -1;
	}
	if (A->cols != B-> cols || A->rows != B->rows){
		printf("division: Incompatible matrix sizes\n"); 
		return -2;
	}
	matrix *div = newMatrix(A->rows, A->cols);
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij = getElement2(A, i, j);
			double B_ij = getElement2(B, i, j);
			setElement(div, i, j, A_ij / B_ij);
		}
	}
	return div;
}

// sums matrix A by scalar x
int sumByScalar(matrix *A, double x, matrix *sum){
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			setElement(sum, i, j, A_ij + x);
		}
	}
}

matrix *sumByScalar2(matrix *A, double x){
	if (!A){
		printf("sumByScalar: input matrix is null\n"); 
		return -1;
	}
	matrix *sum = newMatrix(A->rows, A->cols);
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij = getElement2(A, i, j);
			setElement(sum, i, j, A_ij + x);
		}
	}
	return sum;
}

// subtracts matrix A by scalar x
int subtractionByScalar(matrix *A, double x, matrix *sub){
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			setElement(sub, i, j, A_ij - x);
		}
	}
}

matrix *subtractionByScalar2(matrix *A, double x){
	if (!A){
		printf("sumByScalar: input matrix is null\n"); 
		return -1;
	}
	matrix *sub = newMatrix(A->rows, A->cols);
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij = getElement2(A, i, j);
			setElement(sub, i, j, A_ij - x);
		}
	}
	return sub;
}

// multiplies matrix A by scalar x
int productByScalar(matrix *A, double x, matrix *prod){
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			setElement(prod, i, j, A_ij * x);
		}
	}
}

bool __DELETE_productByScalar2__ = false;
matrix *productByScalar2(matrix *A, double x){
	if (!A){
		printf("productByScalar: input matrix is null\n"); 
		return -1;
	}
	matrix *prod = newMatrix(A->rows, A->cols);
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij = getElement2(A, i, j);
			setElement(prod, i, j, A_ij * x);
		}
	}
	if(__DELETE_productByScalar2__)
		deleteMatrix(A);

	return prod;
}

// divides matrix A by scalar x
int divisionByScalar(matrix *A, double x, matrix *div){
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			setElement(div, i, j, A_ij / x);
		}
	}
}

matrix *divisionByScalar2(matrix *A, double x){
	if (!A){
		printf("divisionByScalar: input matrix is null\n"); 
		return -1;
	}
	matrix *div = newMatrix(A->rows, A->cols);
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			double A_ij = getElement2(A, i, j);
			setElement(div, i, j, A_ij / x);
		}
	}
	return div;
}

int isSymmetric(matrix *A) {
	if (!A){
		printf("isSymmetric: input matrix is null\n"); 
		return -1;
	}

	if (!isSquare(A)) return 0;

	if (areEqual(A, transpose2(A)))
		return 1;
	else
		return 0;
}

// prints the lower triangular part of matrix A
int printLowerTriangular(matrix *A) {
	if (!A){
		printf("printLowerTriangular: input matrix is null\n"); 
		return -1;
	}

	if (!isSquare(A)) return 0;
	int row, col;
	// looks at positions below the diagonal
	for (col = 1; col <= A->cols; col++)
		for (row = 1; row <= col; row++)
			printf("A[%d,%d] = %f\n", row, col, ELEM(A, row, col));
}

int printLowerTriangularCol(matrix *A) {
	if (!A){
		printf("printLowerTriangularCol: input matrix is null\n"); 
		return -1;
	}
  if (!isSquare(A)) return 0;
  int row, col, start;

  start = 0;
  for (col = 1; col <= A->cols; col++){
	  start++;
	  for (row = start; row <= A->rows; row++){
    	    printf("A[%d,%d] = %f\n", row, col, ELEM(A, row, col));
    	  }
  }
}

// Extract the upper triangular from A and stores it in B
// k is the elements on and above  the k-th diagonal of A
// k = 0 is the main diagonal, k > 0 is above the main diagonal
// k <0 is below the main diagonal
int getUpperTriangular(matrix *A, matrix *B, int k) {
  if (!isSquare(A)) return 0;
  int row, col;
  // looks at positions above the diagonal
  for (col = 1; col <= A->cols; col++)
    for (row = 1; row <= col-k; row++){
	    double A_ij;
	    getElement(A, row, col, &A_ij);
	    setElement(B, row, col, A_ij);
//    	printf("A[%d,%d] = %f\n", row, col, ELEM(A, row, col));
    }
}

matrix *getUpperTriangular2(matrix *A, int k) {
  if (!isSquare(A)) return 0;
  matrix *B = newMatrix(A->rows, A->cols);
  int row, col;
  // looks at positions above the diagonal
  for (col = 1; col <= A->cols; col++)
    for (row = 1; row <= col-k; row++){
	    double A_ij = getElement2(A, row, col);
	    setElement(B, row, col, A_ij);
//    	printf("A[%d,%d] = %f\n", row, col, ELEM(A, row, col));
    }
  return B;
}

// Extract the lower triangular from A and stores it in B
// k is the elements on and below the k-th diagonal of A
// k = 0 is the main diagonal, k > 0 is above the main diagonal
// k <0 is below the main diagonal
int getLowerTriangular(matrix *A, matrix *B, int k) {
  if (!isSquare(A)) return 0;
  int row, col;
  // looks at positions below the diagonal
  for (col = 1; col <= A->cols; col++)
    for (row = A->rows; row >= col-k; row--){
	    double A_ij;
	    getElement(A, row, col, &A_ij);
	    setElement(B, row, col, A_ij);
//    	printf("A[%d,%d] = %f\n", row, col, ELEM(A, row, col));
    }
}

matrix *getLowerTriangular2(matrix *A, int k) {
  if (!isSquare(A)) return 0;
  matrix *B = newMatrix(A->rows, A->cols);
  int row, col;
  // looks at positions below the diagonal
  for (col = 1; col <= A->cols; col++)
    for (row = A->rows; row >= col-k; row--){
	    double A_ij = getElement2(A, row, col);
	    setElement(B, row, col, A_ij);
//    	printf("A[%d,%d] = %f\n", row, col, ELEM(A, row, col));
    }
  return B;
}

// NAO É USADA EM LUGAR ALGUM ----------- REMOVER !! --------------------------------------
// NAO É USADA EM LUGAR ALGUM ----------- REMOVER !! --------------------------------------
// NAO É USADA EM LUGAR ALGUM ----------- REMOVER !! --------------------------------------
// NAO É USADA EM LUGAR ALGUM ----------- REMOVER !! --------------------------------------
// NAO É USADA EM LUGAR ALGUM ----------- REMOVER !! --------------------------------------
// NAO É USADA EM LUGAR ALGUM ----------- REMOVER !! --------------------------------------
// NAO É USADA EM LUGAR ALGUM ----------- REMOVER !! --------------------------------------
matrix * addElementMatrix(int rows, int cols, matrix * mtx) {
	mtx->rows = rows;
	mtx->cols = cols;
	mtx->data = (double *) malloc(rows*cols*sizeof(double));
	int i;
	for (i = 0; i < rows*cols; i++)
		mtx->data[i] = 0.0;
	return mtx;
}

// UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
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


// FIM ----> UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// FIM ----> UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// FIM ----> UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// FIM ----> UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// FIM ----> UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO
// FIM ----> UTILIZADO NO symmetricStorageGaxpy - VERIFICAR SE PODE SER REMOVIDO








// swap lines x and y of matrix
int matrixSwapLines(matrix *matrix, int n, int x, int y){

	double tmp, tmp2;
	for (int i=1; i<=n; i++){
		getElement(matrix, x, i, &tmp);
		getElement(matrix, y, i, &tmp2);
		setElement(matrix, x, i, tmp2);
		setElement(matrix, y, i, tmp);
	}
}

int matrixSwapLines2(matrix *A, int x, int y){

	for (int i=1; i<=A->rows; i++){
		double tmp = getElement2(A, x, i);
		setElement(A, x, i, getElement2(A, y, i));
		setElement(A, y, i, tmp);
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

matrix *menorPrincipal2(matrix *A, int n){
	
	matrix *B = newMatrix(n, n);
	for(int i=1; i<=n; i++){
		for(int j=1; j<=n; j++){
			double A_ij = getElement2(A, i, j);
			setElement(B, i, j, A_ij);
		}
	}
	return B;
}

// norma-p vetorial
bool __DELETE_norma__ = false;
double norma(matrix *A, int p){

	double sum = 0;

	for(int i=1; i<=A->rows; i++){
		double A_i;
		getElement(A, i, 1, &A_i);
		sum += pow(fabs(A_i), p);
	}

	if(__DELETE_norma__)
		deleteMatrix(A);

	return pow(sum,(double)1/p);
}

// norma vetorial inferior
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

// norma matricial
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

// norma matricial inferior
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

// norma-1 matricial
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


typedef struct {
  double min;
  int pos;
} _minElemVec;

// retorna o valor e a posição do menor elemento da coluna col da matriz
bool __DELETE_getMinElemVec__ = false;
_minElemVec getMinElemVec(matrix *A, int col){

	double A_i;
	double min;
	int pos = 1;

	getElement(A, col, 1, &min);

	for(int i=1; i<=A->rows; i++){
		getElement(A, i, col, &A_i);
		if (A_i < min){
			min = A_i;
			pos = i;
		}
	}
	_minElemVec str;

	str.min = min;
	str.pos = pos;

	if(__DELETE_getMinElemVec__)
		deleteMatrix(A);

	return str;
}

// retorna os valores absolutos dos elementos da matriz A
void matrixAbs(matrix *A, matrix *B){

	double A_ij;
	for(int i=1; i<=A->rows; i++)
		for(int j=1; j<=A->cols; j++){
			getElement(A, i, j, &A_ij);
			setElement(B, i, j, fabs(A_ij));
		}
}

bool __DELETE_matrixAbs2__ = false;
matrix *matrixAbs2(matrix *A){

	matrix *B = newMatrix(A->rows, A->cols);
	for(int i=1; i<=A->rows; i++)
		for(int j=1; j<=A->cols; j++){
			setElement(B, i, j, fabs(getElement2(A, i, j)));
		}
	if(__DELETE_matrixAbs2__)
		deleteMatrix(A);

	return B;
}

// retorna a linha i da matriz como um vetor
matrix* matrixRowToVector(matrix *A, int i){

	matrix *vector = newMatrix(A->cols, 1);

	for(int j=1; j<=A->cols; j++){
		setElement(vector, j, 1, getElement2(A, i, j));
	}

	return vector;
}

// retorna a coluna j da matriz como um vetor
bool __DELETE_matrixColToVector__ = false;
matrix* matrixColToVector(matrix *A, int j){

	matrix *vector = newMatrix(A->rows, 1);

	for(int i=1; i<=A->rows; i++){
		double A_ij;
		getElement(A, i, j, &A_ij);
		setElement(vector, i, 1, A_ij);
	}

	if(__DELETE_matrixColToVector__)
		deleteMatrix(A);

	return vector;
}

// Copia a submatriz de A em B
void subMatrix(matrix *A, int startRow, int endRow, int startCol, int endCol, matrix *B){

	int b_i = 1;
	for(int i=startRow; i<=endRow; i++){
		int b_j = 1;
		for(int j=startCol; j<=endCol; j++){
			double A_ij;
			getElement(A, i, j, &A_ij);
			setElement(B, b_i, b_j, A_ij);
			b_j++;
		}
//		B->cols = b_j;
		b_i++;
	}
//	B->rows = b_i;
}

matrix *subMatrix2(matrix *A, int startRow, int endRow, int startCol, int endCol){

	matrix *B = newMatrix(endRow - startRow + 1, endCol - startCol + 1);
	int b_i = 1;
	for(int i=startRow; i<=endRow; i++){
		int b_j = 1;
		for(int j=startCol; j<=endCol; j++){
			setElement(B, b_i, b_j, getElement2(A, i, j));
			b_j++;
		}
//		B->cols = b_j;
		b_i++;
	}
//	B->rows = b_i;
	return B;
}

// copia a coluna j da matriz A na matriz B
void copyColumn(matrix *A, matrix *B, int j){

	for(int i=1; i<=A->rows; i++){
		setElement(B, i, 1, getElement2(A, i, j));
	}
}

// copia a coluna j da matriz A, a partir do n-ésimo elemento, na matriz B
void copyNColumn(matrix *A, matrix *B, int n, int j){

	for(int i=n; i<=A->rows; i++){
		//printf("i=%d j=%d = %f\n", i, j, A_ij);
		//printf("c_%d = %f\n", i-n+1, A_ij);
		setElement(B, i-n+1, 1, getElement2(A, i, j));
	}
}

// se os elementos da matriz A forem menores que o valor tol, os igualamos a 0
int tolerance(matrix *A, double tol){
	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			if(fabs(getElement2(A, i, j)) < tol)
				setElement(A, i, j, 0);
		}
	}
}

// copy the diagonal of matrix A into vector v
void diagonalToVector(matrix *A, matrix *v){

	int m = A->rows, n = A->cols;
	int lower = m < n ? m : n;

	for(int i=1; i<=lower; i++)
		for(int j=1; j<=lower; j++)
			if(i==j){
				double A_ij;
				getElement(A, i, j, &A_ij);
				setElement(v, i, 1, A_ij);
			}
}

matrix *diagonalToVector2(matrix *A){

	int m = A->rows, n = A->cols;
	int lower = m < n ? m : n;
	matrix *v = newMatrix(lower, 1);

	for(int i=1; i<=lower; i++)
		for(int j=1; j<=lower; j++)
			if(i==j)
				setElement(v, i, 1, getElement2(A, i, j));
	return v;
			
}


typedef struct {
  double complex x1;
  double complex x2;
} _eqSegGrau;

_eqSegGrau eqSegGrau(double a, double b, double c){
	
	_eqSegGrau x;
	double delta = pow(b, 2) - 4 * a * c;
//	double complex tmp1 = (-b + csqrt(delta))/(2*a);
//	double complex tmp2 = (-b + csqrt(delta))/(2*a);
	x.x1 = (-b + csqrt(delta))/(2*a);
	x.x2 = (-b - csqrt(delta))/(2*a);
//	printf("a = %f, b = %f, c = %f, delta = %f, sqrt(delta) = %f + %fi, tmp1 = %f, tmp2 = %f\n", a, b, c, delta, creal(csqrt(delta)), cimag(csqrt(delta)), tmp1, tmp2);
//	printf("dentro eqSegGrau x1: %f + %fi\n", creal(x.x1), cimag(x.x1));
//	printf("dentro eqSegGrau x2: %f + %fi\n", creal(x.x2), cimag(x.x2));

	return x;
}


// IMPLEMENTAR POSTO
// IMPLEMENTAR NORMA

// IMPLEMENTAR POSTO
// IMPLEMENTAR NORMA

void memoryUsage(void){
	struct rusage r_usage;
	getrusage(RUSAGE_SELF,&r_usage);
	printf("Memory usage: %ld Kbytes\n",r_usage.ru_maxrss);
}
