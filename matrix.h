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
matrix * newMatrix(int rows, int cols);

/* Deletes a matrix.  Returns 0 if successful and -1 if mtx 
 * is NULL.
 */
int deleteMatrix(matrix * mtx);

/* Copies a matrix.  Returns NULL if mtx is NULL.
 */
matrix * copyMatrix(matrix * mtx);

/* Sets the (row, col) element of mtx to val.  Returns 0 if
 * successful, -1 if mtx is NULL, and -2 if row or col are
 * outside of the dimensions of mtx.
 */
int setElement(matrix * mtx, int row, int col, double val);

/* Sets the reference val to the value of the (row, col) 
 * element of mtx.  Returns 0 if successful, -1 if either 
 * mtx or val is NULL, and -2 if row or col are outside of 
 * the dimensions of mtx.
 */
int getElement(matrix * mtx, int row, int col, 
               double * val);

/* Sets the reference n to the number of rows of mtx.
 * Returns 0 if successful and -1 if mtx or n is NULL.
 */
int nRows(matrix * mtx, int * n);

/* Sets the reference n to the number of columns of mtx.
 * Returns 0 if successful and -1 if mtx is NULL.
 */
int nCols(matrix * mtx, int * n);

/* Prints the matrix to stdout.  Returns 0 if successful 
 * and -1 if mtx is NULL.
 */
int printMatrix(matrix * mtx);

/* Writes the transpose of matrix in into matrix out.  
 * Returns 0 if successful, -1 if either in or out is NULL,
 * and -2 if the dimensions of in and out are incompatible.
 */
int transpose(matrix * in, matrix * out);

/* Writes the sum of matrices mtx1 and mtx2 into matrix 
 * sum. Returns 0 if successful, -1 if any of the matrices 
 * are NULL, and -2 if the dimensions of the matrices are
 * incompatible.
 */
int sum(matrix * mtx1, matrix * mtx2, matrix * sum);

/* Writes the product of matrices mtx1 and mtx2 into matrix
 * prod.  Returns 0 if successful, -1 if any of the 
 * matrices are NULL, and -2 if the dimensions of the 
 * matrices are incompatible.
 */
int product(matrix * mtx1, matrix * mtx2, matrix * prod);

int areEqual(matrix * mtx1, matrix * mtx2);

/* Writes the dot product of vectors v1 and v2 into 
 * reference prod.  Returns 0 if successful, -1 if any of
 * v1, v2, or prod are NULL, -2 if either matrix is not a 
 * vector, and -3 if the vectors are of incompatible 
 * dimensions.
 */
int dotProduct(matrix * v1, matrix * v2, double * prod);

int identity(matrix * m);

int isSquare(matrix * mtx);

int isDiagonal(matrix * mtx);

int isUpperTriangular(matrix * mtx);

int diagonal(matrix * v, matrix * mtx);

int isSymmetric(matrix * mtx);

int printLowerTriangular(matrix * mtx);

int printLowerTriangularCol(matrix * mtx);

matrix * addElementMatrix(int rows, int cols, matrix * mtx);

matrix * vec(matrix * mtx);

double vecGetLowerValue(matrix * vec, int i, int j, int n);

double vecGetUpperValue(matrix * vec, int i, int j, int n);

int vecGetLowerIndex(int i, int j, int n);

int vecGetUpperIndex(int i, int j, int n);

double laplace(matrix *A);

double norma(matrix *A, int p);

double normaInf(matrix *A);

double normaMatricial(matrix *A, int p);

double normaMatricialInf(matrix *A);

double normaMatricial1(matrix *A);

int productByScalar(matrix *A, double x, matrix *prod);

int divisionByScalar(matrix *A, double x, matrix *div);

int subtraction(matrix *A, matrix *B, matrix *sub);

int division(matrix *A, matrix *B, matrix *div);

int metodoDasPotencias(matrix *A, matrix *z0, double tol, int maxit, matrix **_z, double *_lambda);

matrix* matrixRowToVector(matrix *A, int i);

matrix* matrixColToVector(matrix *A, int j);

int jacobi(matrix *A, matrix *e, matrix *V);

void householder(matrix *M, matrix *R, matrix *Q);
