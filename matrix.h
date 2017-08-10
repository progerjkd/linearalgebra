#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>


typedef struct {
  int rows;
  int cols;
  double * data;
} matrix;

#define ELEM(mtx, row, col) \
	  mtx->data[(col-1) * mtx->rows + (row-1)]


extern bool __DELETE_norma__;
extern bool __DELETE_matrixColToVector__;
extern bool __DELETE_dotProduct2__;
extern bool __DELETE_product2__;
extern bool __DELETE_subtraction2__;
extern bool __DELETE_productByScalar2__;
extern bool __DELETE_matrixAbs2__;
extern bool __DELETE_getMinElemVec__;
extern bool __DELETE_qr__;
extern bool __DELETE_sum2__;

//bool __MATRIX_DELETE__ = false;

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

void subMatrix(matrix *A, int startRow, int endRow, int startCol, int endCol, matrix *B);

int subtractionByScalar(matrix *A, double x, matrix *sub);

int identity(matrix * m);

int isSquare(matrix * mtx);

int isDiagonal(matrix * mtx);

int isUpperTriangular(matrix * mtx);

int vectorToDiagonal(matrix * v, matrix * mtx);

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

void metodoDasPotencias(matrix *A, matrix *z0, double tol, int maxit, matrix **_z, double *_lambda);

matrix* matrixRowToVector(matrix *A, int i);

matrix* matrixColToVector(matrix *A, int j);

int jacobi(matrix *A, matrix *e, matrix *V);

typedef struct {
	matrix *Q;
	matrix *R;
} _qr;

_qr qr(matrix *A);
void qrold(matrix *M, matrix *Q, matrix *R);
//void qr(matrix *M, matrix *Q, matrix *R);

void copyColumn(matrix *A, matrix *c, int j);

void copyNColumn(matrix *A, matrix *c, int n, int j);

int tolerance(matrix *A, double tol);

void householder(matrix *A, matrix *HH);

void diagonalToVector(matrix *A, matrix *v);

typedef struct {
	double complex x1;
	double complex x2;
} _eqSegGrau;

_eqSegGrau eqSegGrau(double a, double b, double c);

typedef struct {
	double min;
	int pos;
} _minElemVec;

_minElemVec getMinElemVec(matrix *A, int col);

void matrixAbs(matrix *A, matrix *B);

void loadMatrix(matrix **A, char *input);

int getLowerTriangular(matrix *A, matrix *B, int k);

int getUpperTriangular(matrix *A, matrix *B, int k);

matrix *transpose2(matrix * in);

matrix *sum2(matrix * mtx1, matrix * mtx2);

matrix *product2(matrix *A, matrix *B);

double dotProduct2(matrix * v1, matrix * v2);

matrix *identity2(int order);

double getElement2(matrix * mtx, int row, int col);

matrix *subtraction2(matrix *A, matrix *B);

matrix *division2(matrix *A, matrix *B);

matrix *sumByScalar2(matrix *A, double x);

matrix *subtractionByScalar2(matrix *A, double x);

matrix *productByScalar2(matrix *A, double x);

matrix *divisionByScalar2(matrix *A, double x);

matrix *getUpperTriangular2(matrix *A, int k);

matrix *getLowerTriangular2(matrix *A, int k);

int matrixSwapLines(matrix *matrix, int n, int x, int y);

int matrixSwapLines2(matrix *A, int x, int y);

int menorPrincipal(matrix *A, matrix *B, int n);

matrix *menorPrincipal2(matrix *A, int n);

matrix *matrixAbs2(matrix *A);

matrix *subMatrix2(matrix *A, int startRow, int endRow, int startCol, int endCol);

matrix *diagonalToVector2(matrix *A);

matrix *householder2(matrix *A);

int printMatrix2(matrix * mtx, char *name);

matrix *qr2(matrix *A, double epsilon);

void memoryUsage(void);

void gauss(matrix *A, matrix *x, matrix *b, int n);

void gaussP(matrix *A, matrix *x, matrix *b, int n);

void gaussJordan(matrix *A, matrix *x, matrix *b, int n);

int diagonalContainsZero(matrix * mtx);

int retroSubstituicao(matrix *A, matrix *x, matrix *b, int n);

void decomposicaoLU(matrix *A, int n, matrix *L, matrix *U);

int decomposicaoCholesky(matrix *A, matrix *S,  matrix *St, int n);


#endif
