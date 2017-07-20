#include "matrix.h"

matrix *transpose2(matrix * in);
matrix *sum2(matrix * mtx1, matrix * mtx2);
matrix *product2(matrix * mtx1, matrix * mtx2);
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
int matrixSwapLines2(matrix *A, int x, int y);
matrix *menorPrincipal2(matrix *A, int n);
matrix *matrixAbs2(matrix *A);
matrix *subMatrix2(matrix *A, int startRow, int endRow, int startCol, int endCol);

matrix *diagonalToVector2(matrix *A);



int main() {

	matrix *A, *B, *x, *b, *c, *d;
	int n = 4;

	A = newMatrix(n, n);
	B = newMatrix(n, n);
	x = newMatrix(n, 1);
	b = newMatrix(n, 1);
	c = newMatrix(n, 1);
	d = newMatrix(1, n);

	setElement(A, 1, 1,   1);
	setElement(A, 1, 2,   2);
	setElement(A, 1, 3,  -1);
	setElement(A, 1, 4,   0);
	setElement(A, 2, 1,   0);
	setElement(A, 2, 2,  -1);
	setElement(A, 2, 3,   1);
	setElement(A, 2, 4,  -1);
	setElement(A, 3, 1,  -2);
	setElement(A, 3, 2,  -1);
	setElement(A, 3, 3,   4);
	setElement(A, 3, 4,   2);
	setElement(A, 4, 1,   4);
	setElement(A, 4, 2,   3);
	setElement(A, 4, 3,   0);
	setElement(A, 4, 4,   1);



	setElement(B, 1, 1,   1);
	setElement(B, 1, 2,   2);
	setElement(B, 1, 3,  -1);
	setElement(B, 1, 4,   0);
	setElement(B, 2, 1,   0);
	setElement(B, 2, 2,  -1);
	setElement(B, 2, 3,   1);
	setElement(B, 2, 4,  -5);
	setElement(B, 3, 1,  -4);
	setElement(B, 3, 2,  -5);
	setElement(B, 3, 3,   6);
	setElement(B, 3, 4,   7);
	setElement(B, 4, 1,   8);
	setElement(B, 4, 2,   9);
	setElement(B, 4, 3,  10);
	setElement(B, 4, 4,  21);

	setElement(b, 1, 1,  -4);
	setElement(b, 2, 1,   0);
	setElement(b, 3, 1,   7);
	setElement(b, 4, 1, -10);

	setElement(c, 1, 1,  -4);
	setElement(c, 2, 1,   0);
	setElement(c, 3, 1,   7);
	setElement(c, 4, 1, -10);

	setElement(d, 1, 1,  -4);
	setElement(d, 1, 2,   0);
	setElement(d, 1, 3,   7);
	setElement(d, 1, 4, -10);

	printMatrix2(A, "A");
//	matrix **At;
  //      At = transpose2(A);
	printMatrix2(transpose2(A), "At");

	printMatrix2(sum2(b, c), "sum(b, c)");
	printMatrix2(c, "c");
	printMatrix2(d, "d");
	printMatrix2(product2(c, d), "product2(c, d)");

	printf("\ndotProduct(b, c) = %f\n", dotProduct2(b, c));
	printMatrix2(identity2(5), "identity2(5)");
	printf("\ngetElement2(A, 4, 2) = %f\n", getElement2(A, 4, 2));
	double A_31 = getElement2(A, 3, 1);
	printf("\ngetElement2(A, 3, 1) = %f\n", A_31);

	printMatrix2(subtraction2(A, B), "subtraction2(A, B)");
	printMatrix2(division2(A, B), "division2(A, B)");
	printMatrix2(sumByScalar2(A, 10), "sumByScalar2(A, 10)");

	printMatrix2(A, "A");
	printMatrix2(getUpperTriangular2(A, 0), "getUpperTriangular2(A, 0)");
	printMatrix2(getUpperTriangular2(A, 1), "getUpperTriangular2(A, 1)");
	printMatrix2(getUpperTriangular2(A, -1), "getUpperTriangular2(A, -1)");

	printMatrix2(getLowerTriangular2(A, 0), "getLowerTriangular2(A, 0)");
	printMatrix2(getLowerTriangular2(A, 1), "getLowerTriangular2(A, 1)");
	printMatrix2(getLowerTriangular2(A, -1), "getLowerTriangular2(A, -1)");

	printMatrix2(A, "A");
	matrixSwapLines2(A, 2, 3);
	printMatrix2(A, "matrixSwapLines2(A, 2, 3)");

	printMatrix2(menorPrincipal2(A, 1), "menorPrincipal2(A, 1)");
	printMatrix2(menorPrincipal2(A, 2), "menorPrincipal2(A, 2)");
	printMatrix2(menorPrincipal2(A, 3), "menorPrincipal2(A, 3)");
	printMatrix2(menorPrincipal2(A, 4), "menorPrincipal2(A, 4)");

	printf("Norma 2 de norma(menorPrincipal2(A, 4), 2): %f\n", norma(menorPrincipal2(A, 4), 2));

	printMatrix2(A, "A");
	printMatrix2(matrixAbs2(A), "matrixAbs2(A)");
	printMatrix2(matrixRowToVector(A, 2), "matrixRowToVector(A, 2)");
	printMatrix2(matrixColToVector(A, 3), "matrixColToVector(A, 3)");
	printMatrix2(A, "A");
	printMatrix2(subMatrix2(A, 2, 4, 2, 3), "subMatrix2(A, 2, 4, 2, 3)");

	printMatrix2(diagonalToVector2(A), "diagonalToVector2(A)");




}
