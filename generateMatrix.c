#include "matrix.h"
#include <string.h>
#include <stdbool.h>
#include <time.h>


int intRand(int min, int max){

	return(rand() % (max + 1 - min) + min);

}

double doubleRand(int min, int max){

	double range = (max - min);
	double div = RAND_MAX / range;

	return min + (rand() / div);

}

int main(int argc, char *argv[]) {


	if(argc <= 7){
		printf("usage: %s rows cols (integer|double) (symmetric|asymmetric) min max outputFile\n", argv[0]);
		exit(1);
	}

	bool _int = false;
	bool _double = false;
	bool _symmetric = false;
	bool _asymmetric = false;

	int rows = atoi(argv[1]);
	int cols = atoi(argv[2]);
	int min  = atoi(argv[5]);
	int max  = atoi(argv[6]);
	char output[256];
	strcpy(output, argv[7]);
//	string output = argv[6];

	if (strcmp(argv[3], "integer") == 0)
		_int = true;
	else if (strcmp(argv[3], "double") == 0)
		_double = true;
	else {
		printf("usage: %s rows cols (integer|double) (symmetric|asymmetric) min max outputFile\n", argv[0]);
		exit(1);
	}

	if (strcmp(argv[4], "symmetric") == 0)
		_symmetric = true;
	else if (strcmp(argv[4], "asymmetric") == 0)
		_asymmetric = true;
	else {
		printf("usage: %s rows cols (integer|double) (symmetric|asymmetric) min max outputFile\n", argv[0]);
		exit(1);
	}


	//-Wno-pointer-to-int-cast
	srand((unsigned int)**main + (unsigned int)&argc + (unsigned int)time(NULL));
	srand(rand());


	matrix *A = newMatrix(rows, cols);
	FILE *fp;
	fp = fopen(output, "w");

	if(fp == NULL){
		printf("Error writing to output file!\n");
		exit(1);
	}

	for(int i=1; i<=A->rows; i++){
		for(int j=1; j<=A->cols; j++){
			if(_int){
				int intValue = intRand(min,max);
				setElement(A, i, j, intValue);
				fprintf(fp, "%lf", (double)intValue); 
				if(j != A->cols)
					fprintf(fp, " ");
			}
			else{
				double doubleValue = doubleRand(min, max);
				setElement(A, i, j, doubleValue);
				fprintf(fp, "%lf", doubleValue); 
				if(j != A->cols)
					fprintf(fp, " ");
			}
		}
		fprintf(fp, "\n");
	}


	if(_symmetric){

		matrix *triu  = newMatrix(A->rows, A->cols);
		matrix *tril  = newMatrix(A->rows, A->cols);
		matrix *triuT = newMatrix(A->rows, A->cols);
		matrix *tmp1  = newMatrix(A->rows, A->cols);
		matrix *tmp2  = newMatrix(A->rows, A->cols);

		getLowerTriangular(A, tril, -1);
		getUpperTriangular(A, triu, 1);
		transpose(triu, triuT);

		subtraction(A, tril, tmp1);
		sum(tmp1, triuT, tmp2);
		A = tmp2;

		fclose(fp);
		fp = fopen(output, "w");

		if(fp == NULL){
			printf("Error writing to output file!\n");
			exit(1);
		}
		for(int i=1; i<=A->rows; i++){
			for(int j=1; j<=A->cols; j++){
				double A_ij;
				getElement(A, i, j, &A_ij);
				fprintf(fp, "%lf", A_ij); 
				if(j != A->cols)
					fprintf(fp, " ");
				}
			fprintf(fp, "\n");

			}
	}

	printf("\nA =\n");
	printMatrix(A);

	fclose(fp);
}

