#include "matrix.h"
#include <string.h>

int main(int argc, char *argv[]) {

	if(argc <= 1){
		printf("usage: %s inputFile\n", argv[0]);
		exit(1);
	}

	char input[256];
	strcpy(input, argv[1]);
	
	FILE *fp;

	fp = fopen(input, "r");
	if(fp == NULL){
		printf("Error opening input file %s!\n", input);
		exit(1);
	}

	char line[2048];

	// calculating the number of rows
	int numRows = 0;
	while(fgets(line, sizeof(line), fp))
		numRows++;

	// calculating the number of columns
	fgets(line, sizeof(line), fp);

	int numCols = 0;
	char *p;
	p = strtok(line, " ");
	if(p)
		numCols++;
	do{
		p = strtok((char *)'\0', " ");
		if(p)
			numCols++;
	} while(p);



	matrix *A = newMatrix(numRows, numCols);

	rewind(fp);
	for(int i=1; i<=numRows; i++)
		for(int j=1; j<=numCols; j++){
			double A_mn;
			int ret = fscanf(fp, "%lf\n", &A_mn);
			if(ret <=0){
				printf("break return = %d\n", ret);
				break;
			}
			setElement(A, i, j,   A_mn);
		}

	printf("\nA =\n");
	printMatrix(A);

	printf("\nnumRows = %d\n", numRows);
	printf("numCols = %d\n", numCols);

}
