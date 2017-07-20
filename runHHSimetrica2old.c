#include "matrix.h"
#include <sys/time.h>
#include <time.h>

#define VERBOSE

int main(int argc, char *argv[]) {

	if(argc <= 1){
		printf("usage: %s inputFile\n", argv[0]);
		exit(1);
	}

	struct timeval start, stop;
	double duration;

	gettimeofday(&start, NULL);

	matrix *A, *HH, *Q, *R;

	loadMatrix(&A, argv[1]);
	HH = newMatrix(A->rows, A->cols);
	Q  = newMatrix(A->rows, A->rows);
	R  = newMatrix(A->rows, A->cols);

#ifdef VERBOSE
	printMatrix2(A, "A");
#endif	
	isSymmetric(A) ? printf("\nMatrix A is symmetric\n") : printf("\nMatrix A is Asymmetric\n");

	printf("\nHouseholder transformation: (results in a tridiagonal matrix)\n");
	householder(A, HH);
	//HH = householder2(A);

#ifdef VERBOSE
	printf("\nHH =\n");
	printMatrix(HH);
#endif

	printf("\nQR decomposition of HH matrix...\n");

	for (int i=1; i<=50; i++){
//		printf("qr: %d\n", i);
		qrold(HH, Q, R);
		product(R, Q, HH);
	}

#ifdef VERBOSE
	printf("\nHH =\n");
	printMatrix(HH);

	printf("\nQ =\n");
	printMatrix(Q);

	printf("\nEigenvalues diagonal matrix:\n");
	printf("\nR =\n");
	printMatrix(R);

	matrix *v = newMatrix(R->rows, 1);
	diagonalToVector(R, v);
	printf("\nEigenvalues of matrix A:\n\n");
	printMatrix(v);
#endif

	gettimeofday(&stop, NULL);
	duration = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);

	int hour, min, vt;
	int sec;
	hour = duration / 3600;
	vt = (int)duration % 3600;
	min = vt / 60;
	sec = vt % 60;

	int intpart = (int)duration;
	double decpart = duration - intpart;
	double secd = sec + decpart;
	if(hour > 0)
		printf("\nElapsed time: %ldh%ldm%1.3fs\n", hour, min, secd);
	else
		printf("\nElapsed time: %ldm%1.3fs\n", min, secd);

	memoryUsage();

}
