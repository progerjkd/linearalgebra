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

	matrix *A, *HH;
	_qr QR;

	loadMatrix(&A, argv[1]);
	HH = newMatrix(A->rows, A->cols);

#ifdef VERBOSE
	printMatrix2(A, "A");
#endif	
	isSymmetric(A) ? printf("\nMatrix A is symmetric\n") : printf("\nMatrix A is Asymmetric\n");

	printf("\nHouseholder transformation: (results in a tridiagonal matrix)\n");
	//householder(A, HH);
	HH = householder2(A); // TODO

#ifdef VERBOSE
	printMatrix2(HH, "HH");
#endif

	printf("\nQR decomposition of HH matrix...\n");

	__DELETE_product2__ = true;

	for (int i=1; i<=50; i++){
//		printf("i: %d\n", i);
		QR = qr(HH);
		deleteMatrix(HH);

		if(i<50)
			HH = product2(QR.R, QR.Q);
		else{
			__DELETE_product2__ = false;
			HH = product2(QR.R, QR.Q);
		}
	}

#ifdef VERBOSE
	printMatrix2(HH, "HH");

	printMatrix2(QR.Q, "Q");

	printf("\nEigenvalues diagonal matrix:\n");
	printMatrix2(QR.R, "R");

	matrix *v = diagonalToVector2(QR.R);
	printf("\nEigenvalues of matrix A:\n\n");
	printMatrix2(v, "diag(R)");
#endif

	gettimeofday(&stop, NULL);
	duration = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);

	long int hour, min, vt;
	long int sec;
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
