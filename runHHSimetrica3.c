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

	loadMatrix(&A, argv[1]);

#ifdef VERBOSE
	printMatrix2(A, "A");
#endif	
	isSymmetric(A) ? printf("\nMatrix A is symmetric\n") : printf("\nMatrix A is Asymmetric\n");

	printf("\nHouseholder transformation: (results in a tridiagonal matrix)\n");
	HH = householder2(A);
	deleteMatrix(A);

#ifdef VERBOSE
	printMatrix2(HH, "householder(A)");
#endif

	printf("\nQR decomposition of HH matrix...\n");

	double epsilon = 0.001;
	matrix *D = qr2(HH, epsilon);


#ifdef VERBOSE
	printMatrix2(HH, "HH");

	printf("\nEigenvalues of matrix A:\n");
	printMatrix2(D, "D");

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
