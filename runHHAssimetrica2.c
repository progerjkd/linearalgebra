#include "matrix.h"
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define VERBOSE

int main(int argc, char *argv[]) {

	if(argc <= 1){
		printf("usage: %s inputFile\n", argv[0]);
		exit(1);
	}

	struct timeval start, stop;
	double duration;

	gettimeofday(&start, NULL);

	matrix *A, *Hess, *Q, *R;
	loadMatrix(&A, argv[1]);

	Hess = newMatrix(A->rows, A->cols);
	Q  = newMatrix(A->rows, A->rows);
	R  = newMatrix(A->rows, A->cols);

#ifdef VERBOSE
	printf("\nA =\n");
	printMatrix(A);
#endif
	isSymmetric(A) ? printf("\nMatrix A is symmetric\n") : printf("\nMatrix A is Asymmetric\n");

	printf("\nHouseholder transformation: (results in a Hessenberg Matrix)\n");
	householder(A, Hess);

#ifdef VERBOSE
	printf("\nHess =\n");
	printMatrix(Hess);
#endif

	printf("\nQR decomposition of Hess matrix...\n");

	for (int i=1; i<=50; i++){
		qr(Hess, Q, R);
		product(R, Q, Hess);
	}

//	printf("\nMatriz com dentes (ou não) abaixo da diagonal:\n");
#ifdef VERBOSE
	printf("\nHess =\n");
	printMatrix(Hess);
#endif

	double epsilon = 0.001;
	int n = Hess->rows;
	int i = 1;

	printf("\nEigenvalues of matrix A:\n\n");
	while(i <= n){
		double Hesstmp;
		getElement(Hess, i+1, i, &Hesstmp);
		if((i<n) && (fabs(Hesstmp) > epsilon)){
			// seleciona a matriz 2x2 formada com os dentes e calcula os autovalores
			matrix *M2 = newMatrix(2, 2);
			subMatrix(Hess, i, i+1, i, i+1, M2);

			double M2_11, M2_22;
			getElement(M2, 1, 1, &M2_11);
			getElement(M2, 2, 2, &M2_22);


			_eqSegGrau x;
			x = eqSegGrau(1, -(M2_11 + M2_22), laplace(M2));
#ifdef VERBOSE
			printf("%10.4f %s%fi\n", creal(x.x1), (cimag(x.x1) >= 0) ? "+" : "", cimag(x.x1));
			printf("%10.4f %s%fi\n", creal(x.x2), (cimag(x.x2) >= 0) ? "+" : "", cimag(x.x2));
#endif
			i++;
		} else {
			double Hess_ii;
			getElement(Hess, i, i, &Hess_ii);
#ifdef VERBOSE
			printf("%10.4f\n", Hess_ii);
#endif
		}
		i++;
	}

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
	


/*	printf("\nR =\n");
	printMatrix(R);

	printf("\nMatrix diagonal de autovalores:\n");
	printf("Q =\n");
	printMatrix(Q);
*/
}
