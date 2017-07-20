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

	matrix *A, *Hess;
	_qr QR;

	loadMatrix(&A, argv[1]);

//	Hess = newMatrix(A->rows, A->cols);

#ifdef VERBOSE
	printMatrix2(A, "A");
#endif
	isSymmetric(A) ? printf("\nMatrix A is symmetric\n") : printf("\nMatrix A is Asymmetric\n");

	printf("\nHouseholder transformation: (results in a Hessenberg Matrix)\n");
	Hess = householder2(A);
	deleteMatrix(A);

#ifdef VERBOSE
	printMatrix2(Hess, "householder(A)");
#endif

	printf("\nQR decomposition of Hess matrix...\n");

        for (int i=1; i<=50; i++){
//              printf("i: %d\n", i);
		QR = qr(Hess);
		product(QR.R, QR.Q, Hess);
}


//	printf("\nMatriz com dentes (ou não) abaixo da diagonal:\n");
#ifdef VERBOSE
	printMatrix2(Hess, "Hess");
#endif
        double epsilon = 0.001;
        int n = Hess->rows;
        int i = 1;

        printf("\nEigenvalues of matrix A:\n\n");
        while(i <= n){
                if((i<n) && (fabs(getElement2(Hess, i+1, i)) > epsilon)){
                        // seleciona a matriz 2x2 formada com os dentes e calcula os autovalores
                        matrix *M2 = subMatrix2(Hess, i, i+1, i, i+1);

                        _eqSegGrau x;
                        x = eqSegGrau(1, -(getElement2(M2, 1, 1) + getElement2(M2, 2, 2)), laplace(M2));
#ifdef VERBOSE
                        printf("%10.4f %s%fi\n", creal(x.x1), (cimag(x.x1) >= 0) ? "+" : "", cimag(x.x1));
                        printf("%10.4f %s%fi\n", creal(x.x2), (cimag(x.x2) >= 0) ? "+" : "", cimag(x.x2));
#endif
                        i++;
                } else {
#ifdef VERBOSE
                        printf("%10.4f\n", getElement2(Hess, i, i));
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
	
	memoryUsage();


}
