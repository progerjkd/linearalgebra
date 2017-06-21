#include "matrix.h"
#include <math.h>

#define EPS 2.2204E-16
//#define VERBOSE

#ifdef I
#undef I
#endif

double summ(matrix *v, int Pow){

	double _sum = 0;

	for(int i=1; i<=v->rows; i++){
		double v_i;
		getElement(v, i, 1, &v_i);
		_sum += pow(v_i, Pow);
	}

	return _sum;
}

double ss(matrix *A, int jj){

	double _pow;

	matrix *tmp = newMatrix(A->rows - jj, 1);
	copyNColumn(A, tmp, jj+1, jj);
	_pow = summ(tmp, 2);

	return sqrt(_pow);
}

int sign(matrix *A, int i, int j){
	
	double A_ij;
	getElement(A, i, j, &A_ij);
	if(A_ij > 0)
		return 1;
	else if (A_ij < 0)
		return -1;
	else
		return 0;
}

void householder(matrix *A, matrix *HH){

#ifdef VERBOSE
	printf("\nEntrada - Transformação de Householder\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nHH =\n");
	printMatrix(HH);
#endif

	if(A->rows != A->cols){
		printf("A matriz A não é quadrada. Abortando a transformação de Householder...\n");
		return;
	}


	int length = A->rows;
	matrix *v = newMatrix(length, 1);
	matrix *I = newMatrix(length, length);
	identity(I);
	matrix *Aold = copyMatrix(A);

	matrix *Anew;



	for(int jj=1; jj<=length-2; jj++){
		double S;


		// v(1:jj) = 0;
		for(int v_j=1; v_j<=jj; v_j++)
			setElement(v, v_j, 1, 0);

		S = ss(Aold, jj);

		// v(jj+1) = sqrt(.5*(1+abs(Aold(jj+1,jj))/(S+2*EPS)));
		double tmp1, tmp2;
		getElement(Aold, jj+1, jj, &tmp1);
		tmp2 = sqrt(.5 * (1+fabs(tmp1) / (S+2*EPS)));
		setElement(v, jj+1, 1, tmp2);


		// v(jj+2:length) = Aold(jj+2:length,jj)*sign(Aold(jj+1,jj)) /(2*v(jj+1)*S+2*EPS);
		matrix *v1;
		v1 = newMatrix(Aold->rows - jj + 1, 1);
		copyNColumn(Aold, v1, jj+2, jj);

		for(int i=jj+2; i<=length; i++){
			double Atmp, vtmp, tmp;
			getElement(Aold, i, jj, &Atmp);

			getElement(v, jj+1, 1, &vtmp);

			setElement(v, i, 1, Atmp * sign(Aold, jj+1, jj) / (2 * vtmp * S + 2 * EPS));
		}

		
		// P = I-2*v*v';
		matrix *tmp3;
		tmp3 = newMatrix(v->rows, v->cols);
		productByScalar(v, 2, tmp3);

		matrix *vt, *tmp4;
		vt = newMatrix(v->cols, v->rows);
		transpose(v, vt);
		tmp4 = newMatrix(v->rows, vt->cols);
		product(tmp3, vt, tmp4);

		matrix *P;
		P = newMatrix(I->rows, I->cols);
		subtraction(I, tmp4, P);


		//Anew = P*Aold*P
		matrix *tmp5;
		tmp5 = newMatrix(P->rows, Aold->cols);
		product(P, Aold, tmp5);

		Anew = newMatrix(tmp5->rows, P->cols);
		product(tmp5, P, Anew);

		//Aold = Anew;
		Aold = copyMatrix(Anew);
	}

	// Anew(abs(Anew(:))<5e-14)=0; % Tolerence.
	tolerance(Anew, 5E-14);

	//HH = copyMatrix(Anew);
	*HH = *Anew;

#ifdef VERBOSE
	printf("\nSaída - Transformação de Householder\n");

	printf("\nA =\n");
	printMatrix(A);

	printf("\nHH =\n");
	printMatrix(HH);
#endif
}
