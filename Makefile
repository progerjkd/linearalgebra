LIBS=-lm
LDFLAGS=$(LIBS)
CPPFLAGS="-std=c99"

all: runGauss runGaussP runGaussJordan runLU runCholesky runLaplace runMetodoDasPotencias runJacobi runQR runHouseholder ascii runHHSimetrica runHHSimetrica2 runHHAssimetrica runHHAssimetrica2 generateMatrix readMatrix

runGauss: matrix.o gauss.o
	gcc runGauss.c matrix.o  gauss.o -o runGauss $(LDFLAGS)

runGaussP: matrix.o gaussP.o retroSubstituicao.o
	gcc runGaussP.c matrix.o gaussP.o retroSubstituicao.o -o runGaussP $(LDFLAGS)

runGaussJordan: matrix.o gaussJordan.o
	gcc runGaussJordan.c matrix.o gaussJordan.o -o runGaussJordan $(LDFLAGS)

runLU: matrix.o decomposicaoLU.o gaussP.o retroSubstituicao.o
	gcc runLU.c matrix.o decomposicaoLU.o gaussP.o retroSubstituicao.o -o runLU $(LDFLAGS)

runCholesky: matrix.o decomposicaoCholesky.o gaussP.o retroSubstituicao.o laplace.o
	gcc runCholesky.c matrix.o  decomposicaoCholesky.o gaussP.o  retroSubstituicao.o laplace.o -o runCholesky $(LDFLAGS)

runLaplace: matrix.o laplace.o
	gcc runLaplace.c matrix.o laplace.o -o runLaplace $(LDFLAGS)

runMetodoDasPotencias: metodoDasPotencias.o
	gcc runMetodoDasPotencias.c matrix.o metodoDasPotencias.o -o runMetodoDasPotencias $(LDFLAGS)

runJacobi: jacobi.o
	gcc runJacobi.c matrix.o jacobi.o -o runJacobi $(LDFLAGS) $(CPPFLAGS)

runQR: decomposicaoQR.o
	gcc runQR.c matrix.o decomposicaoQR.o -o runQR $(LDFLAGS) $(CPPFLAGS)

runHouseholder: householder.o
	gcc runHouseholder.c matrix.o householder.o -o runHouseholder $(LDFLAGS) $(CPPFLAGS)

runHHSimetrica: householder.o decomposicaoQR.o loadMatrix.o
	gcc runHHSimetrica.c matrix.o householder.o decomposicaoQR.o loadMatrix.o -o runHHSimetrica $(LDFLAGS) $(CPPFLAGS)

runHHSimetrica2: householder.o decomposicaoQR.o loadMatrix.o decomposicaoQR2.o jacobi.o
	gcc runHHSimetrica2.c matrix.o householder.o decomposicaoQR.o decomposicaoQR2.o jacobi.o loadMatrix.o -o runHHSimetrica2 $(LDFLAGS) $(CPPFLAGS)

runHHAssimetrica: householder.o decomposicaoQR.o laplace.o
	gcc runHHAssimetrica.c matrix.o householder.o laplace.o decomposicaoQR.o -o runHHAssimetrica $(LDFLAGS) $(CPPFLAGS)

runHHAssimetrica2: householder.o decomposicaoQR.o laplace.o loadMatrix.o
	gcc runHHAssimetrica2.c matrix.o householder.o laplace.o decomposicaoQR.o loadMatrix.o -o runHHAssimetrica2 $(LDFLAGS) $(CPPFLAGS)

generateMatrix: 
	gcc generateMatrix.c  matrix.o $(LDFLAGS) $(CPPFLAGS) -Wno-pointer-to-int-cast -o generateMatrix

readMatrix: 
	gcc readMatrix.c  matrix.o $(LDFLAGS) $(CPPFLAGS) -Wno-pointer-to-int-cast -o readMatrix 

ascii:
	gcc ascii.c -o ascii




matrix.o:
	gcc -c matrix.c $(CPPFLAGS)

storeByDiagonal:
	gcc storeByDiagonal.c -o storeByDiagonal

symmetricStorageGaxpy:
	gcc symmetricStorageGaxpy.c -o symmetricStorageGaxpy

gaussP.o:
	gcc -c gaussP.c

gauss.o:
	gcc -c gauss.c

retroSubstituicao.o:
	gcc -c retroSubstituicao.c

gaussJordan.o:
	gcc -c gaussJordan.c

decomposicaoLU.o:
	gcc -c decomposicaoLU.c

decomposicaoCholesky.o:
	gcc -c decomposicaoCholesky.c

laplace.o:
	gcc -c laplace.c

metodoDasPotencias.o:
	gcc -c metodoDasPotencias.c $(CPPFLAGS)

jacobi.o:
	gcc -c jacobi.c $(CPPFLAGS)

decomposicaoQR.o:
	gcc -c decomposicaoQR.c $(CPPFLAGS)

decomposicaoQR2.o:
	gcc -c decomposicaoQR2.c $(CPPFLAGS)

householder.o:
	gcc -c householder.c $(CPPFLAGS)

loadMatrix.o:
	gcc -c loadMatrix.c $(CPPFLAGS)

clean:
	-rm -f *.o runGauss runGaussP runGaussJordan runLU runCholesky runLaplace runMetodoDaPotencia runJacobi runQR runHouseholder ascii runHHSimetrica runHHAssimetrica generateMatrix readMatrix runHHSimetrica2 runHHAssimetrica2
