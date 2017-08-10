LIBS=-lm
LDFLAGS=$(LIBS)
CPPFLAGS="-std=c99"

all: runGauss runGaussP runGaussJordan runLU runCholesky runLaplace runMetodoDasPotencias runJacobi runQR runHouseholder ascii runHHSimetrica runHHSimetrica2 runHHSimetrica3 runHHAssimetrica runHHAssimetrica2 runHHAssimetrica3 runHHSimetrica2old generateMatrix readMatrix runNorma storeByDiagonalGaxpy symmetricStorageGaxpy runMatrix

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

runJacobi: jacobi.o matrix.o jacobi.o
	gcc runJacobi.c matrix.o jacobi.o -o runJacobi $(LDFLAGS) $(CPPFLAGS)

runQR: decomposicaoQRold.o matrix.o
	gcc runQR.c matrix.o decomposicaoQRold.o -o runQR $(LDFLAGS) $(CPPFLAGS)

runNorma: matrix.o
	gcc runNorma.c matrix.o -o runNorma $(LDFLAGS) $(CPPFLAGS)

runHouseholder: householder.o
	gcc runHouseholder.c matrix.o householder.o -o runHouseholder $(LDFLAGS) $(CPPFLAGS)

runHHSimetrica: householder.o decomposicaoQRold.o loadMatrix.o
	gcc runHHSimetrica.c matrix.o householder.o decomposicaoQRold.o loadMatrix.o -o runHHSimetrica $(LDFLAGS) $(CPPFLAGS)

runHHSimetrica2: householder.o decomposicaoQR.o loadMatrix.o decomposicaoQR2.o jacobi.o
	gcc runHHSimetrica2.c matrix.o householder.o decomposicaoQR.o decomposicaoQR2.o jacobi.o loadMatrix.o -o runHHSimetrica2 $(LDFLAGS) $(CPPFLAGS)

runHHSimetrica2old: householder.o decomposicaoQR.o loadMatrix.o decomposicaoQRold.o jacobi.o
	gcc runHHSimetrica2old.c matrix.o householder.o decomposicaoQR.o decomposicaoQRold.o jacobi.o loadMatrix.o -o runHHSimetrica2old $(LDFLAGS) $(CPPFLAGS)

runHHSimetrica3: householder.o loadMatrix.o decomposicaoQR.o decomposicaoQR2.o jacobi.o matrix.o
	gcc runHHSimetrica3.c matrix.o householder.o decomposicaoQR.o decomposicaoQR2.o jacobi.o loadMatrix.o -o runHHSimetrica3 $(LDFLAGS) $(CPPFLAGS)

runHHAssimetrica: householder.o decomposicaoQRold.o laplace.o 
	gcc runHHAssimetrica.c matrix.o householder.o laplace.o decomposicaoQRold.o -o runHHAssimetrica $(LDFLAGS) $(CPPFLAGS)

runHHAssimetrica2: householder.o decomposicaoQR.o laplace.o loadMatrix.o
	gcc runHHAssimetrica2.c matrix.o householder.o laplace.o decomposicaoQR.o loadMatrix.o -o runHHAssimetrica2 $(LDFLAGS) $(CPPFLAGS)

runHHAssimetrica3: householder.o decomposicaoQR.o decomposicaoQR2.o laplace.o loadMatrix.o jacobi.o matrix.o
	gcc runHHAssimetrica3.c householder.o decomposicaoQR.o decomposicaoQR2.o laplace.o loadMatrix.o jacobi.o matrix.o -o runHHAssimetrica3 $(LDFLAGS) $(CPPFLAGS)

generateMatrix: matrix.o
	gcc generateMatrix.c  matrix.o $(LDFLAGS) $(CPPFLAGS) -Wno-pointer-to-int-cast -o generateMatrix

readMatrix: matrix.o
	gcc readMatrix.c  matrix.o $(LDFLAGS) $(CPPFLAGS) -Wno-pointer-to-int-cast -o readMatrix 

runMatrix: matrix.o
	gcc runMatrix.c  matrix.o $(LDFLAGS) $(CPPFLAGS) -Wno-pointer-to-int-cast -o runMatrix 

ascii:
	gcc ascii.c -o ascii




matrix.o:
	gcc -c matrix.c $(CPPFLAGS)

storeByDiagonalGaxpy: matrix.o
	gcc storeByDiagonalGaxpy.c matrix.o -o storeByDiagonalGaxpy $(LDFLAGS) $(CPPFLAGS)

symmetricStorageGaxpy: matrix.o
	gcc symmetricStorageGaxpy.c matrix.o -o symmetricStorageGaxpy $(LDFLAGS) $(CPPFLAGS)

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

decomposicaoQRold.o:
	gcc -c decomposicaoQRold.c $(CPPFLAGS)

decomposicaoQR2.o:
	gcc -c decomposicaoQR2.c $(CPPFLAGS)

householder.o:
	gcc -c householder.c $(CPPFLAGS)

loadMatrix.o:
	gcc -c loadMatrix.c $(CPPFLAGS)

clean:
	-rm -f *.o runGauss runGaussP runGaussJordan runLU runCholesky runLaplace runMetodoDaPotencia runJacobi runQR runHouseholder ascii runHHSimetrica runHHAssimetrica generateMatrix readMatrix runHHSimetrica2 runHHAssimetrica2 runHHAssimetrica3 runHHSimetrica2old runHHSimetrica3 runNorma storeByDiagonalGaxpy symmetricStorageGaxpy runMatrix
