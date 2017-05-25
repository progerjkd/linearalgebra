LIBS=-lm
LDFLAGS=$(LIBS)
CPPFLAGS="-std=c99"

all: runGauss runGaussP runGaussJordan runLU runCholesky runLaplace

runGauss: matrix.o gauss.o
	gcc runGauss.c matrix.o  gauss.o -o runGauss

runGaussP: matrix.o gaussP.o retroSubstituicao.o
	gcc runGaussP.c matrix.o gaussP.o retroSubstituicao.o -o runGaussP

runGaussJordan: matrix.o gaussJordan.o
	gcc runGaussJordan.c matrix.o gaussJordan.o -o runGaussJordan

runLU: matrix.o decomposicaoLU.o gaussP.o retroSubstituicao.o
	gcc runLU.c matrix.o decomposicaoLU.o gaussP.o retroSubstituicao.o -o runLU

runCholesky: matrix.o decomposicaoCholesky.o gaussP.o retroSubstituicao.o laplace.o
	gcc runCholesky.c matrix.o  decomposicaoCholesky.o gaussP.o  retroSubstituicao.o laplace.o -o runCholesky $(LDFLAGS)

runLaplace: matrix.o laplace.o
	gcc runLaplace.c matrix.o laplace.o -o runLaplace $(LIBS)




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



clean:
	-rm -f *.o runGauss runGaussP runGaussJordan runLU runCholesky runLaplace
