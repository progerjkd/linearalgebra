#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int main(void){
	time_t hora;

	hora = time(NULL);

	printf("Número de segundos antes: %d\n" , hora);
	printf("Dormindo 5 segundos\n");

	sleep(5);

	time(&hora);

	printf("Número de segundos depois: %d\n", hora);


	printf("\n*** PARTE 2 ***\n");

	time_t hora_sist;
	struct tm *hora_exp;

	hora_sist = time(NULL);
	hora_exp = localtime(&hora_sist);

	printf("Usando ctime(): %s", ctime(&hora_sist));
	printf("Usando asctime(): %s", asctime(hora_exp));
}
