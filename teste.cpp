#include <iostream>
#include <stdlib.h> 
#include <math.h> 
#include <stdio.h> 
#define uni " UNIJUI-Universidade Regional do Noroeste do Estado do Rio Grande do Sul\n"
#define curso " Licenciatura em Matematica\n"
#define alunos " Fulano de Tal e Beltrano de Tal\n\n" 
#define pi 3.1415
int a,b,c,d,e,f;
float x,y;

void main()
{
	cout<<"Este programa foi elaborado pelo Fulano de Tal \n";
	cout<<"para resolver sistemas da forma \n";
	cout<<" ax + by = c \n";
	cout<<" dx + ey = f \n";

	cout<<"\n Entre com o valor de a ";
	cin>>a;

	cout<<"\n Entre com o valor de b ";
	cin>>b;
	cout<<"\n Entre com o valor de c ";
	cin>>c;
	cout<<"\n Entre com o valor de d ";
	cin>>d;
	cout<<"\n Entre com o valor de e ";
	cin>>e;
	cout<<"\n Entre com o valor de f ";
	cin>>f;
	y=(f*a-d*c)/(a*e-d*b);
	x=(c-b*y)/a ;
	cout<<"\n y vale " <<y;
	cout<<"\n x vale " <<x;
	getchar(); 
}
