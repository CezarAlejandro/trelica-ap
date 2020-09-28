/* trelica_ap.h
   Cezar - 02/08/2019 */

#ifndef bar_h
#define bar_h

typedef struct
{
 int id;		//id do n�
 //int gls[2];      
 double x[2];  //coordenadas do n�
 bool free[2];    //n� livre para se deslocar ou n�o - condi��o de contorno essencial
 double p[2];     //carga concentrada aplicada no n�
 double u[2],u0[2];
} Node;

typedef struct
{
 int no[2];           //vetor do �ndice do n� - n� inicial, n� final
 bool on;
 double q;            //carga distribuida no elemento
 double length,E,area,c,s;//comprimento, m�dulo de elasticidade, �rea da se��o transversal, cosseno seno
 double k[4][4];      //matriz de rigidez
 double b[4];         //matriz B
 double epsilon,sigma,sigma0,forca,sigmay;//deforma��o axial, tens�o axial
 int region;          //identifica a regiao no diagrama
} Elem;

#endif
