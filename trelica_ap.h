/* trelica_ap.h
   Cezar - 02/08/2019 */

#ifndef bar_h
#define bar_h

typedef struct
{
 int id;		//id do nó
 //int gls[2];      
 double x[2];  //coordenadas do nó
 bool free[2];    //nó livre para se deslocar ou não - condição de contorno essencial
 double p[2];     //carga concentrada aplicada no nó
 double u[2],u0[2];
} Node;

typedef struct
{
 int no[2];           //vetor do índice do nó - nó inicial, nó final
 bool on;
 double q;            //carga distribuida no elemento
 double length,E,area,c,s;//comprimento, módulo de elasticidade, área da seção transversal, cosseno seno
 double k[4][4];      //matriz de rigidez
 double b[4];         //matriz B
 double epsilon,sigma,sigma0,forca,sigmay;//deformação axial, tensão axial
 int region;          //identifica a regiao no diagrama
} Elem;

#endif
