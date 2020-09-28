#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "trelica_ap.h"

static int col,lin,ne,nn,nr=0,ii,m,nf,qq,g,h,n,flag=0,iso;
static int gl=2;
static double x,lambda,pi=acos(-1.0),alpha,pp;
static double E[3];
//static double mic[4][500];
Elem *elem;                 //ponteiro para vetor de elementos    
FILE *fp=NULL; 
Node *node;                  //ponteiro para vetor de nós
Node* noderesult;            //ponteiro para nó com resultados
static int cont=1;
static int **mic;
static double** kg = NULL;
static double* fg = NULL;
static double* u = NULL;
static int nne=2;
static void L(int id)
{
 //calcula o comprimento do elemento id a partir das coordenadas nodais do elemento
 double dx=node[ elem[id].no[1] ].x[0]-node[ elem[id].no[0] ].x[0];
 double dy=node[ elem[id].no[1] ].x[1]-node[ elem[id].no[0] ].x[1];
 elem[id].length=sqrt(dx*dx+dy*dy); // transformar para milímetros
 elem[id].c = dx/elem[id].length;
 elem[id].s = dy/elem[id].length;
}

static void B(int id)
{
 //define a matriz B do elemento id
//monta a matriz B
 elem[id].b[0]=-elem[id].c/elem[id].length;
 elem[id].b[1]=-elem[id].s/elem[id].length;
 elem[id].b[2]=elem[id].c/elem[id].length;
 elem[id].b[3]=elem[id].s/elem[id].length;
}

static void K(int id)
	{
  double c2, s2, cs;
 	
  // Monta a matriz de rigidez K de cada elemento
  c2 = elem[id].c*elem[id].c;
  s2 = elem[id].s*elem[id].s;
  cs = elem[id].c*elem[id].s;
 	elem[id].k[0][0]=c2;  elem[id].k[0][1]=cs;  elem[id].k[0][2]=-c2; elem[id].k[0][3]=-cs;
	 elem[id].k[1][0]=cs;  elem[id].k[1][1]=s2;	 elem[id].k[1][2]=-cs;	elem[id].k[1][3]=-s2;
	  elem[id].k[2][0]=-c2;	elem[id].k[2][1]=-cs;	elem[id].k[2][2]=c2;		elem[id].k[2][3]=cs;		
		elem[id].k[3][0]=-cs;	elem[id].k[3][1]=-s2;	elem[id].k[3][2]=cs;		elem[id].k[3][3]=s2;
	 
	 double aux= elem[id].E * elem[id].area / elem[id].length;
		for(int lin=0;lin<4;lin++)
		{
			for(int col=0;col<4;col++)
			{
				elem[id].k[lin][col]*=aux;
			}
		}

	}

static void Discretizacao()
{			//montagem da treliça
 int i;
 //Problema Paper
#if 1
 ne = 10;
 nn = 6;
 elem = (Elem*)malloc(ne * sizeof(Elem));
 node = (Node*)malloc(nn * sizeof(Node));
 E[0] = 5E8;
 E[1] = 2.5E8;
 E[2] = 0;     // terceira região
 //definicao do numero de nós
 node[0].x[0] = -1.000; //mm
 node[0].x[1] = 1.000;
 node[0].id = 1;
 node[1].x[0] = 0; //mm
 node[1].x[1] = 1.000;
 node[1].id = 2;
 node[2].x[0] = 1.000; //mm
 node[2].x[1] = 1.000;
 node[2].id = 3;
 node[3].x[0] = -1.000; //mm
 node[3].x[1] = 0;
 node[3].id = 4;
 node[4].x[0] = 0; //mm
 node[4].x[1] = 0;
 node[4].id = 5;
 node[5].x[0] = 1.000; //mm
 node[5].x[1] = 0;
 node[5].id = 6;
 //area otimizada     mm^2
 elem[0].area = 999.9;
 elem[1].area = 983.32;
 elem[2].area = 999.03;
 elem[3].area = 996.86;
 elem[4].area = 101.76;
 elem[5].area = 999.99;
 elem[6].area = 999.13;
 elem[7].area = 100.01;
 elem[8].area = 1011.1;
 elem[9].area = 1004.7;
 //elementos
 elem[0].no[0] = 0;
 elem[0].no[1] = 1;
 elem[1].no[0] = 1;
 elem[1].no[1] = 2;
 elem[2].no[0] = 3;
 elem[2].no[1] = 4;
 elem[3].no[0] = 4;
 elem[3].no[1] = 5;
 elem[4].no[0] = 3;
 elem[4].no[1] = 1;
 elem[5].no[0] = 1;
 elem[5].no[1] = 5;
 elem[6].no[0] = 0;
 elem[6].no[1] = 4;
 elem[7].no[0] = 4;
 elem[7].no[1] = 2;
 elem[8].no[0] = 4;
 elem[8].no[1] = 1;
 elem[9].no[0] = 5;
 elem[9].no[1] = 2;
 for (i = 0; i < nn; i++)
  node[i].u0[0] = node[i].u0[1] = 0.0;
 for (i = 0; i < ne; i++)
 {
	 //elem[i].area *= 1E-6;  //m^2
  elem[i].area = 0.05;
	 elem[i].region = 0;
	 elem[i].E = E[elem[i].region];
	 elem[i].sigmay = 260;  //Pa
	 elem[i].sigma0 = 0.0;
 }

#endif

#if 0 
 //Problema 3 barras
 ne = 3;
 nn = 3;
 elem = (Elem*)malloc(ne * sizeof(Elem));
 node = (Node*)malloc(nn * sizeof(Node));
 E[0] = 200E3;
 E[1] = 0;
 node[0].x[0]=0;
 node[0].x[1]=0;
 node[0].id = 1;	// nó 1
 node[1].x[0]=3; 
 node[1].x[1]=0;
 node[1].id = 2;	// nó 2	
 node[2].x[0]=0;
 node[2].x[1]=4;
 node[2].id = 3;	// nó 3
 
 for (i = 0; i < ne; i++)
 {
 elem[i].area = 0.001164;	//m^2 
 elem[i].E = E[0];
 elem[i].sigmay = 250;  //MPa
 elem[0].no[0] = 0;
 elem[0].no[1] = 1;  // entrada do usuário para montar a treliça
 elem[1].no[0] = 1;
 elem[1].no[1] = 2;
 elem[2].no[0] = 2;
 elem[2].no[1] = 0;
 elem[i].sigma0 = 0.0;
 }

#endif
 //Problema do Gere 2.12.5
#if 0
 ne = 4;
 nn = 5;
 elem = (Elem*)malloc(ne * sizeof(Elem));
 node = (Node*)malloc(nn * sizeof(Node));
 E[0] = 200E9;
 E[1] = 0;
 //coordenadas dos nós

 node[0].x[0]=-1.600;
 node[0].x[1]=1.200;
 node[0].id = 1;	// nó 1
 node[1].x[0]=-0.900; 
 node[1].x[1]=1.200;
 node[1].id = 2;	// nó 2	
 node[2].x[0]=0.900;
 node[2].x[1]=1.200;
 node[2].id = 3;	// nó 3
 node[3].x[0]=1.600;
 node[3].x[1]=1.200;
 node[3].id = 4;	// nó 4
 node[4].x[0]=0.0;
 node[4].x[1]=0.0;
 node[4].id = 5;	// nó 5 
//elementos
  elem[0].no[0] = 0;
  elem[0].no[1] = 4;
  elem[1].no[0] = 1;
  elem[1].no[1] = 4;
  elem[2].no[0] = 2;
  elem[2].no[1] = 4;
  elem[3].no[0] = 3;
  elem[3].no[1] = 4;
  //area
  elem[0].area=0.000200;	//mm^2 
  elem[3].area = 0.000200;	//mm^2
  elem[1].area = 0.000400;	//mm^2
  elem[2].area = 0.000400;	//mm^2
  for (i = 0; i < ne; i++)
  {
   elem[i].region = 0;
	  elem[i].E = E[0];
	  elem[i].sigma0 = 0.0;
	  elem[i].sigmay = 250E6;
  }
  for (i = 0; i < nn; i++)
   node[i].u0[0] = node[i].u0[1] = 0.0;
#endif

//problema Gere 2.12.6
#if 0
  ne = 5;
  nn = 6;
  elem = (Elem*)malloc(ne * sizeof(Elem));
  node = (Node*)malloc(nn * sizeof(Node));
  E[0] = 200E3;
  E[1] = 0;
  //coordenadas dos nos
  node[0].x[0] = -2;
  node[0].x[1] = 2;
  node[0].id = 1;
  node[1].x[0] = -1;
  node[1].x[1] = 2;
  node[1].id = 2;
  node[2].x[0] = 0;
  node[2].x[1] = 2;
  node[2].id = 3;
  node[3].x[0] = 1;
  node[3].x[1] = 2;
  node[3].id = 4;
  node[4].x[0] = 2;
  node[4].x[1] = 2;
  node[4].id = 5;
  node[5].x[0] = 0;
  node[5].x[1] = 0;
  node[5].id = 6;
  //elementos
  elem[0].no[0] = 0;
  elem[0].no[1] = 5;
  elem[1].no[0] = 1;
  elem[1].no[1] = 5;
  elem[2].no[0] = 2;
  elem[2].no[1] = 5;
  elem[3].no[0] = 3;
  elem[3].no[1] = 5;
  elem[4].no[0] = 4;
  elem[4].no[1] = 5;
  double D = 10; //mm^2
  double a = pi * D * D / 4;
  for (i = 0; i < ne; i++)
  {
   elem[i].region = 0;
   elem[i].area = a;
   elem[i].E = E[0];
   elem[i].sigma0 = 0.0;
   elem[i].sigmay = 250;
  }
#endif
  for (i = 0; i < ne; i++)
  { 
   L(i);
   B(i);  	
  }
 }		

static void Ks()
{
	for (int i = 0; i < ne; i++)
	{
		K(i);
	}
}

static void Carregamento()
{
 //define o carregamento no vetor de nós
 int i;
#if 0
 // Gere 2.12.5
 node[4].p[0]=0;
 node[4].p[1]=-1.0;	// entrada do usuario fx=0 e fy=-1
 for(i=1;i<4;i++)
 { 
  node[i].p[0]=node[i].p[1]=0;
 }
 noderesult = &node[4];
#endif
#if 0
  // 3 barras 
 node[0].p[0]=node[0].p[1]=node[2].p[0]=node[2].p[1]=0;
  node[1].p[0]=-2;
  node[1].p[1]=-1;
#endif
#if 1
  //paper
  for (i = 0; i < nn; i++)
  {
	  for (int j = 0; j < 2; j++)
	  {
		  node[i].p[j] = 0;
		  if (i == 5)
		  {
			  node[i].p[1] = 1;
			  node[i].p[0] = 0;
		  }
	  }
  }
  noderesult = &node[5];
#endif
#if 0
 //Gere 2.12.6
  for (i = 0; i < nn; i++)
  {
	  node[i].p[0] = node[i].p[1] = 0;
  }
  node[5].p[1] = -1;
  noderesult = &node[5];
#endif
}

static void Contorno() 
{
 //define as condições de contorno essenciais no vetor de nós
 int i;
#if 0
 //Gere 2.12.5
 for(i=0;i<nn;i++)
	{
		node[i].free[0]=false;
		node[i].free[1]=false;
	}
	node[4].free[0]=node[4].free[1]=true;
#endif
#if 0
// 3 barras
 node[0].free[0]=false;
 node[0].free[1]=false;
 node[1].free[0]=true;
 node[1].free[1]=true;
 node[2].free[0]=false;
 node[2].free[1]=true;*/
#endif
#if 1
 //paper
 for (i = 0; i < nn; i++)
 {
	 for (int j = 0; j < 2; j++)
	 {
		 if (i == 3 || i == 0)
		 {
			 node[i].free[0] = node[i].free[1] = false;
		 }
		 else
			 node[i].free[0] = node[i].free[1] = true;
	 }
 }
#endif
#if 0
 for (i = 0; i < nn; i++)
 {
  node[i].free[0] = node[i].free[1] = false;
 }
 node[5].free[0] = node[5].free[1] = true;
#endif
 for(i=0;i<nn;i++)
  for(int j=0;j<gl;j++)
   if(node[i].free[j]==false)
    nr++;					//calcula o numero de reaçoes 
 iso = nr+ne-2*nn; //grau hiperestático
#if 1 //paper
 iso += 6;
#endif
}

static void IncidenciaCinematica()
{
 int i,j;
//montagem da matriz de incidencia cinematica
 mic=(int**)malloc(nne*gl*sizeof(int*)); // numero de nos por elemento
 for(i=0;i<nne*gl;i++)
 {					 // numero de elemtos   mic[linhas=ne][colunas=nne]
  mic[i]=(int*)malloc(ne*sizeof(int));
 }
	 
  for(i=0;i<ne;i++)
  {
   for(j=0;j<nne;j++)
   {
  	 mic[j*2][i]=node[ elem[i].no[j] ].id*2-1;
    mic[j*2+1][i]=node[ elem[i].no[j] ].id*2;
   }	
  }
  // aplicação das condições de contorno
  for(i=nn-1;i>=0;i--)
  {
   for(j=gl-1;j>=0;j--)
   {
   	if(node[i].free[j]==false)
	   {
	   for(int ii=0;ii<ne;ii++)
	   {
	    for(int jj=0;jj<nne*gl;jj++)
	    {
		    if(mic[jj][ii]==node[i].id*gl-1+j)
		     mic[jj][ii]=0;
		    else if(mic[jj][ii]>node[i].id*gl-1+j)
		     mic[jj][ii]--;   
	     }
	    }
	   }
   } 
  }
}

static void kgFg()
{
	//montas as matrizes de rigidez global e o vetor de forças globais
 int i,j,l,lin,col;
 /* alocação de kg,fg e u */
 lin=nn*gl-nr;		// ta dando lin=2
 if (kg == NULL)
 {
  kg = (double**)malloc(lin * sizeof(double*));
  for (i = 0; i < lin; i++)
   kg[i] = (double*)malloc(lin * sizeof(double)); //kg[2][2]
 }
 if (fg == NULL)
  fg=(double*)malloc(lin*sizeof(double));     //fg[2][1]
 if (u == NULL)
   u=(double*)malloc(lin*sizeof(double));		//u[2][1]
 //inicializa a matriz kg e o vetor fg com zero
 for(i=0;i<lin;i++)
 {
  fg[i]=0.0;
  for(j=0;j<lin;j++)
   kg[i][j]=0.0;
 }
 /* monta a matriz de rigidez global */ 
 for(i=0;i<ne;i++)
 {
  for(j=0;j<2*gl;j++)
  {
   lin=mic[j][i]-1;
   if(lin<0)
    continue;
   for(l=0;l<2*gl;l++)
   {
    col=mic[l][i]-1;
    if(col<0)
     continue;
    kg[lin][col]=kg[lin][col]+elem[i].k[j][l];
   }
  }
  
 }
 /* monta o vetor de forças global */
 lin=0;
 for(i=0;i<nn;i++)
 {
  for(j=0;j<gl;j++)
  {
   if(node[i].free[j]==true)
   fg[lin++]+=node[i].p[j]; 
  }
 }
 
}
static void Solucao()
{
 //inverte a matriz de rigidez e obtém os deslocamentos (solucao)
	int i,j,k;
 	double c,d;
 	int n=nn*gl-nr;	//n=2	reacao de segunda ordem
 	for(i=0;i<n;i++) //i até 2
 	{
  	if(fabs(kg[i][i])>0)
  	{
  	 	c=1.0/kg[i][i];
   		for(k=0;k<n;k++)
   		 kg[i][k]=kg[i][k]*c;
   		kg[i][i]=c;
  	}
  	for(j=0;j<n;j++) //j até 2
   	if(j!=i)
   	{
    	d=kg[j][i];
    	for(k=0;k<n;k++)
     {
     		if(k==i)
      			kg[j][i]=0;
     			kg[j][k]-=d*kg[i][k];
    }
   }
 }		// nao modifiquei nada na rotina de inverção da matriz, devo modificar? nao precisa
 //deslocamentos
  
 	for(i=0;i<n;i++)
 	{
  	u[i]=0;
  	for(j=0;j<n;++j)
	  {	  
    u[i]+=kg[i][j]*fg[j];
	  }
 } 
}

static void Tensoes()
{
 //calcula as deformacoes e tensoes em cada elemento
 int i,j,id;
 for (i = 0; i < ne; i++)
 {

	 elem[i].epsilon = 0;
	 for (j = 0; j < 2 * gl; j++) //loop nos deslocamentos por elemento
	 {
		 id = mic[j][i] - 1;//id do vetor u
		 if (id > -1)
			 elem[i].epsilon += elem[i].b[j] * u[id];  //e=[b][u]
   if (j == 0)
   {
    if (id > -1)
     node[elem[i].no[0]].u[0] = u[id];
    else
     node[elem[i].no[0]].u[0] = 0;
   }
   else if (j == 1)
   {
    if (id > -1)
     node[elem[i].no[0]].u[1] = u[id];
    else
     node[elem[i].no[0]].u[1] = 0;
   }
   else if (j == 2)
   {
    if (id > -1)
     node[elem[i].no[1]].u[0] = u[id];
    else
     node[elem[i].no[1]].u[0] = 0;
   }
   else if (j == 3)
   {
    if (id > -1)
     node[elem[i].no[1]].u[1] = u[id];
    else
     node[elem[i].no[1]].u[1] = 0;
   }
  }
	 elem[i].sigma = elem[i].epsilon * elem[i].E;
	 //printf("\n\n TENSAO DA BARRA %d : %lf \n\n", i + 1, elem[i].sigma);
 }
 //id = 0;
 double min=0,min2;
 for (i = 0; i < ne; i++)
  if (elem[i].region == 0)
  {
   if (fabs(elem[i].sigma) > 0.00001)
    min = (elem[i].sigmay - fabs(elem[i].sigma0)) * 1.001 / fabs(elem[i].sigma);
   id = i;
   break;
  }
 for (i = 0; i < ne; i++)
 {
  if (fabs(elem[i].sigma) > 0.00001)
   min2 = (elem[i].sigmay - fabs(elem[i].sigma0)) / fabs(elem[i].sigma);
  if(min2<min && elem[i].region == 0)
  {
	  min = min2;
  	id=i;
  }
 }
 double lambda0 = lambda;
 lambda = min;// (elem[id].sigmay - fabs(elem[id].sigma0)) / fabs(elem[id].sigma);
 //max*=lambda*0.999;
 int cont=0;
 for(i=0;i<ne;i++)
 {
 	elem[i].sigma*=lambda;
#if 0
	 if(elem[i].sigma>=max)
  {
   elem[i].region++;
	   if (elem[i].region > 1)
		   elem[i].region = 1;
	 	elem[i].E=E[elem[i].region];
	 	//pp=elem[i].sigma/elem[i].area;
   cont++;
  }
#endif
  elem[i].sigma0+=elem[i].sigma;
  if (fabs(elem[i].sigma0) >= elem[i].sigmay * 0.999)
  {
   elem[i].region++;
   if (elem[i].region > 1)
   elem[i].E = E[elem[i].region];
   cont++;
  }
 }
 for (i = 0; i < nn; i++)
 {
  node[i].u0[0] += lambda*node[i].u[0];
  node[i].u0[1] += lambda*node[i].u[1];
 }
 lambda+=lambda0;
 iso-=cont;
}

static void Print()
{
 int flag;
#if 1
 if (fp == NULL)
 {
  fopen_s(&fp, "D:\\Users\\cezaralejandro\\Desktop\\trelica_ap\\saida.out", "w");
  flag = 1;
 }
 else
  flag = 0;
#else
 if (fp == NULL)
 {
  fopen_s(&fp, "D:\\UFSM\\projetos_pesquisa\\cezar\\saida.out", "w");
  flag = 1;
 }
 else
  flag = 0;
#endif
 if (0)
 {
  //saida em kN, MPa, mm
  if(flag==1)
  {
   fprintf(fp, "                kN_P,                mm_u,           kN_abs_P,            mm_abs_u");
   for (int i = 0; i < ne; i++)
    fprintf(fp, ",         MPa_sigma_%d", i + 1);
   fprintf(fp, "\n");
   fprintf(fp, "%20.5f,%20.5f,%20.5f,%20.5f", 0.0, 0.0, 0.0, 0.0);
   for (int i = 0; i < ne; i++)
    fprintf(fp, ",%20.5f", 0.0);
   fprintf(fp, "\n");
  }
  fprintf(fp, "%20.5f,%20.15f,%20.5f,%20.15f", 0.001 * lambda * noderesult->p[1], 1000 * noderesult->u0[1],
   0.001 * fabs(lambda * noderesult->p[1]), 1000 * fabs(noderesult->u0[1]));
  for (int i = 0; i < ne; i++)
   fprintf(fp, ",%20.5f", 1E-6 * elem[i].sigma0);
 }
 else
 {
  //saida em N, Pa, m
  if (flag == 1)
  {
   fprintf(fp, "                 N_P,                 m_u,            N_abs_P,             m_abs_u");
   for (int i = 0; i < ne; i++)
    fprintf(fp, ",          Pa_sigma_%d", i + 1);
   fprintf(fp, "\n");
   fprintf(fp, "%20.5f,%20.5f,%20.5f,%20.5f", 0.0, 0.0, 0.0, 0.0);
   for (int i = 0; i < ne; i++)
    fprintf(fp, ",%20.5f", 0.0);
   fprintf(fp, "\n");
  }
  fprintf(fp, "%20.5f,%20.15f,%20.5f,%20.15f", lambda * noderesult->p[1], noderesult->u0[1],
                                             fabs(lambda * noderesult->p[1]),fabs(noderesult->u0[1]));
  for (int i = 0; i < ne; i++)
   fprintf(fp, ",%20.5f",elem[i].sigma0);
 }
 fprintf(fp, "\n");
}

static void Free()
{
 //libera a memória alocada dinamicamente
 free(elem);
 free(node);
 free(mic);
 for(int i=0;i<nn*gl-nr;i++)
  free(kg[i]);
 free(kg);
 free(fg);
 free(u);
}


int main ()
{
	Discretizacao();
	Carregamento();
	Contorno();
	IncidenciaCinematica();
 lambda = 0;
	while(iso>=0)
	{
	 Ks();
  kgFg();	
	 Solucao();
	 Tensoes();	
	 Print();
 }
 fclose(fp);
	Free();
}
