/*Explicit Method
using boundary condition (II)
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(void)
{
FILE *fp;
double alpha,delt,delx;
int T,N;
fp=fopen("input.dat","r");
//take all other inputs except the initial conc. profile 

fscanf(fp,"%d",&N);//number of nodes 
fscanf(fp,"%le",&delx);
fscanf(fp,"%d",&T);//number of time steps
fscanf(fp,"%le",&delt);
fscanf(fp,"%le",&alpha);
/*
T & N should be chosen such that T*delt=1
and N*deln=1
*/

fclose(fp);

int i,j;
char *filename;
filename=(char*) malloc(50*sizeof(char));

double *c,*cNew;//These vector will hold the value of concentration for all N nodes
c=(double*) malloc(N*sizeof(double));
cNew=(double*) malloc(N*sizeof(double));

fp=fopen("c_0.dat","r");
//inputs the initial conc. profile 

for(i=0;i<N;i++)
{
fscanf(fp,"%le",&c[i]);
}

fclose(fp);

//using boundary condition (II)
double b1,b2,beta;
beta=0;
for(j=1;j<=T;j++)
{
	for(i=0;i<N;i++)
	{
		if(i==0)
		{		
		b1=c[1]-beta*2*delx;		
		cNew[i]=(double)(c[i] + (alpha*delt*(c[i+1]+b1-2*c[i]))/(delx*delx));				
		}
		else if(i==N-1)
		{
		b2=c[N-2]+beta*2*delx;
		cNew[i]=(double)(c[i] + (alpha*delt*(b2+c[i-1]-2*c[i]))/(delx*delx));		
		}		
		else
		{		
		cNew[i]=(double)(c[i] + (alpha*delt*(c[i+1]+c[i-1]-2*c[i]))/(delx*delx));		
		}	
	}
	//transfer the cNew to a file and copy cNew to c
	sprintf(filename,"./output/c_%d.dat",j);	
	fp=fopen(filename,"w");
	
	for(i=0;i<N;i++)
	{
	fprintf(fp,"%le\n",cNew[i]);
	c[i]=cNew[i];
	}
	
fclose(fp);
}

free(filename);
free(c);
free(cNew);
return 0;
}
