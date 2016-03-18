/*making initial file*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(void)
{
FILE *fp;
int N,i;

fp=fopen("input.dat","r");
fscanf(fp,"%d",&N);//number of nodes 
fclose(fp);

double *c;//This vector will hold the value of concentration for all N nodes
c=(double*)malloc(N*sizeof(double));

fp=fopen("c_0.dat","w");
//writes the initial conc. profile 

for(i=0;i<N;i++)
{
	if(i<N/2)
	{c[i]=1;
	}
	else
	{c[i]=0;
	}
	fprintf(fp,"%le\n",c[i]);
}
fclose(fp);

return 0;
}
