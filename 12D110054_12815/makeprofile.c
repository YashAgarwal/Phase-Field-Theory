#include<stdlib.h>
#include<stdio.h>
#include<math.h>

int main(void)
{
FILE *fp;
fp=fopen("input.dat","r");

double N;
int i,j;

fscanf(fp,"%le",&N);
fclose(fp);

fp=fopen("c_0.dat","w");

double *c;
c=(double*)malloc(N*sizeof(double));

for(i=0;i<N;i++)
{
	if(i>N/4 && i<=3*N/4)
	{
		c[i]=1;
	}
	else
	{
		c[i]=0;
	}
	fprintf(fp,"%le\n",c[i]);
}

fclose(fp);

free(c);
return 0;
}

