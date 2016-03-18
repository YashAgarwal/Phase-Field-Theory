#include<stdlib.h>
#include<stdio.h>
#include<math.h>

#include<gsl/gsl_math.h>
#include<gsl/gsl_rng.h>

void getNeighbour(int i1,int i2,int N,int dir,int *k1,int *k2)
{
	if(dir==1)
	{
	*k1=i1;
	*k2=i2-1;
	}
	else if(dir==2)
	{
	*k1=i1;
	*k2=i2+1;
	}
	else if(dir==3)
	{
	*k1=i1+1;
	*k2=i2;
	}
	else if(dir==4)
	{
	*k1=i1-1;
	*k2=i2;
	}
	else if(dir==5)
	{
	*k1=i1-1;
	*k2=i2+1;
	}
	else if(dir==6)
	{
	*k1=i1+1;
	*k2=i2-1;
	}
	else if(dir==7)
	{
	*k1=i1+1;
	*k2=i2+1;
	}
	else if(dir==7)
	{
	*k1=i1-1;
	*k2=i2-1;
	}

	//periodic boundary
	if(*k1==N)*k1=0;
	else if(*k1==-1)*k1=N-1;

	if(*k2==N)*k2=0;
	else if(*k2==-1)*k2=N-1;
}

int main(void)
{

FILE *f;
int N;
double T,Euu,Edd,Edu,conc;
f=fopen("input.dat","r");
fscanf(f,"%d",&N);
fscanf(f,"%le",&conc);
fscanf(f,"%le",&T);
fscanf(f,"%le",&Euu);
fscanf(f,"%le",&Edd);
fscanf(f,"%le",&Edu);

fclose(f);

int *comp;
comp=(int*)malloc(N*N*sizeof(int));

int i,j;

const gsl_rng_type *t;
gsl_rng *r;
gsl_rng_env_setup();
t=gsl_rng_default;
r=gsl_rng_alloc(t);

double u;

//initialisation of the matrix
//print to file
f=fopen("initialProfile.dat","w");
for(i=0;i<N;i++)
{
	for(j=0;j<N;j++)
	{
		u=gsl_rng_uniform(r);
		if(u<conc)
		{
			comp[j+N*i]=1;
		}
		else
		{
			comp[j+N*i]=-1;
		}
		fprintf(f,"%d\t",comp[j+N*i]);
	}
	fprintf(f,"\n");

}

fclose(f);

int k,l,k1,k2;
double Ei,Ef,delE;

//main evolution loop
for(k=0;k<N*N;k++)
{
	//choose random cell
	u=gsl_rng_uniform(r);
	i=(int)(u*1e4)%N;

	u=gsl_rng_uniform(r);
	j=(int)(u*1e4)%N;

	//find Ei
	Ei=0;
	for(l=1;l<=8;l++)
	{
		getNeighbour(i,j,N,l,&k1,&k2);
		if(comp[j+N*i]!=comp[k2+N*k1])
		{
			Ei+=Edu;
		}
		else
		{
			if(comp[j+N*i]==1)Ei+=Euu;
			else Ei+=Edd;
		}
	}

	//flip state
	comp[j+N*i]*=-1;

	//find Ef
	Ef=0;
	for(l=1;l<=8;l++)
	{
		getNeighbour(i,j,N,l,&k1,&k2);
		if(comp[j+N*i]!=comp[k2+N*k1])
		{
			Ef+=Edu;
		}
		else
		{
			if(comp[j+N*i]==1)Ef+=Euu;
			else Ef+=Edd;
		}
	}

	delE=Ef-Ei;

	if(delE>0)
	{
		u=gsl_rng_uniform(r);
		double p=exp(-delE/T);
		if(p<=u)
		{
			//flip state
			comp[j+N*i]*=-1;
		}
	}

}

//print to file
f=fopen("finalProfile.dat","w");
for(i=0;i<N;i++)
{
	for(j=0;j<N;j++)
	{
		fprintf(f,"%d\t",comp[j+N*i]);
	}
	fprintf(f,"\n");
}
fclose(f);

free(comp);
gsl_rng_free(r);
}
