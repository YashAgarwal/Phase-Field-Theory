#include<stdio.h>
#include<math.h>

#include<gsl/gsl_math.h>
#include<gsl/gsl_rng.h>

void getNeighbour(int i1,int i2,int	*k1,int	*k2,int dir,int N)
{
	if(dir==1)
	{
	*k1=i1+1;
	*k2=i2-1;
	}
	else if(dir==2)
	{
	*k1=i1+1;
	*k2=i2;
	}
	else if(dir==3)
	{
	*k1=i1+1;
	*k2=i2+1;
	}
	else if(dir==4)
	{
	*k1=i1;
	*k2=i2+1;
	}
	else if(dir==5)
	{
	*k1=i1-1;
	*k2=i2+1;
	}
	else if(dir==6)
	{
	*k1=i1-1;
	*k2=i2;
	}
	else if(dir==7)
	{
	*k1=i1-1;
	*k2=i2-1;
	}
	else if(dir==8)
	{
	*k1=i1;
	*k2=i2-1;
	}

	if(*k1==N)*k1=0;
	else if(*k1==-1)*k1=N-1;

	if(*k2==N)*k2=0;
	else if(*k2==-1)*k2=N-1;
}

int main(void)
{
FILE *fpw;

int N;
double Eaa,Ebb,Eab,conc,t;

fpw=fopen("input.dat","r");
fscanf(fpw,"%d",&N);
fscanf(fpw,"%le",&conc);//concentration for B
fscanf(fpw,"%le",&Eaa);
fscanf(fpw,"%le",&Ebb);
fscanf(fpw,"%le",&Eab);
fclose(fpw);

int *comp;
comp= (int*)malloc(N*N*sizeof(int));

double u;
const gsl_rng_type *T;
gsl_rng *r;
gsl_rng_env_setup();
T=gsl_rng_default;
r=gsl_rng_alloc(T);


int i1,i2,j,k,dir,k1,k2;
double Ei,Ef,delE;

//Making the initial profile
//B=1
//A=0

for(i1=0;i1<N;i1++)
{
	for(i2=0;i2<N;i2++)
	{
		u=gsl_rng_uniform(r);
		if(u<conc)
		comp[i2+N*i1]=1;
		else
		comp[i2+N*i1]=0;
	}
}


//print
fpw=fopen("initialProfile.dat","w");
for(i1=0;i1<N;i1++)
{
	for(i2=0;i2<N;i2++)
	{
		fprintf(fpw,"%d\t",comp[i2+N*i1]);
	}
	fprintf(fpw,"\n");
}
fclose(fpw);

for(j=0;j<8*N*N;j++)
{
	u=gsl_rng_uniform(r);
	i1=u*1e4;
	i1=((int)i1)%N;

	u=gsl_rng_uniform(r);
	i2=u*1e4;
	i2=((int)i2)%N;

	u=gsl_rng_uniform(r);
	dir=u*1e4;
	dir=((int)dir)%8+1;

	Ei=0;
	for(k=1;k<=8;k++)
	{
		getNeighbour(i1,i2,&k1,&k2,k,N);
		if(comp[i2+N*i1]!=comp[k2+N*k1])
		{
			Ei+=Eab;
		}
		else if(comp[i2+N*i1]==comp[k2+N*k1])
		{
			if(comp[k2+N*k1]==1)
			{
				Ei+=Ebb;
			}
			else
			{
				Ei+=Eaa;
			}
		}
	}

	//swap
	getNeighbour(i1,i2,&k1,&k2,dir,N);
	t=comp[k2+N*k1];
	comp[k2+N*k1]=comp[i2+N*i1];
	comp[i2+N*i1]=t;

	Ef=0;
	for(k=1;k<=8;k++)
	{
		getNeighbour(i1,i2,&k1,&k2,k,N);
		if(comp[i2+N*i1]!=comp[k2+N*k1])
		{
			Ef+=Eab;
		}
		else if(comp[i2+N*i1]==comp[k2+N*k1])
		{
			if(comp[k2+N*k1]==1)
			{
				Ef+=Ebb;
			}
			else
			{
				Ef+=Eaa;
			}
		}
	}

	delE=Ef-Ei;
	if(delE>0)
	{
		//swap
		getNeighbour(i1,i2,&k1,&k2,dir,N);
		double t=comp[k2+N*k1];
		comp[k2+N*k1]=comp[i2+N*i1];
		comp[i2+N*i1]=t;
	}
}

//print
fpw=fopen("finalProfile.dat","w");
for(i1=0;i1<N;i1++)
{
	for(i2=0;i2<N;i2++)
	{
		fprintf(fpw,"%d\t",comp[i2+N*i1]);
	}
	fprintf(fpw,"\n");
}
fclose(fpw);

free(comp);
}
