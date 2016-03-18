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
double T,Eaa,Ebb,Eba,conc;
f=fopen("input.dat","r");
fscanf(f,"%d",&N);
fscanf(f,"%le",&conc);
fscanf(f,"%le",&T);
fscanf(f,"%le",&Eaa);
fscanf(f,"%le",&Ebb);
fscanf(f,"%le",&Eba);
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
			comp[j+N*i]=-1;
		}
		else
		{
			comp[j+N*i]=0;
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
			Ei+=Eba;
		}
		else
		{
			if(comp[j+N*i]==-1)Ei+=Eaa;
			else Ei+=Ebb;
		}
	}

	//flip state
	if(comp[j+N*i]==-1)comp[j+N*i]=0;
	else comp[j+N*i]=-1;

	//find Ef
	Ef=0;
	for(l=1;l<=8;l++)
	{
		getNeighbour(i,j,N,l,&k1,&k2);
		if(comp[j+N*i]!=comp[k2+N*k1])
		{
			Ef+=Eba;
		}
		else
		{
			if(comp[j+N*i]==-1)Ef+=Eaa;
			else Ef+=Ebb;
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
			if(comp[j+N*i]==-1)comp[j+N*i]=0;
			else comp[j+N*i]=-1;
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

//Now that the Monte Carlo Simulation has been done
//We will apply the Hoshen Kopelman Algorithm

int m=1;//label

int *numCluster;//number of atoms in each cluster
numCluster=(int*)calloc(N*N,sizeof(int));//giving it the maximum space //TODO optimise this
int k3,k4;

//first traversal--make the initial groups
for(i=0;i<N;i++)
{
	{
	for(j=0;j<N;j++)
		if(comp[j+N*i]==-1)
		{
					//if cheking the first row dont look on the left side
					if(j==0)
					{
						getNeighbour(i,j,N,4,&k1,&k2);
						if(comp[k2+N*k1])
						{
							comp[j+N*i]=comp[k2+N*k1];
						}
					}
					else
					{
						getNeighbour(i,j,N,4,&k1,&k2);
						getNeighbour(i,j,N,1,&k3,&k4);

						//if both are A
						if(comp[k2+N*k1] && comp[k4+N*k3])
						{
							//find the minimum cluster lable
							if(comp[k2+N*k1] < comp[k4+N*k3])
							{
								comp[j+N*i]=comp[k2+N*k1];
							}
							else
							{
								comp[j+N*i]=comp[k4+N*k3];
							}
						}
						else if(comp[k2+N*k1] && !comp[k4+N*k3])
						{
							comp[j+N*i]=comp[k2+N*k1];
						}
						else if(comp[k4+N*k3] && !comp[k2+N*k1])
						{
							comp[j+N*i]=comp[k4+N*k3];
						}
					}
					if(comp[j+N*i]==-1)
					{
						comp[j+N*i]=m++;
					}
					numCluster[comp[j+N*i]]+=1;
		}
	}
}

//print to file
f=fopen("HKstage1.dat","w");
for(i=0;i<N;i++)
{
	for(j=0;j<N;j++)
	{
		fprintf(f,"%d\t",comp[j+N*i]);
	}
	fprintf(f,"\n");
}
fclose(f);

//second traversal-- make the final groups
//find the maximum in numCluster
int max=0,maxLabel;
for(k1=1;k1<=m;k1++)
{
	if(max<numCluster[k1])
	{
		max=numCluster[k1];
		maxLabel=k1;
	}
}

//now the cluster lable pos wil be made a particular color and all other lables will be merged together
for(i=0;i<N;i++)
{
	for(j=0;j<N;j++)
	{
		if(comp[j+N*i]==maxLabel)
		{
			comp[j+N*i]=200;
		}
		else if(comp[j+N*i])
		{
			comp[j+N*i]=1;
		}
	}
}

//print to file
f=fopen("HKfinal.dat","w");
for(i=0;i<N;i++)
{
	for(j=0;j<N;j++)
	{
		fprintf(f,"%d\t",comp[j+N*i]);
	}
	fprintf(f,"\n");
}
fclose(f);

free(numCluster);

free(comp);
gsl_rng_free(r);
}
