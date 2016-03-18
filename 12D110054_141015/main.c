/************************************************************
DISLOCATION DYNAMICS 2D

gcc main.c -lgsl -lgslcblas -lm

Author: Yash Agarwal
************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h> /*Statistical distributions*/
#include <gsl/gsl_statistics.h> /*Functions like mean and variance*/

int main(void) {

//Initialize the Dislocation Matrix
//we have n dislocations in matrix of size DxD
//n should be even
int n,D;
char name[50];
//nsteps = number of time steps
//dxmax = maximum move per time step
double nsteps,dxmax;

//setting the matrix size
D=1000;

FILE *f,*g;
f=fopen("input.dat","r");
g=fopen("plotAnimation.gp","w");

fprintf(g,"set size square\nset xrange[0:1000]\nset yrange[0:1000]\nunset colorbox\nunset key; unset tics;\n");

fscanf(f,"%d",&n);
fscanf(f,"%d",&D);
fscanf(f,"%le",&nsteps);
fscanf(f,"%le",&dxmax);
fclose(f);

int rc=D/2;
int tWrite=nsteps/50;

double *x, *y, *b;
//const gsl_rng_type *T;
//gsl_rng *r;
//gsl_rng_env_setup();
//T=gsl_rng_default;
//r=gsl_rng_alloc(T);


gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937); /*Define random number t*/
long seed = time(NULL)*getpid(); /*Define a seed for r*/
gsl_rng_set(r,seed); /*Initiate the random number generator with seed*/

//x contains the x coordinate of all the n dislocations
x= (double*)malloc(n*sizeof(double));

int i;
double u;

for(i=0;i<n;i++)
{
	u=gsl_rng_uniform(r);
	x[i]=u*D;
}

//y contains the y coordinate of all the n dislocations
y= (double*)malloc(n*sizeof(double));

for(i=0;i<n;i++)
{
	u=gsl_rng_uniform(r);
	y[i]=u*D;
}
//b contains the sign of the burges vector of all the n dislocations
b= (double*)calloc(n,sizeof(double));

for(i=0;i<n;i++)
{
	if(i<=n/2)
	{
		b[i]=1;
	}
	else
	{
		b[i]=-1;
	}
}

gsl_rng_free(r);

//print the current profile
f=fopen("./output/step0.dat","w");
for(i=0;i<n;i++)
{
	fprintf(f,"%d\t%d\t%d\n",(int)x[i],(int)y[i],(int)b[i]);
}
fclose(f);

fprintf(g,"set title \"Dislocation Dynamics @t=0\"\nplot \"./output/step0.dat\" u 1:2:3 with points palette pt 7 notitle\npause 0.1\n");

//xi=initial x positions,
//x = final positions
//fx = final forces on each dislocation
//xdm = the maximum distance a dislocation has moved in the calculation

//Setting the initial positions
double *xi;
xi= (double*)malloc(n*sizeof(double));
for (i=0;i<n;i++)//note limits
{
	xi[i]=x[i];
}

//start the time steps
int j;
double *fx;
fx=(double*)calloc(n,sizeof(double));

for(j=1;j<=nsteps;j++)
{
		//calculating fx
		//and fMax

		int k;
		for (i=0;i<n;i++)
		{
			fx[i]=0;
		}

		for (i=0;i<n-1;i++)//note limits
			{
				for(k=i+1;k<n;k++)//note limits
				{
				// mimimum image convention
						double dx = x[i] - x[k];
						double dy = y[i] - y[k];
						dx = dx - D*round(dx/D);
						dy = dy - D*round(dy/D);
						double dsq = dx*dx + dy*dy;
						double dist = sqrt(dsq);
						double ffx;
						if (dist <= rc)
						{
							ffx = b[i]*b[k]*dx*(dx*dx-dy*dy)/(dsq*dsq);
							fx[i]+= ffx;
							fx[k]-= ffx; //add -f to sum of force on j
						}
				}
			}

		// calculate maximum value of force
		double fMax=0;
		for(i=0;i<n;i++)
		{
			if(fabs(fx[i])>fMax)
			fMax=fabs(fx[i]);
		}

		//calculate the new positions
		double dt=dxmax/fMax;

		for(i=0;i<n;i++)
		{
			x[i]+=fx[i]*dt;
			if(x[i] > D)
			x[i]-=D;
			if(x[i] < 0)
			x[i]+=D;
		}

		if(j%tWrite==0){
		//print the current profile
		sprintf(name,"./output/step%d.dat",j);
		f=fopen(name,"w");

		for(i=0;i<n;i++)
		{
			fprintf(f,"%d\t%d\t%d\n",(int)x[i],(int)y[i],(int)b[i]);
		}
		fclose(f);
		fprintf(g,"set title \"Dislocation Dynamics @t=%d\"\nplot \"%s\" u 1:2:3 with points palette pt 7 notitle\npause 0.1\n",j,name);
		}
}
//time steps loop ends

//Calculating the Maximum Movement
double *xd;
xd= (double*)malloc(n*sizeof(double));

double xdm=0;
for(i=0;i<n;i++)
{
	//change in position over the run
	xd[i]=x[i]-xi[i];

	// remove movement across periodic boundaries
	xd[i]-=D*round(xi[i]/D);

	// find maximum movement
	if(xdm<fabs(xd[i]))
	xdm=fabs(xd[i]);
}

//print the Maximum Movement
f=fopen("output.dat","w");
fprintf(f,"The Maximum Movement is %le",xdm);
fclose(f);

fclose(g);


free(x);
free(y);
free(b);
free(xi);
free(fx);
free(xd);

}
