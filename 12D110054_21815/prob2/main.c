/************************************************************
use this to compile
gcc main.c -lgsl -lgslcblas -lfftw3 -lm
or just use ./run.sh to compile,run and plot the result 

Cahn hilliard Equation
************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <fftw3.h>

int main(void) {

FILE *fpw;

int i1,i2;
int half_nx;
double kx,delta_kx;
double kx2;
char name[50];
int INDEX;

double A = 1.0;
double alpha = 1.0;
double delta_x = 1.0;
double delta_t = 0.2;
int n_x=128;//number of nodes in the system 
int T,T_write;//(T)number of time steps
//T_write is after how many time steps we write profile in file

//Taking the input parameters from a input file
fpw=fopen("input.dat","r");
fscanf(fpw,"%d",&n_x);
fscanf(fpw,"%le",&delta_x);
fscanf(fpw,"%d",&T);
fscanf(fpw,"%le",&delta_t);
fscanf(fpw,"%d",&T_write);
fclose(fpw);
fftw_complex *comp,*g;
fftw_plan planF,planFg,planB;

comp = fftw_malloc(n_x* sizeof(fftw_complex));
g = fftw_malloc(n_x* sizeof(fftw_complex));

planF = fftw_plan_dft_1d(n_x,comp,comp,FFTW_FORWARD,FFTW_ESTIMATE);
planFg = fftw_plan_dft_1d(n_x,g,g,FFTW_FORWARD,FFTW_ESTIMATE);
planB = fftw_plan_dft_1d(n_x,comp,comp,FFTW_BACKWARD,FFTW_ESTIMATE);

/**Creating Random Numbers uisng  gsl library**/

const gsl_rng_type * T1;
gsl_rng * r;
gsl_rng_env_setup();
T1 = gsl_rng_default;
r = gsl_rng_alloc(T1);
  
//Making Initial Profile 
for(i1=0; i1 < n_x; ++i1)
{
	double u = gsl_rng_uniform(r);
	__real__ comp[i1] = 0.5+ (0.5-u)*1e-4;
	__imag__ comp[i1] = 0.0;
}

half_nx = (int) n_x/2;
delta_kx = (2.0*M_PI)/(n_x*delta_x);

/** Opening a file to write a gnuplot script **/
FILE *gnu;
gnu=fopen("plotAnimation.gp","w");
fprintf(gnu,"set yrange[0:1]\nset xrange[0:%d]\n",n_x);


/** Printing the initial Profile **/
sprintf(name,"./output/c_%d.dat",0);
fpw=fopen(name,"w");
for(i1=0; i1<n_x; ++i1)
{
fprintf(fpw,"%le\n",__real__ comp[i1]);
}
fclose(fpw);
fprintf(gnu,"plot \"%s\" with lines\npause 0.5\n",name);


/** Starting the time loop **/
for(INDEX=1; INDEX<=T; ++INDEX){

	//initialize g
	for(i1=0; i1 < n_x; ++i1){
		g[i1] = 2*A*comp[i1]*(1-comp[i1])*(1-2*comp[i1]);
	}

	/** Let us take comp to the Fourier space **/

	fftw_execute_dft(planF,comp,comp);
	fftw_execute_dft(planFg,g,g);

	/** Evolve composition **/

	for(i1=0; i1 < n_x; ++i1){
		if(i1 < half_nx) kx = i1*delta_kx;
		else kx = (i1-n_x)*delta_kx;
		kx2 = kx*kx;
		comp[i1] = (comp[i1]-g[i1]*delta_t*alpha*kx2)/(1+2*kx2*kx2*delta_t);
	}

	/** Take composition back to real space **/

	fftw_execute_dft(planB,comp,comp);

	for(i1=0; i1<n_x; ++i1){
		comp[i1] = comp[i1]/(n_x);  
	}

	/**Print after a few Time steps**/
	if(INDEX%T_write==0)
	{
		sprintf(name,"./output/c_%d.dat",INDEX);
		fpw=fopen(name,"w");
		for(i1=0; i1<n_x; ++i1)
		{
			fprintf(fpw,"%le\n",__real__ comp[i1]);
			__imag__ comp[i1] = 0.0;
		}
		fclose(fpw);
		fprintf(gnu,"plot \"%s\" with lines\npause 0.5\n",name);
	}
	
}

//Printing the final profile
fprintf(gnu,"set term png\nset output \"finalProfile.png\"\nreplot\nset term x11");

/**Free the memory allocated dynamically**/
gsl_rng_free (r);
fclose(gnu);
fftw_free(comp);
fftw_free(g);
fftw_destroy_plan(planF);
fftw_destroy_plan(planFg);
fftw_destroy_plan(planB);
}
