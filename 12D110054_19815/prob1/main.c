/************************************************************

Compile with
gcc main.c -lgsl -lgslcblas -lfftw3 -lm

or just use ./run.sh to compile,run and plot the result 

************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <fftw3.h>

int main(void) {
FILE *fpw;

int i1,i2;
int half_nx;
double kx, delta_kx;
double kx2,k2, k4;
char name[50];
int INDEX;

double A = 1.0;
double delta_x = 1.0;
double delta_t = 0.2;
int T,T_write;//(T)number of time steps
//T_write is after how many time steps we write profile in file
int n_x;//number of nodes in the system 

n_x = 128;

/**taking the Input from file**/
fpw=fopen("input.dat","r");
fscanf(fpw,"%d",&n_x);
fscanf(fpw,"%le",&delta_x);
fscanf(fpw,"%d",&T);
fscanf(fpw,"%le",&delta_t);
fscanf(fpw,"%d",&T_write);
fclose(fpw);

fftw_complex *comp;
fftw_plan planF, planB;

comp = fftw_malloc(n_x*sizeof(fftw_complex));

planF = fftw_plan_dft_1d(n_x,comp,comp,FFTW_FORWARD,FFTW_ESTIMATE);
planB = fftw_plan_dft_1d(n_x,comp,comp,FFTW_BACKWARD,FFTW_ESTIMATE);
 
//Creating Input Profile
for(i1=0; i1 < n_x; ++i1){
	if( i1>n_x/4 && i1<=(3*n_x/4) ){
	__real__ comp[i1] = 1.0;
	__imag__ comp[i1] = 0.0;
	}
	else
	{
	__real__ comp[i1] = 0.1;
	__imag__ comp[i1] = 0.0;
	}
}

/**Creating a file to make the gnuplot script to view all the profiles generated as an animation**/
FILE *gnu;
gnu=fopen("plotAnimation.gp","w");
fprintf(gnu,"\nset xrange[0:%d]",n_x);
fprintf(gnu,"\nset yrange[0:1]");

/**Adding the initial file to the animation**/	
sprintf(name,"./output/c_%d.dat",0);						
fpw=fopen(name,"w");
for(i1=0; i1<n_x; ++i1)
{
	fprintf(fpw,"%le\n",__real__ comp[i1]);
	__imag__ comp[i1] = 0.0;//Setting the imaginary part to zero every few time steps
}
fprintf(gnu,"\nplot \"%s\" with lines\npause 1\n",name);
fclose(fpw);

half_nx = (int) n_x/2;
delta_kx = (2.0*M_PI)/(n_x*delta_x);

for(INDEX=0; INDEX<T; ++INDEX)
{

	/** Let us take comp to the Fourier space **/

	fftw_execute_dft(planF,comp,comp);

	/** Evolve composition **/

	for(i1=0; i1 < n_x; ++i1){
		if(i1 < half_nx) kx = i1*delta_kx;
		else kx = (i1-n_x)*delta_kx;
		kx2 = kx*kx;
		comp[i1] = comp[i1]/(1+kx2*delta_t);
	}

	/** Take composition back to real space **/

	fftw_execute_dft(planB,comp,comp);

	for(i1=0; i1<n_x; ++i1)
	{
	comp[i1] = comp[i1]/(n_x);  
	}

	if(INDEX % T_write==0)
	{
		//Print The profile to a file
		sprintf(name,"./output/c_%d.dat",INDEX);						
		fpw=fopen(name,"w");
		for(i1=0; i1<n_x; ++i1)
		{
		fprintf(fpw,"%le\n",__real__ comp[i1]);
		__imag__ comp[i1] = 0.0;//Setting the imaginary part to zero a few time steps
		}
		fprintf(gnu,"\nplot \"%s\" with lines\npause 0.3\n",name);
		fclose(fpw);	

	}
}
//Printing the final profile
fprintf(gnu,"set term png\nset output \"finalProfile.png\"\nreplot\nset term x11");


/**Free the memory allocated dynamically**/
fclose(gnu);
fftw_free(comp);
fftw_destroy_plan(planF);
fftw_destroy_plan(planB);
}
