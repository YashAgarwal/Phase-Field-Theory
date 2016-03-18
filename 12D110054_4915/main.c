/************************************************************
gcc main.c -lgsl -lgslcblas -lfftw3 -lm
just run ./runAnalysis.sh to see the complete animation
************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <fftw3.h>

int main(void) {
FILE *fpw,*gnu;
int i1;
int half_nx;
double kx, delta_kx, delta_ky;
double kx2,k2, k4;
char name[50];
int INDEX;

double A = 1.0;
double delta_x = 1.0;
double delta_t = 0.2;
double alpha=1,beta=1;

int n_x;
int T=1000;
int T_write=20;

n_x= 128;

/**Input the parameters**/
fpw=fopen("input.dat","r");
fscanf(fpw,"%d",&n_x);
fscanf(fpw,"%le",&delta_x);
fscanf(fpw,"%d",&T);
fscanf(fpw,"%le",&delta_t);
fscanf(fpw,"%d",&T_write);
fscanf(fpw,"%le",&A);
fscanf(fpw,"%le",&alpha);
fscanf(fpw,"%le",&beta);

fclose(fpw);

fftw_complex *comp,*g;
fftw_plan planF,planFg, planB;

comp = fftw_malloc(n_x* sizeof(fftw_complex));
g = fftw_malloc(n_x* sizeof(fftw_complex));

planF = fftw_plan_dft_1d(n_x,comp,comp,FFTW_FORWARD,FFTW_ESTIMATE);
planFg = fftw_plan_dft_1d(n_x,g,g,FFTW_FORWARD,FFTW_ESTIMATE);

planB = fftw_plan_dft_1d(n_x,comp,comp,FFTW_BACKWARD,FFTW_ESTIMATE);

/**Initial Profile**/
for(i1=0; i1 < n_x; ++i1){
	if(i1<(n_x/4) || i1>(3*n_x/4))
	{
		__real__ comp[i1] = 0.1;
	}
	else
	{
		__real__ comp[i1] = 1;
	}
	__imag__ comp[i1] = 0.0;

}

half_nx = (int) n_x/2;
delta_kx = (2.0*M_PI)/(n_x*delta_x);

/**Make a gnu script**/

gnu=fopen("plotAnimation.gp","w");
fprintf(gnu,"set xrange[0:%d]",n_x);
fprintf(gnu,"\nset yrange[-0.1:1.1]");

/**printing the initial profile**/

sprintf(name,"./output/c_%d.dat",0);
fpw=fopen(name,"w");
for(i1=0; i1<n_x; ++i1){
	fprintf(fpw,"%le\n",__real__ comp[i1]);
	__imag__ comp[i1]=0;
}
fclose(fpw);
fprintf(gnu,"\nplot \"%s\" with lines",name);
fprintf(gnu,"\npause 1");

/**Starting the time loop**/
for(INDEX=1; INDEX<=T; ++INDEX){

	/** calculating g **/
	for(i1=0; i1 < n_x; ++i1){
		g[i1] = 2*A*(comp[i1])*(1-comp[i1])*(1-2*comp[i1]);
	}

	/** Let us take comp to the Fourier space **/
	fftw_execute_dft(planF,comp,comp);
	fftw_execute_dft(planFg,g,g);

	/** Evolve composition **/
	for(i1=0; i1 < n_x; ++i1){
		if(i1 < half_nx) kx = i1*delta_kx;
		else kx = (i1-n_x)*delta_kx;
		kx2 = kx*kx;
		k2 = kx2;
		k4= k2*k2;
		comp[i1] = (comp[i1]-alpha*k2*delta_t*g[i1])/(1+2*beta*k4*delta_t);
	}

	/** Take composition back to real space **/
	fftw_execute_dft(planB,comp,comp);

	for(i1=0; i1<n_x; ++i1){
		comp[i1] = comp[i1]/(n_x);
	}

	/**Printing the Results**/
	if(INDEX%T_write==0)
	{
		sprintf(name,"./output/c_%d.dat",INDEX);
		fpw=fopen(name,"w");
		for(i1=0; i1<n_x; ++i1){
		fprintf(fpw,"%le\n",__real__ comp[i1]);
		__imag__ comp[i1]=0;
		}
		fclose(fpw);
		fprintf(gnu,"\nplot \"%s\" with lines",name);
		fprintf(gnu,"\npause 0.4");
	}
}

//Freeing the dynamically allocated memory
fftw_free(comp);
fftw_free(g);
fftw_destroy_plan(planFg);
fftw_destroy_plan(planF);
fftw_destroy_plan(planB);
fclose(gnu);
}
