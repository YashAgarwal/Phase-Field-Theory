/************************************************************
gcc main.c -lgsl -lgslcblas -lfftw3 -lm

************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <fftw3.h>

int main(void) {
FILE *fpw,*gnu;

int i1,i2;
int half_nx, half_ny;
double kx, ky, delta_kx, delta_ky;
double kx2, ky2, k2, k4;
char name[50];
int INDEX;

double alpha = 1.0;
double A = 1.0;
double delta_x = 1.0;
double delta_y = 1.0;
double delta_t = 0.2;

int n_x, n_y;

int T,T_write;

n_x = n_y = 128;


//Taking input from file
fpw=fopen("input.dat","r");

fscanf(fpw,"%d",&n_x);
fscanf(fpw,"%le",&delta_x);
fscanf(fpw,"%d",&n_y);
fscanf(fpw,"%le",&delta_y);
fscanf(fpw,"%d",&T);
fscanf(fpw,"%le",&delta_t);
fscanf(fpw,"%d",&T_write);//after how many time steps we are writing the output
fscanf(fpw,"%le",&alpha);

fclose(fpw);

fftw_complex *comp,*h;
fftw_plan planF,planFh, planB;

comp = fftw_malloc(n_x*n_y* sizeof(fftw_complex));

h = fftw_malloc(n_x*n_y* sizeof(fftw_complex));

planF = fftw_plan_dft_2d(n_x,n_y,comp,comp,FFTW_FORWARD,FFTW_ESTIMATE);
planFh = fftw_plan_dft_2d(n_x,n_y,comp,comp,FFTW_FORWARD,FFTW_ESTIMATE);

planB = fftw_plan_dft_2d(n_x,n_y,comp,comp,FFTW_BACKWARD,FFTW_ESTIMATE);

//making initial profile
for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
if( ((i1-n_x/2)*(i1-n_x/2) + (i2-n_y/2)*(i2-n_y/2)) < 225){
__real__ comp[i2+n_y*i1] = 0.0;
__imag__ comp[i2+n_y*i1] = 0.0;
}
else
{
__real__ comp[i2+n_y*i1] = 1.0;
__imag__ comp[i2+n_y*i1] = 0.0;
}
}
}

half_nx = (int) n_x/2;
half_ny = (int) n_y/2;

delta_kx = (2.0*M_PI)/(n_x*delta_x);
delta_ky = (2.0*M_PI)/(n_y*delta_y);

//making the gnu script to see the output as an animation
gnu=fopen("plotAnimation.gp","w");
fprintf(gnu,"set cbrange[0:1]\n");
fprintf(gnu,"set xrange[0:%d]\n",n_x);
fprintf(gnu,"set yrange[0:%d]\n",n_y);

//Time loop starts here
for(INDEX=0; INDEX<T; ++INDEX){

/** Let us take comp to the Fourier space **/
for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
h[i2+n_y*i1]=2*A*comp[i2+n_y*i1]*(1-comp[i2+n_y*i1])*(1-2*comp[i2+n_y*i1]);
}
}
	
	fftw_execute_dft(planF,comp,comp);
	fftw_execute_dft(planFh,h,h);
	
	/** Evolve composition **/

	for(i1=0; i1 < n_x; ++i1){
		if(i1 < half_nx) kx = i1*delta_kx;
		else kx = (i1-n_x)*delta_kx;
		kx2 = kx*kx;
	for(i2=0; i2 < n_y; ++i2){
		if(i2 < half_ny) ky = i2*delta_ky;
		else ky = (i2-n_y)*delta_ky;
		ky2 = ky*ky;
		k2 = kx2 + ky2;
		//main equation
		comp[i2+n_y*i1] = (comp[i2+n_y*i1]-h[i2+n_y*i1]*delta_t)/(1+k2*delta_t*alpha);
	}}

	/** Take composition back to real space **/

	fftw_execute_dft(planB,comp,comp);

	for(i1=0; i1<n_x; ++i1){
	for(i2=0; i2<n_y; ++i2){
	comp[i2+n_y*i1] = comp[i2+n_y*i1]/(n_x*n_y);  
	}}

if(INDEX%T_write==0)
{
	sprintf(name,"./output/c_%d.dat",INDEX);
	fpw=fopen(name,"w");
	for(i1=0; i1<n_x; ++i1){
	for(i2=0; i2<n_y; ++i2){
	fprintf(fpw,"%le ",__real__ comp[i2+n_y*i1]);
	__imag__ comp[i2+n_y*i1]=0;
	}
	fprintf(fpw,"\n");
	}
	fprintf(gnu,"plot \"%s\" matrix with image\npause 0.3\n",name);
}
}

fftw_free(comp);
fftw_free(h);
fftw_destroy_plan(planF);
fftw_destroy_plan(planFh);
fftw_destroy_plan(planB);

}
