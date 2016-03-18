/************************************************************
Test function to see if fftw works. Please compile and link with

gcc -LLIBDIR=/usr/local/lib -lgsl -lgslcblas -lfftw3 -lm

Author: M. P. Gururajan

************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>

#include <gsl/gsl_math.h>
#include <fftw3.h>

int main(void) {
FILE *fp;
fp=fopen("input.dat","r");

int i1,i2;
int half_nx, half_ny;
double kx, ky, delta_kx, delta_ky;
double kx2, ky2, k2, k4;
char NAME[50];
char name[50];
int INDEX;

double A = 1.0;
double delta_x = 1.0;
double delta_y = 1.0;
double delta_t = 0.2;

int n_x, n_y;
int T, T_write;


//Input the parameters
fscanf(fp,"%d",&n_x);
fscanf(fp,"%d",&n_y);

fscanf(fp,"%le",&A);
fscanf(fp,"%le",&delta_x);
fscanf(fp,"%le",&delta_y);
fscanf(fp,"%le",&delta_t);
fscanf(fp,"%d",&T);
fscanf(fp,"%d",&T_write);

fclose(fp);

fftw_complex *comp;
fftw_plan planF, planB;

comp = fftw_malloc(n_x*n_y* sizeof(fftw_complex));

planF = fftw_plan_dft_2d(n_x,n_y,comp,comp,FFTW_FORWARD,FFTW_ESTIMATE);
planB = fftw_plan_dft_2d(n_x,n_y,comp,comp,FFTW_BACKWARD,FFTW_ESTIMATE);

//Making Initial Profile
for(i1=0; i1 < n_x; ++i1){
for(i2=0; i2 < n_y; ++i2){
if( ((i1-n_x/2)*(i1-n_x/2) + (i2-n_y/2)*(i2-n_y/2)) < 144){
__real__ comp[i2+n_y*i1] = 1.0;
__imag__ comp[i2+n_y*i1] = 0.0;
}
else
{
__real__ comp[i2+n_y*i1] = 0.1;
__imag__ comp[i2+n_y*i1] = 0.0;
}
}
}

half_nx = (int) n_x/2;
half_ny = (int) n_y/2;

delta_kx = (2.0*M_PI)/(n_x*delta_x);
delta_ky = (2.0*M_PI)/(n_y*delta_y);


/**Creating a file to make the gnuplot script to view all the profiles generated as an animation**/
FILE *gnu;
gnu=fopen("plotAnimation","w");
fprintf(gnu,"set size square");
fprintf(gnu,"\nset xrange[0:%d]",n_x);
fprintf(gnu,"\nset yrange[0:%d]",n_y);
fprintf(gnu,"\nset cbrange [0:1]\nset palette defined ( 0 \"white\", 1 \"black\")");


/**Adding the initial file to the animation**/	
sprintf(name,"./output/c_%d.dat",0);	
		fp=fopen(name,"w");
	
		for(i1=0; i1<n_x; ++i1)
		{
			for(i2=0; i2<n_y; ++i2)
			{
			__imag__ comp[i2+n_y*i1] = 0.0;
			fprintf(fp,"%le	",__real__ comp[i2+n_y*i1]);
			}
			fprintf(fp,"\n");
		}
	fclose(fp);
	//make a entry for the file in the animation
fprintf(gnu,"\nplot \"%s\" matrix with image\npause 1",name);
	

/**Starting the time loop**/
for(INDEX=1; INDEX<=T; ++INDEX){

/** Let us take comp to the Fourier space **/

fftw_execute_dft(planF,comp,comp);

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
	comp[i2+n_y*i1] = comp[i2+n_y*i1]/(1+k2*delta_t);
}}

/** Take composition back to real space **/

fftw_execute_dft(planB,comp,comp);

for(i1=0; i1<n_x; ++i1){
for(i2=0; i2<n_y; ++i2){
comp[i2+n_y*i1] = comp[i2+n_y*i1]/(n_x*n_y);
}}

//print new profile to a file
if(INDEX%T_write==0)//after every T_write time steps
{
	sprintf(name,"./output/c_%d.dat",INDEX);	
		fp=fopen(name,"w");
	
		for(i1=0; i1<n_x; ++i1)
		{
			for(i2=0; i2<n_y; ++i2)
			{
			__imag__ comp[i2+n_y*i1] = 0.0;
			fprintf(fp,"%le	",__real__ comp[i2+n_y*i1]);
			}
			fprintf(fp,"\n");
		}
	fclose(fp);
	//make a entry for the file in the animation
	fprintf(gnu,"\nplot \"%s\" matrix with image\npause 0.5",name);
	
}

}

fclose(gnu);
fftw_free(comp);
fftw_destroy_plan(planF);
fftw_destroy_plan(planB);
}
