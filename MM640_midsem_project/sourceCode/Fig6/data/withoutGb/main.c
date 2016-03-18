/************************************************************
gcc main.c -lgsl -lgslcblas -lfftw3 -lm
************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <fftw3.h>
#include <string.h>
#include <sys/dir.h>

int main(void){

printf("\nComputing the profile for the input parameters....");

FILE *fpw,*gnu;
double max=0, min=1000;
int i,j;
int half_nx, half_ny;
double kx, ky, delta_kx, delta_ky,delc;
double kx2, ky2, k2, k4;
char name[50];
int INDEX;

double A = 1.0,kc= 1.0,M= 1.0,E= 2.0,kn= 1.0,L= 1.0;

double delta_x = 1.0;
double delta_y = 1.0;
double delta_t = 0.2;
int T,T_write;
int n_x, n_y;
n_x = 512;
n_y = 256;

double systemNum=1.1;//by default set to system Ia
//the 4 systems are defined as
//Ia=1.1
//Ib=1.2
//Ic=1.3
//IIa=2.1

/**Take input from file**/
fpw=fopen("input.dat","r");
fscanf(fpw,"%d",&n_x);
fscanf(fpw,"%le",&delta_x);
fscanf(fpw,"%d",&n_y);
fscanf(fpw,"%le",&delta_y);
fscanf(fpw,"%le",&delta_t);
fscanf(fpw,"%d",&T);//number of time steps
fscanf(fpw,"%le",&delc);//initial profile disturbance
fscanf(fpw,"%lf",&systemNum);//system number

fclose(fpw);

T_write=1;

//defining kn value based on system number given 
if(systemNum==1.1)
{
	kn=1.0;
}
else if(systemNum==1.2)
{
	kn=1.0;
}
else if(systemNum==1.3)
{
	kn=0.5;
}
else if(systemNum==2.1)
{
	kn=1.0;
}
else 
{
	printf("\nInvalid System Number\nAborting Program....");
	exit(0);
}

/**Allocating memory to the pointers for FFTW**/
fftw_complex *c,*g;
fftw_plan planF_c, planB_c, planF_g;

c = fftw_malloc(n_x*n_y* sizeof(fftw_complex));
g = fftw_malloc(n_x*n_y* sizeof(fftw_complex));

planF_c = fftw_plan_dft_2d(n_x,n_y,c,c,FFTW_FORWARD,FFTW_ESTIMATE);
planB_c = fftw_plan_dft_2d(n_x,n_y,c,c,FFTW_BACKWARD,FFTW_ESTIMATE);

planF_g = fftw_plan_dft_2d(n_x,n_y,g,g,FFTW_FORWARD,FFTW_ESTIMATE);

fftw_complex *n1,*h1;
fftw_plan planF_n1, planB_n1, planF_h1;

n1 = fftw_malloc(n_x*n_y* sizeof(fftw_complex));
h1 = fftw_malloc(n_x*n_y* sizeof(fftw_complex));

planF_n1 = fftw_plan_dft_2d(n_x,n_y,n1,n1,FFTW_FORWARD,FFTW_ESTIMATE);
planB_n1 = fftw_plan_dft_2d(n_x,n_y,n1,n1,FFTW_BACKWARD,FFTW_ESTIMATE);

planF_h1 = fftw_plan_dft_2d(n_x,n_y,h1,h1,FFTW_FORWARD,FFTW_ESTIMATE);

fftw_complex *n2,*h2;
fftw_plan planF_n2, planB_n2, planF_h2;

n2 = fftw_malloc(n_x*n_y* sizeof(fftw_complex));
h2 = fftw_malloc(n_x*n_y* sizeof(fftw_complex));

planF_n2 = fftw_plan_dft_2d(n_x,n_y,n2,n2,FFTW_FORWARD,FFTW_ESTIMATE);
planB_n2 = fftw_plan_dft_2d(n_x,n_y,n2,n2,FFTW_BACKWARD,FFTW_ESTIMATE);

planF_h2 = fftw_plan_dft_2d(n_x,n_y,h2,h2,FFTW_FORWARD,FFTW_ESTIMATE);

/**Setting up Random Number Generation**/
const gsl_rng_type * R;
gsl_rng * r;
gsl_rng_env_setup();
R = gsl_rng_default;
r = gsl_rng_alloc (R);
double u;

half_nx = (int) n_x/2;
half_ny = (int) n_y/2;

/**Making the initial input**/
for(i=0; i < n_x; ++i)
{
	for(j=0; j < n_y; ++j)
	{
			u=gsl_rng_uniform(r);
			__real__ c[j+n_y*i] = 0.5+((u*2)-1)*delc;
			__imag__ c[j+n_y*i] = 0.0;

			//Setting no GB
			__real__ n1[j+n_y*i] = 1.0;
			__imag__ n1[j+n_y*i] = 0.0;

			__real__ n2[j+n_y*i] = 0.0;
			__imag__ n2[j+n_y*i] = 0.0;			
	}
}

delta_kx = (2.0*M_PI)/(n_x*delta_x);
delta_ky = (2.0*M_PI)/(n_y*delta_y);

/**Making the gnuplot script file**/

char inputFile[50];

if(systemNum==1.1)
{
  sprintf(inputFile,"Ia_%.2lf",delc);
}
else if(systemNum==1.2)
{
  sprintf(inputFile,"Ib_%.2lf",delc);
}
else if(systemNum==1.3)
{
  sprintf(inputFile,"Ic_%.2lf",delc);
}
else if(systemNum==2.1)
{
  sprintf(inputFile,"IIa_%.2lf",delc);
}
else 
{
	printf("\nInvalid System Number\nAborting Program....");
	exit(0);
}

sprintf(name,"./output_%s",inputFile);
mkdir(name);

char *res1 = malloc(strlen("plotAnimation_") + strlen(inputFile) + strlen(".gp")+1);//+1 for the zero-terminator
strcpy(res1,"plotAnimation_");
strcat(res1,inputFile);
strcat(res1,".gp");

gnu=fopen(res1,"w");
fprintf(gnu,"set cbrange[0:1]\n");
fprintf(gnu,"set xrange[0:%d]\n",n_y);
fprintf(gnu,"set yrange[0:%d]\n",n_x);
fprintf(gnu,"set size ratio 2\n");
fprintf(gnu,"set palette gray\n");
fprintf(gnu,"unset colorbox\n");
fprintf(gnu,"unset key; unset tics; unset border\n");
free(res1);

/**Starting time loop**/
for(INDEX=1; INDEX<=T; ++INDEX){

	/** Let us take comp to the Fourier space **/
	for(i=0; i < n_x; ++i)
	{
		for(j=0; j < n_y; ++j)
		{
			fftw_complex V= E * n1[j+n_y*i] * n1[j+n_y*i] * n2[j+n_y*i] * n2[j+n_y*i];
			fftw_complex W= - (n1[j+n_y*i]*n1[j+n_y*i]/2.0) + (pow(n1[j+n_y*i],4)/4) - (n2[j+n_y*i]*n2[j+n_y*i]/2) + (pow(n2[j+n_y*i],4)/4);
			fftw_complex Z=(0.25 + W + V);
			fftw_complex m;
			if(systemNum==1.1)
			{
				m=(1.0 + 0.5*c[j+n_y*i]*c[j+n_y*i]);
			}
			else if(systemNum==1.2)
			{
				m=1.0 + 0.5*c[j+n_y*i]*c[j+n_y*i] - 2.5*c[j+n_y*i]*c[j+n_y*i]*(1-c[j+n_y*i])*(1-c[j+n_y*i]);
			}
			else if(systemNum==1.3)
			{
				m=(2.0 + c[j+n_y*i]*c[j+n_y*i]);
			}
			else if(systemNum==2.1)
			{
				m=(1.0 + 0.1*c[j+n_y*i]*c[j+n_y*i]);
			}
			else 
			{
				printf("\nInvalid System Number\nAborting Program....");
				exit(0);
			}
			
			fftw_complex Y=c[j+n_y*i]*Z;
			fftw_complex X=2.0*A*c[j+n_y*i]*(1-c[j+n_y*i])*(1.0-2*c[j+n_y*i]);
			g[j+n_y*i]=X+Y;

			h1[j+n_y*i]=m*(-n1[j+n_y*i]+n1[j+n_y*i]*n1[j+n_y*i]+2*E*n2[j+n_y*i]*n2[j+n_y*i]*n1[j+n_y*i]);
			h2[j+n_y*i]=m*(-n2[j+n_y*i]+n2[j+n_y*i]*n2[j+n_y*i]+2*E*n1[j+n_y*i]*n1[j+n_y*i]*n2[j+n_y*i]);
		}
	}
	fftw_execute_dft(planF_c,c,c);
	fftw_execute_dft(planF_g,g,g);
	fftw_execute_dft(planF_n1,n1,n1);
	fftw_execute_dft(planF_h1,h1,h1);

	fftw_execute_dft(planF_n2,n2,n2);
	fftw_execute_dft(planF_h2,h2,h2);

	/** Evolve composition **/

	for(i=0; i < n_x; ++i)
	{
		if(i < half_nx) kx = i*delta_kx;
		else kx = (i-n_x)*delta_kx;
		kx2 = kx*kx;
		for(j=0; j < n_y; ++j)
		{
			if(j < half_ny) ky = j*delta_ky;
			else ky = (j-n_y)*delta_ky;
			ky2 = ky*ky;
			k2 = kx2 + ky2;
			k4=k2*k2;
			c[j+n_y*i] = (c[j+n_y*i]-M*k2*delta_t*g[j+n_y*i])/(1.0+2.0*M*kc*k4*delta_t);
			n1[j+n_y*i] = (n1[j+n_y*i]-L*delta_t*h1[j+n_y*i])/(1.0+2.0*L*kn*k2*delta_t);
			n2[j+n_y*i] = (n2[j+n_y*i]-L*delta_t*h2[j+n_y*i])/(1.0+2.0*L*kn*k2*delta_t);
		}
	}

	/** Take composition back to real space **/

	fftw_execute_dft(planB_c,c,c);
	fftw_execute_dft(planB_n1,n1,n1);
	fftw_execute_dft(planB_n2,n2,n2);

	for(i=0; i<n_x; ++i){
	for(j=0; j<n_y; ++j){
		c[j+n_y*i] = c[j+n_y*i]/(n_x*n_y);
		n1[j+n_y*i] = n1[j+n_y*i]/(n_x*n_y);
		n2[j+n_y*i] = n2[j+n_y*i]/(n_x*n_y);
	}}

	/**Print the composition after every T_write time steps**/
	
		if(INDEX%T_write==0)
		{
		sprintf(name,"./output_%s/c_%d.dat",inputFile,INDEX);
		fpw=fopen(name,"w");
		for(i=0; i<n_x; ++i)
		{
			for(j=0; j<n_y; ++j)
			{
				fprintf(fpw,"%le ",__real__ c[j+n_y*i]);
				__imag__ c[j+n_y*i] = 0.0;
				__imag__ n1[j+n_y*i] = 0.0;
				__imag__ n2[j+n_y*i] = 0.0;
				//if(max<__real__ c[j+n_y*i])max=__real__ c[j+n_y*i];
				//if(min>__real__ c[j+n_y*i])min=__real__ c[j+n_y*i];
			}
			fprintf(fpw,"\n");
		}
		fclose(fpw);
		fprintf(gnu,"set title \"t=%d\"\nplot \"%s\" matrix with image\npause 0.1\n",INDEX,name);
		}
	

	}
	fclose(gnu);


//printing final profile
sprintf(name,"finalProfile_%s.dat",inputFile);
fpw=fopen(name,"w");
for(i=0; i<n_x; ++i)
{
	for(j=0; j<n_y; ++j)
	{
		fprintf(fpw,"%le ",__real__ c[j+n_y*i]);
		__imag__ c[j+n_y*i] = 0.0;
		__imag__ n1[j+n_y*i] = 0.0;
		__imag__ n2[j+n_y*i] = 0.0;
		//if(max<__real__ c[j+n_y*i])max=__real__ c[j+n_y*i];
		//if(min>__real__ c[j+n_y*i])min=__real__ c[j+n_y*i];
	}
	fprintf(fpw,"\n");
}
fclose(fpw);


//print plotFig.gp
gnu=fopen("plotFig.gp","w");

fprintf(gnu,"set cbrange[0:1]\n");
fprintf(gnu,"set xrange[0:%d]\n",n_y);
fprintf(gnu,"set yrange[0:%d]\n",n_x);
fprintf(gnu,"set size ratio 2\n");
fprintf(gnu,"set palette gray\n");
fprintf(gnu,"unset colorbox\n");
fprintf(gnu,"unset key; unset tics; unset border\n");
fprintf(gnu,"set title \"Final Profile\"\nset palette defined (0 'black',1 'white')\nset cbrange[0:1]\n");

fprintf(gnu,"plot \"%s\" matrix with image\n",name);

fclose(gnu);

//printf("\n%le:%le\n",min,max);
fftw_free(c);
fftw_free(g);

fftw_free(n1);
fftw_free(h1);

fftw_free(n2);
fftw_free(h2);

fftw_destroy_plan(planF_g);
fftw_destroy_plan(planF_h1);
fftw_destroy_plan(planF_h2);

fftw_destroy_plan(planF_c);
fftw_destroy_plan(planB_c);

fftw_destroy_plan(planF_n1);
fftw_destroy_plan(planB_n1);

fftw_destroy_plan(planF_n2);
fftw_destroy_plan(planB_n2);

gsl_rng_free (r);

printf("\nComputation finished\n\n");

}
