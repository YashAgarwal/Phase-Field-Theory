/************************************************************
gcc main.c -lgsl -lgslcblas -lfftw3 -lm
just use ./run.sh to see the complete animation
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
char name[50];
int INDEX;
int half_nx;
double kx, delta_kx, delta_ky;
double kx2,k2, k4;


double A = 1.0;
double delta_x = 1.0;
double delta_t = 0.2;
double alpha=1,beta=1;

int n_x;
int T=1000;
int T_write=20;

n_x= 128;
int flag=0;
/**Input the parameters**/
fpw=fopen("input.dat","r");
fscanf(fpw,"%d",&flag);//to check for the case
fscanf(fpw,"%d",&n_x);
fscanf(fpw,"%le",&delta_x);
fscanf(fpw,"%d",&T);
fscanf(fpw,"%d",&T_write);
fscanf(fpw,"%le",&delta_t);
fclose(fpw);

//checking which case is used and assigning the value to A and beta
if(flag==1)
{
	A=1;
	beta=1;
}
else if(flag==2)
{
	A=1;
	beta=4;
}
else if(flag==3)
{
	A=4;
	beta=1;
}
else
{
	printf("\nInvalid input\nExiting program......");
	exit(0);
}

fftw_complex *comp,*compDel,*g;
fftw_plan planF,planFg, planB;
fftw_plan planF_Del,planB_Del;

comp = fftw_malloc(n_x* sizeof(fftw_complex));
compDel = fftw_malloc(n_x* sizeof(fftw_complex));
g = fftw_malloc(n_x* sizeof(fftw_complex));

planF = fftw_plan_dft_1d(n_x,comp,comp,FFTW_FORWARD,FFTW_ESTIMATE);
planF_Del = fftw_plan_dft_1d(n_x,compDel,compDel,FFTW_FORWARD,FFTW_ESTIMATE);

planFg = fftw_plan_dft_1d(n_x,g,g,FFTW_FORWARD,FFTW_ESTIMATE);

planB = fftw_plan_dft_1d(n_x,comp,comp,FFTW_BACKWARD,FFTW_ESTIMATE);
planB_Del = fftw_plan_dft_1d(n_x,compDel,compDel,FFTW_BACKWARD,FFTW_ESTIMATE);

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
		//main equation implementation take place here
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

sprintf(name,"output%d.dat",flag);
fpw=fopen(name,"w");
/**calculating the interface width for the final profile**/

//first find the point where the composition value just increases 0.5
int pos=1;
for(i1=0; i1<n_x/2; ++i1){
			if(__real__ comp[i1]>=0.5){
				pos=i1;
				break;
			}
		}


//now find the slope by taking a point 10 points ahead of this and 10 points before this
double m=(__real__ comp[pos]-__real__ comp[pos-1]);
//the interface width will be 1/slope
double interfaceWidth=1/m;

fprintf(fpw,"Interface Width = %le\n",interfaceWidth);

/**calculating the interfacial energy**/
//fisrt find dc/dx

for(i1=0; i1 < n_x; ++i1){
	compDel[i1] = comp[i1];
	}

/** Let us take comp to the Fourier space **/
fftw_execute_dft(planF_Del,compDel,compDel);

/** Evolve composition **/
for(i1=0; i1 < n_x; ++i1){
	if(i1 < half_nx) kx = i1*delta_kx;
	else kx = (i1-n_x)*delta_kx;
	kx2 = kx*kx;

	//main equation implementation take place here
	__real__ compDel[i1] = -kx * __imag__ compDel[i1];
	__imag__ compDel[i1] = __real__ compDel[i1] * kx;
	}

/** Take composition back to real space **/
fftw_execute_dft(planB_Del,compDel,compDel);

for(i1=0; i1<n_x; ++i1){
	compDel[i1] = compDel[i1]/(n_x);
	}

double energy=0;
for(i1=0; i1<n_x; ++i1){
		energy+=(A*comp[i1]*comp[i1]*(1-comp[i1])*(1-comp[i1]) + beta*compDel[i1]*compDel[i1])*delta_x;
		}

fprintf(fpw,"Interfacial Energy = %le\n",energy);
fclose(fpw);

/**Freeing the dynamically allocated memory**/
fftw_destroy_plan(planFg);
fftw_free(comp);
fftw_free(g);
fftw_destroy_plan(planF);
fftw_destroy_plan(planB);
fclose(gnu);
}
