/************************************************************
gcc analysis.c -lgsl -lgslcblas -lfftw3 -lm
Analysis of the main data
just run ./runAnalysis.sh to see the final Analysis Report
************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <fftw3.h>

int main(void) {
FILE *fpw1,*fpw2,*fpw3,*gnu;
int i1,j;
char name[50];

double delta_x = 1.0;
double delta_t = 0.2;

int n_x=128;
int T=1000;
int T_write=20;

/**Input the parameters**/
fpw1=fopen("input.dat","r");
fscanf(fpw1,"%d",&n_x);
fscanf(fpw1,"%le",&delta_x);
fscanf(fpw1,"%d",&T);
fscanf(fpw1,"%le",&delta_t);
fscanf(fpw1,"%d",&T_write);

fclose(fpw1);


/**finding the length of the ppt**/
/**declare a array to store the length values**/
//no. of values will come by T/T_write + 1
int numSample=T/T_write +1;
double *pptLength = malloc(numSample* sizeof(double));
double *interfacePosition = malloc(numSample* sizeof(double));
double *velocity = malloc((numSample-1)* sizeof(double));

/**declare a array to store the comp values**/
double *comp = malloc(n_x* sizeof(double));

int filenum=0;//counter to keep the number of file

//Read each file and find the lenght of the precipitate
for(j=0; j<numSample; j++){

	/**Read the output file one by one**/
	sprintf(name,"./output/c_%d.dat",j*T_write);
	fpw1=fopen(name,"r");
	for(i1=0; i1<n_x; ++i1){
	fscanf(fpw1,"%le",&comp[i1]);
	}
	fclose(fpw1);
	/**now analyse the data**/
	/**need to find the interface in the first half of the data**/
	int pos;
	for(i1=0; i1<n_x/2; ++i1){
			if(comp[i1]>=0.5){
				pos=i1;
				break;
			}
		}

	/**now find the point just below the 0.5 line and join the 2 points to form a line**/
	/**take the intersection of that line with the 0.5 line**/
	double m=(comp[pos]-comp[pos-1]);
	interfacePosition[j]=delta_x*(0.5+m*pos-comp[pos])/m;
	pptLength[j]=delta_x*(n_x - 2*(0.5+m*pos-comp[pos])/m);

	//now find the velocity of the interface
  if(j)
	{
	velocity[j-1]=fabs((interfacePosition[j]-interfacePosition[j-1])/T_write);
	}
}

//Write the value of precipitate length, interface position and interface velocity with time in the data files
fpw1=fopen("./interfacePosition.dat","w");
fpw2=fopen("./pptLength.dat","w");
fpw3=fopen("./velocity.dat","w");

for(i1=0; i1<numSample; ++i1){
fprintf(fpw1,"%le\t%d\n",interfacePosition[i1],i1*T_write);
fprintf(fpw2,"%le\t%d\n",pptLength[i1],i1*T_write);
if(i1)
	{
		fprintf(fpw3,"%le\t%d\n",velocity[i1-1],i1*T_write);
	}
}

fclose(fpw1);
fclose(fpw2);
fclose(fpw3);

//making the gnu script to plot the final data in a more readable format
gnu=fopen("plotAnalysis.gp","w");

fprintf(gnu,"set xtics 0,%d,%d\nset tics font \"Times-Roman,8\"\nset multiplot layout 2, 2 title \"Analysis Report\" font \"Times-Roman,14\"\nset tmargin 2\nset title \"Interface Position Vs Time\"\nunset key\nplot \"./interfacePosition.dat\" using 2:1 with line\nset title \"Precipitate length Vs Time\"\nunset key\nplot \"./pptLength.dat\" using 2:1 with line\nset title \"Interface Velocity Vs Time\"\nunset key\nplot \"./velocity.dat\" using 2:1 with line\nunset multiplot",T/4,T);

fclose(gnu);
}
