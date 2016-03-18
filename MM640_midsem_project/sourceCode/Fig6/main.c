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

int main(void)
{

printf("\nStarting the Analysis of the data....\n");

int i,j;
int n_x=512,n_y=256;

FILE *fp;

double *c;
c = calloc(n_x*n_y,sizeof(double));

int k=0;
int t=0;//counter for f 
char name[50];
	
double f[6][200];

int spoint=253; //where the first band starts
//all band widths are 6

int bandNum;
double Vb;
//calculation of fD with GB

for(t=1;t<=200;t++)//going through time steps 1-200
{
		sprintf(name,"./data/withGb/output_Ia_0.00/c_%d.dat",t);
		
		if((fp=fopen(name,"r"))==NULL)
		{
			printf("Can't open file\nAborting the Program......\n");
			exit(0);
		}

		for(i=0; i<n_x; ++i)
		{
			for(j=0; j<n_y; ++j)
			{
				fscanf(fp,"%le",&c[j+n_y*i]);				
			}
		}

		fclose(fp);
		
		//calculate fD for all bands 1 to 6
		//for time t
		//store it in array

		
		for(bandNum=0;bandNum<6;bandNum++)
		{
			f[bandNum][t-1]=0;
			
			for(i=spoint+bandNum*6; i<=spoint+(bandNum+1)*6; ++i)
			{
				for(j=0; j<n_y; ++j)
				{
					f[bandNum][t-1] += (c[j+n_y*i]-0.5)*(c[j+n_y*i]-0.5);				
				}
			}
			
			Vb= 6*n_y;
			f[bandNum][t-1]=f[bandNum][t-1]/(Vb*0.5*0.5);
			//both delta_x and delta_y are 1, hence not written in the multiplication		
		}
}

//write array to file
for(bandNum=0;bandNum<6;bandNum++){
	
		sprintf(name,"./fD_%d.dat",bandNum+1);

		if((fp=fopen(name,"w"))==NULL)
		{
			printf("Can't open file\nAborting the Program......\n");
			exit(0);
		}

		for(i=0; i<200; ++i)
		{
			fprintf(fp,"%le\n",f[bandNum][i]);
		}

		fclose(fp);
}

//calculation of fD without GB

double SS1[200];
double SS4[200];

for(t=1;t<=200;t++)//going through time steps 1-200
{
		
		//for delc=0.01
		sprintf(name,"./data/withoutGb/output_Ia_0.01/c_%d.dat",t);

		if((fp=fopen(name,"r"))==NULL)
		{
			printf("Can't open file\nAborting the Program......\n");
			exit(0);
		}

		for(i=0; i<n_x; ++i)
		{
			for(j=0; j<n_y; ++j)
			{
				fscanf(fp,"%le",&c[j+n_y*i]);				
			}
		}

		fclose(fp);
		
		//calculate fD for all bands 1 to 6
		//for time t
		//store it in array

		SS1[t-1]=0;
			
		for(i=0; i<n_x; ++i)
		{
			for(j=0; j<n_y; ++j)
			{
				SS1[t-1] += (c[j+n_y*i]-0.5)*(c[j+n_y*i]-0.5);				
			}
		}
		
		Vb= n_x*n_y;
		SS1[t-1]=SS1[t-1]/(Vb*0.5*0.5);
		//both delta_x and delta_y are 1, hence not written in the multiplication		
		
		//for delc=0.04
		sprintf(name,"./data/withoutGb/output_Ia_0.04/c_%d.dat",t);

		if((fp=fopen(name,"r"))==NULL)
		{
			printf("Can't open file\nAborting the Program......\n");
			exit(0);
		}

		for(i=0; i<n_x; ++i)
		{
			for(j=0; j<n_y; ++j)
			{
				fscanf(fp,"%le",&c[j+n_y*i]);				
			}
		}

		fclose(fp);
		
		//calculate fD for all bands 1 to 6
		//for time t
		//store it in array

		SS4[t-1]=0;
			
		for(i=0; i<n_x; ++i)
		{
			for(j=0; j<n_y; ++j)
			{
				SS4[t-1] += (c[j+n_y*i]-0.5)*(c[j+n_y*i]-0.5);				
			}
		}
		
		Vb= n_x*n_y;
		SS4[t-1]=SS4[t-1]/(Vb*0.5*0.5);
		//both delta_x and delta_y are 1, hence not written in the multiplication		
		
}

//write array to file

if((fp=fopen("./SS_1.dat","w"))==NULL)
{
	printf("Can't open file\nAborting the Program......\n");
	exit(0);
}

for(i=0; i<200; ++i)
{
	fprintf(fp,"%le\n",SS1[i]);
}

fclose(fp);


if((fp=fopen("./SS_4.dat","w"))==NULL)
{
	printf("Can't open file\nAborting the Program......\n");
	exit(0);
}

for(i=0; i<200; ++i)
{
	fprintf(fp,"%le\n",SS4[i]);
}

fclose(fp);
free(c);

}
