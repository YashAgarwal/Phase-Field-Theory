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
int t[3]={10,20,50};
char name[50];

int numLayer=100;
int layerWidth=1;

int s= numLayer * layerWidth;
double *avg;
avg=calloc(numLayer,sizeof(double));
int L=n_x/2 - s/2;
int U=n_x/2 + s/2;

	
for(k=0;k<3;k++)
{
	sprintf(name,"../mainData/output_Ia_0.04/c_%d.dat",t[k]);
	fp=fopen(name,"r");
	
	for(i=0; i<n_x; ++i)
	{
		for(j=0; j<n_y; ++j)
		{
			fscanf(fp,"%le",&c[j+n_y*i]);
		}
	}

	fclose(fp);
	//calculate average of a row and store in array 
	//repeat for all rows

	int h,l;
	i=L;
	for(h=0;h<numLayer;h++)
	{
		avg[h]=0;
		for(l=0;l<layerWidth;l++)
		{
			for(j=0; j<n_y; ++j)
			{
				avg[h]+=c[j+n_y*i];
			}
			i++;
		}
		avg[h]/=(double)(layerWidth*n_y);
		avg[h]-=0.5;
	}

	//write array to file
	sprintf(name,"./analysis_%d.dat",t[k]);
	//printf("%s\n",name);
	
	if((fp=fopen(name,"w"))==NULL)
	{
		printf("Can't open file\nAborting the Program......\n");
		exit(0);
	}

	for(i=0; i<numLayer; ++i)
	{
		fprintf(fp,"%le\t%d\n",avg[i],i-(numLayer/2));
	}

	fclose(fp);
}

free(c);
free(avg);
}
