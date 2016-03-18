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
int t=0;//counter for R 
char name[50];

int numLayer=100;
int layerWidth=1;


//printf("%s\n", "dvdv");

int s= numLayer * layerWidth;
double *avg;
avg=calloc(numLayer,sizeof(double));
int L=n_x/2 - s/2;
int U=n_x/2 + s/2;
	
double R[100];

char counter;
//printf("%s\n", "dvdv");
for(counter='a';counter<='c';counter++)
{
	for(t=1;t<=100;t++)//going through time steps 1-100
	{
		sprintf(name,"../mainData/output_I%c_0.04/c_%d.dat",counter,t);
		//fp=fopen(name,"r");
//		printf("%s\n", name);
		
		if((fp=fopen(name,"r"))==NULL)
			{
				printf("Can't open file\nAborting the Program......\n");
				exit(0);
			}

		for(i=0; i<n_x; ++i)
		{
			for(j=0; j<n_y; ++j)
			{

//		printf("%s\n", name);
		
				fscanf(fp,"%le",&c[j+n_y*i]);
		
//		printf("%s\n", name);
		
			}
		}

		fclose(fp);
//		printf("%s\n", name);
		
		//calculate average of a row and store in array 
		//repeat for all rows

		//shortening the avg vector as we just need to the see the point it intersects the x axis starting from the orign
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

		//find the first zero ----R(t)
		for(h=51;h<60;h++){
			if(avg[h]<0)
				continue;
			else
				break;
		}

		if(avg[h]==0)
			R[t-1]=h-50.5;
		else
			R[t-1]=h - avg[h]/(avg[h]-avg[h-1])-50.5;
	
	}

	//write array to file
	sprintf(name,"./R(t)_I%c.dat",counter);
	
	if((fp=fopen(name,"w"))==NULL)
	{
		printf("Can't open file\nAborting the Program......\n");
		exit(0);
	}

	for(i=0; i<100; ++i)
	{
		fprintf(fp,"%le\n",R[i]);
	}

	fclose(fp);
}

free(c);
free(avg);
}
