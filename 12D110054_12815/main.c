/*
Precipitate diffusion with periodic boundary condition
*/
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

int main(void)
{

/**Opening the file for input parameters**/
FILE *fp;
fp=fopen("input.dat","r");

int i,j;
double N,T,delx,delt,alpha,Tol;

fscanf(fp,"%le",&N);
fscanf(fp,"%le",&T);
fscanf(fp,"%le",&alpha);
fscanf(fp,"%le",&Tol);
fclose(fp);

/**Calculating the values for delx and delt according the input parameters**/
delx=1/N;
delt=1/T;

/**Setting the time step for saving the profiles**/
int tStep=T/400;

/**Allocating the dynamic memory according the input N**/
double *c,*cNew,*b;
c=(double*)malloc(N*sizeof(double));
cNew=(double*)malloc(N*sizeof(double));
b=(double*)malloc(N*sizeof(double));

/**Allocating memory for a variable to store dynamically generated filenames**/
char *filename;
filename=(char*)malloc(50*sizeof(char));

//input the Initial Profile
fp=fopen("c_0.dat","r");
for(i=0;i<N;i++)
{
	fscanf(fp,"%le",&c[i]);
}

fclose(fp);

/**The matrix can be completely defiend with only phi,gamma and 0**/
double sig,phi,gamma,error=0,errormax=0;
phi=alpha*delt/(delx*delx);
gamma=(2*phi+1);

/**Creating a file to make the octave script to view all the profiles generated as an animation**/
FILE *oct;
oct=fopen("plotAnimation.oct","w");

/**Adding the initial file to the animation**/
sprintf(filename,"./output/c_%d.dat",0);	
fp=fopen(filename,"w");
for(i=0;i<N;i++)
{
fprintf(fp,"%le\n",c[i]);
}
fclose(fp);
fprintf(oct,"\ndata=load(\"%s\");\nplot(data(:,1)','LineWidth',4);\nylim([0,1]);\nhold off\nsleep(0.5);",filename);

/**Starting the time loop**/
for(j=1;j<=T;j++)
{
	//make b matrix
	//b matrix is equal to c matrix
	for(i=0;i<N;i++)
	{
		b[i]=c[i];
	}

	//apply Gauss Seidel
	//initial guess for GS
	// take the previous step solution as the guess for the next to reduce the convergence steps
	for(i=0;i<N;i++)
	{
		cNew[i]=c[i];
	}

	do
	{
	//Storing previous Value to check for convergence later	
		for(i=0;i<N;i++)
		{
			c[i]=cNew[i];
		}	
		for(i=0;i<N;i++)
		{
		sig=0;
		if(i==0)	
			{
			sig+=-1*phi*cNew[1];
			sig+=-1*phi*cNew[(int)N-1];		
			}
		else if(i==N-1)
			{
			sig+=-1*phi*cNew[0];
			sig+=-1*phi*cNew[(int)N-2];		
			}
		else
			{
			sig+=-1*phi*cNew[i-1];
			sig+=-1*phi*cNew[i+1];		
			}
		cNew[i]=(b[i]-sig)/gamma;
		}
		//cheking for convergence
		errormax=0;	
		for(i=0;i<N;i++)
		{
			error=fabs(c[i]-cNew[i]);
			if(errormax<error)errormax=error;
		}
	}while(!(errormax<Tol));

		//Now cNew is copied in c
		for(i=0;i<N;i++)
		{
			c[i]=cNew[i];
		}

		//print new profile to a file
		if(j%tStep==0)//after every tStep time steps
		{
		sprintf(filename,"./output/c_%d.dat",j);	
		fp=fopen(filename,"w");
		for(i=0;i<N;i++)
		{
		fprintf(fp,"%le\n",cNew[i]);
		}
		fclose(fp);
		//make a entry for the file in the animation
		fprintf(oct,"\ndata=load(\"%s\");\nplot(data(:,1)','LineWidth',4);\nylim([0,1]);\nhold off\nsleep(0.5);",filename);
		}
}
//freeing up memory and closing the necessary files
fclose(oct);
free(c);
free(cNew);
free(b);
free(filename);

return 0;
}

