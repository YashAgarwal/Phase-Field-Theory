#include<stdlib.h>
#include<stdio.h>
#include<math.h>


int main(void)
{
int n,i,j;
double sig,tol;//Tolerance value	
//Open file for input
FILE *fp;
fp=fopen("input.dat","r");
/*
Format of input file

n

tol

A[0][0] A[0][1] ....
A[1][0]  
A[2][0]  
.
.
.
A[n-1][0].....A[n-1][n-1]

b[0]....b[n-1]

*/
//take inputs
fscanf(fp,"%d",&n);//matrix size
fscanf(fp,"%le",&tol);//tolerance value

//dynamic declaration of matrix
double *A;
double *b;
double *phi;
double *phiOld;
double *error;

A=(double *) malloc(n*n*sizeof(double));
b=(double *) malloc(n*sizeof(double));
phi=(double *) malloc(n*sizeof(double));
phiOld=(double *) malloc(n*sizeof(double));
error=(double *) malloc(n*sizeof(double));

//take values of matrix A
for(i=0;i<n;i++)
{
	for(j=0;j<n;j++)
	fscanf(fp,"%le",&A[i*n+j]);
}

//take values of vector b
for(i=0;i<n;i++)
{
	fscanf(fp,"%le",&b[i]);
}

//close input file
fclose(fp);

//open output file
fp=fopen("output.dat","w");

/*output test

for(i=0;i<n;i++)
{
	for(j=0;j<n;j++)
	{	
	fprintf(fp,"%lf  ",A[i*n+j]);
	}	
	fprintf(fp,"\n");
}

fprintf(fp,"\n\n\n");

for(i=0;i<n;i++)
{
	fprintf(fp,"%lf  ",b[i]);
}
*/

//Gauss Seidel Algorithm
//initial guess for phi... [1 1 1 1 ... ]
for(i=0;i<n;i++)
{
phi[i]=1;
}
int convergenceNotReached=1; // check for convergence, is 1 is convergence is not reached 
while(convergenceNotReached)//till convergence means until convergence=1
{
	//copying the current value of phi in phiOld to check for convergence later 
	for(i=0;i<n;i++)
	{
		phiOld[i]=phi[i];
	}

	for(i=0;i<n;i++)
	{
		sig=0;
		for(j=0;j<n;j++)
		{		
		if(j!=i)
		{
		sig=sig+A[n*i+j]*phi[j];
		}
		}
		phi[i]=(b[i]-sig)/A[n*i+i];		
	}

	//check convergence
	/*
	finding the maximum absolute error 
	*/	
	for(i=0;i<n;i++)//calculating the error vector
	{
		error[i]=fabs(phiOld[i]-phi[i]);
	}
	double max=error[0];
	for(i=0;i<n;i++)//finding max in error vector
	{
		if(max<error[i])max=error[i];
	}
	if(max<=tol)
	{
		convergenceNotReached=0;//means convergence is reached
	}
}

//ouput result
//printing solution vector
for(i=0;i<n;i++)
{
	fprintf(fp,"%le  ",phi[i]);
}

/*
printing residual vector 
It is calculated by 
[A]*[phi]-[b] or known as Ax-b
*/
fprintf(fp,"\n\n\n");

for(i=0;i<n;i++)
{
	sig=0;
	for(j=0;j<n;j++)
	{
		sig+=A[i*n+j]*phi[j]; 
	}
	sig-=b[i];
	fprintf(fp,"%le  ",sig);
} 

free(phiOld);
free(error);
free(phi);
free(b);
free(A);
fclose(fp);
}
