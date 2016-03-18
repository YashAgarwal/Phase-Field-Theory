/*
Implicit Method
for general boundary condition case
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(void)
{
FILE *fp;
double alpha,delt,delx;
int T,N;
fp=fopen("input.dat","r");
//take all other inputs except the initial conc. profile 

fscanf(fp,"%d",&N);//number of nodes 
fscanf(fp,"%le",&delx);
fscanf(fp,"%d",&T);//number of time steps
fscanf(fp,"%le",&delt);
fscanf(fp,"%le",&alpha);
/*
T & N should be chosen such that T*delt=1
and N*deln=1
*/
int left,right;//for storing which boundary condition will be used for which side : 1,2 or 3 

//input the boundary condition on left side or right side
fscanf(fp,"%d",&left);
fscanf(fp,"%d",&right);

fclose(fp);

int i,j;

//Allocatin space for the dynamic filename for the output files
char *filename;
filename=(char*) malloc(50*sizeof(char));

double *c,*cNew;//These vector will hold the value of concentration for all N nodes
c=(double*) malloc(N*sizeof(double));
cNew=(double*) malloc(N*sizeof(double));


double *b;//For Gauss Seidel Method
b=(double*) malloc(N*sizeof(double));

//inputs the initial conc. profile 
//Make sure you have executed the ./make(creates the input profile) file before this 
fp=fopen("c_0.dat","r");
for(i=0;i<N;i++)
{
fscanf(fp,"%le",&c[i]);
}

fclose(fp);

//using general boundary condition
double b1,b2,beta,gamma,phi;
b1=c[0];
b2=c[N-1];
beta=0.01;
phi=(alpha*delt)/(delx*delx);//saving some recalculation
gamma=-1*(2*phi+1);//saving some recalculation

/*
To Apply the Gauss Seidel Method we need
A and b
here A is a constant matrix

for Boundary Condition I
A =	[ 0     gamma   phi]
	[phi	gamma	phi]
	[phi	gamma	phi]
	:......
	:......
	.
	[phi	gamma	 0 ]

and  
 b= [-c[i]-phi*b1]
 	[   -c[i]  	]
 	[   -c[i]   ]
	:
	:
	[-c[i]-phi*b2]

for Boundary Condition II and III
A =	[ 0     gamma  2phi]
	[phi	gamma	phi]
	[phi	gamma	phi]
	:......
	:......
	.
	[2phi	gamma	 0 ]

and  
 b= [-c[i]+2phi*beta*delx]
 	[       -c[i]		]
 	[       -c[i]		]
	:
	:
	[-c[i]-2phi*beta*delx]

Thus we don't need to declare A matrix no matter what N is 
*/


FILE *oct;
oct=fopen("plotAnimation.oct","w");
for(j=1;j<=T;j++)
{	
	
	//This for loop is the values in vector b
	for(i=0;i<N;i++)
	{	
		if(i==0)
		{		
			if(left==1)//Boundary Condition I
			{			
			b[i]=-1*c[i]-phi*b1;
			}
			else if(left==2)//Boundary Condition II
			{		
			b[i]=-1*c[i];
			}
			else//Boundary Condition III
			{		
			b[i]=-1*c[i]+2*0.01*phi*delx;
			}
		}
		else if(i==N-1)
		{
			if(right==1)//Boundary Condition I
			{								
			b[i]=-1*c[i]-phi*b2;
			}
			else if(right==2)//Boundary Condition II
			{		
			b[i]=-1*c[i];
			}
			else//Boundary Condition III
			{		
			b[i]=-1*c[i]-2*0.01*phi*delx;
			}
		}		
		else
		{		
		b[i]=-1*c[i];
		}
			
	}
	
	//We apply Gauss Seidel on the Matrices
			
	int convergenceReached=0;//A counter to Check if Convergence is reached
	double sig,errorMax,error,Tol=1e-6;
	//Set complete cNew to zero for the initial solution
	for(i=0;i<N;i++)
		{
			cNew[i]=0;
		}
	while(!convergenceReached)
	{
		for(i=0;i<N;i++)
		{
			c[i]=cNew[i];
		}		
		errorMax=0;
		for(i=0;i<N;i++)
		{
			sig=0;
			
			if(i==0)
			{
				if(left==1)//Boundary Condition I
				{
					sig+=phi*cNew[1];
				}
				else//Boundary Condition II & III
				{
					sig+=2*phi*cNew[1];
				}
			}	
			else if(i==N-1)
			{
				if(left==1)//Boundary Condition I
				{
					sig+=phi*cNew[N-2];
				}
				else//Boundary Condition II & III
				{
					sig+=2*phi*cNew[N-2];
				}
			}
			else
			{
				sig+=phi*(cNew[i-1]+cNew[i+1]) ;
			}	
			
			cNew[i]=(b[i]-sig)/gamma;				
			error=fabs(cNew[i]-c[i]);
			if(error>errorMax)
			{
				errorMax=error;
			}
		}
		if(errorMax<Tol)
			convergenceReached=1;
	}
	//copy cNew to c
	for(i=0;i<N;i++)
		{
		c[i]=cNew[i];
		}
	//transfer the cNew to a file
	if(j%50==0)//after every 50 time steps
	{
		sprintf(filename,"./output/c_%d.dat",j);	
		fp=fopen(filename,"w");
		for(i=0;i<N;i++)
		{
		fprintf(fp,"%le\n",cNew[i]);
		}
		fclose(fp);
		fprintf(oct,"\ndata=load(\"%s\");\nplot(data(:,1)','LineWidth',4);\nylim([0,1]);\nhold off\nsleep(0.5);",filename);
	}
}
//freeing up memory and closing the necessary files
fclose(oct);
free(filename);
free(c);
free(b);
free(cNew);
return 0;
}
