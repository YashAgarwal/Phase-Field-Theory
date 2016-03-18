/*
Explicit Method
for general boundary condition case
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(void)
{
FILE *fp,*oct;
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
filename=(char*) malloc(50*sizeof(char));//assuming filename will not exceed 50

double *c,*cNew;//These vector will hold the value of concentration for all N nodes
c=(double*) malloc(N*sizeof(double));
cNew=(double*) malloc(N*sizeof(double));


//inputs the initial conc. profile 
//Make sure you have executed the ./make(creates the input profile) file before this 
fp=fopen("c_0.dat","r");
for(i=0;i<N;i++)
{
fscanf(fp,"%le",&c[i]);
}

fclose(fp);

//creating a octave script side by side and writing the output file names directly into the script
oct=fopen("plotAnimation.oct","w");

//using general boundary condition
double b1,b2,beta,phi;
b1=c[0];
b2=c[N-1];
beta=0.01;
phi=alpha*delt/(delx*delx);//just saving some recalculation

for(j=1;j<=T;j++)
{
	for(i=0;i<N;i++)
	{
		if(i==0)
		{		
			if(left==1)//Boundary Condition I
			{					
			cNew[i]=c[i] + phi*(c[i+1]+b1-2*c[i]);				
			}
			else if(left==2)//Boundary Condition II
			{		
			beta=0;			
			b1=c[1]-beta*2*delx;		
			cNew[i]=c[i] + phi*(c[i+1]+b1-2*c[i]);				
			}
			else//Boundary Condition III
			{		
			beta=0.01;			
			b1=c[1]-beta*2*delx;		
			cNew[i]=c[i] + phi*(c[i+1]+b1-2*c[i]);				
			}		
		}
		else if(i==N-1)
		{
			if(right==1)//Boundary Condition I
			{					
			cNew[i]=c[i] + phi*(b2+c[i-1]-2*c[i]);		
			}
			else if(right==2)//Boundary Condition II
			{		
			beta=0;			
			b2=c[N-2]+beta*2*delx;
			cNew[i]=c[i] + phi*(b2+c[i-1]-2*c[i]);		
			}
			else//Boundary Condition III
			{		
			beta=0.01;			
			b2=c[N-2]+beta*2*delx;
			cNew[i]=c[i] + phi*(b2+c[i-1]-2*c[i]);		
			}
		}		
		else
		{		
		cNew[i]=c[i] + phi*(c[i+1]+c[i-1]-2*c[i]);		
		}	
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
		fprintf(oct,"\ndata=load(\"%s\");\nplot(data(:,1)','LineWidth',4);\nylim([0,1]);\nhold off\nsleep(0.1);",filename);
	}
	
}
//freeing up memory and closing the necessary files
fclose(oct);
free(filename);
free(c);
free(cNew);
return 0;
}
