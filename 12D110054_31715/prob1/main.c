/*

Lecture 3
Problem  1
Improved Euler Algorithm 
for y'+0.1*y=0, y(0)=2
find y(1)

Author:Yash Agarwal

*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double f(double x,double y)//y'=f(x,y)
{
return -0.1*y; 
}

int main(void)
{
FILE *fp,*fp1;
double h;//time step
double x;
double y;
int N;//number of steps
int i;
double k1,k2,k3,k4;
double ans,real_ans;

//take all the inputs
/*
Format of input file

x0
y0
h
N

*/
//Opening Output file
if((fp1=fopen("output.dat","w"))==NULL)
{
	printf("\nUnable to create \"output.dat\"\n");
	exit(0);
}

//Opening Input file
if((fp=fopen("input.dat","r"))==NULL)
{
	fprintf(fp1,"\nUnable to find \"input.dat\"\n");
	exit(0);
}
fscanf(fp,"%le",&x);
fscanf(fp,"%le",&y);
fscanf(fp,"%le",&h);
fscanf(fp,"%d",&N);
fclose(fp);

//main Algorithm
for(i=0;i<N;i++)
{
	k1=h*f(x,y);
	k2=h*f(x+h,y+k1);
	x+=h;
	y+=(k1+k2)/2;	
	//output x,y
	fprintf(fp1,"%le  ",x);
	fprintf(fp1,"%le\n",y);
	ans=y;
}

//our solution
fprintf(fp1,"\n\nThe answer is %le\n",ans);

//analytical answer is y(1)=2*exp(-0.1)
real_ans=2*exp(-0.1);
fprintf(fp1,"\nThe analytical answer is  %le\n",real_ans);
//printing error
fprintf(fp1,"\n The error is  %le\n",real_ans-ans);

fclose(fp1);
return 0;
}
