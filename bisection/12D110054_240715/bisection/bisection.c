#include<stdio.h>
#include<math.h>

int main(void)
{

FILE *fout;
fout=fopen("outputData.dat","w");
double N=0,Nmax=1e10;
while(N++<Nmax)//finding the root
{
fprintf(fout,"garbage\n");
}
fclose(fout);
return 0;
}
