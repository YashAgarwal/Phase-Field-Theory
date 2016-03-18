#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(void)
{
	int *M;							//Crystal matrix data
	int Q;							//crystallographic orientation
	int i,j,i1,j1;						//iteration variables 
	
	int x_index;
	int y_index;

	FILE *fp,*fp2,*fp3; 						//file pointer 
	double N_MCS,k;
	int e_old;
	
	M=(int *) malloc(40*40*sizeof(int));
	
	if((fp=fopen("crystal_data_input.dat","r"))==NULL)
	{
		printf("Unable to open file\n\n");				//file check
		exit(0);
	}
	else 
	{
		fp=fopen("crystal_data_input.dat","r");		//file opening 
	}
	
/*Taking Data from crystal input file*/	
	for(i=0;i<40;i++)
	{
		for(j=0;j<40;j++)
		{
			fscanf(fp,"%d",&M[40*i+j]);
		}
	}
/*-------------------------------------*/

/*----monte-carlo simulation-----------*/
	N_MCS=1600*1e6;
	while(k<N_MCS)
	{
		int l=rand()%40;				//random_x	
		int m=rand()%40;				//random_y
		int likes_old=0;				//like spins in before 
		int likes_new=0;				//like spins after
		int swap_x;
		int swap_y;
		int temp;
		int e_new=9;
		int delta_E;					//change in energy
		double probablity;
		double temperature=1000;
		double random_gen;
		
/*		while(1)
		{
			delta_l=rand()%3-1;
			delta_m=rand()%3-1;
			if(delta_l==0 && delta_m==0)
				{continue;}
			else if(l+delta_l==0 || m+delta_m==0 || l+delta_l==101 || m+delta_m==0)
				{continue;}
			else
				{break;}
		}*/
		k++;
		
		/*_____ENERGY IN CURRENT STATE______*/
		for(i=-1;i<2;i++)
		{
			for(j=-1;j<2;j++)
			{
				x_index=(l+j);
				y_index=(m+i);
				
				if(x_index==-1)x_index=39;
				if(x_index==40)x_index=0;
				if(y_index==-1)y_index=39;
				if(y_index==40)y_index=0;
				
				if(M[x_index*40+y_index]==M[l*40+m])
				{
					likes_old++;
				}
			}
		}
		likes_old--;
		/*__________________________________*/
		if(likes_old==8)
			continue;
		

		
		/*______ENERGY IN SWITCHED STATE_____*/
		for(i=-1;i<2;i++)
		{
			for(j=-1;j<2;j++)
			{
				x_index=(l+j);
				y_index=(m+i);
				
				if(x_index==-1)x_index=39;
				if(x_index==40)x_index=0;
				if(y_index==-1)y_index=39;
				if(y_index==40)y_index=0;
				
				printf("%d\t%d\n",x_index,y_index);
				if(M[l*40+m]==M[x_index*40+y_index])continue;
				else
				{
					M[l*40+m]=-M[l*40+m];
					M[x_index*40+y_index]=-M[x_index*40+y_index];
				}
				
				likes_new=0;
				for(i1=-1;i1<2;i1++)
				{
					for(j1=-1;j1<2;j1++)
					{
						int x_index_1=(l+j1);
						int y_index_1=(m+i1);
				
						if(x_index_1==-1)x_index_1=39;
						if(x_index_1==40)x_index_1=0;
						if(y_index_1==-1)y_index_1=39;
						if(y_index_1==40)y_index_1=0;
				
						if(M[x_index_1*40+y_index_1]==M[l*40+m])
						{
							likes_new++;
						}
					}
				}
				likes_new--;
				printf("%d\t%d\n",x_index,y_index);
				M[l*40+m]=-M[l*40+m];
				M[x_index*40+y_index]=-M[x_index*40+y_index];
				
				if(likes_new<e_old)
				{
					swap_x=x_index;
					swap_y=y_index;
					printf("%d\t%d\n",swap_x,swap_y);
					e_old=likes_new;
				}
			}
		}
		
		/*__________________________________*/
		
		delta_E=likes_old-likes_new;
		printf("%d\n",delta_E);
		if(delta_E<0)
		{
			//swap accepted
			printf("%d\t%d\n",swap_x,swap_y);
			M[40*swap_x+swap_y]=-1*M[40*swap_x+swap_y];
			M[l*40+m]=-M[l*40+m];
		}
		else
		{
			//swap not accepted
		}
		/*
		else 
		{
			probablity=exp(((double)delta_E)/2);
			random_gen=(rand()%10000)/10000;
			if(random_gen > probablity)
			{
				//spin swapping is reverted
				M[102*l+m]=M[102*l+m]-1;
			}
			else
			{
				//spin swap is accepted 
			}
		}*/
		if(((int)k)%20000==0)printf("\r%.2le",k);
	}
/*-------------------------------------*/	
	
	fp2=fopen("SimulationResult.dat","w");
	
	//----------writing crystal matrix data in file-------------////////
	for(i=0;i<40;i++)
	{
		for(j=0;j<40;j++)
		{
			fprintf(fp2,"%d ",M[40*i+j]);
		}
		fprintf(fp2,"\n");
	}



	return 0;
}
