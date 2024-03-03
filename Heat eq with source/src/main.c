#include <stdio.h>
#include <math.h>
#include <string.h>
#define _USE_MATH_DEFINES // for C
#include <math.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

int main()
{
	int N = 99;
	float l = 1;
	float Tc = 0.56*(36.5+273.15)/(pow(0.02,2)*944000);
	float ta=0.025;

	float dx = l/N;

	float dt1=0.51*dx*dx;
	float dt2=0.49*dx*dx;
	float dt3=0.25*dx*dx;
	float dtsol=0.13*dx*dx;
	
	int len_t1=480+1;
	int len_t2=500;
	int len_t3=980;
	int len_sol=1536;

	double T1[len_t1][N];
	double T2[len_t2][N];
	double T3[len_t3][N];
	double TSol[len_sol][N];

	for (int i = 0; i < len_t1; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			T1[i][j] = Tc;
		}
	}
	for (int i = 0; i <= len_t2; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			T2[i][j] = Tc;
		}

	}
	for (int i = 0; i <= len_t3; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			T3[i][j] = Tc;
		}

	}
	for (int i = 0; i <= len_sol; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			TSol[i][j] = Tc;
		}

	}
	

	for (int i = 0; i < len_t1-1; i++)
	{
		for (int j = 1; j < N; j++)
		{
			T1[i+1][j]= dt1/(dx*dx)*(T1[i][j+1]-2*T1[i][j]+T1[i][j-1])+dt1+T1[i][j];
		}

	}
	
	for (int i = 0; i < len_t2; i++)
	{
		for (int j = 1; j < N; j++)
		{
			T2[i+1][j]= dt2/(dx*dx)*(T2[i][j+1]-2*T2[i][j]+T2[i][j-1])+dt2+T2[i][j];
		}

	}

	
	for (int i = 0; i < len_t3; i++)
	{
		for (int j = 1; j < N; j++)
		{
			T3[i+1][j]= dt3/(dx*dx)*(T3[i][j+1]-2*T3[i][j]+T3[i][j-1])+dt3+T3[i][j];
		}

	}

	for (int i = 0; i < len_sol; i++)
	{
		for (int j = 1; j < N; j++)
		{
			TSol[i+1][j]= dtsol/(dx*dx)*(TSol[i][j+1]-2*TSol[i][j]+TSol[i][j-1])+dtsol+TSol[i][j];
		}

	}

FILE *f1= fopen("T1.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		
		fprintf(f1,"%f\t", 1.0/0.56*(T1[len_t1-1][i])*(pow(0.02,2)*944000)-273.15);	
	}
	
	fclose(f1);

FILE *f2 = fopen("T2.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		fprintf(f2,"%f\t",1.0/0.56*(T2[len_t2][i])*(pow(0.02,2)*944000)-273.15);	
	}
	
	fclose(f2);

FILE *fptr;
	fptr = fopen("T3.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		fprintf(fptr,"%f\t",1.0/0.56*(T3[len_t3][i])*(pow(0.02,2)*944000)-273.15);	
	}
	
	fclose(fptr);

	fptr = fopen("TSol.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		fprintf(fptr,"%f\t",1.0/0.56*(TSol[len_sol][i])*(pow(0.02,2)*944000)-273.15);	
	}
	
	fclose(fptr);

	fptr = fopen("TSol_matrix.txt", "wb");
	for (int i = 0; i <= len_sol; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			fprintf(fptr,"%f,",1.0/0.56*(TSol[i][j])*(pow(0.02,2)*944000)-273.15);	
		}
		fprintf(fptr,"\n");	
	}
	
	fclose(fptr);

///////

	double dti = dx;
	double dti2 = 0.5*dx;
	double dti3 = 0.125*dx;

	double gamma1 = dti/(dx*dx);
	double gamma2 = dti2/(dx*dx);
	double gamma3 = dti3/(dx*dx);
	double gammaSol = dt3/(dx*dx);

	double Ti1[3][N];
	double Ti2[5][N];
	double Ti3[20][N];
	double TiSol[len_sol][N]; //usem el mallat dt=0.25(dx)^2


	for (int i = 0; i <= 3; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			Ti1[i][j] = Tc;
		}

	}

	for (int i = 0; i <= 5; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			Ti2[i][j] = Tc;
		}

	}

	for (int i = 0; i <= 20; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			Ti3[i][j] = Tc;
		}

	}

	for (int i = 0; i <= len_sol; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			TiSol[i][j] = Tc;
		}

	}

	for (int i = 1; i <= 3; i++)
	{
		double Tans[N];
		double Tgs[N];

		for (int k = 0; k <= N; k++)
		{
			Tans[k] = Ti1[i-1][k];
		}

		int N_GS = 0;
		while (N_GS < 100000)
		{
			for (int j = 1; j < N; j++)
			{
				Tgs[j] = (gamma1*(Tans[j+1]+Tans[j-1])+Ti1[i-1][j]+dti)/(1.0+2.0*gamma1);
			}
			
			for (int k = 1; k < N; k++)
			{
				Tans[k] = Tgs[k];
			}
			N_GS++;
		}
		for (int z = 1; z < N; z++)
		{
			Ti1[i][z]=Tgs[z];
		}
		

	}

	for (int i = 1; i <= 5; i++)
	{
		double Tans[N];
		double Tgs[N];

		for (int k = 0; k <= N; k++)
		{
			Tans[k] = Ti2[i-1][k];
		}

		int N_GS = 0;
		while (N_GS < 100000)
		{
			for (int j = 1; j < N; j++)
			{
				Tgs[j] = (gamma2*(Tans[j+1]+Tans[j-1])+Ti2[i-1][j]+dti2)/(1.0+2.0*gamma2);
			}
			
			for (int k = 1; k < N; k++)
			{
				Tans[k] = Tgs[k];
			}
			N_GS++;
		}
		for (int z = 1; z < N; z++)
		{
			Ti2[i][z]=Tgs[z];
		}

	}

	for (int i = 1; i <= 20; i++)
	{
		double Tans[N];
		double Tgs[N];

		for (int k = 0; k <= N; k++)
		{
			Tans[k] = Ti3[i-1][k];
		}

		int N_GS = 0;
		while (N_GS < 100000)
		{
			for (int j = 1; j < N; j++)
			{
				Tgs[j] = (gamma3*(Tans[j+1]+Tans[j-1])+Ti3[i-1][j]+dti3)/(1.0+2.0*gamma3);
			}
			
			for (int k = 1; k < N; k++)
			{
				Tans[k] = Tgs[k];
			}
			N_GS++;
		}
		for (int z = 1; z < N; z++)
		{
			Ti3[i][z]=Tgs[z];
		}

	}

	for (int i = 1; i <= len_sol; i++)
	{
		double Tans[N];
		double Tgs[N];

		for (int k = 0; k <= N; k++)
		{
			Tans[k] = TiSol[i-1][k];
		}

		int N_GS = 0;
		while (N_GS < 100000)
		{
			for (int j = 1; j < N; j++)
			{
				Tgs[j] = (gammaSol*(Tans[j+1]+Tans[j-1])+TiSol[i-1][j]+dt3)/(1.0+2.0*gammaSol);
			}
			
			for (int k = 1; k < N; k++)
			{
				Tans[k] = Tgs[k];
			}
			N_GS++;
		}
		for (int z = 1; z < N; z++)
		{
			TiSol[i][z]=Tgs[z];
		}

	}

	

	fptr = fopen("Ti1.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		fprintf(fptr,"%f\t",1.0/0.56*(Ti1[3][i])*(pow(0.02,2)*944000)-273.15);	
	}
	
	fclose(fptr);

	fptr = fopen("Ti2.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		fprintf(fptr,"%f\t",1.0/0.56*(Ti2[5][i])*(pow(0.02,2)*944000)-273.15);	
	}
	
	fclose(fptr);

	fptr = fopen("Ti3.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		fprintf(fptr,"%f\t",1.0/0.56*(Ti3[20][i])*(pow(0.02,2)*944000)-273.15);	
	}
	
	fclose(fptr);

	fptr = fopen("TiSol.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		fprintf(fptr,"%f\t",1.0/0.56*(TiSol[len_sol][i])*(pow(0.02,2)*944000)-273.15);	
	}
	
	fclose(fptr);
	


//en sol analitica posem el temps on ens quedem més propers en la discritització tempopal a tc
//si en graph error num oscilla => la serie infinita encara no conv => hem d'augmentar la N del sumatori inf
//la sol inicial de gs és el vector de CI
//T=T/P_joule
//bn= 0, n=2m; 4/(n*pi) n=2m+1

//sol analitica
int N_sum = 100000;

double T1_an[N];
double T2_an[N];
double T3_an[N];

double Ti1_an[N];
double Ti2_an[N];

double T_sol_an[N];

double t1 = (len_t1-1.0)*dt1;
double t2 = (len_t2)*dt2;
double t3 = (len_t3)*dt3;

double ti1 = 3.0*dti;
double ti2 = 5.0*dti2;

for (int i = 0; i <= N; i++)
{
	double sum = 0;
	for (int n = 0; n < N_sum; n++)
	{
		sum += sin((2.0*n+1.0)*M_PI*i*dx)*(1.0-exp(-M_PI*M_PI*t1*pow((2.0*n+1.0),2)))/(pow(M_PI*(2*n+1),3));
	}
	T1_an[i]=4.0*sum+Tc;
}

for (int i = 0; i <= N; i++)
{
	double sum = 0;
	for (int n = 0; n < N_sum; n++)
	{
		sum += sin((2.0*n+1.0)*M_PI*i*dx)*(1.0-exp(-M_PI*M_PI*t2*pow((2.0*n+1.0),2)))/(pow(M_PI*(2*n+1),3));
	}
	T2_an[i]=4.0*sum+Tc;
}

for (int i = 0; i <= N; i++)
{
	double sum = 0;
	for (int n = 0; n < N_sum; n++)
	{
		sum += sin((2.0*n+1.0)*M_PI*i*dx)*(1.0-exp(-M_PI*M_PI*t3*pow((2.0*n+1.0),2)))/(pow(M_PI*(2.0*n+1.0),3));
	}
	T3_an[i]=4.0*sum+Tc;
}

for (int i = 0; i <= N; i++)
{
	double sum = 0;
	for (int n = 0; n < N_sum; n++)
	{
		sum += sin((2.0*n+1.0)*M_PI*i*dx)*(1.0-exp(-M_PI*M_PI*ti1*pow((2.0*n+1.0),2)))/(pow(M_PI*(2.0*n+1.0),3));
	}
	Ti1_an[i]=4.0*sum+Tc;
}

for (int i = 0; i <= N; i++)
{
	double sum = 0;
	for (int n = 0; n < N_sum; n++)
	{
		sum += sin((2*n+1)*M_PI*i*dx)*(1-exp(-M_PI*M_PI*ti2*pow((2*n+1),2)))/(pow(M_PI*(2*n+1),3));
	}
	Ti2_an[i]=4*sum+Tc;
}
double t = 0.02036;
for (int i = 0; i <= N; i++)
{
	double sum = 0;
	for (int n = 0; n < N_sum; n++)
	{
		sum += sin((2*n+1)*M_PI*i*dx)*(1-exp(-M_PI*M_PI*t*pow((2*n+1),2)))/(pow(M_PI*(2*n+1),3));
	}
	T_sol_an[i]=4*sum+Tc;
}

	fptr = fopen("T1_an.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		fprintf(fptr,"%f\t",1.0/0.56*(T1_an[i])*(pow(0.02,2)*944000)-273.15);	
	}
	
	fclose(fptr);

	fptr = fopen("T2_an.txt", "wb");
	//printf("T2 \n");
	for (int i = 0; i <= N; i++)
	{
		//printf("%f ",T2[len_t2][i]-T2_an[i]);
		fprintf(fptr,"%f\t",1.0/0.56*(T2_an[i])*(pow(0.02,2)*944000)-273.15);
	}
	 
	
	fclose(fptr);

	//printf("\n\nT3 \n");
	fptr = fopen("T3_an.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		//printf("%f ",T3[len_t3][i]-T3_an[i]);
		fprintf(fptr,"%f\t",1.0/0.56*(T3_an[i])*(pow(0.02,2)*944000)-273.15);	
	}
	fclose(fptr);


	//printf("\n\nTi1 \n");
	fptr = fopen("Ti1_an.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		//printf("%f ",Ti1[3][i]-Ti1_an[i]);
		fprintf(fptr,"%f\t",1.0/0.56*(Ti1_an[i])*(pow(0.02,2)*944000)-273.15);	
	}

	fclose(fptr);

	//printf("\n\nTi2 \n");
	fptr = fopen("Ti2_an.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		//printf("%f ",Ti2[5][i]-Ti2_an[i]);
		fprintf(fptr,"%f\t",1.0/0.56*(Ti2_an[i])*(pow(0.02,2)*944000)-273.15);	
	}

	fclose(fptr);
	//printf("\n\nTsol \n");
	fptr = fopen("T_sol_an.txt", "wb");
	for (int i = 0; i <= N; i++)
	{
		if(i==38){
			printf("%f", 1.0/0.56*(T_sol_an[i])*(pow(0.02,2)*944000)-273.15);
		}
		//printf("%f ",TSol[len_sol][i]-T_sol_an[i]);
		fprintf(fptr,"%f\t",1.0/0.56*(T_sol_an[i])*(pow(0.02,2)*944000)-273.15);	
	}

	fclose(fptr);

}
