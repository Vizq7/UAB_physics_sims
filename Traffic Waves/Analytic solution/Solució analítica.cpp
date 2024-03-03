#include <stdio.h>
#include <math.h>
#include <string.h>

int main(){
	int C = 20000;
	double m = 1500;
	double rho_max = 0.25;
	double rho_cri = 0.99;
	double D = 1;
	int L = 4;
	double v_max = 33.3;
	int d = 26;
	int tc = 1;
	int ta = 1;
	int vegades = 30;
	int cops = 26;
	int num_posicions = 5;
	int N = 4;
	double posicions[num_posicions];
	double rho[num_posicions];
	double rho2[vegades];
	double velocitat[num_posicions];
	double velocitat2[vegades];
	double flux[num_posicions];
	double flux2[vegades];
	int i;
	int j;
	int simulacions = 34;
	posicions[0]= 4;
	posicions[1]= 34;
	posicions[2]= 64;
	posicions[3]= 94;
	posicions[4]= 124;
	double relacio[vegades];
	
	
//	for (i=1; i<=N; i++){
//		rho[i-1] = 1/fabs((posicions[i]-posicions[i-1]));
//		velocitat[i-1] = (v_max*log(rho_max/rho[i-1]))/log(rho_max/rho_cri);
//		flux[i-1] = (v_max*log(rho_max/rho[i-1])*rho[i-1])/log(rho_max/rho_cri);
//	}
//	
//	for (i=1; i<=N; i++){
//		velocitat[i-1] = (v_max*log(rho_max/rho[i-1]))/log(rho_max/rho_cri);
//	}
//	
//	
//	for (i=1; i<=N; i++){
//		flux[i-1] = (v_max*log(rho_max/rho[i-1])*rho[i-1])/log(rho_max/rho_cri);
//	}

//	FILE*analitica; 
//	
//	analitica = fopen("Analítica.txt", "wb");
//
//	for (i=0; i<N; i++){
//		
//		fprintf(analitica,"%f\t %f\t %f\n ", rho[i], velocitat[i], flux[i]);	
//	}
//
//	fclose(analitica);
//	
	
	for (j=4; j<=vegades; j++){
		rho2[31] = 0.0000001;
		rho2[32] = 0.02;
		rho2[33] = 0.03;
		rho2[34] = 0.031;
		rho2[j] = 1/fabs(j);
	}

	for (j=4; j<=vegades; j++){
		velocitat2[31] = fabs((v_max*log(rho_max/rho2[j]))/log(rho_max/rho_cri));
		velocitat2[32] = fabs((v_max*log(rho_max/rho2[j]))/log(rho_max/rho_cri));
		velocitat2[33] = fabs((v_max*log(rho_max/rho2[j]))/log(rho_max/rho_cri));
		velocitat2[34] = fabs((v_max*log(rho_max/rho2[j]))/log(rho_max/rho_cri));
		velocitat2[j] = fabs((v_max*log(rho_max/rho2[j]))/log(rho_max/rho_cri));
	}
	
	for (j=4; j<=vegades; j++){
		
		flux2[j] = fabs((v_max*log(rho_max/rho2[j])*rho2[j])/log(rho_max/rho_cri));
		flux2[31] = fabs((v_max*log(rho_max/rho2[31])*rho2[31])/log(rho_max/rho_cri));
		flux2[32] = fabs((v_max*log(rho_max/rho2[32])*rho2[32])/log(rho_max/rho_cri));
		flux2[33] = fabs((v_max*log(rho_max/rho2[33])*rho2[33])/log(rho_max/rho_cri));
		flux2[34] = fabs((v_max*log(rho_max/rho2[34])*rho2[34])/log(rho_max/rho_cri));
	}
	
	for (j=4; j<=simulacions; j++){
		
		relacio[j] = fabs(flux2[j]/rho_cri);
	}
	
	FILE*output; 
	
	output = fopen("Analítica_2.txt", "wb");

	for (j=4; j<=simulacions; j++){
		
		fprintf(output,"%f\t %f\t %f\t %f\n ", rho2[j], velocitat2[j], flux2[j], relacio[j]);	
	}

	fclose(output);


}



