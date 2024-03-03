#include <stdio.h>
#include <math.h>

double adimx(double x){
	return x*1500.0/(1.0*20000.0);
}
double adimv(double v){
	return v*1500.0/20000.0;
}

double SIx(double x){
	return (1*20000.0)*x/1500.0;
}
double SIv(double v){
	return 20000.0*v/1500.0;
}

double SIa(double a){
	return 20000.0*a/(1.0*1500.0);
}

double kmhToms(double v){
	return v/3.6;
}
double msTokmh(double v){
	return 3.6*v;
}

int main()
{
	int tc=1;
	int m=1500;
	double L=4.0;
	int c=20000;
	double d=26.0;
	double v_eq=120.0;


	double dt=0.001;
	int tfrenat = 20;
	int tf=tc+tfrenat;

	int Nt=(tf/dt);

	double veqAdim = adimv(kmhToms(v_eq));
	
	double carsx[5][Nt+1];
	double carsv[5][Nt+1];	
	double carserr[5][(int)(tfrenat/dt)];

	double LAdim = adimx(L);
	double dAdim = adimx(d);

	double K0 = ((1.0/2.0*m*SIv(veqAdim)*SIv(veqAdim)));

	carsx[0][0]=5.0*LAdim+4.0*dAdim;
	carsx[1][0]=4.0*LAdim+3.0*dAdim;
	carsx[2][0]=3.0*LAdim+2.0*dAdim;
	carsx[3][0]=2.0*LAdim+1.0*dAdim;
	carsx[4][0]=1.0*LAdim+0.0*dAdim;

	carsv[0][0]=veqAdim;
	carsv[1][0]=veqAdim;
	carsv[2][0]=veqAdim;
	carsv[3][0]=veqAdim;
	carsv[4][0]=veqAdim;

	//euler method per les tc/dc posicions que tenim MRU
	
	for (int i = 0; i < tc/dt; i++)
	{
		carsx[0][i+1]=carsx[0][i]+veqAdim*dt;
		carsx[1][i+1]=carsx[1][i]+veqAdim*dt;
		carsx[2][i+1]=carsx[2][i]+veqAdim*dt;
		carsx[3][i+1]=carsx[3][i]+veqAdim*dt;
		carsx[4][i+1]=carsx[4][i]+veqAdim*dt;

		carsv[0][i+1]=veqAdim;
		carsv[1][i+1]=veqAdim;
		carsv[2][i+1]=veqAdim;
		carsv[3][i+1]=veqAdim;
		carsv[4][i+1]=veqAdim;
	}
	
	//tau = 0

	//---cas 1----

	//euler cotxe 0
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		carsx[0][i+1]=carsx[0][i]+0.2*veqAdim*dt;

		carsv[0][i+1]=0.2*veqAdim;
	}
	int crash = 0;
	int tBreak;
	double wi1T = 0;
	double wi2T = 0;
	double wi3T = 0;
	double wi4T = 0;
	//RK4
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		double K10 = 0;//K pel cotxe 0 a ordre 1
		double K20 = 0+dt/2.0*K10;
		double K30 = 0+dt/2.0*K20;
		
		double L10 = 0.2*veqAdim;//L pel cotxe 0 a ordre 1
		double L20 = 0.2*veqAdim+dt/2*K10;
		double L30 = 0.2*veqAdim+dt/2*K20;

		double K11 = -(carsv[1][i]-carsv[0][i])/fabs(carsx[1][i]-carsx[0][i]);
		double L11 = carsv[1][i];
		double K21 = - ((carsv[1][i]+dt/2.0*K11)-(carsv[0][i]+dt/2.0*K10))/fabs((carsx[1][i]+dt/2.0*L11)-(carsx[0][i]+dt/2.0*L10));
		double L21 = carsv[1][i]+dt/2*K11;
		double K31 = - ((carsv[1][i]+dt/2*K21)-(carsv[0][i]+dt/2*K20))/fabs((carsx[1][i]+dt/2*L21)-(carsx[0][i]+dt/2*L20));
		double L31 = carsv[1][i]+dt/2*K21;
		double K41 = - ((carsv[1][i]+dt*K31)-(carsv[0][i]+dt*K30))/fabs((carsx[1][i]+dt*L31)-(carsx[0][i]+dt*L30));
		double L41 = carsv[1][i]+dt*K31;

		carsx[1][i+1] = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if((carsx[1][i+1] >= carsx[0][i]-LAdim)){
			 crash = 1;	
		}
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);

		double Ki = 1.0/2.0*m*((SIv(carsv[1][i+1]))*(SIv(carsv[1][i+1]))); 
		double wi1 = m*SIa(K11)*SIx(dt/6*(L11+2*L21+2*L31+L41));
		wi1T += wi1;
		carserr[1][(int)(i-tc/dt)]=fabs(K0-(Ki-wi1T))/K0*100;
		//carserr[1][(int)(i-tc/dt)]=K0-(Ki-wi);
		if(i==tc/dt+1){
		//printf("\n1: %f %f %f", Ki,wi,K0-(Ki-wi));
		}
		double K12 = -(carsv[2][i]-carsv[1][i])/fabs(carsx[2][i]-carsx[1][i]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][i]+dt/2*K12)-(carsv[1][i]+dt/2*K11))/fabs((carsx[2][i]+dt/2*L12)-(carsx[1][i]+dt/2*L11));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][i]+dt/2*K21))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][i]+dt/2*L21));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][i]+dt*K31))/fabs((carsx[2][i]+dt*L32)-(carsx[1][i]+dt*L31));
		double L42 = carsv[2][i]+dt*K32;

		carsx[2][i+1] = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if((carsx[2][i+1] >= carsx[1][i]-LAdim)){
			 crash = 1;
		}
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);
		
		double wi2 = m*SIa(K12)*SIx(dt/6*(L12+2*L22+2*L32+L42));
		Ki = 1.0/2.0*m*((SIv(carsv[2][i+1]))*(SIv(carsv[2][i+1])));
		wi2T += wi2;
		carserr[2][(int)(i-tc/dt)]=fabs(K0-(Ki-wi2T))/K0*100;
		//carserr[2][(int)(i-tc/dt)]=K0-(Ki-wi);
		if(i==tc/dt+1){
		//printf("\n2: %f %f %f", Ki,wi,K0-(Ki-wi));
		}
		double K13 = -(carsv[3][i]-carsv[2][i])/fabs(carsx[3][i]-carsx[2][i]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][i]+dt/2*K12))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][i]+dt/2*L12));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][i]+dt/2*K22))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][i]+dt/2*L22));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][i]+dt*K32))/fabs((carsx[3][i]+dt*L33)-(carsx[2][i]+dt*L32));
		double L43 = carsv[3][i]+dt*K33;

		carsx[3][i+1] = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if((carsx[3][i+1] >= carsx[2][i]-LAdim)){
			 crash = 1;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		double wi3 = m*SIa(K13)*SIx(dt/6*(L13+2*L23+2*L33+L43));
		Ki = 1.0/2.0*m*((SIv(carsv[3][i+1]))*(SIv(carsv[3][i+1]))); 
		wi3T += wi3;
		carserr[3][(int)(i-tc/dt)]=fabs(K0-(Ki-wi3T))/K0*100;
		//carserr[3][(int)(i-tc/dt)]=K0-(Ki-wi);
		if(i==tc/dt+1){
		//printf("\n3: %f %f %f", Ki,wi,K0-(Ki-wi));
		}
		double K14 = -(carsv[4][i]-carsv[3][i])/fabs(carsx[4][i]-carsx[3][i]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][i]+dt/2*K13))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][i]+dt/2*L13));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][i]+dt/2*K23))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][i]+dt/2*L23));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][i]+dt*K33))/fabs((carsx[4][i]+dt*L34)-(carsx[3][i]+dt*L33));
		double L44 = carsv[4][i]+dt*K34;

		carsx[4][i+1] = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if ((carsx[4][i+1] >= carsx[3][i]-LAdim))
		{
			 crash = 1;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		double wi4 = m*SIa(K14)*SIx(dt/6*(L14+2*L24+2*L34+L44));
		wi4T += wi4;
		Ki = 1.0/2.0*m*((SIv(carsv[4][i+1]))*(SIv(carsv[4][i+1]))); 
		carserr[4][(int)(i-tc/dt)]=fabs(K0-(Ki-wi4T))/K0*100;
		//carserr[4][(int)(i-tc/dt)]=K0-(Ki-wi);
		if(i==tc/dt+1){
		//printf("\n4: %f %f %f", Ki,wi,K0-(Ki-wi));
		}
		//carserr[(int)(i-tc/dt)]=fabs(K0-(Ki-(wi1+wi2+wi3+wi4)))/K0*100;
		if (crash == 1)
		{	
			tBreak = i+1;
			break;
		}
		
	}
	int idx = Nt+1;
	if(crash==1){
		idx = tBreak+1;
	}
	FILE *fptr;
	fptr = fopen("C1_tau0_x.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C1_tau0_v.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",msTokmh(SIv(carsv[j][i])));
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);
	idx = (int)(tfrenat/dt);
	if(crash == 1){
		idx = tBreak-(tc/dt);
	}
	fptr = fopen("C1_tau0_err.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 1; j < 5; j++)
		{
			fprintf(fptr,"%f, ",carserr[j][i]);
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);

	//---cas 2----

	//euler cotxe 0
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		carsx[0][i+1]=carsx[0][i]+dt*(veqAdim*(1-(i*dt-tc)*exp(1-(i*dt-tc))));

		carsv[0][i+1]=veqAdim*(1-(i*dt-tc)*exp(1-(i*dt-tc)));
	}
	
	//RK4
	crash = 0;
	double wiT= 0;
	FILE *fptr2;
	fptr2 = fopen("c2_tau0_err_total.txt", "wb");
	fptr = fopen("c2_tau0_err_ctx1.txt", "wb");
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		double K10 = -veqAdim*exp(1+tc-i*dt)*(1+tc-i*dt);//K pel cotxe 0 a ordre 1
		double K20 = -veqAdim*exp(1+tc-(i*dt+dt/2))*(1+tc-(i*dt+dt/2));
		double K30 = -veqAdim*exp(1+tc-(i*dt+dt/2))*(1+tc-(i*dt+dt/2));
		
		double L10 = veqAdim*(1-(i*dt-tc)*exp(1-(i*dt-tc)));//L pel cotxe 0 a ordre 1
		double L20 = veqAdim*(1-(i*dt+dt/2-tc)*exp(1-(i*dt+dt/2-tc)));
		double L30 = veqAdim*(1-(i*dt+dt/2-tc)*exp(1-(i*dt+dt/2-tc)));
		double L40 = veqAdim*(1-(i*dt+dt-tc)*exp(1-(i*dt+dt-tc)));

		double wi = (m*SIa(K10)*SIx(dt/6*(L10+2*L20+2*L30+L40)));
		double Ki = 1.0/2.0*m*((SIv(carsv[0][i+1]))*(SIv(carsv[0][i+1]))); 

		wiT += wi;
		double KiT = Ki;

		carserr[0][(int)(i-tc/dt)]=(K0-(Ki-wi))/K0*100;

		double K11 = -(carsv[1][i]-carsv[0][i])/fabs(carsx[1][i]-carsx[0][i]);
		double L11 = carsv[1][i];
		double K21 = - ((carsv[1][i]+dt/2*K11)-(carsv[0][i]+dt/2*K10))/fabs((carsx[1][i]+dt/2*L11)-(carsx[0][i]+dt/2*L10));
		double L21 = carsv[1][i]+dt/2*K11;
		double K31 = - ((carsv[1][i]+dt/2*K21)-(carsv[0][i]+dt/2*K20))/fabs((carsx[1][i]+dt/2*L21)-(carsx[0][i]+dt/2*L20));
		double L31 = carsv[1][i]+dt/2*K21;
		double K41 = - ((carsv[1][i]+dt*K31)-(carsv[0][i]+dt*K30))/fabs((carsx[1][i]+dt*L31)-(carsx[0][i]+dt*L30));
		double L41 = carsv[1][i]+dt*K31;

		carsx[1][i+1] = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if ((carsx[1][i+1] >= carsx[0][i]-LAdim))
		{
			crash = 1;
		}
		
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);
		
		wi = (m*SIa(K11)*SIx(dt/6*(L11+2*L21+2*L31+L41)));
		Ki = 1.0/2.0*m*((SIv(carsv[1][i+1]))*(SIv(carsv[1][i+1]))); 
		wiT += wi;
		KiT += Ki;
		//printf("K0: %f, ki: %f, wi: %f,Ei: %f\n\n\n",K0,Ki,wi, Ki-wi);
		fprintf(fptr,"K0: %f, ki: %f, wi: %f,Ei: %f\n",K0,Ki,wi, Ki-wi);
		carserr[1][(int)(i-tc/dt)]=(K0-(Ki-wi))/K0*100;

		double K12 = -(carsv[2][i]-carsv[1][i])/fabs(carsx[2][i]-carsx[1][i]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][i]+dt/2*K12)-(carsv[1][i]+dt/2*K11))/fabs((carsx[2][i]+dt/2*L12)-(carsx[1][i]+dt/2*L11));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][i]+dt/2*K21))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][i]+dt/2*L21));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][i]+dt*K31))/fabs((carsx[2][i]+dt*L32)-(carsx[1][i]+dt*L31));
		double L42 = carsv[2][i]+dt*K32;

		carsx[2][i+1] = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if ((carsx[2][i+1] >= carsx[1][i]-LAdim))
		{
			 crash = 1;
		}
		
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);
		
		wi = (m*SIa(K12)*SIx(dt/6*(L12+2*L22+2*L32+L42)));
		Ki = 1.0/2.0*m*((SIv(carsv[2][i+1]))*(SIv(carsv[2][i+1]))); 
		wiT += wi;
		KiT += Ki;

		carserr[2][(int)(i-tc/dt)]=(K0-(Ki-wi))/K0*100;

		double K13 = -(carsv[3][i]-carsv[2][i])/fabs(carsx[3][i]-carsx[2][i]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][i]+dt/2*K12))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][i]+dt/2*L12));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][i]+dt/2*K22))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][i]+dt/2*L22));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][i]+dt*K32))/fabs((carsx[3][i]+dt*L33)-(carsx[2][i]+dt*L32));
		double L43 = carsv[3][i]+dt*K33;

		carsx[3][i+1]=carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if((carsx[3][i+1]>=carsx[2][i]-LAdim)){
			 crash = 1;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		
		wi = (m*SIa(K13)*SIx(dt/6*(L13+2*L23+2*L33+L43)));
		Ki = 1.0/2.0*m*((SIv(carsv[3][i+1]))*(SIv(carsv[3][i+1]))); 
		wiT += wi;
		KiT += Ki;

		carserr[3][(int)(i-tc/dt)]=(K0-(Ki-wi))/K0*100;

		double K14 = -(carsv[4][i]-carsv[3][i])/fabs(carsx[4][i]-carsx[3][i]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][i]+dt/2*K13))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][i]+dt/2*L13));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][i]+dt/2*K23))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][i]+dt/2*L23));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][i]+dt*K33))/fabs((carsx[4][i]+dt*L34)-(carsx[3][i]+dt*L33));
		double L44 = carsv[4][i]+dt*K34;


		carsx[4][i+1] = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if((carsx[4][i+1]>=carsx[3][i]-LAdim)){
			 crash = 1;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		
		wi = (m*SIa(K14)*SIx(dt/6*(L14+2*L24+2*L34+L44)));
		Ki = 1.0/2.0*m*((SIv(carsv[4][i+1]))*(SIv(carsv[4][i+1]))); 
		wiT += wi;
		KiT += Ki;

		fprintf(fptr2,"%f\n",fabs(5*K0-KiT+wiT)/(K0*5)*100);

		carserr[4][(int)(i-tc/dt)]=(K0-(Ki-wi))/K0*100;


		if(crash == 1){
			tBreak = i+1;
			break;
		}
	}
	fclose(fptr);
	fclose(fptr2);
	//printf("\n\n c2 tau 0 %d \n\n",tBreak);
	//printf("\n\n pos car 2 break %f \n\n",carsx[1][tBreak]);
	idx = Nt+1;
	if (crash == 1)
	{	
		idx = tBreak+1;
	}
	
	fptr = fopen("C2_tau0_x.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C2_tau0_v.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",msTokmh(SIv(carsv[j][i])));
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);
	idx = (int)(tfrenat/dt)+1;
	if(crash == 1){
		idx = tBreak-(tc/dt)+1;
	}
	fptr = fopen("C2_tau0_err.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",carserr[j][i]);
		}
		
		fprintf(fptr, "\n");
	}

	fclose(fptr);

	//---cas 3----
	double omega = 0.2;
	//euler cotxe 0
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		carsx[0][i+1]=carsx[0][i]+dt*(veqAdim*(1-0.8*sin(omega*(i*dt-tc)*sin(omega*(i*dt-tc)))));

		carsv[0][i+1]=veqAdim*(1-0.8*(sin(omega*(i*dt-tc))*sin(omega*(i*dt-tc))));
	}
	
	//RK4
	crash = 0;
	wiT= 0;
	fptr2 = fopen("c3_tau0_err_total.txt", "wb");
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		double K10 = 0.8*veqAdim*omega*sin(2*omega*(tc-dt*i));//K pel cotxe 0 a ordre 1
		double K20 = 0.8*veqAdim*omega*sin(2*omega*(tc-(dt*i+dt/2)));
		double K30 = K20;
		
		double L10 = veqAdim*(1-0.8*sin(omega*(i*dt-tc)*sin(omega*(i*dt-tc))));//L pel cotxe 0 a ordre 1
		double L20 = veqAdim*(1-0.8*sin(omega*(i*dt+dt/2-tc)*sin(omega*(i*dt+dt/2-tc))));
		double L30 = L20;
		double L40 = veqAdim*(1-0.8*sin(omega*(i*dt+dt-tc)*sin(omega*(i*dt+dt-tc))));

		double wi = (m*SIa(K10)*SIx(dt/6*(L10+2*L20+2*L30+L40)));
		double Ki = 1.0/2.0*m*((SIv(carsv[0][i+1]))*(SIv(carsv[0][i+1]))); 

		wiT += wi;
		double KiT = Ki;

		double K11 = -(carsv[1][i]-carsv[0][i])/fabs(carsx[1][i]-carsx[0][i]);
		double L11 = carsv[1][i];
		double K21 = - ((carsv[1][i]+dt/2*K11)-(carsv[0][i]+dt/2*K10))/fabs((carsx[1][i]+dt/2*L11)-(carsx[0][i]+dt/2*L10));
		double L21 = carsv[1][i]+dt/2*K11;
		double K31 = - ((carsv[1][i]+dt/2*K21)-(carsv[0][i]+dt/2*K20))/fabs((carsx[1][i]+dt/2*L21)-(carsx[0][i]+dt/2*L20));
		double L31 = carsv[1][i]+dt/2*K21;
		double K41 = - ((carsv[1][i]+dt*K31)-(carsv[0][i]+dt*K30))/fabs((carsx[1][i]+dt*L31)-(carsx[0][i]+dt*L30));
		double L41 = carsv[1][i]+dt*K31;

		carsx[1][i+1] = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if(carsx[1][i+1] >= carsx[0][i]-LAdim){
			crash = 1;
		}

		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);
		wi = (m*SIa(K11)*SIx(dt/6*(L11+2*L21+2*L31+L41)));
		Ki = 1.0/2.0*m*((SIv(carsv[1][i+1]))*(SIv(carsv[1][i+1]))); 

		wiT += wi;
		KiT += Ki;

		double K12 = -(carsv[2][i]-carsv[1][i])/fabs(carsx[2][i]-carsx[1][i]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][i]+dt/2*K12)-(carsv[1][i]+dt/2*K11))/fabs((carsx[2][i]+dt/2*L12)-(carsx[1][i]+dt/2*L11));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][i]+dt/2*K21))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][i]+dt/2*L21));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][i]+dt*K31))/fabs((carsx[2][i]+dt*L32)-(carsx[1][i]+dt*L31));
		double L42 = carsv[2][i]+dt*K32;

		carsx[2][i+1] = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if(carsx[2][i+1] >= carsx[1][i]-LAdim){
			 crash = 1;
		}
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);
		wi = (m*SIa(K12)*SIx(dt/6*(L12+2*L22+2*L32+L42)));
		Ki = 1.0/2.0*m*((SIv(carsv[2][i+1]))*(SIv(carsv[2][i+1]))); 

		wiT += wi;
		KiT += Ki;

		double K13 = -(carsv[3][i]-carsv[2][i])/fabs(carsx[3][i]-carsx[2][i]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][i]+dt/2*K12))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][i]+dt/2*L12));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][i]+dt/2*K22))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][i]+dt/2*L22));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][i]+dt*K32))/fabs((carsx[3][i]+dt*L33)-(carsx[2][i]+dt*L32));
		double L43 = carsv[3][i]+dt*K33;

		carsx[3][i+1] = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if(carsx[3][i+1] >= carsx[2][i]-LAdim){
			 crash = 1;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		wi = (m*SIa(K13)*SIx(dt/6*(L13+2*L23+2*L33+L43)));
		Ki = 1.0/2.0*m*((SIv(carsv[3][i+1]))*(SIv(carsv[3][i+1]))); 

		wiT += wi;
		KiT += Ki;

		double K14 = -(carsv[4][i]-carsv[3][i])/fabs(carsx[4][i]-carsx[3][i]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][i]+dt/2*K13))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][i]+dt/2*L13));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][i]+dt/2*K23))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][i]+dt/2*L23));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][i]+dt*K33))/fabs((carsx[4][i]+dt*L34)-(carsx[3][i]+dt*L33));
		double L44 = carsv[4][i]+dt*K34;

		carsx[4][i+1] = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if(carsx[4][i+1] >= carsx[3][i]-LAdim){
			crash = 1;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		wi = (m*SIa(K14)*SIx(dt/6*(L14+2*L24+2*L34+L44)));
		Ki = 1.0/2.0*m*((SIv(carsv[4][i+1]))*(SIv(carsv[4][i+1]))); 

		wiT += wi;
		KiT += Ki;

		fprintf(fptr2,"%f\n",fabs(5*K0-KiT+wiT)/(K0*5)*100);

		if(crash == 1){
			tBreak = i+1;
			break;
		}
	}
	fclose(fptr2);
	idx = Nt+1;
	if (crash == 1)
	{	
		idx = tBreak+1;
	}
	fptr = fopen("C3_tau0_x.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C3_tau0_v.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",msTokmh(SIv(carsv[j][i])));
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);

	//tau = 0.5	
	double tau = 0.3;

	//---cas 1----

	//euler cotxe 0
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		carsx[0][i+1]=carsx[0][i]+0.2*veqAdim*dt;

		carsv[0][i+1]=0.2*veqAdim;
	}
	
	//RK4
	crash = 0;
	wi1T = 0;
	wi2T = 0;
	wi3T = 0;
	wi4T = 0;
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		double K11 = -(carsv[1][i]-carsv[0][(int)((i*dt-tau)/dt)])/fabs(carsx[1][i]-carsx[0][(int)((i*dt-tau)/dt)]);
		double L11 = carsv[1][i];
		double K21 = - ((carsv[1][i]+dt/2*K11)-(carsv[0][(int)((i*dt-tau)/dt)]))/fabs((carsx[1][i]+dt/2*L11)-(carsx[0][(int)((i*dt-tau)/dt)]));
		double L21 = carsv[1][i]+dt/2*K11;
		double K31 = - ((carsv[1][i]+dt/2*K21)-(carsv[0][(int)((i*dt-tau)/dt)]))/fabs((carsx[1][i]+dt/2*L21)-(carsx[0][(int)((i*dt-tau)/dt)]));
		double L31 = carsv[1][i]+dt/2*K21;
		double K41 = - ((carsv[1][i]+dt*K31)-(carsv[0][(int)((i*dt-tau)/dt)]))/fabs((carsx[1][i]+dt*L31)-(carsx[0][(int)((i*dt-tau)/dt)]));
		double L41 = carsv[1][i]+dt*K31;

		carsx[1][i+1] = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if(carsx[1][i+1] >= carsx[0][i]-LAdim){
			 crash = 1;
		}
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);
		double Ki = 1.0/2.0*m*((SIv(carsv[1][i+1]))*(SIv(carsv[1][i+1]))); 
		double wi1 = m*SIa(K11)*SIx(dt/6*(L11+2*L21+2*L31+L41));
		wi1T += wi1;
		carserr[1][(int)(i-tc/dt)]=fabs(K0-(Ki-wi1T))/K0*100;


		double K12 = -(carsv[2][i]-carsv[1][(int)((i*dt-tau)/dt)])/fabs(carsx[2][i]-carsx[1][(int)((i*dt-tau)/dt)]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][(int)((i*dt-tau)/dt)])-(carsv[1][i]+dt/2*K11))/fabs((carsx[2][(int)((i*dt-tau)/dt)])-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt*L32)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L42 = carsv[2][i]+dt*K32;

		carsx[2][i+1] = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if(carsx[2][i+1] >= carsx[1][i]-LAdim){
			 crash = 1;
		}
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);
		double wi2 = m*SIa(K12)*SIx(dt/6*(L12+2*L22+2*L32+L42));
		Ki = 1.0/2.0*m*((SIv(carsv[2][i+1]))*(SIv(carsv[2][i+1]))); 
		wi2T += wi2;
		carserr[2][(int)(i-tc/dt)]=fabs(K0-(Ki-wi2T))/K0*100;

		double K13 = -(carsv[3][i]-carsv[2][(int)((i*dt-tau)/dt)])/fabs(carsx[3][i]-carsx[2][(int)((i*dt-tau)/dt)]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt*L33)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L43 = carsv[3][i]+dt*K33;

		carsx[3][i+1] = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if(carsx[3][i+1] >= carsx[2][i]-LAdim){
			 crash = 1;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		double wi3 = m*SIa(K13)*SIx(dt/6*(L13+2*L23+2*L33+L43));
		Ki = 1.0/2.0*m*((SIv(carsv[3][i+1]))*(SIv(carsv[3][i+1]))); 
		wi3T += wi3;
		carserr[3][(int)(i-tc/dt)]=fabs(K0-(Ki-wi3T))/K0*100;

		double K14 = -(carsv[4][i]-carsv[3][(int)((i*dt-tau)/dt)])/fabs(carsx[4][i]-carsx[3][(int)((i*dt-tau)/dt)]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt*L34)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L44 = carsv[4][i]+dt*K34;

		carsx[4][i+1] = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if(carsx[4][i+1] >= carsx[3][i]-LAdim){
			 crash = 1;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		double wi4 = m*SIa(K14)*SIx(dt/6*(L14+2*L24+2*L34+L44));
		Ki = 1.0/2.0*m*((SIv(carsv[4][i+1]))*(SIv(carsv[4][i+1]))); 
		wi4T += wi4;
		carserr[4][(int)(i-tc/dt)]=fabs(K0-(Ki-wi4T))/K0*100;

		if(crash == 1){
			tBreak = i+1;
			break;
		}
	}
	idx = Nt+1;
	if (crash == 1)
	{	
		idx = tBreak+1;
	}
	fptr = fopen("C1_tau_x.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C1_tau_v.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",msTokmh(SIv(carsv[j][i])));
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);

	idx = (int)(tfrenat/dt);
	if(crash == 1){
		idx = tBreak-(tc/dt);
	}
	fptr = fopen("C1_tau_err.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 1; j < 5; j++)
		{
			fprintf(fptr,"%f, ",carserr[j][i]);
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);

	//---cas 2----

	//euler cotxe 0
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		carsx[0][i+1]=carsx[0][i]+dt*(veqAdim*(1-(i*dt-tc)*exp(1-(i*dt-tc))));

		carsv[0][i+1]=veqAdim*(1-(i*dt-tc)*exp(1-(i*dt-tc)));
	}
	
	//RK4
	crash = 0;
	wiT= 0;
	fptr2 = fopen("c2_tau_err_total.txt", "wb");
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		double K10 = -veqAdim*exp(1+tc-i*dt)*(1+tc-i*dt);//K pel cotxe 0 a ordre 1
		double K20 = -veqAdim*exp(1+tc-(i*dt+dt/2))*(1+tc-(i*dt+dt/2));
		double K30 = -veqAdim*exp(1+tc-(i*dt+dt/2))*(1+tc-(i*dt+dt/2));

		double L10 = veqAdim*(1-(i*dt-tc)*exp(1-(i*dt-tc)));//L pel cotxe 0 a ordre 1
		double L20 = veqAdim*(1-(i*dt+dt/2-tc)*exp(1-(i*dt+dt/2-tc)));
		double L30 = veqAdim*(1-(i*dt+dt/2-tc)*exp(1-(i*dt+dt/2-tc)));
		double L40 = veqAdim*(1-(i*dt+dt-tc)*exp(1-(i*dt+dt-tc)));

		double wi = (m*SIa(K10)*SIx(dt/6*(L10+2*L20+2*L30+L40)));
		double Ki = 1.0/2.0*m*((SIv(carsv[0][i+1]))*(SIv(carsv[0][i+1]))); 

		wiT += wi;
		double KiT = Ki;

		double K11 = -(carsv[1][i]-carsv[0][(int)((i*dt-tau)/dt)])/fabs(carsx[1][i]-carsx[0][(int)((i*dt-tau)/dt)]);
		double L11 = carsv[1][i];
		double K21 = - ((carsv[1][i]+dt/2*K11)-(carsv[0][(int)((i*dt-tau)/dt)]))/fabs((carsx[1][i]+dt/2*L11)-(carsx[0][(int)((i*dt-tau)/dt)]));
		double L21 = carsv[1][i]+dt/2*K11;
		double K31 = - ((carsv[1][i]+dt/2*K21)-(carsv[0][(int)((i*dt-tau)/dt)]))/fabs((carsx[1][i]+dt/2*L21)-(carsx[0][(int)((i*dt-tau)/dt)]));
		double L31 = carsv[1][i]+dt/2*K21;
		double K41 = - ((carsv[1][i]+dt*K31)-(carsv[0][(int)((i*dt-tau)/dt)]))/fabs((carsx[1][i]+dt*L31)-(carsx[0][(int)((i*dt-tau)/dt)]));
		double L41 = carsv[1][i]+dt*K31;

		carsx[1][i+1] = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if(carsx[1][i+1]>=carsx[0][i]-LAdim){
			 crash = 1;	
		}
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);
		
		wi = (m*SIa(K11)*SIx(dt/6*(L11+2*L21+2*L31+L41)));
		Ki = 1.0/2.0*m*((SIv(carsv[1][i+1]))*(SIv(carsv[1][i+1]))); 
		wiT += wi;
		KiT += Ki;


		double K12 = -(carsv[2][i]-carsv[1][(int)((i*dt-tau)/dt)])/fabs(carsx[2][i]-carsx[1][(int)((i*dt-tau)/dt)]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][i]+dt/2*K12)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt/2*L12)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt*L32)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L42 = carsv[2][i]+dt*K32;

		carsx[2][i+1] = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if(carsx[2][i+1] >= carsx[1][i]-LAdim){
			 crash = 1;	
		}
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);

		wi = (m*SIa(K12)*SIx(dt/6*(L12+2*L22+2*L32+L42)));
		Ki = 1.0/2.0*m*((SIv(carsv[2][i+1]))*(SIv(carsv[2][i+1]))); 
		wiT += wi;
		KiT += Ki;


		double K13 = -(carsv[3][i]-carsv[2][(int)((i*dt-tau)/dt)])/fabs(carsx[3][i]-carsx[2][(int)((i*dt-tau)/dt)]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt*L33)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L43 = carsv[3][i]+dt*K33;

		carsx[3][i+1] = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if(carsx[3][i+1] >= carsx[2][i]-LAdim){
			 crash = 1;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		
		wi = (m*SIa(K13)*SIx(dt/6*(L13+2*L23+2*L33+L43)));
		Ki = 1.0/2.0*m*((SIv(carsv[3][i+1]))*(SIv(carsv[3][i+1]))); 
		wiT += wi;
		KiT += Ki;

		double K14 = -(carsv[4][i]-carsv[3][(int)((i*dt-tau)/dt)])/fabs(carsx[4][i]-carsx[3][(int)((i*dt-tau)/dt)]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt*L34)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L44 = carsv[4][i]+dt*K34;

		carsx[4][i+1] = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if(carsx[4][i+1] >= carsx[3][i]-LAdim){
			crash = 1;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);

		wi = (m*SIa(K14)*SIx(dt/6*(L14+2*L24+2*L34+L44)));
		Ki = 1.0/2.0*m*((SIv(carsv[4][i+1]))*(SIv(carsv[4][i+1]))); 
		wiT += wi;
		KiT += Ki;
		
		fprintf(fptr2,"%f\n",fabs(5*K0-KiT+wiT)/(K0*5)*100);

		if(crash == 1){
			tBreak = i+1;
			break;
		}
	}
	fclose(fptr2);

	idx = Nt+1;
	if (crash == 1)
	{	
		idx = tBreak+1;
	}
	fptr = fopen("C2_tau_x.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C2_tau_v.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",msTokmh(SIv(carsv[j][i])));
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);

	//---cas 3----
	omega = 0.2;
	//euler cotxe 0
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		carsx[0][i+1]=carsx[0][i]+dt*(veqAdim*(1-0.8*sin(omega*(i*dt-tc)*sin(omega*(i*dt-tc)))));

		carsv[0][i+1]=veqAdim*(1-0.8*(sin(omega*(i*dt-tc))*sin(omega*(i*dt-tc))));
	}
	
	//RK4
	crash = 0;
	wiT= 0;
	fptr2 = fopen("c3_tau_err_total.txt", "wb");
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		double K10 = 0.8*veqAdim*omega*sin(2*omega*(tc-dt*i));//K pel cotxe 0 a ordre 1
		double K20 = 0.8*veqAdim*omega*sin(2*omega*(tc-(dt*i+dt/2)));
		double K30 = K20;
		
		double L10 = veqAdim*(1-0.8*sin(omega*(i*dt-tc)*sin(omega*(i*dt-tc))));//L pel cotxe 0 a ordre 1
		double L20 = veqAdim*(1-0.8*sin(omega*(i*dt+dt/2-tc)*sin(omega*(i*dt+dt/2-tc))));
		double L30 = L20;
		double L40 = veqAdim*(1-0.8*sin(omega*(i*dt+dt-tc)*sin(omega*(i*dt+dt-tc))));

		double wi = (m*SIa(K10)*SIx(dt/6*(L10+2*L20+2*L30+L40)));
		double Ki = 1.0/2.0*m*((SIv(carsv[0][i+1]))*(SIv(carsv[0][i+1]))); 

		wiT += wi;
		double KiT = Ki;

		double K11 = -(carsv[1][i]-carsv[0][(int)((i*dt-tau)/dt)])/fabs(carsx[1][i]-carsx[0][(int)((i*dt-tau)/dt)]);
		double L11 = carsv[1][i];
		double K21 = - ((carsv[1][i]+dt/2*K11)-(carsv[0][(int)((i*dt-tau)/dt)]))/fabs((carsx[1][i]+dt/2*L11)-(carsx[0][(int)((i*dt-tau)/dt)]));
		double L21 = carsv[1][i]+dt/2*K11;
		double K31 = - ((carsv[1][i]+dt/2*K21)-(carsv[0][(int)((i*dt-tau)/dt)]))/fabs((carsx[1][i]+dt/2*L21)-(carsx[0][(int)((i*dt-tau)/dt)]));
		double L31 = carsv[1][i]+dt/2*K21;
		double K41 = - ((carsv[1][i]+dt*K31)-(carsv[0][(int)((i*dt-tau)/dt)]))/fabs((carsx[1][i]+dt*L31)-(carsx[0][(int)((i*dt-tau)/dt)]));
		double L41 = carsv[1][i]+dt*K31;

		carsx[1][i+1] = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if(carsx[1][i+1] >= carsx[0][i]-LAdim){
			 crash = 1;
		}
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);
		
		wi = (m*SIa(K11)*SIx(dt/6*(L11+2*L21+2*L31+L41)));
		Ki = 1.0/2.0*m*((SIv(carsv[1][i+1]))*(SIv(carsv[1][i+1]))); 

		wiT += wi;
		KiT += Ki;


		double K12 = -(carsv[2][i]-carsv[1][(int)((i*dt-tau)/dt)])/fabs(carsx[2][i]-carsx[1][(int)((i*dt-tau)/dt)]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][i]+dt/2*K12)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt/2*L12)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt*L32)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L42 = carsv[2][i]+dt*K32;

		carsx[2][i+1] = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if(carsx[2][i+1] >= carsx[1][i]-LAdim){
			 crash = 1;
		}
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);
		
		wi = (m*SIa(K12)*SIx(dt/6*(L12+2*L22+2*L32+L42)));
		Ki = 1.0/2.0*m*((SIv(carsv[2][i+1]))*(SIv(carsv[2][i+1]))); 

		wiT += wi;
		KiT += Ki;


		double K13 = -(carsv[3][i]-carsv[2][(int)((i*dt-tau)/dt)])/fabs(carsx[3][i]-carsx[2][(int)((i*dt-tau)/dt)]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt*L33)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L43 = carsv[3][i]+dt*K33;

		carsx[3][i+1] = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if(carsx[3][i+1] >= carsx[2][i]-LAdim){
			 crash = 1;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		
		wi = (m*SIa(K13)*SIx(dt/6*(L13+2*L23+2*L33+L43)));
		Ki = 1.0/2.0*m*((SIv(carsv[3][i+1]))*(SIv(carsv[3][i+1]))); 

		wiT += wi;
		KiT += Ki;

		double K14 = -(carsv[4][i]-carsv[3][(int)((i*dt-tau)/dt)])/fabs(carsx[4][i]-carsx[3][(int)((i*dt-tau)/dt)]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt*L34)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L44 = carsv[4][i]+dt*K34;

		carsx[4][i+1] = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if(carsx[4][i+1] >= carsx[3][i]-LAdim){
			 crash = 1;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		
		wi = (m*SIa(K14)*SIx(dt/6*(L14+2*L24+2*L34+L44)));
		Ki = 1.0/2.0*m*((SIv(carsv[4][i+1]))*(SIv(carsv[4][i+1]))); 

		wiT += wi;
		KiT += Ki;

		fprintf(fptr2,"%f\n",fabs(5*K0-KiT+wiT)/(K0*5)*100);

		if(crash == 1){
			tBreak = i+1;
			break;
		}
	}
	fclose(fptr2);
	
	idx = Nt+1;
	if (crash == 1)
	{	
		idx = tBreak+1;
	}
	fptr = fopen("C3_tau_x.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C3_tau_v.txt", "wb");
	for (int i = 0; i < idx; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",msTokmh(SIv(carsv[j][i])));
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);


	// ----- IDM -----


	double v0 = kmhToms(120);
	double T = 1.5;
	double a = 0.73;
	double b = 1.67;
	double delta = 4;
	double s0 = 5;

	LAdim = L/s0;
	dAdim = d/s0;
	veqAdim = kmhToms(v_eq)/v0;

	omega = 0.4;

	carsx[0][0]=5.0*LAdim+4.0*dAdim;
	carsx[1][0]=4.0*LAdim+3.0*dAdim;
	carsx[2][0]=3.0*LAdim+2.0*dAdim;
	carsx[3][0]=2.0*LAdim+1.0*dAdim;
	carsx[4][0]=1.0*LAdim+0.0*dAdim;

	carsv[0][0]=veqAdim;
	carsv[1][0]=veqAdim;
	carsv[2][0]=veqAdim;
	carsv[3][0]=veqAdim;
	carsv[4][0]=veqAdim;

	//euler method per les tc/dc posicions que tenim MRU
	
	for (int i = 0; i < tc/dt; i++)
	{
		carsx[0][i+1]=carsx[0][i]+veqAdim*dt;
		carsx[1][i+1]=carsx[1][i]+veqAdim*dt;
		carsx[2][i+1]=carsx[2][i]+veqAdim*dt;
		carsx[3][i+1]=carsx[3][i]+veqAdim*dt;
		carsx[4][i+1]=carsx[4][i]+veqAdim*dt;

		carsv[0][i+1]=veqAdim;
		carsv[1][i+1]=veqAdim;
		carsv[2][i+1]=veqAdim;
		carsv[3][i+1]=veqAdim;
		carsv[4][i+1]=veqAdim;
	}
	//euler cotxe 0
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		carsx[0][i+1]=carsx[0][i]+dt*(veqAdim*(1-0.8*sin(omega*(i*dt-tc)*sin(omega*(i*dt-tc)))));

		carsv[0][i+1]=veqAdim*(1-0.8*(sin(omega*(i*dt-tc))*sin(omega*(i*dt-tc))));
	}
	
	//RK4
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		double K10 = 0.8*veqAdim*omega*sin(2*omega*(tc-dt*i));//K pel cotxe 0 a ordre 1
		double K20 = 0.8*veqAdim*omega*sin(2*omega*(tc-(dt*i+dt/2)));
		double K30 = K20;
		
		double L10 = veqAdim*(1-0.8*sin(omega*(i*dt-tc)*sin(omega*(i*dt-tc))));//L pel cotxe 0 a ordre 1
		double L20 = veqAdim*(1-0.8*sin(omega*(i*dt+dt/2-tc)*sin(omega*(i*dt+dt/2-tc))));
		double L30 = L20;

		double K11 = a*T*T/s0*(1-pow(carsv[1][i],delta)-pow(((1+v0*T/s0*carsv[1][i]+v0*v0/s0*(carsv[1][i]*(carsv[1][i]-carsv[0][i]))/(2*sqrt(a*b)))/(carsx[0][i]-carsx[1][i]-L/s0)),2));
		double L11 = carsv[1][i];
		double K21 = a*T*T/s0*(1-pow(carsv[1][i]+dt/2*K11,delta)-pow(((1+v0*T/s0*(carsv[1][i]+dt/2*K11)+v0*v0/s0*((carsv[1][i]+dt/2*K11)*(carsv[1][i]+dt/2*K11-(carsv[0][i]+dt/2*K10)))/(2*sqrt(a*b)))/(carsx[0][i]+dt/2*L10-(carsx[1][i]+dt/2*L11)-L/s0)),2));
		double L21 = carsv[1][i]+dt/2*K11;
		double K31 = a*T*T/s0*(1-pow(carsv[1][i]+dt/2*K21,delta)-pow(((1+v0*T/s0*(carsv[1][i]+dt/2*K21)+v0*v0/s0*((carsv[1][i]+dt/2*K21)*(carsv[1][i]+dt/2*K21-(carsv[0][i]+dt/2*K20)))/(2*sqrt(a*b)))/(carsx[0][i]+dt/2*L20-(carsx[1][i]+dt/2*L21)-L/s0)),2));
		double L31 = carsv[1][i]+dt/2*K21;
		double K41 = a*T*T/s0*(1-pow(carsv[1][i]+dt/2*K11,delta)-pow(((1+v0*T/s0*(carsv[1][i]+dt*K31)+v0*v0/s0*((carsv[1][i]+dt*K31)*(carsv[1][i]+dt*K31-(carsv[0][i]+dt*K30)))/(2*sqrt(a*b)))/(carsx[0][i]+dt*L30-(carsx[1][i]+dt*L31)-L/s0)),2));
		double L41 = carsv[1][i]+dt*K31;

		carsx[1][i+1] = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);


		double K12 = a*T*T/s0*(1-pow(carsv[2][i],delta)-pow(((1+v0*T/s0*carsv[2][i]+v0*v0/s0*(carsv[2][i]*(carsv[2][i]-carsv[1][i]))/(2*sqrt(a*b)))/(carsx[1][i]-carsx[2][i]-L/s0)),2));
		double L12 = carsv[2][i];
		double K22 = a*T*T/s0*(1-pow(carsv[2][i]+dt/2*K12,delta)-pow(((1+v0*T/s0*(carsv[2][i]+dt/2*K12)+v0*v0/s0*((carsv[2][i]+dt/2*K12)*(carsv[2][i]+dt/2*K12-(carsv[1][i]+dt/2*K11)))/(2*sqrt(a*b)))/(carsx[1][i]+dt/2*L11-(carsx[2][i]+dt/2*L12)-L/s0)),2));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = a*T*T/s0*(1-pow(carsv[2][i]+dt/2*K22,delta)-pow(((1+v0*T/s0*(carsv[2][i]+dt/2*K22)+v0*v0/s0*((carsv[2][i]+dt/2*K22)*(carsv[2][i]+dt/2*K22-(carsv[1][i]+dt/2*K21)))/(2*sqrt(a*b)))/(carsx[1][i]+dt/2*L21-(carsx[2][i]+dt/2*L22)-L/s0)),2));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = a*T*T/s0*(1-pow(carsv[2][i]+dt/2*K32,delta)-pow(((1+v0*T/s0*(carsv[2][i]+dt*K32)+v0*v0/s0*((carsv[2][i]+dt*K32)*(carsv[2][i]+dt*K32-(carsv[1][i]+dt*K31)))/(2*sqrt(a*b)))/(carsx[1][i]+dt*L31-(carsx[2][i]+dt*L32)-L/s0)),2));
		double L42 = carsv[2][i]+dt*K32;

		carsx[2][i+1] = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);


		double K13 = a*T*T/s0*(1-pow(carsv[3][i],delta)-pow(((1+v0*T/s0*carsv[3][i]+v0*v0/s0*(carsv[3][i]*(carsv[3][i]-carsv[2][i]))/(2*sqrt(a*b)))/(carsx[2][i]-carsx[3][i]-L/s0)),2));
		double L13 = carsv[3][i];
		double K23 = a*T*T/s0*(1-pow(carsv[3][i]+dt/2*K13,delta)-pow(((1+v0*T/s0*(carsv[3][i]+dt/2*K13)+v0*v0/s0*((carsv[3][i]+dt/2*K13)*(carsv[3][i]+dt/2*K13-(carsv[2][i]+dt/2*K12)))/(2*sqrt(a*b)))/(carsx[2][i]+dt/2*L12-(carsx[3][i]+dt/2*L13)-L/s0)),2));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = a*T*T/s0*(1-pow(carsv[3][i]+dt/2*K23,delta)-pow(((1+v0*T/s0*(carsv[3][i]+dt/2*K23)+v0*v0/s0*((carsv[3][i]+dt/2*K23)*(carsv[3][i]+dt/2*K23-(carsv[2][i]+dt/2*K22)))/(2*sqrt(a*b)))/(carsx[2][i]+dt/2*L22-(carsx[3][i]+dt/2*L23)-L/s0)),2));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = a*T*T/s0*(1-pow(carsv[3][i]+dt/2*K33,delta)-pow(((1+v0*T/s0*(carsv[3][i]+dt*K33)+v0*v0/s0*((carsv[3][i]+dt*K33)*(carsv[3][i]+dt*K33-(carsv[2][i]+dt*K32)))/(2*sqrt(a*b)))/(carsx[2][i]+dt*L32-(carsx[3][i]+dt*L33)-L/s0)),2));
		double L43 = carsv[3][i]+dt*K33;

		carsx[3][i+1] = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);

		double K14 = a*T*T/s0*(1-pow(carsv[4][i],delta)-pow(((1+v0*T/s0*carsv[4][i]+v0*v0/s0*(carsv[4][i]*(carsv[4][i]-carsv[3][i]))/(2*sqrt(a*b)))/(carsx[3][i]-carsx[4][i]-L/s0)),2));
		double L14 = carsv[4][i];
		double K24 = a*T*T/s0*(1-pow(carsv[4][i]+dt/2*K14,delta)-pow(((1+v0*T/s0*(carsv[4][i]+dt/2*K14)+v0*v0/s0*((carsv[4][i]+dt/2*K14)*(carsv[4][i]+dt/2*K14-(carsv[3][i]+dt/2*K13)))/(2*sqrt(a*b)))/(carsx[3][i]+dt/2*L13-(carsx[4][i]+dt/2*L14)-L/s0)),2));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = a*T*T/s0*(1-pow(carsv[4][i]+dt/2*K24,delta)-pow(((1+v0*T/s0*(carsv[4][i]+dt/2*K24)+v0*v0/s0*((carsv[4][i]+dt/2*K24)*(carsv[4][i]+dt/2*K24-(carsv[3][i]+dt/2*K23)))/(2*sqrt(a*b)))/(carsx[3][i]+dt/2*L23-(carsx[4][i]+dt/2*L24)-L/s0)),2));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = a*T*T/s0*(1-pow(carsv[4][i]+dt/2*K34,delta)-pow(((1+v0*T/s0*(carsv[4][i]+dt*K34)+v0*v0/s0*((carsv[4][i]+dt*K34)*(carsv[4][i]+dt*K34-(carsv[3][i]+dt*K33)))/(2*sqrt(a*b)))/(carsx[3][i]+dt*L33-(carsx[4][i]+dt*L34)-L/s0)),2));
		double L44 = carsv[4][i]+dt*K34;

		carsx[4][i+1] = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		
	}

	fptr = fopen("x_IDM.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",s0*(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("v_IDM.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",msTokmh(v0*(carsv[j][i])));
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);


}