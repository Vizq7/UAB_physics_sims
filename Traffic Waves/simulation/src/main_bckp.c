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
/*
void RK4tau0Cas1(double xi, double xant, double vi, double vant, double dt, double vect[]){
	double K10 = ;//K pel cotxe 0 a ordre 1
	double K20 = ;
	double K30 = ;
	double K40 = ;
	
	double K1 = -(vi-vant)/fabs(xi-xant);
	double L1 = vi;
	double K2 = - ((vi+dt/2*K1)-(vant+dt/2*aaaaa))/fabs((xi+dt/2*L1)-(xant+dt/2*aaaaa));
	double L2 = vi+dt/2*K1;
	double K3 = - ((vi+dt/2*K2)-(vant+dt/2*aaaaa))/fabs((xi+dt/2*L2)-(xant+dt/2*aaaaa));
	double L3 = vi+dt/2*K2;
	double K4 = - ((vi+dt*K3)-(vant+dt*aaaaa))/fabs((xi+dt*L3)-(xant+dt*aaaaa));
	double L4 = vi+dt*K3;

	double vtNext = vi+dt/6*(K1+2*K2+2*K3+K4);
	double xtNext = xi+dt/6*(L1+2*L2+2*L3+L4);

	vect = [xtNext, vtNext]
}
*/
int main()
{
	int tc=1;
	int m=1500;
	double L=4;
	int c=20000;
	double d=26;
	double v_eq=120;


	double dt=0.001;
	int tfrenat = 20;
	int tf=tc+tfrenat;

	int Nt=(tf/dt);

	double veqAdim = adimv(v_eq);
	double carsx[5][Nt+1];
	double carsv[5][Nt+1];	
	double carserr[5][(int)(tfrenat/dt)-1];

	double LAdim = adimx(L);
	double dAdim = adimx(d);

	double K0 = 5.0*(1.0/2.0*m*v_eq*v_eq);
	printf("%f",K0);


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

		double newx = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if(newx < carsx[0][i]-LAdim){
			carsx[1][i+1] = newx;	
		}else{
			carsx[1][i+1] = carsx[0][i]-LAdim;
		}
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);

		carserr[0][(int)(i-tc/dt)]=0;

		double wi = m*SIa(K11)*SIx(dt/6*(L11+2*L21+2*L31+L41));
		double Ki = 0.5*m*(SIv(carsv[0][i])*SIv(carsv[0][i])+SIv(carsv[1][i])*SIv(carsv[1][i])+SIv(carsv[2][i])*SIv(carsv[2][i])+SIv(carsv[3][i])*SIv(carsv[3][i])+SIv(carsv[4][i])*SIv(carsv[4][i]));
		
		carserr[1][(int)(i-tc/dt)]=K0-(Ki-wi);

		double K12 = -(carsv[2][i]-carsv[1][i])/fabs(carsx[2][i]-carsx[1][i]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][i]+dt/2*K12)-(carsv[1][i]+dt/2*K11))/fabs((carsx[2][i]+dt/2*L12)-(carsx[1][i]+dt/2*L11));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][i]+dt/2*K21))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][i]+dt/2*L21));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][i]+dt*K31))/fabs((carsx[2][i]+dt*L32)-(carsx[1][i]+dt*L31));
		double L42 = carsv[2][i]+dt*K32;

		newx = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if(newx < carsx[1][i]-LAdim){
			carsx[2][i+1] = newx;
		}else{
			carsx[2][i+1]= carsx[1][i]-LAdim;
		}
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);
		
		wi = m*SIa(K12)*SIx(dt/6*(L12+2*L22+2*L32+L42));
		carserr[2][(int)(i-tc/dt)]=K0-(Ki-wi);

		double K13 = -(carsv[3][i]-carsv[2][i])/fabs(carsx[3][i]-carsx[2][i]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][i]+dt/2*K12))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][i]+dt/2*L12));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][i]+dt/2*K22))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][i]+dt/2*L22));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][i]+dt*K32))/fabs((carsx[3][i]+dt*L33)-(carsx[2][i]+dt*L32));
		double L43 = carsv[3][i]+dt*K33;

		newx = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if(newx < carsx[2][i]-LAdim){
			carsx[3][i+1] = newx;
		}else{
			carsx[3][i+1] = carsx[2][i]-LAdim;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		wi = m*SIa(K13)*SIx(dt/6*(L13+2*L23+2*L33+L43));
		carserr[3][(int)(i-tc/dt)]=K0-(Ki-wi);

		double K14 = -(carsv[4][i]-carsv[3][i])/fabs(carsx[4][i]-carsx[3][i]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][i]+dt/2*K13))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][i]+dt/2*L13));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][i]+dt/2*K23))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][i]+dt/2*L23));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][i]+dt*K33))/fabs((carsx[4][i]+dt*L34)-(carsx[3][i]+dt*L33));
		double L44 = carsv[4][i]+dt*K34;

		newx = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if (newx < carsx[3][i]-LAdim)
		{
			carsx[4][i+1] = newx;	
		}else{
			carsx[4][i+1] = carsx[3][i]-LAdim;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		wi = m*SIa(K14)*SIx(dt/6*(L14+2*L24+2*L34+L44));
		carserr[4][(int)(i-tc/dt)]=K0-(Ki-wi);
		
	}

	FILE *fptr;
	fptr = fopen("C1_tau0_x.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C1_tau0_v.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIv(carsv[j][i]));
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);

	fptr = fopen("C1_tau0_err.txt", "wb");
	for (int i = 0; i < (int)(tfrenat/dt); i++)
	{
		for (int j = 0; j < 5; j++)
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
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		double K10 = -veqAdim*exp(1+tc-i*dt)*(1+tc-i*dt);//K pel cotxe 0 a ordre 1
		double K20 = -veqAdim*exp(1+tc-(i*dt+dt/2))*(1+tc-(i*dt+dt/2));
		double K30 = -veqAdim*exp(1+tc-(i*dt+dt/2))*(1+tc-(i*dt+dt/2));
		
		double L10 = veqAdim*(1-(i*dt-tc)*exp(1-(i*dt-tc)));//L pel cotxe 0 a ordre 1
		double L20 = veqAdim*(1-(i*dt+dt/2-tc)*exp(1-(i*dt+dt/2-tc)));
		double L30 = veqAdim*(1-(i*dt+dt/2-tc)*exp(1-(i*dt+dt/2-tc)));
		double L40 = veqAdim*(1-(i*dt+dt-tc)*exp(1-(i*dt+dt-tc)));

		double wi = K10*dt/6*(L10+2*L20+2*L30+L40);
		double Ki = 1/2*m*(SIv(carsv[0][i])*SIv(carsv[0][i])+SIv(carsv[1][i])*SIv(carsv[1][i])+SIv(carsv[2][i])*SIv(carsv[2][i])+SIv(carsv[3][i])*SIv(carsv[3][i])+SIv(carsv[4][i])*SIv(carsv[4][i]));

		carserr[1][(int)(i-tc/dt)]=K0-(Ki-wi);

		double K11 = -(carsv[1][i]-carsv[0][i])/fabs(carsx[1][i]-carsx[0][i]);
		double L11 = carsv[1][i];
		double K21 = - ((carsv[1][i]+dt/2*K11)-(carsv[0][i]+dt/2*K10))/fabs((carsx[1][i]+dt/2*L11)-(carsx[0][i]+dt/2*L10));
		double L21 = carsv[1][i]+dt/2*K11;
		double K31 = - ((carsv[1][i]+dt/2*K21)-(carsv[0][i]+dt/2*K20))/fabs((carsx[1][i]+dt/2*L21)-(carsx[0][i]+dt/2*L20));
		double L31 = carsv[1][i]+dt/2*K21;
		double K41 = - ((carsv[1][i]+dt*K31)-(carsv[0][i]+dt*K30))/fabs((carsx[1][i]+dt*L31)-(carsx[0][i]+dt*L30));
		double L41 = carsv[1][i]+dt*K31;

		double newx = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if (newx < carsx[0][i]-LAdim)
		{
			carsx[1][i+1] = newx;
		}else{
			carsx[1][i+1] = carsx[0][i]-LAdim;
		}
		
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);
		 


		double K12 = -(carsv[2][i]-carsv[1][i])/fabs(carsx[2][i]-carsx[1][i]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][i]+dt/2*K12)-(carsv[1][i]+dt/2*K11))/fabs((carsx[2][i]+dt/2*L12)-(carsx[1][i]+dt/2*L11));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][i]+dt/2*K21))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][i]+dt/2*L21));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][i]+dt*K31))/fabs((carsx[2][i]+dt*L32)-(carsx[1][i]+dt*L31));
		double L42 = carsv[2][i]+dt*K32;

		newx = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if (newx < carsx[1][i]-LAdim)
		{
			carsx[2][i+1] = newx;
		}else{
			carsx[2][i+1] = carsx[1][i]-LAdim;
		}
		
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);


		double K13 = -(carsv[3][i]-carsv[2][i])/fabs(carsx[3][i]-carsx[2][i]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][i]+dt/2*K12))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][i]+dt/2*L12));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][i]+dt/2*K22))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][i]+dt/2*L22));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][i]+dt*K32))/fabs((carsx[3][i]+dt*L33)-(carsx[2][i]+dt*L32));
		double L43 = carsv[3][i]+dt*K33;

		newx=carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if(newx<carsx[2][i]-LAdim){
			carsx[3][i+1] = newx;
		}else{
			carsx[3][i+1] = carsx[2][i]-LAdim;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		

		double K14 = -(carsv[4][i]-carsv[3][i])/fabs(carsx[4][i]-carsx[3][i]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][i]+dt/2*K13))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][i]+dt/2*L13));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][i]+dt/2*K23))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][i]+dt/2*L23));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][i]+dt*K33))/fabs((carsx[4][i]+dt*L34)-(carsx[3][i]+dt*L33));
		double L44 = carsv[4][i]+dt*K34;


		newx = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if(newx<carsx[3][i]-LAdim){
			carsx[4][i+1] = newx;
		}else{
			carsx[4][i+1] = carsx[3][i]-LAdim;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		
	}

	fptr = fopen("C2_tau0_x.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C2_tau0_v.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIv(carsv[j][i]));
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
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		double K10 = 0.8*veqAdim*omega*sin(2*omega*(tc-dt*i));//K pel cotxe 0 a ordre 1
		double K20 = 0.8*veqAdim*omega*sin(2*omega*(tc-(dt*i+dt/2)));
		double K30 = K20;
		
		double L10 = veqAdim*(1-0.8*sin(omega*(i*dt-tc)*sin(omega*(i*dt-tc))));//L pel cotxe 0 a ordre 1
		double L20 = veqAdim*(1-0.8*sin(omega*(i*dt+dt/2-tc)*sin(omega*(i*dt+dt/2-tc))));
		double L30 = L20;

		double K11 = -(carsv[1][i]-carsv[0][i])/fabs(carsx[1][i]-carsx[0][i]);
		double L11 = carsv[1][i];
		double K21 = - ((carsv[1][i]+dt/2*K11)-(carsv[0][i]+dt/2*K10))/fabs((carsx[1][i]+dt/2*L11)-(carsx[0][i]+dt/2*L10));
		double L21 = carsv[1][i]+dt/2*K11;
		double K31 = - ((carsv[1][i]+dt/2*K21)-(carsv[0][i]+dt/2*K20))/fabs((carsx[1][i]+dt/2*L21)-(carsx[0][i]+dt/2*L20));
		double L31 = carsv[1][i]+dt/2*K21;
		double K41 = - ((carsv[1][i]+dt*K31)-(carsv[0][i]+dt*K30))/fabs((carsx[1][i]+dt*L31)-(carsx[0][i]+dt*L30));
		double L41 = carsv[1][i]+dt*K31;

		double newx = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if(newx < carsx[0][i]-LAdim){
			carsx[1][i+1] = newx;
		}else{
			carsx[1][i+1] = carsx[0][i]-LAdim;
		}

		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);


		double K12 = -(carsv[2][i]-carsv[1][i])/fabs(carsx[2][i]-carsx[1][i]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][i]+dt/2*K12)-(carsv[1][i]+dt/2*K11))/fabs((carsx[2][i]+dt/2*L12)-(carsx[1][i]+dt/2*L11));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][i]+dt/2*K21))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][i]+dt/2*L21));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][i]+dt*K31))/fabs((carsx[2][i]+dt*L32)-(carsx[1][i]+dt*L31));
		double L42 = carsv[2][i]+dt*K32;

		newx = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if(newx < carsx[1][i]-LAdim){
			carsx[2][i+1] = newx;
		}else{
			carsx[2][i+1] = carsx[1][i]-LAdim;
		}
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);


		double K13 = -(carsv[3][i]-carsv[2][i])/fabs(carsx[3][i]-carsx[2][i]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][i]+dt/2*K12))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][i]+dt/2*L12));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][i]+dt/2*K22))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][i]+dt/2*L22));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][i]+dt*K32))/fabs((carsx[3][i]+dt*L33)-(carsx[2][i]+dt*L32));
		double L43 = carsv[3][i]+dt*K33;

		newx = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if(newx < carsx[2][i]-LAdim){
			carsx[3][i+1] = newx;
		}else{
			carsx[3][i+1] = carsx[2][i]-LAdim;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);

		double K14 = -(carsv[4][i]-carsv[3][i])/fabs(carsx[4][i]-carsx[3][i]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][i]+dt/2*K13))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][i]+dt/2*L13));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][i]+dt/2*K23))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][i]+dt/2*L23));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][i]+dt*K33))/fabs((carsx[4][i]+dt*L34)-(carsx[3][i]+dt*L33));
		double L44 = carsv[4][i]+dt*K34;

		newx = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if(newx<carsx[3][i]-LAdim){
			carsx[4][i+1] = newx;
		}else{
			carsx[4][i+1] = carsx[3][i]-LAdim;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		
	}

	fptr = fopen("C3_tau0_x.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C3_tau0_v.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIv(carsv[j][i]));
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);

	//tau = 0.5	
	double tau = 0.6;

	//---cas 1----

	//euler cotxe 0
	
	for (int i = (tc/dt); i < tf/dt; i++)
	{
		carsx[0][i+1]=carsx[0][i]+0.2*veqAdim*dt;

		carsv[0][i+1]=0.2*veqAdim;
	}
	
	//RK4
	
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

		double newx = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if(newx < carsx[0][i]-LAdim){
			carsx[1][i+1] = newx;
		}else{
			carsx[1][i+1] = carsx[0][i]-LAdim;
		}
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);
		


		double K12 = -(carsv[2][i]-carsv[1][(int)((i*dt-tau)/dt)])/fabs(carsx[2][i]-carsx[1][(int)((i*dt-tau)/dt)]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][(int)((i*dt-tau)/dt)])-(carsv[1][i]+dt/2*K11))/fabs((carsx[2][(int)((i*dt-tau)/dt)])-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt*L32)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L42 = carsv[2][i]+dt*K32;

		newx = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if(newx < carsx[1][i]-LAdim){
			carsx[2][i+1] = newx;
		}else{
			carsx[2][i+1] = carsx[1][i]-LAdim;
		}
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);


		double K13 = -(carsv[3][i]-carsv[2][(int)((i*dt-tau)/dt)])/fabs(carsx[3][i]-carsx[2][(int)((i*dt-tau)/dt)]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt*L33)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L43 = carsv[3][i]+dt*K33;

		newx = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if(newx < carsx[2][i]-LAdim){
			carsx[3][i+1] = newx;
		}else{
			carsx[3][i+1] = carsx[2][i]-LAdim;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		

		double K14 = -(carsv[4][i]-carsv[3][(int)((i*dt-tau)/dt)])/fabs(carsx[4][i]-carsx[3][(int)((i*dt-tau)/dt)]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt*L34)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L44 = carsv[4][i]+dt*K34;

		newx = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if(newx < carsx[3][i]-LAdim){
			carsx[4][i+1] = newx;
		}else{
			carsx[4][i+1] = carsx[3][i]-LAdim;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		
	}

	fptr = fopen("C1_tau_x.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C1_tau_v.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIv(carsv[j][i]));
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

		double newx = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if(newx<carsx[0][i]-LAdim){
			carsx[1][i+1] = newx;	
		}else{
			carsx[1][i+1] = carsx[0][i]-LAdim;
		}
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);
		


		double K12 = -(carsv[2][i]-carsv[1][(int)((i*dt-tau)/dt)])/fabs(carsx[2][i]-carsx[1][(int)((i*dt-tau)/dt)]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][i]+dt/2*K12)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt/2*L12)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt*L32)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L42 = carsv[2][i]+dt*K32;

		newx = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if(newx < carsx[1][i]-LAdim){
			carsx[2][i+1] = newx;	
		}else{
			carsx[2][i+1] = carsx[1][i]-LAdim;
		}
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);
		


		double K13 = -(carsv[3][i]-carsv[2][(int)((i*dt-tau)/dt)])/fabs(carsx[3][i]-carsx[2][(int)((i*dt-tau)/dt)]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt*L33)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L43 = carsv[3][i]+dt*K33;

		newx = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if(newx < carsx[2][i]-LAdim){
			carsx[3][i+1] = newx;
		}else{
			carsx[3][i+1] = carsx[2][i]-LAdim;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		

		double K14 = -(carsv[4][i]-carsv[3][(int)((i*dt-tau)/dt)])/fabs(carsx[4][i]-carsx[3][(int)((i*dt-tau)/dt)]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt*L34)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L44 = carsv[4][i]+dt*K34;

		newx = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if(newx < carsx[3][i]-LAdim){
			carsx[4][i+1] = newx;
		}else{
			carsx[4][i+1] = carsx[3][i]-LAdim;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		
	}

	fptr = fopen("C2_tau_x.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C2_tau_v.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIv(carsv[j][i]));
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

		double newx = carsx[1][i]+dt/6*(L11+2*L21+2*L31+L41);
		if(newx < carsx[0][i]-LAdim){
			carsx[1][i+1] = newx;
		}else{
			carsx[1][i+1] = carsx[0][i]-LAdim;
		}
		carsv[1][i+1] = carsv[1][i]+dt/6*(K11+2*K21+2*K31+K41);
		


		double K12 = -(carsv[2][i]-carsv[1][(int)((i*dt-tau)/dt)])/fabs(carsx[2][i]-carsx[1][(int)((i*dt-tau)/dt)]);
		double L12 = carsv[2][i];
		double K22 = - ((carsv[2][i]+dt/2*K12)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt/2*L12)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L22 = carsv[2][i]+dt/2*K12;
		double K32 = - ((carsv[2][i]+dt/2*K22)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt/2*L22)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L32 = carsv[2][i]+dt/2*K22;
		double K42 = - ((carsv[2][i]+dt*K32)-(carsv[1][(int)((i*dt-tau)/dt)]))/fabs((carsx[2][i]+dt*L32)-(carsx[1][(int)((i*dt-tau)/dt)]));
		double L42 = carsv[2][i]+dt*K32;

		newx = carsx[2][i]+dt/6*(L12+2*L22+2*L32+L42);
		if(newx < carsx[1][i]-LAdim){
			carsx[2][i+1] = newx;
		}else{
			carsx[2][i+1] = carsx[1][i]-LAdim;
		}
		carsv[2][i+1] = carsv[2][i]+dt/6*(K12+2*K22+2*K32+K42);
		


		double K13 = -(carsv[3][i]-carsv[2][(int)((i*dt-tau)/dt)])/fabs(carsx[3][i]-carsx[2][(int)((i*dt-tau)/dt)]);
		double L13 = carsv[3][i];
		double K23 = - ((carsv[3][i]+dt/2*K13)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L13)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L23 = carsv[3][i]+dt/2*K13;
		double K33 = - ((carsv[3][i]+dt/2*K23)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt/2*L23)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L33 = carsv[3][i]+dt/2*K23;
		double K43 = - ((carsv[3][i]+dt*K33)-(carsv[2][(int)((i*dt-tau)/dt)]))/fabs((carsx[3][i]+dt*L33)-(carsx[2][(int)((i*dt-tau)/dt)]));
		double L43 = carsv[3][i]+dt*K33;

		newx = carsx[3][i]+dt/6*(L13+2*L23+2*L33+L43);
		if(newx < carsx[2][i]-LAdim){
			carsx[3][i+1] = newx;
		}else{
			carsx[3][i+1] = carsx[2][i]-LAdim;
		}
		carsv[3][i+1] = carsv[3][i]+dt/6*(K13+2*K23+2*K33+K43);
		

		double K14 = -(carsv[4][i]-carsv[3][(int)((i*dt-tau)/dt)])/fabs(carsx[4][i]-carsx[3][(int)((i*dt-tau)/dt)]);
		double L14 = carsv[4][i];
		double K24 = - ((carsv[4][i]+dt/2*K14)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L14)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L24 = carsv[4][i]+dt/2*K14;
		double K34 = - ((carsv[4][i]+dt/2*K24)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt/2*L24)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L34 = carsv[4][i]+dt/2*K24;
		double K44 = - ((carsv[4][i]+dt*K34)-(carsv[3][(int)((i*dt-tau)/dt)]))/fabs((carsx[4][i]+dt*L34)-(carsx[3][(int)((i*dt-tau)/dt)]));
		double L44 = carsv[4][i]+dt*K34;

		newx = carsx[4][i]+dt/6*(L14+2*L24+2*L34+L44);
		if(newx < carsx[3][i]-LAdim){
			carsx[4][i+1] = newx;
		}else{
			carsx[4][i+1] = carsx[3][i]-LAdim;
		}
		carsv[4][i+1] = carsv[4][i]+dt/6*(K14+2*K24+2*K34+K44);
		
	}

	fptr = fopen("C3_tau_x.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIx(carsx[j][i]));
		}
		fprintf(fptr, "\n");
	}
	
	fclose(fptr);

	fptr = fopen("C3_tau_v.txt", "wb");
	for (int i = 0; i < Nt+1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			fprintf(fptr,"%f, ",SIv(carsv[j][i]));
		}
		fprintf(fptr, "\n");
	}

	fclose(fptr);


}