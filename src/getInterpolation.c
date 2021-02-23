#include<stdio.h>
#include<stdlib.h>
#include<math.h>
extern double Pn_Interpolation(int order,double *x,double *y,double target);
extern double nPolExtrapolation(double *x,double *y, int N,double target);
double doInterpolation(double *xin,double *yin){

	double value;
	double x[7],y[7];	
	 x[0]=xin[0];       y[0]=yin[1];
	 x[1]=xin[1];       y[1]=yin[2];
	 x[2]=xin[2];       y[2]=yin[3];
	 x[3]=xin[3];       y[3]=yin[4];
	 x[4]=xin[4];       y[4]=yin[5];
	 x[5]=xin[5];       y[5]=yin[6];
	 x[6]=xin[6];       y[6]=yin[7];
	 value = nPolExtrapolation(x,y,7,0.000000);


	return value;
}

void getHartreePotentialRealValues(double *vhpotential,double *femnodes,double *uhpotential,int nodes){

	double factor,screen;
	double lambda = 1.0;
	int fembasis  = nodes-2;
	for(int i=0; i<fembasis; i++){
		factor = femnodes[i];
		//screen = exp(-factor/lambda);
		vhpotential[i+1] = uhpotential[i]/femnodes[i];
		//vhpotential[i+1] = vhpotential[i+1]*screen;
	
	}

	vhpotential[0] = doInterpolation(femnodes,vhpotential);
	vhpotential[nodes-1] = 0.0;

	/*for(int i=0; i<fembasis; i++){
		screen = exp(-femnodes[i]/lambda);
		vhpotential[i+1] = vhpotential[i+1]*screen;
	}*/


}
