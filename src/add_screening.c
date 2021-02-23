#include<stdio.h>
#include<math.h>




void addScreening(double *potnew,double *potold, double *rnodes,int fembasis){

	double lambda = 0.80;
	double factor;
	for(int i=0; i<fembasis; i++){
		factor = exp(-rnodes[i]/lambda);	
		potnew[i] = potold[i]*factor;
	}
}
