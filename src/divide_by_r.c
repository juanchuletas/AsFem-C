#include<stdio.h>
#include<stdlib.h>


void divideByR(double *newrho, double *rho,double *femnodes,int nfembasis){


	double factor;
	int l=0;
	for(int i=0; i<nfembasis; i++){
		
		newrho[i] = rho[i]/femnodes[i];
	
	}



}
