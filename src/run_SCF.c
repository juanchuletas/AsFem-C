#include<stdio.h>
#include<stdlib.h>



extern void PoissonSolver(double *hp,double *lij,double *sij,double *rho,double *bc, double boundVal,int fembasis);

extern void print_matrix(char *name,int m,int n,double *matA);

void startSCF(double *sij,double *lij,double *rho_r,int Ne, int order,int phase,double *bound_vec,double bcValue){

	printf("***** SCF RUNNING *****\n");
	int nodes = Ne*order + 1;
	int fembasis = nodes-2;
	double hpot[fembasis];
	//for(int i =0; i<fembasis*fembasis; i++)
	//printf("lij[%d] = %lf\n",i,lij[i]);
	//print_matrix("*** L Matrix ***",fembasis,fembasis,lij);
	//print_matrix("*** S Matrix ***",fembasis,fembasis,sij);

	PoissonSolver(hpot,lij,sij,rho_r,bound_vec,bcValue,fembasis);




}
