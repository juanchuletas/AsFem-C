#include<stdio.h>
#include<stdlib.h>
#include "data_structures.h"
#include "wfn_operations.h"
#include <math.h>
extern void addScreening(double *rhonew,double *rhoold, double *rnodes,int fembasis);
extern void FillZeroMat(double *mat,int M,int N);
extern void GetHmatrix(double *h,double *v,double *k,int M,int N);
extern void getFemNodes(double *femnodes, struct Element *e, struct Vertex *n, int Ne, int order);
extern void diag (int n, double *h, double *s, double *e, double *v);
extern void ComputeDensity(double *rho,double *wfn,int Ne, int order, int atomicN);
extern void NormWfn(double *wfn_vec,int *link_mat,struct Element *e,int order,int Ne);
extern void GetHartreePotential(double *uh_pot,double *sij,double *kij,double *rho_r,int Ne,int order,int phase,double *b,double hp);
extern void divideByR(double *vhpotential,double *uhpotential,double *femnodes,int nfembasis);
extern void startSCF(double *sij,double *lij,double *rho_r,int Ne, int order,int phase,double *bound_vec,double bcValue);
//extern void PoissonSolver(double *hp,double *lij,double *sij,double *rho,double *bc, double boundVal,int fembasis);
extern void getHartreePotentialRealValues(double *vhpotential, double*femnodes,double *uhpotential,int nodes);
extern void IntegrateHartreePot(int Ne,int order,int *link_mat,double *vhij,struct Element *e,double *v);
extern void getDensityMatrix(double *densMat,double *wfn, int fembasis,int atomicN);
extern void SumMatrices(double *A, double *B, double *C, int size);
extern double getHFenergy(double *hcore, double *fmat, double *densMat,int totSize);
void startProgram(double *hij,double *sij,double *kij,double *vij,double *lij,double *wfn,double *vec_ei,int *link_mat,struct Element *e,double *femRnodes,int Ne,int order,double *boundC,double *a_vec,int atomicN){

	printf(" *********** RUNNING PROGRAM **************** \n");
	printf("Atomic Number: %d\n",atomicN);
	int nodes = Ne*order + 1;
	int fembasis = nodes-2;
	int totSize = fembasis*fembasis;
	double *f_mat = (double *)malloc(sizeof(double)*(totSize));
	double *densMat = (double *)malloc(sizeof(double)*(totSize));
	double *vhij_mat = (double *)malloc(sizeof(double)*(totSize));

	double vhpotential[nodes];

	int phase = 1;


	//DENS ARRAY
	
	double rho[fembasis];
	double rho_r[fembasis];
	double uh_pot[fembasis];


	getFemNodes(femRnodes,e,e->n,Ne,order);

	FillZeroMat(hij,fembasis,fembasis);
	FillZeroMat(f_mat,fembasis,fembasis);
	FillZeroMat(densMat,fembasis,fembasis);

	GetHmatrix(hij,kij,vij,fembasis,fembasis);

	

	//FIRST DIAGONALIZATION TO OBTAIN THE WFN
	//
	//nite differences
	printf("*** FIRST CALCULATION OF EIGENVALUES ***\n");
	diag(fembasis,hij,sij,vec_ei,wfn);
	double first_orb_energy = 0.5*vec_ei[0];
	printf("First step orbital energy: %lf    hartrees\n",first_orb_energy);
	
	//GET THE WFN AND COMPUTE THE DENSITY
	NormWfn(wfn,link_mat,e,order,Ne);
	GetWfnPhase(fembasis,0,&phase,wfn);
	getDensityMatrix(densMat,wfn,fembasis,atomicN);
	double energy0 = getHFenergy(hij,f_mat,densMat,totSize);

	printf("FIRST HATREE-FOCK ENERGY: %lf\n",energy0);

	ComputeDensity(rho,wfn,Ne,order,atomicN);
	divideByR(rho_r,rho,femRnodes,fembasis);
	//addScreening(rhoScreen,rho_r,femRnodes,fembasis);

	//FIRST APPROXIMATION OF THE HARTREE POTENTIAL
	double lambda = 1.0;
	double screen = exp(-e[Ne-1].n[order].x*lambda);
	double hp = 2.0*screen;

	printf("Boundary conditon Hartree-Pot: %lf\n",hp);
	
	
	GetHartreePotential(uh_pot,sij,lij,rho_r,Ne,order,phase,boundC,hp);
	double carga = 2.0;
	getHartreePotentialRealValues(vhpotential,femRnodes,uh_pot,nodes);
        vhpotential[nodes-1] = ((carga*screen)/e[Ne-1].n[order].x);
	
	//for(int i = 0; i<nodes; i++){printf("VH(r) = %lf\n",vhpotential[i]);}
	//printf("V(0) = %lf\n",vhpotential[0]);
        //printf("V(Rc) = %lf\n",vhpotential[nodes-1]);

	
	printf("-----------------------------------------------------------\n");
        printf("-------------------SCF MODULE------------------------------\n");


	int iter=0;
	double deltaEnergy = 10000000.0;
	double new_HF = energy0;
	double old_HF,orb_energy;
	while(deltaEnergy>0.0000001){
		old_HF = new_HF;
		IntegrateHartreePot(Ne,order,link_mat,vhij_mat,e,vhpotential);
		SumMatrices(hij,vhij_mat,f_mat,totSize);
		diag(fembasis,f_mat,sij,vec_ei,wfn);
		orb_energy = 0.5*vec_ei[0];
		//printf("new orbital energy: %lf    hartrees at step: %d\n",orb_energy,iter+1);
		NormWfn(wfn,link_mat,e,order,Ne);
		GetWfnPhase(fembasis,0,&phase,wfn);
		getDensityMatrix(densMat,wfn,fembasis,atomicN);
		new_HF = getHFenergy(hij,f_mat,densMat,totSize);
		printf("Hartree-Fock energy=%lf  orb=%lf   iteration %d\n",new_HF,orb_energy,iter);
		ComputeDensity(rho,wfn,Ne,order,atomicN);
		divideByR(rho_r,rho,femRnodes,fembasis);
		//addScreening(rhoScreen,rho_r,femRnodes,fembasis);
		GetHartreePotential(uh_pot,sij,lij,rho_r,Ne,order,phase,boundC,hp);
		getHartreePotentialRealValues(vhpotential,femRnodes,uh_pot,nodes);
        	//0vhpotential[nodes-1] = hp/e[Ne-1].n[order].x;
        	vhpotential[nodes-1] = ((carga*screen)/e[Ne-1].n[order].x);
		
		deltaEnergy = fabs(new_HF-old_HF);
		iter++;
	
	}
	printf("-----------------------------------------------------------\n");
	printf("-----------------------------------------------------------\n");
	printf("--------------Hartree-Fock Final Results-------------------\n");
	printf("-----------------------------------------------------------\n");
	printf("Hatree-Fock Final Energy=%20.10lf\n",new_HF);
	printf("-----------------------------------------------------------\n");
	printf("Occupied Orbital Energy = %lf\n",orb_energy);
	printf("-----------------------------------------------------------\n");
	printf("-----------------------------------------------------------\n");



	
	free(f_mat);
	free(vhij_mat);
	free(densMat);

}


