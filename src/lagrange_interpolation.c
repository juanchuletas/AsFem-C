#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_structures.h"
#include <math.h>
extern int AtomicNumber(char target[3]);
extern void PrintWaveFunction(int Ne,int order,double *femnodes,double *wfn,double *vecE,double *hpotential);
extern void GetHartreePotential(double *mat_vh,double *sij,double *kij,double *wfn,int Ne,int order,int phase,double *b,double hp);
extern void EnergyResults(double *eigenval,int N,int r_nodes);
extern void PrintPotentialData(int order,int Ne,struct Element *e,struct Vertex *n,struct Vertex *pot);
extern void GenerateNodes(double *x,int nodes,double r0,double rN,char mesh[180],char atom[3]);
extern void SetElements(struct Element *e,struct Vertex *n,int Ne, double *x,int order,double rN);
extern void EvaluatePotential(char kindpot[180],char atom[3],int angular,int Ne,int order,struct Element *e,struct Vertex *n,struct Vertex *pot);
extern void print_matrix( char* desc, int m, int n, double* a);
extern void GetLinkMatrix(int *link_mat,int Ne,int order);
extern void AssambleGlobalMatrices(int Ne,int order,int *link_mat,double *sij,double *kij,double *vij,double *uij,double *lij,struct Element *e,struct Vertex *pot,double *b_vec,double *avec);
extern void startProgram(double *hij,double *sij,double *kij,double *vij,double *lij,double *wfn,double *vec_ei,int *link_mat,struct Element *e,double *femnodes,int Ne,int order,double *b_vec,double *a_vec,int atomicN);
void LagrangeInterpolation(int Ne,double r0,double rN,char kindpot[180],char mesh[180],char atom[3],int order,int angular)
{
	int nodes = Ne*order + 1;
	int fembasis = nodes-2;
	int simplyNodes = Ne+1;

	int *link_mat;
	struct Element *e;
	double *x;

	//---------- MEMORY ALLOCATION FOR THE NECCESARY MATRICES ----------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//******************Structures for nodes********************************
        e = (struct Element *)malloc(sizeof(struct Element)*Ne);
	//***************** Nodes in x ****************************************
	x = (double *)malloc(sizeof(double)*(simplyNodes));
	x[simplyNodes-1]=rN;
	//***************** Link Matrix ****************************************
	link_mat = (int *)malloc(sizeof(int)*(Ne*(order+1)));
	
	// **** Memory allocation: Reduced matrices for eigenvalue problem
	double *sij,*kij,*vij,*hij,*vh_mat,*uij,*lij,*ci,*ei,*vhij;
	int K = nodes-2;
	sij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	kij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	vij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	hij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	vhij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	//***************** Eigenvectos and eigenvalues ************************
	ci = (double *)malloc(sizeof(double)*(K*K));
        ei = (double *)malloc(sizeof(double)*fembasis);
	// **** Memory allocation: Reduced matrices for the Poisson problem
	lij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	uij = (double *)malloc(sizeof(double)*(fembasis*fembasis));
	double b_vec[nodes-2];
	double a_vec[nodes-2];
	double femnodes[fembasis];
	
	int atomicN = AtomicNumber(atom);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------



	//------------ GENERATION OF NODES AND ELEMENTS------------------------
	GenerateNodes(x,Ne,r0,rN,mesh,atom);
	SetElements(e,e->n,Ne,x,order,rN);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//------------------- EVALUATE THE POTENTIAL ---------------------------
	//
	EvaluatePotential(kindpot,atom,angular,Ne,order,e,e->n,e->pot);
	//
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------

	//------------------- GET THE ELEMENTAL MATRICES AND -------------------
	//--------------------------- THE LINK MATRIX ---------------------------
	GetLinkMatrix(link_mat,Ne,order);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//------------------- GET THE GLOBAL MATRICES ---------------------------
	AssambleGlobalMatrices(Ne,order,link_mat,sij,kij,vij,uij,lij,e,e->pot,b_vec,a_vec);
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//print_matrix("*** Overlap Matrix ***",K,K,sij);
	//print_matrix("*** Kinect Matrix ***",K,K,kij);
        /*print_matrix("*** Potential Matrix ***",K,K,vij);
        print_matrix("*** U Poisson Matrix ***",fembasis_poisson,fembasis_poisson,uij);
        print_matrix("*** L Poisson Matrix ***",fembasis_poisson,fembasis_poisson,lij);*/
	//PrintPotentialData(order,Ne,e,e->n,e->pot);
        startProgram(hij,sij,kij,vij,lij,ci,ei,link_mat,e,femnodes,Ne,order,b_vec,a_vec,atomicN);
	//EnergyResults(ei,10,fembasis);
	//print_matrix("Hartree-Potential",K,K,mat_vh);
	//PrintWaveFunction(Ne,order,femnodes,ci,ei,hpotential);

	free(link_mat);
	free(hij);
	free(vhij);
	free(sij);
	free(kij);
	free(uij);
	free(lij);
	free(vij);
	free(e);
	free(x);
	free(ei);
	free(ci);

}
