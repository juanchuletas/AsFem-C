#include<stdio.h>
#include<stdlib.h>
#include "data_structures.h"
#include "Elemental_Matrices.h"
extern double *ElementalOverlapMatrix(int order);
extern void MatrixProduct(double *A,double *B,double *C,int N, int M, int P);
extern void ScalarXMatrix(double coeff,double *mat_A,double *mat_Res,int N,int M);
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );


void PoissonSolver(double *hp,double *lij,double *sij,double *rho,double *bc, double boundVal,int fembasis)
{
	printf("POISSON SOLVER\n");
	int NRHS = 1;
	int N=fembasis;
	int LDA=N;
	int LDB=N;
	int ipiv[N];
	int Norb=1;
	int info;
	const double PI = 3.14159265358979323846;
	const double coeff = 4.0*PI;

	
        double aux_vec[fembasis];
        double f_bc[fembasis];
        double f_vec[fembasis];

	//THIS SHIT ALMOST 
	double aux_mat[fembasis*fembasis];

	ScalarXMatrix(1.0,lij,aux_mat,fembasis,fembasis);

	MatrixProduct(sij,rho,f_vec,fembasis,fembasis,NRHS);

	for(int i=0; i<fembasis; i++){
		
		f_bc[i] = f_vec[i] - boundVal*bc[i];
	
	}

	dgesv_(&N,&NRHS,aux_mat,&LDA,ipiv,f_bc,&LDB,&info);

	if( info > 0 )
        {
        printf( "The diagonal element of the triangular factor of A,\n" );
        printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
        printf( "the solution could not be computed.\n" );
        exit(1);
        }

	for(int i=0; i<fembasis; i++)
	{
		hp[i] = f_bc[i];
		printf("H_Pot[%d] = %lf\n",i,hp[i]);

	}
	//printf("HartreePot[%d] = %lf\n",nodes-1,hartree_vec[nodes-1]);

	//free(right_vec);
	//free(aux_mat);


}
