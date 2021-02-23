#include <stdio.h>
#include <stdlib.h>
extern void MatrixProduct(double *A,double *B,double *C,int N, int M, int P);
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
extern void MatXvec(double *res,double *matA,double *vec,int SIZE);
extern void ScalarXMatrix(double coeff,double *mat_A,double *mat_Res,int N,int M);
extern void SumMatrices(double *A, double *B, double *C, int size);
void GetHartreePotential(double *hartree_vec,double *sij,double *lij,double *rho,int Ne,int order,int phase,double *b,double hp)
{
	int nodes = Ne*order + 1;
	int r_nodes = nodes-2;
	int totSize = r_nodes*r_nodes;
	int N=r_nodes;
	int NRHS=1;
	int LDA=N;
	int LDB=N;
	int ipiv[N];
	int info;
	int Norb=1;
	double sum=0.0;
	double coeffwfn,coeffrho;
	const double PI = 3.14159265358979323846;

	//LOCAL ARRAYS:
	double right_vec[r_nodes];
	double aux_vec[r_nodes];
	
	double *aux_mat = (double *)malloc(sizeof(double)*(totSize));



	ScalarXMatrix(1.0,lij,aux_mat,r_nodes,r_nodes);
		

	MatrixProduct(sij,rho,right_vec,r_nodes,r_nodes,NRHS);
	
	//double qtot = 0.0;
	for(int i=0; i<r_nodes; i++){
		aux_vec[i] = right_vec[i]  - b[i]*hp;
		//qtot  = qtot + right_vec[i];
		//printf("f = %lf\n",rho[i]);
	}
	//printf("THE SUM IS: %lf\n",qtot);


	dgesv_(&N,&NRHS,aux_mat,&LDA,ipiv,aux_vec,&LDB,&info);
	
	

	if( info > 0 ) 
	{
        	printf( "The diagonal element of the triangular factor of A,\n" );
       		 printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
        	printf( "the solution could not be computed.\n" );
        	exit(1);
	}

	for(int i=0; i<r_nodes; i++)
	{
		hartree_vec[i] = aux_vec[i];
	        //printf("HP[%d] = % lf\n",i,aux_vec[i]);
	}

	free(aux_mat);

}
