#include<stdio.h>
#include <stdlib.h>
void ReduceMatrices(double *sij,double *kij,double *vij,double *uij,double *lij,double *s_mat,double *k_mat,double *v_mat,double *u_mat,double *l_mat,int nodes)
{
	int l=0;
	int n = nodes;
	for(int i=1;i<n-1;i++)
	{
		for(int j=1;j<n-1;j++)
		{
			sij[l] = s_mat[i*n+j];
		       	vij[l] = v_mat[i*n+j];
			kij[l] = k_mat[i*n+j];
			uij[l] = u_mat[i*n+j];
			lij[l] = l_mat[i*n+j];
			l++;
		}
    	}

}
