#include<stdio.h>
#include<stdlib.h>


double getHFenergy(double *hcore, double *fmat, double *densMat,int totSize)
{
	double energy = 0.0;

	for(int i=0; i<totSize; i++)
	{
		energy = energy + 0.5*(hcore[i]*densMat[i]);
		energy = energy + 0.5*(fmat[i]*densMat[i]);
	}
	return 0.5*energy;


}
