#include <stdio.h>
#include <stdlib.h>




void getDensityMatrix(double *densMat,double *wfn, int fembasis,int atomicN)
{
	int k=0;
	int occ = 1;
	for(int i=0; i<fembasis; i++)
	{
		for(int j=0; j<fembasis; j++)
		{
			for(int orb=0; orb<occ; orb++)
			{
				densMat[k] =  2.0*(wfn[i + orb*fembasis]*wfn[j + orb*fembasis]);
			}
			k++;
		}
	}
}
