#include <stdio.h>
#include <stdlib.h>
void ComputeDensity(double *rho, double *wfn,int Ne,int order,int atomicN)
{
	int nodes = Ne*order + 1;
	int fembasis = nodes-2;
	int Occ = 1;
	for(int i=0; i<fembasis; i++)
	{
		rho[i]=0.0;
		for(int orb=0; orb<Occ; orb++)
		{
			rho[i] += wfn[i + orb*fembasis]*wfn[i + orb*fembasis];
			//printf("%lf\n",rho[i]);
			rho[i] = 2.0*rho[i];
		//	printf("%lf\n",rho[i]);
		}
	}
}
