#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"
#include "Elemental_Matrices.h"
extern double *ElementalOverlapMatrix(int order);
void GetDens(double *rho,double *wfn,int nodes)
{
	int fembasis = nodes-2;
	int Occ = 1;
	for(int i=0; i<fembasis; i++)
        {
                rho[i] = 0.0;
                for(int orb=0; orb<Occ; orb++)
                {
                        rho[i] += wfn[i + orb*fembasis]*wfn[i + orb*fembasis];
                       	//rho[i] = 2.0*rho[i];
			//printf("Density = %lf\n",rho[i]);

                }
        }

}


double IntegrateDens(double *wfn,struct Element *e,int *link_mat,int Ne, int order)
{
	int nodes = Ne*order+1;
	int fembasis = nodes-2;
	int p=order+1;
	double dens[fembasis];
	double rho[nodes];
	rho[0] = 0.0;
	rho[nodes-1] = 0.0;
	double integ[p];
	integ[0] = 1.f;
	integ[1] = 1.f;
	//integ[2] = 1.3333;
	double *eMatS;
	eMatS = ElementalOverlapMatrix(order);

	GetDens(dens,wfn,nodes);

	for(int i=1; i<nodes-1; i++)
	{
		rho[i] = dens[i-1];

	}
	double qtot = 0.0;
	for(int ei=0; ei<Ne; ei++)
	{
		double coeff = 0.5*e[ei].h;
		for(int mu=0; mu<p; mu++)
		{
			int l = link_mat[p*ei + mu];
			qtot = qtot + coeff*rho[l]*eMatS[p*mu + mu];

		}
	}
	printf("MODULE CHARGE: %lf\n",qtot);

	return qtot;

}
