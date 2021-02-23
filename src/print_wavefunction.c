#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"
#include "wfn_operations.h"
#include <math.h>
void PrintWaveFunction(int Ne,int order,double *femnodes,double *wfn,double *vecE,double *hpotential)
{
	
	FILE *wfn_data,*dens_data,*hp_data,*analitic,*dens_analitic;
	double cfwfn,coefrho;
	int orb = 0;
	int nodes = Ne*order+1;
	int r_nodes = nodes-2;
	int phase=1;
	double rhoOriginal[r_nodes];
	wfn_data = fopen("wfn.dat","w");
	dens_data = fopen("dens.dat","w");
	hp_data = fopen("hp.dat","w");
	analitic = fopen("hp_analitic.dat","w");
	dens_analitic = fopen("dens_analitic.dat","w");
	/*if(wfn_data==NULL);
	{
		printf("ERROR IN CREATING THE WFN\n");
		exit(1);
	}*/
	printf("Wave Function in: wfn.dat\n");
	fprintf(wfn_data,"#Node\tPhi(r_i)\tRho(r_i)\n");
	GetWfnPhase(r_nodes,orb,&phase,wfn);
	printf("Wave function phase: %d\n",phase);
	for(int i=0; i<r_nodes; i++)
	{
		cfwfn = wfn[i + r_nodes*orb]*phase;
		coefrho = cfwfn*cfwfn;
		rhoOriginal[i] = coefrho;
		//cfwfn += 0.5*vecE[orb];
		//coefrho += 0.5*vecE[orb];
		fprintf(wfn_data,"%d\t% lf\t% lf\n",i,cfwfn,coefrho);

		/*for(int j=0;j<order;j++)
		{
			cwfn = wfn[(k)]
			fprintf(wfn_data,"%d\t%lf\t%lf\n",k,e[i].n[j].x,wfn[k]);
			k++;

		}*/
	}
	double rhodiv[r_nodes];
	printf("******************* Values at the boundary  *****************\n");
	printf("                    Vh(0) = 0  Vh(Rc) = qtot                 \n");
	printf("******************* Values at the grid points ***************\n");
	/*int l=1;
	for(int i=0; i<Ne; i++)
	{
		for(int j=1; j<(order+1); j++)
		{
			factor = (e[i].n[j].x);
			rhodiv[l-1] = rhoOriginal[l-1]/factor;
			double dens = 4.0*(e[i].n[j].x)*exp(-2.0*e[i].n[j].x);
			//double dens = 4.0*exp(-2.0*e[i].n[j].x);
			double pot = -1.0*(e[i].n[j].x + 1)*exp(-2.0*e[i].n[j].x) + 1;
			fprintf(dens_data,"%lf\t %lf\n",e[i].n[j].x,rhoOriginal[l-1]);
			fprintf(hp_data,"%lf\t %lf\n",e[i].n[j].x,hpotential[l-1]);
printf("%d  r = % lf    Vh = % lf    Vh_exact = % lf     rho_analitic = % lf      rho = %lf   rho_div = %lf\n",l-1,e[i].n[j].x,hpotential[l-1],pot,dens,rhoOriginal[l-1],rhodiv[l-1]);
			l++;
		}
	}*/
	double factor;
	for(int i=0; i<r_nodes; i++){
		factor = femnodes[i];
		rhodiv[i] = rhoOriginal[i]/factor;
		double dens = 4.0*(factor)*exp(-2.0*factor);
		double pot = -1.0*(factor + 1)*exp(-2.0*factor) + 1;
		fprintf(dens_data,"%lf\t %lf\n",factor,rhodiv[i]);
		fprintf(hp_data,"%lf\t %lf\n",factor,hpotential[i]);
		fprintf(analitic,"%lf\t %lf\n",factor,pot);
		fprintf(dens_analitic,"%lf\t %lf\n",factor,dens);
		printf("%d  r = % lf    Vh = % lf    Vh_exact = % lf     rho_analitic = % lf      rho = %lf   rho_div = %lf\n",i,factor,hpotential[i],pot,dens,rhoOriginal[i],rhodiv[i]);
	}
	printf("**************************************************************\n");
	//double dens = 4.0*(e[Ne-1].n[order].x)*exp(-2.0*e[Ne-1].n[order].x);
        //double pot = -1.0*(e[Ne-1].n[order].x + 1)*exp(-2.0*e[Ne-1].n[order].x) + 1;
	//printf("r = %lf        dens = %lf         hpot = %lf\n",e[Ne-1].n[order].x,dens,pot);
	fclose(wfn_data);
	fclose(dens_data);
	fclose(hp_data);
	fclose(analitic);
	fclose(dens_analitic);
}
