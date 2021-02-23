#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"


void getFemNodes(double *femnodes, struct Element *e, struct Vertex *n, int Ne, int order){


	int nodes = Ne*order + 1;
	int fembasis = nodes-2;
	double aux[nodes]; 

	aux[0] = e[0].n[0].x;
	int l=1;
	double factor;
        for(int i=0; i<Ne; i++)
        {
                for(int j=1; j<(order+1); j++)
                {
                        factor = (e[i].n[j].x);
                        aux[l] =  factor;
                        l++;
                }
        }
	for(int i=0; i<fembasis; i++){
		femnodes[i] = aux[i+1];
	
	}

}
