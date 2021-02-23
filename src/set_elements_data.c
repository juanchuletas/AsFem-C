#include <stdio.h>
#include <stdlib.h>
#include "data_structures.h"

void SetElements(struct Element *e,struct Vertex *n,int Ne,double *x,int order,double rN)
{
	int p;
	p=order+1;
	switch(order)
	{
		case 1:
			for(int i=0; i<Ne; i++)
			{
				e[i].n = (struct Vertex *)malloc(sizeof(struct Vertex)*2);	
		       		e[i].n[0].x = x[i];
		       		e[i].n[1].x = x[i+1];
	       			e[i].h = e[i].n[1].x - e[i].n[0].x;
				//printf("e[%d] = %lf\n",i,e[i].h);
        		}
			/*for(int i=0; i<Ne; i++)
			{
				for(int j=0; j<p-1; j++){


					printf("e[%d].n[%d].x = %lf\n",i,j,e[i].n[j].x);
				
				}
			}
					printf("e[%d].n[%d].x = %lf\n",Ne,1,e[Ne-1].n[1].x);*/
	                break;
		case 2:
        		for(int i=0; i<Ne; i++)
       			{       
				e[i].n = (struct Vertex *)malloc(sizeof(struct Vertex)*3);	
				e[i].n[0].x = x[i];
			       	e[i].n[1].x = (x[i]+x[i+1])/2.0;
			       	e[i].n[2].x = x[i+1];
			      	e[i].h = e[i].n[2].x - e[i].n[0].x;
			//	printf("e[%d] = %lf\n",i,e[i].h);
        		}
			/*for(int i=0; i<Ne; i++)
			{
				for(int j=0; j<p-1; j++){


					printf("e[%d].n[%d].x = %lf\n",i,j,e[i].n[j].x);
				
				}
			}
					printf("e[%d].n[%d].x = %lf\n",Ne-1,order,e[Ne-1].n[order].x);*/
			break;
		case 3:
        		for(int i=0; i<Ne; i++)
       			{       
				e[i].n = (struct Vertex *)malloc(sizeof(struct Vertex)*p);	
				e[i].n[0].x = x[i];
			       	e[i].n[1].x = x[i] + (x[i+1]-x[i])/3.0;
			       	e[i].n[2].x = x[i] + (x[i+1]-x[i])*(2.0/3.0);
			       	e[i].n[3].x = x[i+1];
			      	e[i].h = e[i].n[3].x - e[i].n[0].x;
        		}
			/*for(int i=0; i<Ne; i++)
			{
				for(int j=0; j<p-1; j++){


					printf("e[%d].n[%d].x = %lf\n",i,j,e[i].n[j].x);
				
				}
			}
					printf("e[%d].n[%d].x = %lf\n",Ne-1,order,e[Ne-1].n[order].x);*/
			break;



	};
}
