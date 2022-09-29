//compile with: g++ -lm -O3 *.c -o Konvergenz -lm
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <c++/4.8/iostream>
#include <c++/4.8/fstream>
#include "functions.h"

 #define NDOM 2

using namespace std;

int main(int argc, char **argv)
{
	long double **pot[2], **pot_ref[2], diff[NDOM], wert;
//	long double sh = 1.5; //gesetzter scheinbarer horizont
	FILE *fp, *fp_ref;
	char string[255], *name, *name_ref; //ein Dummy, der Dateiname der Felder und der Referenz
	int n_A_ref[2], n_B_ref, n_A[2], n_B, nA, nB, n_check, i_aux, dom, ipot, i, j;
	int n_points_sigma, n_points_mu; //Anzahl der zu pruefenden Punkte

	n_check = 1000;

	if(3==argc)
	{
		name=argv[1];
		name_ref=argv[2];
		//cout << name <<"\t\t" << name_ref << endl;
	}
	else
	{
		printf("Fehler: Argumente in der Form 'Konvergenz Feld Referenzfeld' übergeben.\n");
		return -1;
	}

	fp = fopen(name_ref,"r");
	for(int i=1; i<=12; i++)
	{
		i_aux = fscanf(fp,"%Lf %s\n",&wert,string);
	}
	i_aux = fscanf(fp,"%d %s\n",&n_A_ref[0],string);
	i_aux = fscanf(fp,"%d %s\n",&n_A_ref[1],string);
	i_aux = fscanf(fp,"\n");
	i_aux = fscanf(fp,"%d %s\n",&n_B_ref,string);

	//cout << "n_A_ref[0]=" <<n_A_ref[0] <<"\tn_A_ref[1]="<<n_A_ref[1] << "\tn_B_ref="<<n_B_ref<< endl;

	pot_ref[0]=ldmatrix(0,n_A_ref[0],0,n_B_ref);
	pot_ref[1]=ldmatrix(0,n_A_ref[1],0,n_B_ref);	

	for(int k=1; k <= (n_A_ref[0]+n_A_ref[1])*n_B_ref*NDOM; k++)
	{
		i_aux = fscanf(fp,"%Lf \t (idom,ipot,i,j)=(%d,%d,%d,%d)\n",&wert, &dom, &ipot, &i, &j);
		pot_ref[dom][i][j]=wert;
	}
	fclose(fp);

	fp = fopen(name,"r");
	for(int i=1; i<=12; i++)
	{
		i_aux = fscanf(fp,"%Lf %s\n",&wert,string);
	}
	i_aux = fscanf(fp,"%d %s\n",&n_A[0],string);
	i_aux = fscanf(fp,"%d %s\n",&n_A[1],string);
	i_aux = fscanf(fp,"\n");
	i_aux = fscanf(fp,"%d %s\n",&n_B,string);

	//cout << "n_A[0]=" <<n_A[0] <<"\tn_A[1]="<<n_A[1] << "\tn_B="<<n_B<< endl;

	pot[0]=ldmatrix(0,n_A[0],0,n_B);
	pot[1]=ldmatrix(0,n_A[1],0,n_B);

	for(int k=1; k <= (n_A[0]+n_A[1])*n_B*NDOM; k++)
	{
		i_aux = fscanf(fp,"%Lf \t (idom,ipot,i,j)=(%d,%d,%d,%d)\n",&wert, &dom, &ipot, &i, &j);
		pot[dom][i][j]=wert;
	}
	fclose(fp);

	chebft_Extremes_2D(pot_ref[0], n_A_ref[0], n_B_ref);
	chebft_Extremes_2D(pot_ref[1], n_A_ref[1], n_B_ref);
	chebft_Extremes_2D(pot[0], n_A[0], n_B);
	chebft_Extremes_2D(pot[1], n_A[1], n_B);

	//Vergleich der Lösungen auf einem n_check x n_check Gitter
	// for(int idom = 0; idom<NDOM; idom++)
	// {
	// 	diff[idom]=-1;
	// 	for (long double x = 0; x <= 1; x = x + (1 /((long double) n_check)))
	// 	{
	// 		for (long double y = -1; y <= 1; y = y + (2 / ((long double) n_check)))
	// 		{
	// 			if (fabsl(chebev_2D(0, 1, -1, 1, pot[idom], n_A[idom], n_B, x, y)
	// 					- chebev_2D(0, 1, -1, 1, pot_ref[idom], n_A_ref[idom], n_B_ref, x, y)) > diff[idom]);
	// 			{
	// 				diff[idom] = fabsl(chebev_2D(0, 1, -1, 1, pot[idom], n_A[idom], n_B, x, y)
	// 								 - chebev_2D(0, 1, -1, 1, pot_ref[idom], n_A_ref[idom], n_B_ref, x, y));
	// 			}
	// 		}
	// 	}
	// }


	//Vergleich der Koeffizienten
	for(int idom = 0; idom<NDOM; idom++)
	{
		diff[idom]=-1;

		nA = (n_A[idom]<n_A_ref[idom])? n_A[idom] : n_A_ref[idom];
		nB = (n_B<n_B_ref)? n_B : n_B;

		for (int i=0; i<nA; i++)
		{
			for (int j=0; j<nB; j++)
			{
				if(fabsl(pot_ref[idom][i][j] - pot[idom][i][j])>diff[idom])
				{
					diff[idom] = fabsl(pot_ref[idom][i][j] - pot[idom][i][j]);	
				}
			}
		}
	}

	printf("%d \t %Le \t %Le\n",n_B,diff[0],diff[1]);

	return 0;
}
