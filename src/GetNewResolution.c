#include "main.h"


// Calculating an approximate solution corresponding to new resolution orders through
// spectral interpolation

void v1_To_v2(parameters par1, derivs_2D v1, parameters par2, derivs_2D v2)
{
	int ipot, idom, i_1, j_1, i_2, j_2, indx;
	ftype p[NMAX], *v21;
	parameters par21 = par1;

	par21.nst[0] = 0;
	par21.nt = par1.nt;
	for (idom = 0; idom < NDOM; idom++)
	{
		par21.ns[idom] = par2.ns[idom];
		par21.nst[idom + 1] = par21.nst[idom] + par21.ns[idom] * par21.nt;
	}
	v21 = new ftype[NPOT * par21.nst[NDOM]];

	for (idom = 0; idom < NDOM; idom++)
	{
		int ns_1 = par1.ns[idom], nt_1 = par1.nt, ns_2 = par2.ns[idom], nt_2 =
				par2.nt, Ns_2 = ns_2 - 1, Nt_2 = nt_2 - 1;
		for (ipot = 0; ipot < NPOT; ipot++)
		{
			for (j_1 = 0; j_1 < nt_1; j_1++)
			{
				for (i_1 = 0; i_1 < ns_1; i_1++)
				{
					indx = Index(par1, idom, ipot, i_1, j_1);
					p[i_1] = v1.d0[indx];
				}
				chebft_Extremes(p, ns_1, 0);
				for (i_2 = 0; i_2 < ns_2; i_2++)
				{
					ftype a = Pih * i_2 / Ns_2, sa = sin(a), s = sa * sa;
					indx = Index(par21, idom, ipot, i_2, j_1);
					v21[indx] = chebev(0., 1., p, ns_1, s);
				}
			}
			for (i_2 = 0; i_2 < ns_2; i_2++)
			{
				for (j_1 = 0; j_1 < nt_1; j_1++)
				{
					indx = Index(par21, idom, ipot, i_2, j_1);
					p[j_1] = v21[indx];
				}
				chebft_Extremes(p, nt_1, 0);
				for (j_2 = 0; j_2 < nt_2; j_2++)
				{
					ftype b = Pih * j_2 / Nt_2, sb = sin(b), t = sb * sb;
					indx = Index(par2, idom, ipot, i_2, j_2);
					v2.d0[indx] = chebev(0., 1., p, nt_1, t);
				}
			}
		}
	}
	delete[] v21;
}

void GetNewResolution(parameters par0, parameters par, ftype *X0, ftype *X)
{
	int n_v;
	derivs_2D v0, v;

	allocate_derivs_2D(&v0, par0.num_v);
	Get_v_From_X(par0, X0, v0);

	allocate_derivs_2D(&v, par.num_v);
	v1_To_v2(par0, v0, par, v);

	for (int n = 0; n < par.n_2D; n++)
	{
		//n_v = par.n_v_of_n[n];
		n_v = n;
		X[n] = v.d0[n_v];
	}

	free_derivs_2D(&v);
	free_derivs_2D(&v0);
}
