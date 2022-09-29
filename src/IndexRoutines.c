#include "main.h"


int Index_i(parameters par, int idom, int i)
{
	int Ns = par.ns[idom] - 1, i1 = i;

	if (i1 < 0)
		i1 = -i1;
	if (i1 > Ns)
		i1 = 2 * Ns - i1;

	return i1;
}

int Index_j(parameters par, int idom, int j)
{
	int Nt = par.nt - 1, j1 = j;

	if (j1 < 0)
		j1 = -j1;
	if (j1 > Nt)
		j1 = 2 * Nt - j1;

	return j1;
}

int Index(parameters par, int idom, int ipot, int i, int j)
{// nst[idom] = Sum[ ns[jdom]*nt[jdom] , {jdom, 0, idom-1} ]; nst[0] = 0; see page 6
	int ns = par.ns[idom], nst = par.nst[idom], ii = Index_i(par, idom, i), jj =
			Index_j(par, idom, j);

	return ipot + NPOT * (ii + ns * jj + nst);
}

void Get_Arrays_n_And_n_v_domain_0(parameters *par, int *n)
{	 //Inside the Colony
	int idom = 0, ns = (*par).ns[idom], Ns = ns - 1, nt = (*par).nt, Nt = nt
			- 1, i, j, ipot, n_v;

	for (i = 0; i <= Ns; i++)
	{     // Indices of v^{0;ipot}_{ij}
		for (j = 0; j <= Nt; j++)
		{
			for (ipot = 0; ipot < NPOT; ipot++)
			{
				n_v = Index(*par, idom, ipot, i, j);
				(*par).n_v_of_n[*n] = n_v;
				(*par).n_of_n_v[n_v] = *n;
				*n += 1;
			}
		}
	}
}

void Get_Arrays_n_And_n_v_domain_1(parameters *par, int *n)
{	// Outside the Colony
	int idom = 1, ns = (*par).ns[idom], Ns = ns - 1, nt = (*par).nt, Nt = nt
			- 1, i, j, ipot, n_v;

	for (i = Ns; i >= 0; i--)
	{  // Indices of v^{1;ipot}_{ij}
		for (j = 0; j <= Nt; j++)
		{
			for (ipot = 0; ipot < NPOT; ipot++)
			{
				n_v = Index(*par, idom, ipot, i, j);
				if (i == Ns)
					(*par).n_of_n_v[n_v] = -1;
				else
				{
					(*par).n_v_of_n[*n] = n_v;
					(*par).n_of_n_v[n_v] = *n;
					*n += 1;
				}
			}
		}
	}
}

void Get_Arrays_n_And_n_v(parameters *par)
{
	int n;

	(*par).num_v = NPOT * (*par).nst[NDOM]; // #(all Potential-values v)
	(*par).n_of_n_v = ivector(0, (*par).num_v - 1);
	(*par).n_v_of_n = ivector(0, (*par).num_v - 1);

	n = 0;
	Get_Arrays_n_And_n_v_domain_0(par, &n);
	Get_Arrays_n_And_n_v_domain_1(par, &n);

	(*par).n_2D = n;
	(*par).ntotal = (*par).n_2D;
}

void Get_Indices_From_n_v(parameters par, int n_v, int *idom, int *ipot, int *i,
		int *j)
{
	int aux1, aux2, ns;

	*idom = NDOM;
	while (n_v < NPOT * par.nst[*idom])
		*idom -= 1;

	if (*idom < NDOM)
	{
		aux1 = n_v - NPOT * par.nst[*idom];
		ns = par.ns[*idom];
		*j = aux1 / (NPOT * ns);
		aux2 = aux1 - NPOT * ns * (*j);
		*i = aux2 / NPOT;
		*ipot = aux2 - NPOT * (*i);
	}
}

void Free_Arrays_n_And_n_v(parameters *par)
{
	free_ivector((*par).n_of_n_v, 0, (*par).num_v - 1);
	free_ivector((*par).n_v_of_n, 0, (*par).num_v - 1);
}

