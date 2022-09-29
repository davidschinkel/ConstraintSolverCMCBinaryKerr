#include "main.h"

void Derivatives_2D(parameters par, derivs_2D v)
{
	int ns, nt, indx[NMAX];
	ftype *psi, *psid1, *c, *dc, *ddc, facA, fac2A, facB, fac2B;
	facA = -2;
	fac2A = sqr(facA);
	facB = -1;
	fac2B = 1;

	for (int idom = 0; idom < NDOM; idom++)
	{
		ns = par.ns[idom];
		nt = par.nt;

		for (int ipot = 0; ipot < NPOT; ipot++)
		{
			psi = new ftype[ns];
			c = new ftype[ns];
			dc = new ftype[ns];
			ddc = new ftype[ns];

			for (int j = 0; j < nt; j++)
			{		// Calculation of Derivatives w.r.t. s-Direction
				for (int i = 0; i < ns; i++)
				{	// (Chebyshev_Extremes)
					indx[i] = Index(par, idom, ipot, i, j);
					psi[i] = v.d0[indx[i]];
				}
				fast_Chebyshev_Coefficients_Lobatto(psi, c, ns - 1, par.FFTWplanLobatto[idom][0]);
				Chebyshev_Coefficients_Derivative(c, dc, ns - 1); //1. Ableitung
				Chebyshev_Coefficients_Derivative(dc, ddc, ns - 1); //2. Ableitung

				fast_Chebyshev_Collocations_Lobatto(psi, dc, ns - 1, par.FFTWplanLobatto[idom][0]);
				for (int i = 0; i < ns; i++)
				{
					v.d1[indx[i]] = psi[i]* facA;
				}

				fast_Chebyshev_Collocations_Lobatto(psi, ddc, ns - 1, par.FFTWplanLobatto[idom][0]);
				for (int i = 0; i < ns; i++)
				{
					v.d11[indx[i]] = psi[i]* fac2A;
				}
			}
			delete[] psi;
			delete[] c;
			delete[] dc;
			delete[] ddc;
			//derivative with respect to B
			psi = new ftype[nt];
			psid1 = new ftype[nt];
			c = new ftype[nt];
			dc = new ftype[nt];
			ddc = new ftype[nt];

			for (int i = 0; i < ns; i++)
			{// Calculation of Derivatives w.r.t. t-Direction
				for (int j = 0; j < nt; j++)
				{	// (Chebyshev_Extremes)
					indx[j] = Index(par, idom, ipot, i, j);
					psi[j] = v.d0[indx[j]];
					psid1[j] = v.d1[indx[j]];
				}
				fast_Chebyshev_Coefficients_Lobatto(psi, c, nt - 1, par.FFTWplanLobatto[idom][1]);
				Chebyshev_Coefficients_Derivative(c, dc, nt - 1); //1. Ableitung
				Chebyshev_Coefficients_Derivative(dc, ddc, nt - 1); //2. Ableitung

				fast_Chebyshev_Collocations_Lobatto(psi, dc, nt - 1, par.FFTWplanLobatto[idom][1]);
				for (int j = 0; j < nt; j++)
				{
					v.d2[indx[j]] = psi[j]*facB;
				}
				fast_Chebyshev_Collocations_Lobatto(psi, ddc, nt - 1, par.FFTWplanLobatto[idom][1]);
				for (int j = 0; j < nt; j++)
				{
					v.d22[indx[j]] = psi[j]*fac2B;
				}
				//gemischte Ableitungen
				fast_Chebyshev_Coefficients_Lobatto(psid1, c, nt - 1, par.FFTWplanLobatto[idom][1]);
				Chebyshev_Coefficients_Derivative(c, dc, nt - 1); //d12
				fast_Chebyshev_Collocations_Lobatto(psi, dc, nt - 1, par.FFTWplanLobatto[idom][1]);
				for (int j = 0; j < nt; j++)
				{
					v.d12[indx[j]] = psi[j]*facB;
				}
			}
			delete[] psi;
			delete[] psid1;
			delete[] c;
			delete[] dc;
			delete[] ddc;
		}
	}
}

void Derivatives_1D(parameters par, derivs_1D v, ftype fac, int n)
{
	ftype *psi, *c, *dc, *ddc, fac2;
	fac2 = sqr(fac);

	psi = new ftype[n];
	c = new ftype[n];
	dc = new ftype[n];
	ddc = new ftype[n];

	for (int i = 0; i < n; i++)
	{	// (Chebyshev_Extremes)
		psi[i] = v.d0[i];
	}
	Chebyshev_Coefficients_Lobatto(psi, c, n - 1);
	Chebyshev_Coefficients_Derivative(c, dc, n - 1); //1. Ableitung
	Chebyshev_Coefficients_Derivative(dc, ddc, n - 1); //2. Ableitung

	Chebyshev_Collocations_Lobatto(psi, dc, n - 1);
	for (int i = 0; i < n; i++)
	{
		v.d1[i] = psi[i]* fac;
	}
	Chebyshev_Collocations_Lobatto(psi, ddc, n - 1);
	for (int i = 0; i < n; i++)
	{
		v.d11[i] = psi[i]* fac2;
	}

	delete[] psi;
	delete[] c;
	delete[] dc;
	delete[] ddc;
}

void Derivatives_2D_3(parameters par, derivs_2D v)
{//Get the 3rd derivatives - assuming the first & second are already available
	int ns, nt, indx[NMAX];
	ftype *psi, *psid11, *psid12, *psid22, *c, *dc, facA, facB;
	facA = -2;
	facB = -1;

	for (int idom = 0; idom < NDOM; idom++)
	{
		ns = par.ns[idom];
		nt = par.nt;

		for (int ipot = 0; ipot < NPOT; ipot++)
		{
			psi = new ftype[ns];
			psid11 = new ftype[ns];
			c = new ftype[ns];
			dc = new ftype[ns];

			for (int j = 0; j < nt; j++)
			{		// Calculation of Derivatives w.r.t. s-Direction
				for (int i = 0; i < ns; i++)
				{	// (Chebyshev_Extremes)
					indx[i] = Index(par, idom, ipot, i, j);
					psid11[i] = v.d11[indx[i]];
				}
				fast_Chebyshev_Coefficients_Lobatto(psid11, c, ns - 1, par.FFTWplanLobatto[idom][0]);
				Chebyshev_Coefficients_Derivative(c, dc, ns - 1); //1. Ableitung

				fast_Chebyshev_Collocations_Lobatto(psi, dc, ns - 1, par.FFTWplanLobatto[idom][0]);
				for (int i = 0; i < ns; i++)
				{
					v.d111[indx[i]] = psi[i]* facA;
				}
			}
			delete[] psi;
			delete[] psid11;
			delete[] c;
			delete[] dc;
			//derivative with respect to B
			psi = new ftype[nt];
			psid11 = new ftype[nt];
			psid12 = new ftype[nt];
			psid22 = new ftype[nt];
			c = new ftype[nt];
			dc = new ftype[nt];

			for (int i = 0; i < ns; i++)
			{// Calculation of Derivatives w.r.t. t-Direction
				for (int j = 0; j < nt; j++)
				{	// (Chebyshev_Extremes)
					indx[j] = Index(par, idom, ipot, i, j);
					psid11[j] = v.d11[indx[j]];
					psid12[j] = v.d12[indx[j]];
					psid22[j] = v.d22[indx[j]];
				}
				fast_Chebyshev_Coefficients_Lobatto(psid11, c, nt - 1, par.FFTWplanLobatto[idom][1]);
				Chebyshev_Coefficients_Derivative(c, dc, nt - 1); //1. Ableitung
				fast_Chebyshev_Collocations_Lobatto(psi, dc, nt - 1, par.FFTWplanLobatto[idom][1]);
				for (int j = 0; j < nt; j++)
				{
					v.d112[indx[j]] = psi[j]*facB;
				}


				fast_Chebyshev_Coefficients_Lobatto(psid12, c, nt - 1, par.FFTWplanLobatto[idom][1]);
				Chebyshev_Coefficients_Derivative(c, dc, nt - 1); //1. Ableitung
				fast_Chebyshev_Collocations_Lobatto(psi, dc, nt - 1, par.FFTWplanLobatto[idom][1]);
				for (int j = 0; j < nt; j++)
				{
					v.d122[indx[j]] = psi[j]*facB;
				}

				fast_Chebyshev_Coefficients_Lobatto(psid22, c, nt - 1, par.FFTWplanLobatto[idom][1]);
				Chebyshev_Coefficients_Derivative(c, dc, nt - 1); //1. Ableitung
				fast_Chebyshev_Collocations_Lobatto(psi, dc, nt - 1, par.FFTWplanLobatto[idom][1]);
				for (int j = 0; j < nt; j++)
				{
					v.d222[indx[j]] = psi[j]*facB;
				}
			}
			delete[] psi;
			delete[] psid11;
			delete[] psid12;
			delete[] psid22;
			delete[] c;
			delete[] dc;
		}
	}
}

void Get_v_From_X(parameters par, ftype *X, derivs_2D v)
{
	copy_dvector(v.d0,X,0,par.ntotal);
	Derivatives_2D(par, v);
}

void F_of_X(parameters par, ftype *X, derivs_2D v, ftype *F)
{
	int ns, nt;

	Get_v_From_X(par, X, v);

	for (int idom = 0; idom < NDOM; idom++)
	{
		ns = par.ns[idom];
		nt = par.nt;
		for (int j = 0; j < nt; j++)
		{
			for (int i = 0; i < ns; i++)
			{
				NonLinFieldEqns(par, idom, i, j, v, F);
			}
		}
	}
}

void J_times_DX(parameters par, ftype *X, ftype *DX, ftype *JDX)
{
	int ntotal = par.ntotal, ns, nt;
	derivs_2D v, Dv;

	allocate_derivs_2D(&v, ntotal);
	allocate_derivs_2D(&Dv, ntotal);

	Get_v_From_X(par, X, v);
	Get_v_From_X(par, DX, Dv);

	for (int idom = 0; idom < NDOM; idom++)
	{
		ns = par.ns[idom];
		nt = par.nt;
		for (int j = 0; j < nt; j++)
		{
			for (int i = 0; i < ns; i++)
			{
				LinFieldEqns(par, idom, i, j, v, Dv, JDX);
			}
		}
	}
	free_derivs_2D(&v);
	free_derivs_2D(&Dv);
}
