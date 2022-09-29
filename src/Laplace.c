#include "main.h"

ftype LaplaceSolution(parameters par,int idom, int i,int j)
{
	ftype rho = par.rho.d0[Index(par,idom,0,i,j)],
		  z = par.z.d0[Index(par,idom,0,i,j)];
	return jn(0,0.1*rho)*cosh(0.1*z);
}

void Check_Laplace(parameters par, ftype *X)
{
	ftype diff;
	derivs_2D v;

	allocate_derivs_2D(&v, par.ntotal);
	Get_v_From_X(par,X,v);
	for(int idom = 0; idom<NDOM; idom++)
	{
		ftype **psi0 = new ftype*[par.ns[idom]];
		ftype **psi1 = new ftype*[par.ns[idom]];
		ftype **c0 = new ftype*[par.ns[idom]];
		ftype **c1 = new ftype*[par.ns[idom]];
		for (int i = 0; i < par.ns[idom]; i++)
		{
			psi0[i] = new ftype[par.nt];
			psi1[i] = new ftype[par.nt];
			c0[i] = new ftype[par.nt];
			c1[i] = new ftype[par.nt];
		}

		for(int i=0; i<par.ns[idom]; i++)
		{
			for(int j=0; j<par.nt; j++)
			{
				psi0[i][j]=LaplaceSolution(par,idom,i,j);
				psi1[i][j]=v.d0[Index(par,idom,0,i,j)];
			}
		}

		Chebyshev_Coefficients_Lobatto_2D(psi0, c0, par.ns[idom], par.nt);
		Chebyshev_Coefficients_Lobatto_2D(psi1, c1, par.ns[idom], par.nt);

		diff=-1;

		for(int i=0; i<par.ns[idom]; i++)
		{
			for(int j=0; j<par.nt; j++)
			{
				diff = (fabs(c0[i][j]-c1[i][j])>diff) ? fabs(c0[i][j]-c1[i][j]) : diff;
			}
		}
		cout << "Abweichung in den Chebyshevkoeffizienten im Gebiet" << idom << " von " << scientific << diff << fixed << "." << endl;

		for (int i = 0; i < par.ns[idom]; i++)
		{
			delete[] psi0[i];
			delete[] psi1[i];
			delete[] c0[i];
			delete[] c1[i];
		}
		delete[] psi0;
		delete[] psi1;
		delete[] c0;
		delete[] c1;
	}
}
