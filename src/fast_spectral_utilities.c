#include "fast_spectral_utilities.h"

#define NDOM 2 //TODO Die Konstanten auslagern, NDOM wird zweimal definiert

void fast_Chebyshev_Coefficients_Lobatto(ftype *psi, ftype *c, int N, fftw_plan plan)
{//c output, psi input, N size; plan: FFTW_REDFT00

	//create & execute plan
	fftw_execute_r2r(plan, psi, c);
	//normalize result
	double norm = 1/ftype(N);
	for(int i=0; i<=N; i++)
	{
	  c[i] = c[i]*norm;
	}
	//the last coefficient is different
	c[N]=c[N]/ftype(2);
}

void fast_Chebyshev_Coefficients_Gauss(double *psi, double *c, int N,fftw_plan plan) 
{//c output, psi input, N size; plan: FFTW_REDFT10

	//create & execute plan
	fftw_execute_r2r(plan, psi, c);

	double norm = 1./(N+1);

	for(int i=0; i<=N; i++)
	{
		c[i] = c[i]*norm;
	}
}

void fast_Chebyshev_Collocations_Lobatto(double *psi, double *c, int N,fftw_plan plan)
{//c input, psi output, N size

	//create & execute plan
	fftw_execute_r2r(plan, c, psi);

	for(int i=0; i<=N; i++)
	{
		psi[i] = psi[i]/ftype(2);
	}
}

void fast_Chebyshev_Collocations_Gauss(ftype *psi, ftype *c, int N,fftw_plan plan)
{//c input, psi output, N size

	//create & execute plan
	fftw_execute_r2r(plan, c, psi);

	for(int i=0; i<=N; i++)
	{
		psi[i] = psi[i]/ftype(2);
	}
}
