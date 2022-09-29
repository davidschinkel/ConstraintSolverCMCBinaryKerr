#include "main.h"

void allocate_derivs_2D(derivs_2D *v, int n)
{
	(*v).d0 = new ftype[n];
	(*v).d1 = new ftype[n];
	(*v).d2 = new ftype[n];
	(*v).d11 = new ftype[n];
	(*v).d12 = new ftype[n];
	(*v).d22 = new ftype[n];
}

void fill0_derivs_2D(derivs_2D v, int n)
{
	fill0_dvector(v.d0, 0, n);
	fill0_dvector(v.d1, 0, n);
	fill0_dvector(v.d2, 0, n);
	fill0_dvector(v.d11, 0, n);
	fill0_dvector(v.d12, 0, n);
	fill0_dvector(v.d22, 0, n);
}

void free_derivs_2D(derivs_2D *v)
{
	delete[] (*v).d0;
	delete[] (*v).d1;
	delete[] (*v).d2;
	delete[] (*v).d11;
	delete[] (*v).d12;
	delete[] (*v).d22;
}

void allocate_derivs_2D_3rd_derivatives(derivs_2D *v, int n)
{
	(*v).d111 = new ftype[n];
	(*v).d112 = new ftype[n];
	(*v).d122 = new ftype[n];
	(*v).d222 = new ftype[n];
}

void free_derivs_2D_3rd_derivatives(derivs_2D *v)
{
	delete[] (*v).d111;
	delete[] (*v).d112;
	delete[] (*v).d122;
	delete[] (*v).d222;
}
