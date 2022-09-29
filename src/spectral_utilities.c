#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "utilities.h"
#include "spectral_utilities.h"

//---------------------------------------------------------------
// ********* Clenshaw Algorithms ********************************
//---------------------------------------------------------------
ftype Clenshaw_cos(ftype *alpha, int N, ftype phi)
{
	ftype bk = 0, bkp1 = 0, bkp2, cos_phi = cos(phi), c2 = 2 * cos_phi;
	int k;

	for (k = N; k >= 0; k--)
	{
		bkp2 = bkp1;
		bkp1 = bk;
		bk = alpha[k] + c2 * bkp1 - bkp2;
	}
	return (bk - bkp1 * cos_phi - 0.5 * alpha[0]);
}
//---------------------------------------------------------------
ftype Clenshaw_sin(ftype *bpsi, int N, ftype phi)
{
	ftype bk = 0, bkp1 = 0, bkp2, cos_phi = cos(phi), sin_phi = sin(phi), c2 = 2
			* cos_phi;
	int k;

	for (k = N; k >= 1; k--)
	{
		bkp2 = bkp1;
		bkp1 = bk;
		bk = bpsi[k] + c2 * bkp1 - bkp2;
	}
	return (bk * sin_phi);
}
//---------------------------------------------------------------
ftype Clenshaw_Fourier(ftype *alpha, ftype *bpsi, int N, ftype phi)
{
	ftype bck = 0, bckp1 = 0, bckp2, bsk = 0, bskp1 = 0, bskp2, cos_phi = cos(
			phi), sin_phi = sin(phi), c2 = 2 * cos_phi;
	int k;

	for (k = N; k >= 0; k--)
	{
		bckp2 = bckp1;
		bckp1 = bck;
		bck = alpha[k] + c2 * bckp1 - bckp2;
		if (k > 0)
		{
			bskp2 = bskp1;
			bskp1 = bsk;
			bsk = bpsi[k] + c2 * bskp1 - bskp2;
		}
	}
	return (bck - bckp1 * cos_phi - 0.5 * alpha[0] + bsk * sin_phi);
}
//---------------------------------------------------------------
ftype Clenshaw_Chebyshev(ftype *c, int N, ftype x)
{
	ftype bk = 0, bkp1 = 0, bkp2 = 0, x2 = 2 * x;
	int k;

	for (k = N; k >= 0; k--)
	{
		bkp2 = bkp1;
		bkp1 = bk;
		bk = c[k] + x2 * bkp1 - bkp2;
	}
	return (0.5 * (bk - bkp2));
}
//---------------------------------------------------------------
// ********* Chebyshev-Coefficients  ****************************
//---------------------------------------------------------------
void Chebyshev_Coefficients_Radau_RHS(ftype *psi, ftype *c, int N)
{
	int m;
	ftype fac1 = 1 / ((ftype) (2 * N + 1)), delta_phi = 2 * Pi* fac1, fac3 = 4 * fac1;

	for (m = 0; m <= N; m++)
		c[m] = Clenshaw_cos(psi, N, delta_phi * m) * fac3; // phi_m = 2*Pi*m/(2N+1)
}
//---------------------------------------------------------------
void Chebyshev_Coefficients_Radau_LHS(ftype *psi, ftype *c, int N)
{
	int m, sign = 1;
	ftype fac1 = 1 / ((ftype) (2 * N + 1)), delta_phi = 2 * Pi* fac1, fac3 = 4 * fac1;

	for (m = 0; m <= N; m++)
	{
		c[m] = Clenshaw_cos(psi, N, delta_phi * m) * fac3 * sign; // tilde_phi_m = 2*Pi*m/(2N+1)
		sign = -sign;
	}
}
//---------------------------------------------------------------
void Chebyshev_Coefficients_Gauss(ftype *psi, ftype *c, int N)
{
	int m;
	ftype fac1 = 1 / ((ftype) (N + 1)), delta_phi = Pi* fac1, fac3 =
	2 * fac1, aux = 0.5 * psi[0];

	for (m = 0; m <= N; m++)
	{
		ftype tilde_phi_m = delta_phi * m,  // tilde_phi_m = Pi*m/(N+1)
		arg = 0.5 * tilde_phi_m, fac_cos = fac3 * cos(arg), fac_sin = fac3
				* sin(arg);

		c[m] = fac_cos * (Clenshaw_cos(psi, N, tilde_phi_m) + aux)
				- fac_sin * Clenshaw_sin(psi, N, tilde_phi_m);
	}
}
//---------------------------------------------------------------
void Chebyshev_Coefficients_Lobatto(ftype *psi, ftype *c, int N)
{
	int m, sign = 1;
	ftype fac1 = 1 / ((ftype) N), delta_phi = Pi* fac1, fac3 = 2
	* fac1, aux = 0.5 * psi[N];

	for (m = 0; m <= N; m++)
	{
		c[m] = fac3 * (Clenshaw_cos(psi, N, delta_phi * m) - aux * sign); // phi_m = Pi*m/N
		sign = -sign;
	}
	c[N] *= 0.5;
}
//---------------------------------------------------------------
void Chebyshev_Coefficients_Gauss_2D(ftype** psi, ftype** c, int ns, int nt)
{
	int j, k;
	ftype *p, *p_coeff;
	p=new ftype[ns];
	p_coeff=new ftype[ns];

	for (k = 0; k < nt; k++)
	{
		for (j = 0; j < ns; j++)
			p[j] = psi[j][k];
		//chebft_Zeros(p, ns, 0);
		Chebyshev_Coefficients_Gauss(p,p_coeff,ns-1);
		for (j = 0; j < ns; j++)
			c[j][k] = p[j];
	}
	for (j = 0; j < ns; j++)
	{
		//chebft_Zeros(u[j], nt, 0);
		Chebyshev_Coefficients_Gauss(psi[j],c[j],nt-1);
	}
	delete[] p;
	delete[] p_coeff;
}

void Chebyshev_Coefficients_Lobatto_2D(ftype** psi, ftype** c, int ns, int nt)
{
	int j, k;
	ftype *p, *p_coeff;
	p=new ftype[ns];
	p_coeff=new ftype[ns];

	for (k = 0; k < nt; k++)
	{
		for (j = 0; j < ns; j++)
			p[j] = psi[j][k];
		Chebyshev_Coefficients_Lobatto(p,p_coeff,ns-1);
		for (j = 0; j < ns; j++)
			c[j][k] = p[j];
	}
	for (j = 0; j < ns; j++)
	{
		Chebyshev_Coefficients_Lobatto(psi[j],c[j],nt-1);
	}
	delete[] p;
	delete[] p_coeff;
}
//---------------------------------------------------------------
// ********* Chebyshev- Collocation values from Coefficients ****
//---------------------------------------------------------------
void Chebyshev_Collocations_Radau_RHS(ftype *psi, ftype *c, int N)
{
	int j;
	ftype h = 2 * Pi/ ((ftype) (2 * N + 1)); // h: Stepsize in phi=arccos(x)
	for (j = 0; j <= N; j++)
		psi[j] = Clenshaw_Chebyshev(c, N, cos(h * j));
}
//---------------------------------------------------------------
void Chebyshev_Collocations_Radau_LHS(ftype *psi, ftype *c, int N)
{
	int j;
	ftype h = 2 * Pi/ ((ftype) (2 * N + 1)); // h: Stepsize in phi=arccos(x)
	for (j = 0; j <= N; j++)
		psi[j] = Clenshaw_Chebyshev(c, N, -cos(h * j));
}
//---------------------------------------------------------------
void Chebyshev_Collocations_Gauss(ftype *psi, ftype *c, int N)
{
	int j;
	ftype h = Pi/ ((ftype) (N + 1)), hh = 0.5 * h; // h: Stepsize in phi=arccos(x)
	for (j = 0; j <= N; j++)
		psi[j] = Clenshaw_Chebyshev(c, N, cos(h * j + hh));
}
//---------------------------------------------------------------
void Chebyshev_Collocations_Lobatto(ftype *psi, ftype *c, int N)
{
	int j;
	ftype h = Pi/ ((ftype) N); // h: Stepsize in phi=arccos(x)
	for (j = 0; j <= N; j++)
		psi[j] = Clenshaw_Chebyshev(c, N, cos(h * j));
}
//---------------------------------------------------------------
// ********* Chebyshev-Coefficients of the derivative ***********
//---------------------------------------------------------------
void Chebyshev_Coefficients_Derivative(ftype *c, ftype *dc, int N)
{
	int k;

	dc[N] = 0;
	dc[N - 1] = 2 * N * c[N];
	for (k = N - 1; k >= 1; k--)
		dc[k - 1] = 2 * k * c[k] + dc[k + 1];
}
//---------------------------------------------------------------
// ********* Chebyshev-Coefficients of the integral *************
//---------------------------------------------------------------
void Chebyshev_Coefficients_Integral(ftype *c, ftype *C, int N)
{ // undetermined value C[0] (put to 0); C is indexed 0..N+1
	int k;

	for (k = 1; k <= N; k++)
		C[k] = 0.5 * (c[k - 1] - c[k + 1]) / ((ftype) k);

	C[N + 1] = 0.5 * c[N] / ((ftype) (N + 1));
	C[0] = 0;
}
//---------------------------------------------------------------
// ********* Chebyshev- definite integrals **********************
//---------------------------------------------------------------
ftype Chebyshev_Definite_Integral_Coefficients(ftype *c, int N)
{
	int k, kmax = 2 * floor(0.5 * N);
	ftype sum = 0;

	for (k = kmax; k >= 2; k -= 2)
		sum += c[k] / ((ftype) (-1 + k * k));

	return (c[0] - 2 * sum);
}
//---------------------------------------------------------------
void Chebyshev_Integration_Vector_Gauss(int N, ftype *I)
{
	int j;
	ftype *c, h = Pi/ ((ftype) (N + 1)), hh = 0.5 * h, aux = 4 / ((ftype) (N + 1));

	c = new ftype[N];

	for (j = 0; j <= N; j++)
	{
		if (j % 2 == 0)
			c[j] = aux / ((ftype) (1 - j * j));
		else
			c[j] = 0;
	}

	for (j = 0; j <= N; j++)
		I[j] = Clenshaw_Chebyshev(c, N, cos(h * j + hh));

	delete[] c;
}
//---------------------------------------------------------------
ftype Chebyshev_Definite_Integral_Collocations(ftype *psi, ftype *I, int N)
{
	int j;
	ftype s = 0;

	for (j = 0; j <= N; j++)
		s += I[j] * psi[j];

	return (s);
}
//---------------------------------------------------------------
// ********* Chebyshev- Differentiation matrices ****************
//---------------------------------------------------------------
void Chebyshev_Differentiation_Matrices_Radau_RHS(int N, ftype **D1, ftype **D2)
{
	ftype h = 2 * Pi/ ((ftype) (2 * N + 1)), m1m = -1; // h: Stepsize in phi=arccos(x)
	int m, j;

	for (m = 0; m <= N; m++)
	{
		ftype argm = h * m, argm2 = 0.5 * argm, xm = cos(argm), // arg2m = argj/m
		root_1pxm = Sqrt_2 * cos(argm2), xmp1 = sqr(root_1pxm),// root_1pxm = Sqrt[1+xm], xmp1 = xm+1
		root_1mxm = Sqrt_2 * sin(argm2), xmm1 = -sqr(root_1mxm),// root_1mxm = Sqrt[1-xm], xmm1 = xm-1
		xmm1_2 = xmp1 * xmm1, m1j = -1;// xmm1_2    = xm^2-1
		m1m *= -1;// m1m = (-1)^m
		for (j = 0; j <= N; j++)
		{
			ftype argj = h * j, argj2 = 0.5 * argj, xj = cos(argj), // arg2j = argj/2
			root_1pxj = Sqrt_2 * cos(argj2),// root_1pxj = Sqrt[1+xj],
			root_1mxj = Sqrt_2 * sin(argj2), xjm1 = -sqr(root_1mxj),// root_1mxj = Sqrt[1-xj], xjm1 = xj-1
			xmmxj = -2 * sin(argm2 + argj2) * sin(argm2 - argj2);// xmmxj = xm-xj = cos(argm)-cos(argj)= -2.*sin(argm2+argj2)*sin(argm2-argj2)
			m1j *= -1;// m1j = (-1)^j
			if (m == 0 && j == 0)
			{
				D1[m][j] = Third * N * (N + 1);
				D2[m][j] = Third * ftype(0.2) * (N - 1) * N * (N + 1) * (N + 2);
			}
			if (m == 0 && j != 0)
			{
				D1[m][j] = -m1j * Sqrt_2 * root_1pxj / xjm1;
				D2[m][j] = m1j * 2 * Sqrt_2 * Third * root_1pxj
				* (-N * (N + 1) * xjm1 - 3) / sqr(xjm1);
			}
			if (m != 0 && j == 0)
			{
				D1[m][j] = m1m / (Sqrt_2 * xmm1 * root_1pxm);
				D2[m][j] = -m1m * (2 * xm + 1)
				/ (Sqrt_2 * sqr(xmm1) * xmp1 * root_1pxm);
			}
			if (m != 0 && j == m)
			{
				D1[m][j] = 0.5 / xmm1_2;
				D2[m][j] = Third * N * (N + 1) / xmm1_2 - xm / sqr(xmm1_2);
			}
			if (m != 0 && j != 0 && m != j)
			{
				D1[m][j] = m1j * m1m * root_1pxj / (xmmxj * root_1pxm);
				D2[m][j] = -m1j * m1m * (2 * sqr(xm) - xm + xj - 2) * root_1pxj
				/ (sqr(xmmxj) * xmm1_2 * root_1pxm);
			}
		}
	}
}
//---------------------------------------------------------------
void Chebyshev_Differentiation_Matrices_Radau_LHS(int N, ftype **D1, ftype **D2)
{
	int m, j;
	Chebyshev_Differentiation_Matrices_Radau_RHS(N, D1, D2);
	for (m = 0; m <= N; m++)
		for (j = 0; j <= N; j++)
			D1[m][j] *= -1;
}
//---------------------------------------------------------------
void Chebyshev_Differentiation_Matrices_Gauss(int N, ftype **D1, ftype **D2)
{
	ftype h = Pi/ ((ftype) (N + 1)), hh = 0.5 * h, m1m = -1.; // h: Stepsize in phi=arccos(x)
	int m, j;

	for (m = 0; m <= N; m++)
	{
		ftype argm = h * m + hh, argm2 = 0.5 * argm, xm = cos(argm), // arg2m = argj/m
		root_1pxm = Sqrt_2 * cos(argm2), xmp1 = sqr(root_1pxm),// root_1pxm = Sqrt[1+xm], xmp1 = xm+1
		root_1mxm = Sqrt_2 * sin(argm2), xmm1 = -sqr(root_1mxm),// root_1mxm = Sqrt[1-xm], xmm1 = xm-1
		xmm1_2 = xmp1 * xmm1, m1j = -1;// xmm1_2    = xm^2-1
		m1m *= -1;// m1m = (-1)^m
		for (j = 0; j <= N; j++)
		{
			ftype argj = h * j + hh, argj2 = 0.5 * argj,  // arg2j = argj/2
			root_1pxj = Sqrt_2 * cos(argj2),// root_1pxj = Sqrt[1+xj]
			root_1mxj = Sqrt_2 * sin(argj2),// root_1mxj = Sqrt[1-xj]
			xmmxj = -2 * sin(argm2 + argj2) * sin(argm2 - argj2);// xmmxj = xm-xj = cos(argm)-cos(argj)= -2.*sin(argm2+argj2)*sin(argm2-argj2)
			m1j *= -1;// m1j = (-1)^j
			if (m != j)
			{
				ftype big_root = root_1pxj * root_1mxj
				/ (root_1pxm * root_1mxm);
				D1[m][j] = m1m * m1j * big_root / xmmxj;
				D2[m][j] = -D1[m][j] * (xm / xmm1_2 + 2 / xmmxj);
			}
			if (m == j)
			{
				D1[m][j] = -0.5 * xm / xmm1_2;
				D2[m][j] = sqr(xm / xmm1_2) + Third * N * (N + 2) / xmm1_2;
			}
		}
	}
}
//---------------------------------------------------------------
void Chebyshev_Differentiation_Matrices_Lobatto(int N, ftype **D1, ftype **D2)
{
	ftype h = Pi/ ((ftype) N), m1m = -1, m1N; // h: Stepsize in phi=arccos(x)
	int m, j;

	if (N % 2 == 0)
	m1N = 1.;// m1N = (-1)^N
	else
	m1N = -1.;

	for (m = 0; m <= N; m++)
	{
		ftype argm = h * m, argm2 = 0.5 * argm, xm = cos(argm), // arg2m = argj/m
		root_1pxm = Sqrt_2 * cos(argm2), xmp1 = sqr(root_1pxm),// root_1pxm = Sqrt[1+xm], xmp1 = xm+1
		root_1mxm = Sqrt_2 * sin(argm2), xmm1 = -sqr(root_1mxm),// root_1mxm = Sqrt[1-xm], xmm1 = xm-1
		xmm1_2 = xmp1 * xmm1, m1j = -1, ka_m, ka_j;// xmm1_2    = xm^2-1
		if (m == 0 || m == N)
		ka_m = 2.;
		else
		ka_m = 1.;
		m1m *= -1;// m1m = (-1)^m
		for (j = 0; j <= N; j++)
		{
			ftype argj = h * j, argj2 = 0.5 * argj, xj = cos(argj), // arg2j = argj/2
			root_1pxj = Sqrt_2 * cos(argj2), xjp1 = sqr(root_1pxj),// root_1pxj = Sqrt[1+xj], xjp1 = xj+1
			root_1mxj = Sqrt_2 * sin(argj2), xjm1 = -sqr(root_1mxj),// root_1mxj = Sqrt[1-xj], xjm1 = xj-1
			xjm1_2 = xjp1 * xjm1,// xjm1_2    = xj^2-1
			xmmxj = -2* sin(argm2 + argj2) * sin(argm2 - argj2);// xmmxj = xm-xj = cos(argm)-cos(argj)= -2.*sin(argm2+argj2)*sin(argm2-argj2)
			if (j == 0 || j == N)
			ka_j = 2;
			else
			ka_j = 1;
			m1j *= -1;// m1j = (-1)^j
			if (m == N && j == N)
			{
				D1[m][j] = -0.5 * Third * (2 * N * N + 1);
				D2[m][j] = ftype(0.2) * Third * (N * N * N * N - 1);
			}
			if (m == 0 && j == 0)
			{
				D1[m][j] = 0.5 * Third * (2 * N * N + 1);
				D2[m][j] = ftype(0.2) * Third * (N * N * N * N - 1);
			}
			if (m > 0 && m < N && m == j)
			{
				D1[m][j] = 0.5 * xj / xjm1_2;
				D2[m][j] = Third * (N * N - 1) / xmm1_2 - 1 / sqr(xmm1_2);
			}
			if (m != j)
			{
				D1[m][j] = ka_m * m1m * m1j / (ka_j * xmmxj);
				if (m == 0)
				D2[m][j] = -2 * Third * m1j / ka_j
				* ((2 * N * N + 1) * xjm1 + 6) / sqr(xjm1);
				if (m == N)
				D2[m][j] = 2. * Third * m1j / ka_j
				* ((2 * N * N + 1) * xjp1 - 6) / sqr(xjp1) * m1N;
				if (m > 0 && m < N)
				D2[m][j] = -m1m * m1j / ka_j * (sqr(xm) + xm * xj - 2)
				/ (xmm1_2 * sqr(xmmxj));
			}
		}
	}
}
//---------------------------------------------------------------
void Chebyshev_Collocations_Derivatives(ftype *psi, ftype *d1psi, ftype *d2psi,
		ftype **D1, ftype **D2, int N)
{
	int m, j;
	ftype s1, s2, s3, s4;

	for (m = 0; m <= N; m++)
	{
		s1 = s2 = s3 = s4 = 0;
		for (j = 0; j <= m; j++)
		{
			s1 += D1[m][j] * psi[j];
			s2 += D2[m][j] * psi[j];
		}
		for (j = N; j > m; j--)
		{
			s3 += D1[m][j] * psi[j];
			s4 += D2[m][j] * psi[j];
		}
		d1psi[m] = s1 + s3;
		d2psi[m] = s2 + s4;
	}
}
//---------------------------------------------------------------
