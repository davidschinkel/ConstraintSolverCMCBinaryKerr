#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "utilities.h"


void nrerror(const char error_text[])
/* Numerical Recipes standard error handler */
{
	cout<<"Numerical Recipes run-time error...\n"<<endl;
	cout<<error_text<<endl;
	cout<<"...now exiting to system...\n"<<endl;

	exit(1);
}

int *ivector(int nl, int nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v = (int *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(int)));
	if (!v)
		nrerror("allocation failure in ivector()");
	return v - nl + NR_END;
}

ftype *dvector(int nl, int nh)
/* allocate a ftype vector with subscript range v[nl..nh] */
{
	ftype *v;

	v = (ftype *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(ftype)));
	if (!v)
		nrerror("allocation failure in dvector()");
	return v - nl + NR_END;
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	int **m;

	/* allocate pointers to rows */
	m = (int **) malloc((size_t) ((nrow + NR_END) * sizeof(int*)));
	if (!m)
		nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (int *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(int)));
	if (!m[nrl])
		nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++)
		m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

ftype **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a ftype matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	ftype **m;

	/* allocate pointers to rows */
	m = (ftype **) malloc((size_t) ((nrow + NR_END) * sizeof(ftype*)));
	if (!m)
		nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (ftype *) malloc(
			(size_t) ((nrow * ncol + NR_END) * sizeof(ftype)));
	if (!m[nrl])
		nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++)
		m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

ftype ***d3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
/* allocate a ftype 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	int i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
	ftype ***t;

	/* allocate pointers to pointers to rows */
	t = (ftype ***) malloc((size_t) ((nrow + NR_END) * sizeof(ftype**)));
	if (!t)
		nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl] = (ftype **) malloc(
			(size_t) ((nrow * ncol + NR_END) * sizeof(ftype*)));
	if (!t[nrl])
		nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl] = (ftype *) malloc(
			(size_t) ((nrow * ncol * ndep + NR_END) * sizeof(ftype)));
	if (!t[nrl][ncl])
		nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for (j = ncl + 1; j <= nch; j++)
		t[nrl][j] = t[nrl][j - 1] + ndep;
	for (i = nrl + 1; i <= nrh; i++)
	{
		t[i] = t[i - 1] + ncol;
		t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
		for (j = ncl + 1; j <= nch; j++)
			t[i][j] = t[i][j - 1] + ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_ivector(int *v, int nl, int nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v + nl - NR_END));
}

void free_dvector(ftype *v, int nl, int nh)
/* free a ftype vector allocated with dvector() */
{
	free((FREE_ARG) (v + nl - NR_END));
}

void free_dpvector(ftype **v, int nl, int nh)
/* free a ftype pointer vector allocated with dpvector() */
{
	free((FREE_ARG) (v + nl - NR_END));
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl] + ncl - NR_END));
	free((FREE_ARG) (m + nrl - NR_END));
}

void free_dmatrix(ftype **m, int nrl, int nrh, int ncl, int nch)
/* free a ftype matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl] + ncl - NR_END));
	free((FREE_ARG) (m + nrl - NR_END));
}

void free_d3tensor(ftype ***t, int nrl, int nrh, int ncl, int nch, int ndl,
		int ndh)
/* free a ftype f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
	free((FREE_ARG) (t[nrl] + ncl - NR_END));
	free((FREE_ARG) (t + nrl - NR_END));
}

void fill0_dvector(ftype *X, int n0, int n1)
{
	for (int n = n0; n < n1; n++)
		X[n] = 0.;
}

void fill0_ivector(int *X, int n0, int n1)
{
	for (int n = n0; n < n1; n++)
		X[n] = 0;
}

void fill0_dmatrix(ftype **X, int m0, int m1, int n0, int n1)
{
	for (int m = m0; m < m1; m++)
		for (int n = n0; n < n1; n++)
			X[m][n] = 0.;
}

void fill0_imatrix(int **X, int m0, int m1, int n0, int n1)
{
	for (int m = m0; m < m1; m++)
		for (int n = n0; n < n1; n++)
			X[m][n] = 0;
}

void copy_dvector(ftype *aout, ftype *ain, int n0, int n1)
{
	for (int n = n0; n < n1; n++)
		aout[n] = ain[n];
}

ftype norm1(ftype *v, int n)
{
	ftype result = -1;

	for (int i = 0; i < n; i++)
		if (fabs(v[i]) > result)
			result = fabs(v[i]);

	return result;
}

ftype norm2(ftype *v, int n)
{
	ftype result = 0;

	for (int i = 0; i < n; i++)
	{
		result += v[i] * v[i];
	}

	return sqrt(result);
}

ftype normalized_norm(ftype *v, int n)
{
	return norm2(v,n)/ftype(n);
}

ftype scalarproduct(ftype *v, ftype *w, int n)
{
	int i;
	ftype result = 0;

	for (i = 0; i < n; i++)
		result += v[i] * w[i];

	return result;
}





int minimum2(int i, int j)
{
	int result = i;
	if (j < result)
		result = j;
	return result;
}

int minimum3(int i, int j, int k)
{
	int result = i;
	if (j < result)
		result = j;
	if (k < result)
		result = k;
	return result;
}

int maximum2(int i, int j)
{
	int result = i;
	if (j > result)
		result = j;
	return result;
}

ftype dmaximum2(ftype a, ftype b)
{
	ftype result = a;
	if (b > result)
		result = b;
	return result;
}

int maximum3(int i, int j, int k)
{
	int result = i;
	if (j > result)
		result = j;
	if (k > result)
		result = k;
	return result;
}

int pow_int(int mantisse, int exponent)
{
	int result = 1;

	for (int i = 1; i <= exponent; i++)
		result *= mantisse;

	return result;
}

ftype sinch(ftype x)
{
	ftype result = 1.;

	if (fabs(x) > TINY)
		result = sinh(x) / x;

	return result;
}

ftype Sqrt(ftype x)
{
	return sqrt(fabs(x));
}

ftype sqr(ftype x)
{
	return x * x;
}

dcomplex Cadd(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r = a.r + b.r;
	c.i = a.i + b.i;
	return c;
}

dcomplex Csub(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r = a.r - b.r;
	c.i = a.i - b.i;
	return c;
}

dcomplex Cmul(dcomplex a, dcomplex b)
{
	dcomplex c;
	c.r = a.r * b.r - a.i * b.i;
	c.i = a.i * b.r + a.r * b.i;
	return c;
}

dcomplex RCmul(ftype x, dcomplex a)
{
	dcomplex c;
	c.r = x * a.r;
	c.i = x * a.i;
	return c;
}

dcomplex Cdiv(dcomplex a, dcomplex b)
{
	dcomplex c;
	ftype r, den;
	if (fabs(b.r) >= fabs(b.i))
	{
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = (a.r + r * a.i) / den;
		c.i = (a.i - r * a.r) / den;
	}
	else
	{
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = (a.r * r + a.i) / den;
		c.i = (a.i * r - a.r) / den;
	}
	return c;
}

dcomplex Complex(ftype re, ftype im)
{
	dcomplex c;
	c.r = re;
	c.i = im;
	return c;
}

dcomplex Conjg(dcomplex z)
{
	dcomplex c;
	c.r = z.r;
	c.i = -z.i;
	return c;
}

ftype Cabs(dcomplex z)
{
	ftype x, y, ans, temp;
	x = fabs(z.r);
	y = fabs(z.i);
	if (x == 0.0)
		ans = y;
	else if (y == 0.0)
		ans = x;
	else if (x > y)
	{
		temp = y / x;
		ans = x * sqrt(1.0 + temp * temp);
	}
	else
	{
		temp = x / y;
		ans = y * sqrt(1.0 + temp * temp);
	}
	return ans;
}

dcomplex Csqrt(dcomplex z)
{
	dcomplex c;
	ftype x, y, w, r;
	if ((z.r == 0.0) && (z.i == 0.0))
	{
		c.r = 0.0;
		c.i = 0.0;
		return c;
	}
	else
	{
		x = fabs(z.r);
		y = fabs(z.i);
		if (x >= y)
		{
			r = y / x;
			w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
		}
		else
		{
			r = x / y;
			w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
		}
		if (z.r >= 0.0)
		{
			c.r = w;
			c.i = z.i / (2.0 * w);
		}
		else
		{
			c.i = (z.i >= 0) ? w : -w;
			c.r = z.i / (2.0 * c.i);
		}
		return c;
	}
}

dcomplex Cexp(dcomplex z)
{
	dcomplex c;
	ftype exp_r = exp(z.r);

	c.r = exp_r * cos(z.i);
	c.i = exp_r * sin(z.i);

	return c;
}

dcomplex Clog(dcomplex z)
{
	dcomplex c;

	c.r = 0.5 * log(z.r * z.r + z.i * z.i);
	c.i = atan2(z.i, z.r);

	return c;
}

dcomplex Csin(dcomplex z)
{
	dcomplex c;

	c.r = sin(z.r) * cosh(z.i);
	c.i = cos(z.r) * sinh(z.i);

	return c;
}
dcomplex Ccos(dcomplex z)
{
	dcomplex c;

	c.r = cos(z.r) * cosh(z.i);
	c.i = -sin(z.r) * sinh(z.i);

	return c;
}

dcomplex Ctan(dcomplex z)
{
	return Cdiv(Csin(z), Ccos(z));
}

dcomplex Ccot(dcomplex z)
{
	return Cdiv(Ccos(z), Csin(z));
}

dcomplex Csinh(dcomplex z)
{
	dcomplex c;

	c.r = sinh(z.r) * cos(z.i);
	c.i = cosh(z.r) * sin(z.i);

	return c;
}
dcomplex Ccosh(dcomplex z)
{
	dcomplex c;

	c.r = cosh(z.r) * cos(z.i);
	c.i = sinh(z.r) * sin(z.i);

	return c;
}

dcomplex Ctanh(dcomplex z)
{
	return Cdiv(Csinh(z), Ccosh(z));
}

dcomplex Ccoth(dcomplex z)
{
	return Cdiv(Ccosh(z), Csinh(z));
}

void chebft_Zeros(ftype u[], int n, int inv)
{
	int isignum;
	ftype fac, sum, Pion, *c;

	c = new ftype[n + 1];
	Pion = Pi / n;
	if (inv == 0)
	{
		fac = 2.0 / n;
		isignum = 1;
		for (int j = 0; j < n; j++)
		{
			sum = 0.0;
			for (int k = 0; k < n; k++)
				sum += u[k] * cos(Pion * j * (k + 0.5));
			c[j] = fac * sum * isignum;
			isignum = -isignum;
		}
	}
	else
	{
		for (int j = 0; j < n; j++)
		{
			sum = -0.5 * u[0];
			isignum = 1;
			for (int k = 0; k < n; k++)
			{
				sum += u[k] * cos(Pion * (j + 0.5) * k) * isignum;
				isignum = -isignum;
			}
			c[j] = sum;
		}
	}
	for (int j = 0; j < n; j++)
		u[j] = c[j];
	delete[] c;
}

void chebft_Extremes(ftype u[], int n, int inv)
{
	int isignum, N = n - 1;
	ftype fac, sum, PioN, *c;

	c = new ftype[n];
	PioN = Pi / N;
	if (inv == 0)
	{
		fac = 2.0 / N;
		isignum = 1;
		for (int j = 0; j < n; j++)
		{
			sum = 0.5 * (u[0] + u[N] * isignum);
			for (int k = 1; k < N; k++)
				sum += u[k] * cos(PioN * j * k);
			c[j] = fac * sum * isignum;
			isignum = -isignum;
		}
		c[N] = 0.5 * c[N];
	}
	else
	{
		for (int j = 0; j < n; j++)
		{
			sum = -0.5 * u[0];
			isignum = 1;
			for (int k = 0; k < n; k++)
			{
				sum += u[k] * cos(PioN * j * k) * isignum;
				isignum = -isignum;
			}
			c[j] = sum;
		}
	}
	for (int j = 0; j < n; j++)
		u[j] = c[j];
	delete[] c;
}

void chebft_Extremes_2D(ftype** u, int ns, int nt)
{
	ftype p[NMAX];

	for (int k = 0; k < nt; k++)
	{
		for (int j = 0; j < ns; j++)
			p[j] = u[j][k];
		chebft_Extremes(p, ns, 0);
		for (int j = 0; j < ns; j++)
			u[j][k] = p[j];
	}
	for (int j = 0; j < ns; j++)
		chebft_Extremes(u[j], nt, 0);
}

void chebft_Zeros_2D(ftype** u, int ns, int nt)
{
	ftype p[NMAX];

	for (int k = 0; k < nt; k++)
	{
		for (int j = 0; j < ns; j++)
			p[j] = u[j][k];
		chebft_Zeros(p, ns, 0);
		for (int j = 0; j < ns; j++)
			u[j][k] = p[j];
	}
	for (int j = 0; j < ns; j++)
		chebft_Zeros(u[j], nt, 0);
}

ftype chebevy(ftype a, ftype b, ftype **c, int k, int m, ftype y)
{
	ftype d = 0., dd = 0., xi = (2. * y - a - b) / (b - a), xi2 = 2. * xi, sv;

	for (int j = m - 1; j >= 1; j--)
	{
		sv = d;
		d = xi2 * d - dd + c[k][j];
		dd = sv;
	}
	return xi * d - dd + 0.5 * c[k][0];
}

ftype chebevxy(ftype ax, ftype bx, ftype ay, ftype by, ftype **c, int mx,
		int my, ftype x, ftype y)
{
	ftype d = 0., dd = 0., xi = (2. * x - ax - bx) / (bx - ax), xi2 = 2. * xi,
			sv;

	for (int k = mx - 1; k >= 1; k--)
	{
		sv = d;
		d = xi2 * d - dd + chebevy(ay, by, c, k, my - 1, y);
		dd = sv;
	}

	return xi * d - dd + 0.5 * chebevy(ay, by, c, 0., my - 1, y);
}

void chder(ftype a, ftype b, ftype c[], ftype cder[], int n)
{
	ftype con;

	cder[n - 1] = 0.0;
	cder[n - 2] = 2 * (n - 1) * c[n - 1];
	for (int j = n - 3; j >= 0; j--)
		cder[j] = cder[j + 2] + 2 * (j + 1) * c[j + 1];
	con = 2.0 / (b - a);
	for (int j = 0; j < n; j++)
		cder[j] *= con;
}

ftype chebev(ftype a, ftype b, ftype c[], int m, ftype x)
{
	ftype d = 0.0, dd = 0.0, sv, y, y2;

	y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a));
	for (int j = m - 1; j >= 1; j--)
	{
		sv = d;
		d = y2 * d - dd + c[j];
		dd = sv;
	}
	return y * d - dd + 0.5 * c[0];
}

void chint(ftype a, ftype b, ftype c[], ftype cint[], int n)
{
	ftype sum = 0.0, fac = 1.0, con;

	con = 0.25 * (b - a);
	for (int j = 1; j <= n - 2; j++)
	{
		cint[j] = con * (c[j - 1] - c[j + 1]) / j;
		sum += fac * cint[j];
		fac = -fac;
	}
	cint[n - 1] = con * c[n - 2] / (n - 1);
	sum += fac * cint[n - 1];
	cint[0] = 2.0 * sum;
}

void fourft(ftype *u, int n, int inv)
{
	int N = n - 1, iy, M, m_max;
	ftype phi_1, phi_k, arg, fac, *a, *b;

	if (n % 2 == 0)
		M = n / 2;
	else
		M = N / 2;

	a = new ftype[M+1];
	b = new ftype[M+1];

	fac = 2. / n;
	phi_1 = Pi * fac;

	if (inv == 0)
	{
		for (int m = 0; m <= M; m++)
		{
			a[m] = b[m] = 0.;
			for (int k = 0; k < n; k++)
			{
				phi_k = phi_1 * k;
				arg = phi_k * m;
				a[m] += fac * u[k] * cos(arg);
				b[m] += fac * u[k] * sin(arg);
			}
		}

		u[0] = a[0];
		u[M] = a[M];
		for (int m = 1; m < M; m++)
		{
			u[m] = a[m];
			u[m + M] = b[m];
		}
		if (n % 2 != 0)
			u[N] = b[M];
	}
	else
	{
		a[0] = u[0];
		a[M] = u[M];
		for (int m = 1; m < M; m++)
		{
			a[m] = u[m];
			b[m] = u[M + m];
		}
		if (n % 2 != 0)
			b[M] = u[N];

		iy = 1;
		for (int k = 0; k < n; k++)
		{
			u[k] = 0.5 * a[0];
			if (n % 2 == 0)
			{
				u[k] += 0.5 * a[M] * iy;
				m_max = M - 1;
			}
			else
				m_max = M;
			phi_k = phi_1 * k;
			for (int m = 1; m <= m_max; m++)
			{
				arg = phi_k * m;
				u[k] += a[m] * cos(arg) + b[m] * sin(arg);
			}
			iy = -iy;
		}
	}
	delete[] a;
	delete[] b;
}

void fourder(ftype u[], ftype du[], int n)
{
	int N = n - 1, M, mpM, m_max;

	du[0] = 0.;

	if (n % 2 == 0)
	{
		M = n / 2;
		m_max = M - 1;
		du[M] = 0.;
	}
	else
	{
		M = N / 2;
		m_max = M;
	}

	for (int m = 1; m <= m_max; m++)
	{
		mpM = m + M;
		du[m] = u[mpM] * m;
		du[mpM] = -u[m] * m;
	}
}

void fourder2(ftype u[], ftype d2u[], int n)
{
	int N = n - 1, m2, M, mpM;

	d2u[0] = 0.;
	if (n % 2 == 0)
		M = n / 2;
	else
		M = N / 2;

	for (int m = 1; m <= M; m++)
	{
		m2 = m * m;
		d2u[m] = -u[m] * m2;
		mpM = m + M;
		if (m < M || n % 2 != 0)
			d2u[mpM] = -u[mpM] * m2;
	}
}

ftype fourev(ftype *u, int n, ftype phi)
{
	int N = n - 1, m_max, M;
	ftype arg, result;

	result = 0.5 * u[0];

	if (n % 2 == 0)
	{
		M = n / 2;
		m_max = M - 1;
		result += 0.5 * u[M] * cos(phi * M);
	}
	else
	{
		M = N / 2;
		m_max = M;
	}

	for (int m = 1; m <= m_max; m++)
	{
		arg = phi * m;
		result += u[m] * cos(arg) + u[M + m] * sin(arg);
	}
	return result;
}

void ludcmp(ftype **a, int n, int *indx, ftype *d)
{   // Version of 'ludcmp' of the numerical recipes for
// matrices a[0:n-1][0:n-1]
	int imax = 0;
	ftype big, dum, sum, temp;
	ftype *vv;

	vv = new ftype[n+1];
	*d = 1.0;
	for (int i = 0; i < n; i++)
	{
		big = 0.0;
		for (int j = 0; j < n; j++)
		{
			if ((temp = fabs(a[i][j])) > big)
			{
				big = temp;
			}
		}
		if (big == 0.0)
		{
			 cout<<"ludcmp: Row i="<<i<<" is identical 0"<<endl;
		}
		if (big == 0.0)
			//nrerror("Singular matrix in routine ludcmp");
		vv[i] = 1.0 / big;
	}
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < j; i++)
		{
			sum = a[i][j];
			for (int k = 0; k < i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (int i = j; i < n; i++)
		{
			sum = a[i][j];
			for (int k = 0; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big)
			{
				big = dum;
				imax = i;
			}
		}
		if (j != imax)
		{
			for (int k = 0; k < n; k++)
			{
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0)
			a[j][j] = TINY;
		if (j != n - 1)
		{
			dum = 1.0 / (a[j][j]);
			for (int i = j + 1; i < n; i++)
				a[i][j] *= dum;
		}
	}
	delete[] vv;
}

void ludcmp_1(ftype **a, int n, int *indx, ftype *d)
{  // Version of 'lubksb' of the numerical recipes for
// matrices a[1:n][1:n] and vectors b[1:n]
	int imax = 0;
	ftype big, dum, sum, temp;
	ftype *vv;

	vv = new ftype[n+1];
	*d = 1.0;
	for (int i = 1; i <= n; i++)
	{
		big = 0.0;
		for (int j = 1; j <= n; j++)
			if ((temp = fabs(a[i][j])) > big)
				big = temp;
		if (big == 0.0)
			cout<<"ludcmp_1: Row i="<<i<<" is identical 0"<<endl;
		if (big == 0.0)
			nrerror("Singular matrix in routine ludcmp_1");
		vv[i] = 1.0 / big;
	}
	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i < j; i++)
		{
			sum = a[i][j];
			for (int k = 1; k < i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (int i = j; i <= n; i++)
		{
			sum = a[i][j];
			for (int k = 1; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big)
			{
				big = dum;
				imax = i;
			}
		}
		if (j != imax)
		{
			for (int k = 1; k <= n; k++)
			{
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0)
			a[j][j] = TINY;
		if (j != n)
		{
			dum = 1.0 / (a[j][j]);
			for (int i = j + 1; i <= n; i++)
				a[i][j] *= dum;
		}
	}
	delete[] vv;
}

void lubksb(ftype **a, int n, int *indx, ftype b[])
{   // Version of 'lubksb' of the numerical recipes for
// matrices a[0:n-1][0:n-1] and vectors b[0:n-1]

	int ii = 0, ip;
	ftype sum;

	for (int i = 0; i < n; i++)
	{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
		{
			for (int j = ii; j <= i - 1; j++)
			{
				sum -= a[i][j] * b[j];
			}
		}
		else if (sum)
		{
			ii = i;
		}
		b[i] = sum;
		//cout << "i="<<i<<"\tsum="<<sum<<"\tb="<<b[i]<<endl;
	}
	for (int i = n - 1; i >= 0; i--)
	{
		sum = b[i];
		//cout << "i="<<i<<"\tsum="<<sum<<"\tb="<<b[i]<<endl;
		for (int j = i + 1; j < n; j++)
		{
			sum -= a[i][j] * b[j];
//			cout <<"i="<<i<<"\tb="<<b[j]<<"\tsum="<<sum<<"\ta="<<a[i][j]<<endl;
//			if(boost::math::isnan(sum))
//			{
//				exit(1);
//			}
		}
		b[i] = sum / a[i][i];
		//cout <<"i="<<i<<"\tb="<<b[i]<<"\tsum="<<sum<<"\ta="<<a[i][i]<<endl;
	}
}

void lubksb_1(ftype **a, int n, int *indx, ftype b[])
{
	int ii = 0, ip;
	ftype sum;

	for (int i = 1; i <= n; i++)
	{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
			for (int j = ii; j <= i - 1; j++)
				sum -= a[i][j] * b[j];
		else if (sum)
			ii = i;
		b[i] = sum;
	}
	for (int i = n; i >= 1; i--)
	{
		sum = b[i];
		for (int j = i + 1; j <= n; j++)
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}

void bandec(ftype **a, int n, int m1, int m2, ftype **al, int indx[], ftype *d)
{
	int i, j, k, l, mm;
	ftype temp;

	mm = m1 + m2 + 1;
	l = m1;
	for (i = 1; i <= m1; i++)
	{
		for (j = m1 + 2 - i; j <= mm; j++)
			a[i][j - l] = a[i][j];
		l--;
		for (j = mm - l; j <= mm; j++)
			a[i][j] = 0.0;
	}
	*d = 1.0;
	l = m1;
	for (k = 1; k <= n; k++)
	{
		temp = a[k][1];
		i = k;
		if (l < n)
			l++;
		for (j = k + 1; j <= l; j++)
		{
			if (fabs(a[j][1]) > fabs(temp))
			{
				temp = a[j][1];
				i = j;
			}
		}
		indx[k] = i;
		if (temp == 0.0)
			a[k][1] = TINY;
		if (i != k)
		{
			*d = -(*d);
			for (j = 1; j <= mm; j++)
				SWAP(a[k][j], a[i][j])
		}
		for (i = k + 1; i <= l; i++)
		{
			temp = a[i][1] / a[k][1];
			al[k][i - k] = temp;
			for (j = 2; j <= mm; j++)
				a[i][j - 1] = a[i][j] - temp * a[k][j];
			a[i][mm] = 0.0;
		}
	}
}

void banbks(ftype **a, int n, int m1, int m2, ftype **al, int indx[], ftype b[])
{
	int i, k, l, mm;
	ftype temp;
	mm = m1 + m2 + 1;
	l = m1;
	for (k = 1; k <= n; k++)
	{
		i = indx[k];
		if (i != k)
			SWAP(b[k], b[i])
		if (l < n)
			l++;
		for (i = k + 1; i <= l; i++)
			b[i] -= al[k][i - k] * b[k];
	}
	l = 1;
	for (i = n; i >= 1; i--)
	{
		temp = b[i];
		for (k = 2; k <= l; k++)
			temp -= a[i][k] * b[k + i - 1];
		b[i] = temp / a[i][1];
		if (l < mm)
			l++;
	}
}

