/*
 * functions.h
 *
 *  Created on: 11.08.2012
 *      Author: david
 */

#define PI       3.14159265358979323846264338328
#define NR_END 1
#define FREE_ARG char*

long double *ldvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	long double *v;

	v = (long double *) malloc(
			(size_t) ((nh - nl + 1 + NR_END) * sizeof(long double)));
	//      if (!v) nrerror("allocation failure in ldvector()");
	return v - nl + NR_END;
}

void free_ldvector(long double *v, int nl, int nh)
/* free a double vector allocated with ldvector() */
{
	free((FREE_ARG) (v + nl - NR_END));
}

void chebft_Extremes(long double u[], int n, int inv) {
	int k, j, isignum, N = n - 1;
	long double fac, sum, PioN, *c;

	c = ldvector(0, N);
	PioN = PI / N;
	if (inv == 0) {
		fac = 2.0 / N;
		isignum = 1;
		for (j = 0; j < n; j++) {
			sum = 0.5 * (u[0] + u[N] * isignum);
			for (k = 1; k < N; k++)
				sum += u[k] * cos(PioN * j * k);
			c[j] = fac * sum * isignum;
			isignum = -isignum;
		}
		c[N] = 0.5 * c[N];
	} else {
		for (j = 0; j < n; j++) {
			sum = -0.5 * u[0];
			isignum = 1;
			for (k = 0; k < n; k++) {
				sum += u[k] * cos(PioN * j * k) * isignum;
				isignum = -isignum;
			}
			c[j] = sum;
		}
	}
	for (j = 0; j < n; j++)
		u[j] = c[j];
	free_ldvector(c, 0, N);
}

void chebft_Extremes_2D(long double** u, int ns, int nt) {
	int j, k;
	long double p[ns];

	for (k = 0; k < nt; k++) {
		for (j = 0; j < ns; j++)
			p[j] = u[j][k];
		chebft_Extremes(p, ns, 0);
		for (j = 0; j < ns; j++)
			u[j][k] = p[j];
	}
	for (j = 0; j < ns; j++)
		chebft_Extremes(u[j], nt, 0);
}

long double chebevy(long double a, long double b, long double **c, int k, int m,
		long double y) {
	int j;
	long double d = 0., dd = 0., xi = (2. * y - a - b) / (b - a), xi2 = 2. * xi,
			sv;

	for (j = m - 1; j >= 1; j--) {
		sv = d;
		d = xi2 * d - dd + c[k][j];
		dd = sv;
	}
	return xi * d - dd + 0.5 * c[k][0];
}

long double chebev_2D(long double ax, long double bx, long double ay,
		long double by, long double **c, int mx, int my, long double x,
		long double y) {
	int k;
	long double d = 0., dd = 0., xi = (2. * x - ax - bx) / (bx - ax), xi2 = 2.
			* xi, sv;

	for (k = mx - 1; k >= 1; k--) {
		sv = d;
		d = xi2 * d - dd + chebevy(ay, by, c, k, my - 1, y);
		dd = sv;
	}

	return xi * d - dd + 0.5 * chebevy(ay, by, c, 0., my - 1, y);
}

long double **ldmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	long double **m;

	/* allocate pointers to rows */
	m = (long double **) malloc(
			(size_t) ((nrow + NR_END) * sizeof(long double*)));
	//if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (long double *) malloc(
			(size_t) ((nrow * ncol + NR_END) * sizeof(long double)));
	//if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++) {
		m[i] = m[i - 1] + ncol;
	}

	/* return pointer to array of pointers to rows */
	return m;
}

void free_ldmatrix(long double **m, int nrl, int nrh, int ncl, int nch)
/* free a double matrix allocated by ldmatrix() */
{
	free((FREE_ARG) (m[nrl] + ncl - NR_END));
	free((FREE_ARG) (m + nrl - NR_END));
}
