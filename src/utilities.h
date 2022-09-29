#include "ftype.h"

#define Pi       pi<ftype>()
#define Pih      half_pi<ftype>()
//#define Piq      M_PI_4q // Pi/4
#define Third    third<ftype>()
#define TwoThird twothirds<ftype>()

#define TINY 1.0e-20
#define BIG  1.0e+10
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define NR_END 1
#define FREE_ARG char*
#define NMAX 200

typedef struct DCOMPLEX
{
	ftype r, i;
} dcomplex;

void nrerror(const char error_text[]);
int *ivector(int nl, int nh);
ftype *dvector(int nl, int nh);
ftype **dpvector(int nl, int nh);
int **imatrix(int nrl, int nrh, int ncl, int nch);
ftype **dmatrix(int nrl, int nrh, int ncl, int nch);
ftype ***d3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh);
void free_ivector(int *v, int nl, int nh);
void free_dvector(ftype *v, int nl, int nh);
void free_dpvector(ftype **v, int nl, int nh);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(ftype **m, int nrl, int nrh, int ncl, int nch);
void free_d3tensor(ftype ***t, int nrl, int nrh, int ncl, int nch, int ndl,
		int ndh);
void fill0_dvector(ftype *X, int n0, int ntotal);
void fill0_ivector(int *X, int n0, int ntotal);
void fill0_dmatrix(ftype **X, int m0, int mtotal, int n0, int ntotal);
void fill0_imatrix(int **X, int m0, int mtotal, int n0, int ntotal);
void copy_dvector(ftype *aout, ftype *ain, int n0, int ntotal);
ftype norm1(ftype *v, int n);
ftype norm2(ftype *v, int n);
ftype normalized_norm(ftype *v, int n);
ftype scalarproduct(ftype *v, ftype *w, int n);

int minimum2(int i, int j);
int minimum3(int i, int j, int k);
int maximum2(int i, int j);
ftype dmaximum2(ftype a, ftype b);
int maximum3(int i, int j, int k);
int pow_int(int mantisse, int exponent);
ftype sinch(ftype x);
ftype Sqrt(ftype x);
ftype sqr(ftype x);

::dcomplex Cadd(::dcomplex a, ::dcomplex b);
::dcomplex Csub(::dcomplex a, ::dcomplex b);
::dcomplex Cmul(::dcomplex a, ::dcomplex b);
::dcomplex RCmul(ftype x, ::dcomplex a);
::dcomplex Cdiv(::dcomplex a, ::dcomplex b);
::dcomplex Complex(ftype re, ftype im);
::dcomplex Conjg(::dcomplex z);
ftype Cabs(::dcomplex z);

::dcomplex Csqrt(::dcomplex z);
::dcomplex Cexp(::dcomplex z);
::dcomplex Clog(::dcomplex z);
::dcomplex Csin(::dcomplex z);
::dcomplex Ccos(::dcomplex z);
::dcomplex Ctan(::dcomplex z);
::dcomplex Ccot(::dcomplex z);
::dcomplex Csinh(::dcomplex z);
::dcomplex Ccosh(::dcomplex z);
::dcomplex Ctanh(::dcomplex z);
::dcomplex Ccoth(::dcomplex z);

void chebft_Zeros(ftype u[], int n, int inv);
void chebft_Extremes(ftype u[], int n, int inv);
void chebft_Extremes_2D(ftype** u, int ns, int nt);
void chebft_Zeros_2D(ftype** u, int ns, int nt);
ftype chebevy(ftype a, ftype b, ftype **c, int k, int m, ftype y);
ftype chebevxy(ftype ax, ftype bx, ftype ay, ftype by, ftype **c, int mx,
		int my, ftype x, ftype y);
void chder(ftype a, ftype b, ftype c[], ftype cder[], int n);
ftype chebev(ftype a, ftype b, ftype c[], int m, ftype x);
void chint(ftype a, ftype b, ftype c[], ftype cint[], int n);
void fourft(ftype *u, int N, int inv);
void fourder(ftype u[], ftype du[], int N);
void fourder2(ftype u[], ftype d2u[], int N);
ftype fourev(ftype *u, int N, ftype x);

void ludcmp(ftype **a, int n, int *indx, ftype *d);
void ludcmp_1(ftype **a, int n, int *indx, ftype *d);
void lubksb(ftype **a, int n, int *indx, ftype b[]);
void lubksb_1(ftype **a, int n, int *indx, ftype b[]);
void bandec(ftype **a, int n, int m1, int m2, ftype **al, int indx[], ftype *d);
void banbks(ftype **a, int n, int m1, int m2, ftype **al, int indx[],
		ftype b[]);

