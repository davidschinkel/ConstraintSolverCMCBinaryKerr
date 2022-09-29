#include "ftype.h"

#define Sqrt_2 root_two<ftype>()

ftype Clenshaw_cos(ftype *alpha, int N, ftype phi);
ftype Clenshaw_sin(ftype *bpsi, int N, ftype phi);
ftype Clenshaw_Fourier(ftype *alpha, ftype *bpsi, int N,
		ftype phi);
ftype Clenshaw_Chebyshev(ftype *c, int N, ftype x);

void Chebyshev_Coefficients_Radau_RHS(ftype *psi, ftype *c, int N);
void Chebyshev_Coefficients_Radau_LHS(ftype *psi, ftype *c, int N);
void Chebyshev_Coefficients_Gauss(ftype *psi, ftype *c, int N);
void Chebyshev_Coefficients_Lobatto(ftype *psi, ftype *c, int N);
void Chebyshev_Coefficients_Gauss_2D(ftype** psi, ftype** c, int ns, int nt);
void Chebyshev_Coefficients_Lobatto_2D(ftype** psi, ftype** c, int ns, int nt);

void Chebyshev_Collocations_Radau_RHS(ftype *psi, ftype *c, int N);
void Chebyshev_Collocations_Radau_LHS(ftype *psi, ftype *c, int N);
void Chebyshev_Collocations_Gauss(ftype *psi, ftype *c, int N);
void Chebyshev_Collocations_Lobatto(ftype *psi, ftype *c, int N);

void Chebyshev_Coefficients_Derivative(ftype *c, ftype *dc, int N);
void Chebyshev_Coefficients_Integral(ftype *c, ftype *C, int N);
ftype Chebyshev_Definite_Integral_Coefficients(ftype *c, int N);

void Chebyshev_Integration_Vector_Gauss(int N, ftype *I);
ftype Chebyshev_Definite_Integral_Collocations(ftype *psi,
		ftype *I, int N);

void Chebyshev_Differentiation_Matrices_Radau_RHS(int N, ftype **D1,
		ftype **D2);
void Chebyshev_Differentiation_Matrices_Radau_LHS(int N, ftype **D1,
		ftype **D2);
void Chebyshev_Differentiation_Matrices_Gauss(int N, ftype **D1,
		ftype **D2);
void Chebyshev_Differentiation_Matrices_Lobatto(int N, ftype **D1,
		ftype **D2);

void Chebyshev_Collocations_Derivatives(ftype *psi, ftype *d1psi,
		ftype *d2psi, ftype **D1, ftype **D2, int N);

