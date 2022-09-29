#include "ftype.h"
#include <fftw3.h>

void fast_Chebyshev_Coefficients_Lobatto(ftype *psi, ftype *c, int N, fftw_plan plan);
void fast_Chebyshev_Coefficients_Gauss(ftype *psi, ftype *c, int N, fftw_plan plan);
void fast_Chebyshev_Collocations_Lobatto(ftype *psi, ftype *c, int N, fftw_plan plan);
void fast_Chebyshev_Collocations_Gauss(ftype *psi, ftype *c, int N,fftw_plan plan);
