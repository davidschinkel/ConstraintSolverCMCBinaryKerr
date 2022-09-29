#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include "spectral_utilities.h"
#include "fast_spectral_utilities.h"
#include "utilities.h"
#include "cast.h"

#define NDOM 		2	//number of domains
#define NPOT 		1	//number of potentials

#define bool int

#define Hmax 9			//number of Hamilton-Constraint coefficients
#define hmax 3			//number of apparent-horizon coefficients

//summation limits for the indices in the coefficient functions, as this is
//axial-symmetric derivatives are only summed up via [1,2] while in the whole
//index range is [1,3]
#define n1_1 3

#define n1_h1 3

#define n1_h2 2
#define n2_h2 3

typedef struct DERIVS_2D
{
	ftype *d0, *d1, *d2, *d11, *d12, *d22,
		  *d111, *d112, *d122, *d222,
		  **cheb[NDOM], **chebd1[NDOM], **chebd2[NDOM],
		  **chebd11[NDOM], **chebd12[NDOM], **chebd22[NDOM],
		  **chebd111[NDOM], **chebd112[NDOM], **chebd122[NDOM], **chebd222[NDOM];
} derivs_2D;

typedef struct DERIVS_1D
{
	ftype *d0, *d1, *d11;
} derivs_1D;

typedef struct HORIZON_CHEB
{
	ftype **h1[NDOM], **h2[NDOM], **h3[NDOM], **phi[NDOM], **phid1[NDOM], **phid2[NDOM],
		  **s1[NDOM], **s2[NDOM];
} horizon_cheb;

typedef struct JFD_COMPONENTS
{
	ftype **J, **K, **Kl;
	int *ncols_J, **cols_J, *iK, m1, m2;
} JFD_Components;

typedef struct PARAMETERS
{
	ftype eta1, eta2, c1, c2, a0, Rp, Ktr, **A, *B, rH1, rH2, d1, d2, Sz1, Sz2, Pz1, Pz2, coefficientswitch,
			*H, *s,	*a1h, *a2h, *a3h, d, MBondi, J_common, J_up, J_down,
			Newton_tol, Newton_tol_factor_steps, Newton_tol_finish,
			*r, *dr, *ddr, nu, *X_aux, h0_old, *MOTS_common, *MOTS_up, *MOTS_down, *MITS_common, *MITS_up, *MITS_down,
			Area_common, Area_down, Area_up, Mass_C_common, Mass_C_down, Mass_C_up, Mass_irr_common, Mass_irr_down, Mass_irr_up,
			proper_circumference_common, proper_circumference_up, proper_circumference_down,
			proper_distance;
	derivs_2D f, kappa, psi, omega, z, rho;
	horizon_cheb h_cheb;
	int ntotal, num_v, n_2D, ns_dat, nt_dat, ns[NDOM], nt, nst[NDOM + 1], nH, trapping_surface, horizon_finder_surface,
			*n_of_n_v, *n_v_of_n, n_Seq, idom_test, i_test, j_test,
			Newton_itmin, Newton_itmax, Newton_verb,
			Scalarlength, verbose, num_threads,
			MOTS_up_exist, MOTS_down_exist, MOTS_common_exist;
	fftw_plan FFTWplanLobatto[NDOM][2]; //[NDOM][1] : ns	[NDOM][2] : nt
	time_t timestamp;
} parameters;

// Subroutines in "ACMC_Data.c"
int main(int argc, char **argv);

// Subroutines in "FuncAndJacobian.c"
void Derivatives_2D(parameters par, derivs_2D v);
void Derivatives_1D(parameters par, derivs_1D v, ftype fac, int n);
void Derivatives_2D_3(parameters par, derivs_2D v);
void Get_v_From_X(parameters par, ftype *X, derivs_2D v);
void F_of_X(parameters par, ftype *X, derivs_2D v, ftype *F);
void J_times_DX(parameters par, ftype *X, ftype *DX, ftype *JDX);

// Subroutines in "IndexRoutines.c"
int Index_i(parameters par, int idom, int i);
int Index_j(parameters par, int idom, int j);
int Index(parameters par, int idom, int ipot, int i, int j);
void Get_Arrays_n_And_n_v_domain_0(parameters *par, int *n);
void Get_Arrays_n_And_n_v_domain_1(parameters *par, int *n);
void Get_Arrays_n_And_n_v(parameters *par);
void Get_Indices_From_n_v(parameters par, int n_v, int *idom, int *ipot, int *i,
		int *j);
void Free_Arrays_n_And_n_v(parameters *par);

// Subroutines in "FieldAndBoundEqns.c"
void NonLinFieldEqns(parameters par, int idom, int i, int j, derivs_2D v,
		ftype *F);
void LinFieldEqns(parameters par, int idom, int i, int j, derivs_2D v,
		derivs_2D Dv, ftype *JDX);

// Subroutines in "ScanAndPrint.c"
void ScanConfig(parameters *par);
void PrintToFile(parameters par, ftype *X, char const *filename);
void CreateInitialData(parameters par, ftype *X);
void TestData(parameters par, ftype *X);
void PrintCheb(ftype *c, int length, const char *name);
void Plot_v(parameters par,derivs_2D v, const char *name);
void Plot_data(parameters par,derivs_2D v, const char *name);
void ScanFile(parameters *par, const char *name);
void PrintHorizon(parameters par, ftype *h, const char *name);
void Print_physical_data(parameters par, const char *name);
void PrintCoordinateLines(parameters par);

// Subroutines in "GetNewResolution.c"
void v1_To_v2(parameters par1, derivs_2D v1, parameters par2, derivs_2D v2);
void GetNewResolution(parameters par0, parameters par, ftype *X0, ftype *X);

// Subroutines in "allocate_and_free_derivs.c"
void allocate_derivs_2D(derivs_2D *v, int n);
void fill0_derivs_2D(derivs_2D v, int n);
void free_derivs_2D(derivs_2D *v);
void allocate_derivs_2D_3rd_derivatives(derivs_2D *v, int n);
void free_derivs_2D_3rd_derivatives(derivs_2D *v);

// Subroutines in "newton_LU.c"
void Jacobian(parameters par, ftype *X, ftype **DF);
void DF_of_X(parameters par, ftype *X, ftype **DF);
void newton_LU(parameters par, ftype *X);

// Subroutines in "newton_Eigen.c"
void newton_Eigen(parameters par, ftype *X);

// Subroutines in "Get_Coefficients.c"
int indx_H(parameters par,int i1,int idom,int i,int j);
int indx_ah0(parameters par, int i, int B);
int indx_ah1(parameters par, int i1, int i, int B);
void Get_Coefficients(parameters *par);
void Free_Coefficients(parameters *par);

//Subroutines in "functions.c"
void GetGrid(parameters *par);
void InvertIndex(parameters par, int index);
void ChebyshevTestH(parameters par);
void ChebyshevTestv(parameters par, derivs_2D v, int step);
void max_vector(ftype *out, ftype *in1, ftype *in2, int length);
void Get_psi(parameters *par);
void Get_kappa(parameters *par);
void Get_z(parameters *par);
void Get_rho(parameters *par);
ftype Get_dphi_drho(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j);
ftype Get_dphi_dz(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j);
ftype Get_d2phi_drho2(parameters par, derivs_2D phi_A_B, ftype dphidrho, ftype dphidz, int idom, int ipot, int i, int j);
ftype Get_d2phi_drho_dz(parameters par, derivs_2D phi_A_B, ftype dphidrho, ftype dphidz, int idom, int ipot, int i, int j);
ftype Get_d2phi_dz2(parameters par, derivs_2D phi_A_B, ftype dphidrho, ftype dphidz, int idom, int ipot, int i, int j);
ftype Get_d3phi_drho3(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j);
ftype Get_d3phi_drho2dz(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j);
ftype Get_d3phi_drhodz2(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j);
ftype Get_d3phi_dz3(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j);
void CreateFFTWplansLobatto(parameters *par);
void DestroyFFTWplansLobatto(parameters *par);
ftype vector_transition(parameters par, int i, int j, int component);
void ChangeResolution_X(parameters parI, parameters parG, ftype *XI, ftype *XG);
void Average_at_scri(parameters par, derivs_2D v);

//Subroutines in "Get_r_B.c"
void Getr(parameters *par);
ftype implicit_r(parameters par, ftype r, ftype B);
ftype implicit_r_dr(parameters par, ftype r, ftype B);

//Subroutines in "math_functions.c"
ftype csch(ftype x);
ftype csc(ftype x);
ftype sech(ftype x);
ftype cot(ftype x);
ftype ArcTan(ftype e1, ftype e2);
ftype maximum(ftype e1, ftype e2);
ftype minimum(ftype e1, ftype e2);
int sgn(ftype x);

//Laplace.c
ftype LaplaceSolution(parameters par,int idom, int i,int j);
void Check_Laplace(parameters par, ftype *X);

//horizonfinder.c
int Get_A_B_from_rho_z(parameters par, int &dom, ftype &A, ftype &B, ftype rho, ftype z);
void Get_Add(parameters par, ftype &A11, ftype &A12, ftype &A21, ftype &A22, ftype R, ftype mu);
void Get_Chebyshev_Coefficients_v(parameters par, derivs_2D *v, int level);
void free_Chebyshev_Coefficients_v(parameters par, derivs_2D *v, int level);
int FieldEquations_HorizonFinder(parameters par, ftype *F, derivs_1D hh, derivs_2D field, ftype *E);
int Check_Domain(parameters par, ftype *h);
int guess_common_horizon(parameters par, ftype &h0, derivs_2D field);
int Jacobian_HorizonFinder(parameters par, derivs_1D h, derivs_2D v, ftype **DF, ftype *E);
int newton_HorizonFinder(parameters par, derivs_1D h, derivs_2D v, ftype *E);
int HorizonFinder(parameters *par, derivs_2D v, ftype h0, ftype &h_old);
void Check_dTheta(parameters par, derivs_1D hh, derivs_2D field);
void MITS_Finder(parameters par, derivs_2D v);
int Check_horizon_normal_derivative_sign(ftype *dTheta, int n);

//BondiMass.c
ftype Get_phi3_R_mu(parameters par, derivs_2D field, ftype R, ftype mu);
void BondiMass(parameters *par, derivs_2D v);
void Horizon_Area(parameters *par, derivs_2D v, derivs_1D h);
void Local_Mass(parameters *par);

//AngularMomentum.c
void AngularMomentum(parameters *par, derivs_2D v, derivs_1D h);

//solution_test.c
void solution_test(parameters par, derivs_2D v, int idom, int i, int j);

//proper_sizes.c
//void proper_circumference(parameters *par, derivs_2D v, derivs_1D h);
void proper_distance(parameters *par, derivs_2D v);
