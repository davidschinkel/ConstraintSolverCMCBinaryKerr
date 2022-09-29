#include "main.h"

#define SEPDOM 1

void GetGrid(parameters *par)
{
	int s, t, ns, Ns, nt = (*par).nt, Nt=nt-1;
	(*par).A = new ftype*[NDOM];
	for(int i=0;i<NDOM;i++)
	{
		(*par).A[i] = new ftype[(*par).ns[i]];
	}
	(*par).B = new ftype[nt];
	for(int idom=0;idom<NDOM;idom++)
	{
		ns=(*par).ns[idom];
		Ns=ns-1;
		for (s = 0; s < ns; s++)
		{
			(*par).A[idom][s] =sqr(sin(Pih * s / cast<ftype>(Ns)));
		}
	}
	for (t = 0; t < nt; t++)
	{
		(*par).B[t] = -cos(Pi * t / cast<ftype>(Nt));
	}
}

void InvertIndex(parameters par, int index)
{
	for(int idom=0;idom<NDOM;idom++)
	{
		for(int ipot=0;ipot<NPOT;ipot++)
		{
			for(int i=0;i<par.ns[idom];i++)
			{
				for(int j=0;j<par.nt;j++)
				{
					if(Index(par,idom,ipot,i,j)==index)
					{
						cout<< "index=" << index << "\t(idom,ipot,i,j)=("<<idom<<","<<ipot<<","<<i<<","<<j<<")"<<endl;
					}
				}
			}
		}
	}
}

void ChebyshevTestH(parameters par)
{
	ftype *psi, *c, *cmax_ns_dom, *cmax_nt_dom;
	char filename[NMAX];

	for(int i1=1; i1<=Hmax;i1++)
	{
		for(int idom=0; idom<NDOM; idom++)
		{
			//maximale Koeffizienten je Gebiet
			cmax_ns_dom = new ftype[par.ns[idom]];
			cmax_nt_dom = new ftype[par.nt];
			fill0_dvector(cmax_ns_dom,0,par.ns[idom]);
			fill0_dvector(cmax_nt_dom,0,par.nt);

			psi = new ftype[par.ns[idom]];
			c = new ftype[par.ns[idom]];
			for (int t = 0; t < par.nt; t++)
			{
				for (int s = 0; s < par.ns[idom]; s++)
				{
					psi[s] = par.H[indx_H(par,i1,idom,s,t)];
				}
				Chebyshev_Coefficients_Lobatto(psi, c, par.ns[idom] - 1);
				max_vector(cmax_ns_dom, cmax_ns_dom, c, par.ns[idom]);
			}
			delete[] psi;
			delete[] c;
			psi = new ftype[par.nt];
			c = new ftype[par.nt];
			for (int s = 0; s < par.ns[idom]; s++)
			{
				for (int t = 0; t < par.nt; t++)
				{
					psi[t] = par.H[indx_H(par,i1,idom,s,t)];
				}
				Chebyshev_Coefficients_Lobatto(psi, c, par.nt - 1);
				max_vector(cmax_nt_dom, cmax_nt_dom, c, par.nt);
			}
			delete[] psi;
			delete[] c;
			snprintf(filename, NMAX, "%s%d_A_dom%d", "../run/test/H", i1, idom);
			PrintCheb(cmax_ns_dom,par.ns[idom],filename);
			snprintf(filename, NMAX, "%s%d_B_dom%d", "../run/test/H", i1, idom);
			PrintCheb(cmax_nt_dom,par.nt,filename);
			delete[] cmax_ns_dom;
			delete[] cmax_nt_dom;
		}

	}
}

void ChebyshevTestv(parameters par, derivs_2D v, int step)
{
	ftype *psi, *c, *cmax_ns_dom, *cmax_nt_dom;
	char filename[NMAX];

	for(int ipot=0; ipot<NPOT; ipot++)
	{
		for(int idom=0; idom<NDOM; idom++)
		{
			//maximale Koeffizienten je Gebiet
			cmax_ns_dom = new ftype[par.ns[idom]];
			cmax_nt_dom = new ftype[par.nt];
			fill0_dvector(cmax_ns_dom,0,par.ns[idom]);
			fill0_dvector(cmax_nt_dom,0,par.nt);

			psi = new ftype[par.ns[idom]];
			c = new ftype[par.ns[idom]];
			for (int t = 0; t < par.nt; t++)
			{
				for (int s = 0; s < par.ns[idom]; s++)
				{
					psi[s] = v.d0[Index(par,idom,ipot,s,t)];
				}
				Chebyshev_Coefficients_Lobatto(psi, c, par.ns[idom] - 1);
				max_vector(cmax_ns_dom, cmax_ns_dom, c, par.ns[idom]);
			}
			delete[] psi;
			delete[] c;
			psi = new ftype[par.nt];
			c = new ftype[par.nt];

			for (int s = 0; s < par.ns[idom]; s++)
			{
				for (int t = 0; t < par.nt; t++)
				{
					psi[t] = v.d0[Index(par,idom,ipot,s,t)];
				}
				Chebyshev_Coefficients_Lobatto(psi, c, par.nt - 1);
				max_vector(cmax_nt_dom, cmax_nt_dom, c, par.nt);
			}
			delete[] psi;
			delete[] c;
			snprintf(filename, NMAX, "%s%d_A_dom%d_step%d", "../run/test/v", ipot, idom, step);
			PrintCheb(cmax_ns_dom,par.ns[idom],filename);
			snprintf(filename, NMAX, "%s%d_B_dom%d_step%d", "../run/test/v", ipot, idom, step);
			PrintCheb(cmax_nt_dom,par.nt,filename);
			delete[] cmax_ns_dom;
			delete[] cmax_nt_dom;
		}
	}
}

void max_vector(ftype *out, ftype *in1, ftype *in2, int length)
{
	int i;
	for (i = 0; i < length; i++)
	{
		if (fabs(in1[i]) > fabs(in2[i]))
		{
			out[i] = fabs(in1[i]);
		}
		else
		{
			out[i] = fabs(in2[i]);
		}
	}
}

void Get_psi(parameters *par)
{
	int idom, ns, Ns, nt = (*par).nt, Nt=nt-1;
	ftype eta1=(*par).eta1, eta2=(*par).eta2, A, B, r;

	(*par).psi.d0=new ftype[(*par).ntotal];

	//domain0
	idom=0;
	ns = (*par).ns[idom];
	Ns= ns-1;
	for(int s=0;s<ns;s++)
	{
		for(int t=0;t<nt;t++)
		{
			A = sqr(sin(Pih * s / cast<ftype>(Ns))),
			B = -cos(Pi * t / cast<ftype>(Nt)),
			r = (*par).r[t];
			(*par).psi.d0[Index(*par, idom, 0, s,t)]=(1-A)*r*cos(Pi*(B+1)/ftype(4)) + SEPDOM*A*(1-B)*sqrt(eta1)/ftype(2);
		}
	}
	//domain1
	idom=1;
	ns = (*par).ns[idom];
	Ns= ns-1;
	for(int s=0;s<ns;s++)
	{
		for(int t=0;t<nt;t++)
		{
			A = sqr(sin(Pih * s / cast<ftype>(Ns))),
			B = -cos(Pi * t / cast<ftype>(Nt)),
			r = (*par).r[t];
			(*par).psi.d0[Index(*par, idom, 0, s,t)]=(A*(-1+B)*cos(atan2(Pi,eta1)/ftype(2))*pow(pow(eta1,2)+pow(Pi,2),ftype(0.25))-A*cos(atan2(Pi,eta2)/ftype(2))*pow(pow(eta2,2)+pow(Pi,2),ftype(0.25))-A*B*cos(atan2(Pi,eta2)/ftype(2))*pow(pow(eta2,2)+pow(Pi,2),ftype(0.25))+A*cos(atan2(2*Pi,eta1-B*eta1+eta2+B*eta2)/ftype(2))*pow(2,0.5)*pow(pow(eta1-B*eta1+eta2+B*eta2,2)+4*pow(Pi,2),ftype(0.25))+cos(atan2(A*Pi,eta1)/ftype(2))*pow(pow(eta1,2)+pow(A,2)*pow(Pi,2),ftype(0.25))-B*cos(atan2(A*Pi,eta1)/ftype(2))*pow(pow(eta1,2)+pow(A,2)*pow(Pi,2),ftype(0.25))+cos(atan2(A*Pi,eta2)/ftype(2))*pow(pow(eta2,2)+pow(A,2)*pow(Pi,2),ftype(0.25))+B*cos(atan2(A*Pi,eta2)/ftype(2))*pow(pow(eta2,2)+pow(A,2)*pow(Pi,2),ftype(0.25)))/ftype(2);
		}
	}
	//Derivatives_2D(*par, (*par).psi);


	for(int idom=0; idom<NDOM; idom++)
	{
		(*par).psi.cheb[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		for(int i=0; i<(*par).ns[idom];i++)
		{
			for(int j=0;j<(*par).nt;j++)
			{
				(*par).psi.cheb[idom][i][j] = (*par).psi.d0[Index(*par,idom,0,i,j)];
			}
		}
		chebft_Extremes_2D((*par).psi.cheb[idom], (*par).ns[idom], (*par).nt);
	}

	if(1==(*par).verbose)
	{
		printf("Get_psi\n");
	}
}

void Get_kappa(parameters *par)
{
	int idom, ns, Ns, nt = (*par).nt, Nt=nt-1;
	ftype eta1=(*par).eta1, eta2=(*par).eta2, A, B, r;

	(*par).kappa.d0=new ftype[(*par).ntotal];

	//domain0
	idom=0;
	ns = (*par).ns[idom];
	Ns= ns-1;
	for(int s=0;s<ns;s++)
	{
		for(int t=0;t<nt;t++)
		{
			A = sqr(sin(Pih * s / cast<ftype>(Ns))),
			B = -cos(Pi * t / cast<ftype>(Nt)),
			r = (*par).r[t];
			(*par).kappa.d0[Index(*par, idom, 0, s,t)]=(1-A)*r*sin(Pi*(B+1)/ftype(4)) + SEPDOM*A*(1+B)*sqrt(-eta2)/ftype(2);
		}
	}
	//domain1
	idom=1;
	ns = (*par).ns[idom];
	Ns= ns-1;
	for(int s=0;s<ns;s++)
	{
		for(int t=0;t<nt;t++)
		{
			A = sqr(sin(Pih * s / cast<ftype>(Ns))),
			B = -cos(Pi * t / cast<ftype>(Nt)),
			r = (*par).r[t];
			(*par).kappa.d0[Index(*par, idom, 0, s,t)]=(A*(-1+B)*pow(pow(eta1,2)+pow(Pi,2),ftype(0.25))*sin(atan2(Pi,eta1)/ftype(2))-A*pow(pow(eta2,2)+pow(Pi,2),ftype(0.25))*sin(atan2(Pi,eta2)/ftype(2))-A*B*pow(pow(eta2,2)+pow(Pi,2),ftype(0.25))*sin(atan2(Pi,eta2)/ftype(2))+A*pow(2,0.5)*pow(pow(eta1-B*eta1+eta2+B*eta2,2)+4*pow(Pi,2),ftype(0.25))*sin(atan2(2*Pi,eta1-B*eta1+eta2+B*eta2)/ftype(2))+pow(pow(eta1,2)+pow(A,2)*pow(Pi,2),ftype(0.25))*sin(atan2(A*Pi,eta1)/ftype(2))-B*pow(pow(eta1,2)+pow(A,2)*pow(Pi,2),ftype(0.25))*sin(atan2(A*Pi,eta1)/ftype(2))+pow(pow(eta2,2)+pow(A,2)*pow(Pi,2),ftype(0.25))*sin(atan2(A*Pi,eta2)/ftype(2))+B*pow(pow(eta2,2)+pow(A,2)*pow(Pi,2),ftype(0.25))*sin(atan2(A*Pi,eta2)/ftype(2)))/ftype(2);
		}
	}
	//Derivatives_2D(*par,(*par).kappa);

	for(int idom=0; idom<NDOM; idom++)
	{
		(*par).kappa.cheb[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		for(int i=0; i<(*par).ns[idom];i++)
		{
			for(int j=0;j<(*par).nt;j++)
			{
				(*par).kappa.cheb[idom][i][j] = (*par).kappa.d0[Index(*par,idom,0,i,j)];
			}
		}
		chebft_Extremes_2D((*par).kappa.cheb[idom], (*par).ns[idom], (*par).nt);
	}

	if(1==(*par).verbose)
	{
		printf("Getkappa\n");
	}
}

void Get_rho(parameters *par)
{
	int ns, nt = (*par).nt;
	ftype a0=(*par).a0, kappa, psi;

	allocate_derivs_2D(&(*par).rho,(*par).ntotal);
	allocate_derivs_2D_3rd_derivatives(&(*par).rho,(*par).ntotal);

	for(int idom=0; idom<NDOM; idom++)
	{
		ns=(*par).ns[idom];
		for(int i=0; i<ns; i++)
		{
			for(int j=0; j<nt; j++)
			{
				kappa=(*par).kappa.d0[Index(*par,idom,0,i,j)];
				psi=(*par).psi.d0[Index(*par,idom,0,i,j)];
				(*par).rho.d0[Index(*par,idom,0,i,j)]=-(a0*pow(cos(2*kappa*psi) - cosh(pow(kappa,2) - pow(psi,2)),-1)*sin(2*kappa*psi));
			}
		}
	}
	Derivatives_2D(*par,(*par).rho);
	Derivatives_2D_3(*par, (*par).rho);

	for(int idom=0; idom<NDOM; idom++)
	{
		(*par).rho.cheb[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).rho.chebd1[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).rho.chebd2[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).rho.chebd11[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).rho.chebd12[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).rho.chebd22[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).rho.chebd111[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).rho.chebd112[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).rho.chebd122[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).rho.chebd222[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		for(int i=0; i<(*par).ns[idom];i++)
		{
			for(int j=0;j<(*par).nt;j++)
			{
				(*par).rho.cheb[idom][i][j] = (*par).rho.d0[Index(*par,idom,0,i,j)];
				(*par).rho.chebd1[idom][i][j] = (*par).rho.d1[Index(*par,idom,0,i,j)];
				(*par).rho.chebd2[idom][i][j] = (*par).rho.d2[Index(*par,idom,0,i,j)];
				(*par).rho.chebd11[idom][i][j] = (*par).rho.d11[Index(*par,idom,0,i,j)];
				(*par).rho.chebd12[idom][i][j] = (*par).rho.d12[Index(*par,idom,0,i,j)];
				(*par).rho.chebd22[idom][i][j] = (*par).rho.d22[Index(*par,idom,0,i,j)];
				(*par).rho.chebd111[idom][i][j] = (*par).rho.d111[Index(*par,idom,0,i,j)];
				(*par).rho.chebd112[idom][i][j] = (*par).rho.d112[Index(*par,idom,0,i,j)];
				(*par).rho.chebd122[idom][i][j] = (*par).rho.d122[Index(*par,idom,0,i,j)];
				(*par).rho.chebd222[idom][i][j] = (*par).rho.d222[Index(*par,idom,0,i,j)];
			}
		}
		chebft_Extremes_2D((*par).rho.cheb[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).rho.chebd1[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).rho.chebd2[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).rho.chebd11[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).rho.chebd12[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).rho.chebd22[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).rho.chebd111[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).rho.chebd112[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).rho.chebd122[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).rho.chebd222[idom], (*par).ns[idom], (*par).nt);
	}

	if(1==(*par).verbose)
	{
		printf("Get_rho\n");
	}
}

void Get_z(parameters *par)
{
	int ns, nt = (*par).nt;
	ftype a0=(*par).a0, kappa, psi;

	allocate_derivs_2D(&(*par).z,(*par).ntotal);
	allocate_derivs_2D_3rd_derivatives(&(*par).z,(*par).ntotal);

	for(int idom=0; idom<NDOM; idom++)
	{
		ns=(*par).ns[idom];
		for(int i=0; i<ns; i++)
		{
			for(int j=0; j<nt; j++)
			{
				kappa=(*par).kappa.d0[Index(*par,idom,0,i,j)];
				psi=(*par).psi.d0[Index(*par,idom,0,i,j)];
				(*par).z.d0[Index(*par,idom,0,i,j)]=a0*pow(cos(2*kappa*psi) - cosh(pow(kappa,2) - pow(psi,2)),-1)*sinh(pow(kappa,2) - pow(psi,2));
			}
		}
	}
	Derivatives_2D(*par,(*par).z);
	Derivatives_2D_3(*par, (*par).z);

	for(int idom=0; idom<NDOM; idom++)
	{
		(*par).z.cheb[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).z.chebd1[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).z.chebd2[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).z.chebd11[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).z.chebd12[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).z.chebd22[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).z.chebd111[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).z.chebd112[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).z.chebd122[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		(*par).z.chebd222[idom] = dmatrix(0,(*par).ns[idom],0,(*par).nt);
		for(int i=0; i<(*par).ns[idom];i++)
		{
			for(int j=0;j<(*par).nt;j++)
			{
				(*par).z.cheb[idom][i][j] = (*par).z.d0[Index(*par,idom,0,i,j)];
				(*par).z.chebd1[idom][i][j] = (*par).z.d1[Index(*par,idom,0,i,j)];
				(*par).z.chebd2[idom][i][j] = (*par).z.d2[Index(*par,idom,0,i,j)];
				(*par).z.chebd11[idom][i][j] = (*par).z.d11[Index(*par,idom,0,i,j)];
				(*par).z.chebd12[idom][i][j] = (*par).z.d12[Index(*par,idom,0,i,j)];
				(*par).z.chebd22[idom][i][j] = (*par).z.d22[Index(*par,idom,0,i,j)];
				(*par).z.chebd111[idom][i][j] = (*par).z.d111[Index(*par,idom,0,i,j)];
				(*par).z.chebd112[idom][i][j] = (*par).z.d112[Index(*par,idom,0,i,j)];
				(*par).z.chebd122[idom][i][j] = (*par).z.d122[Index(*par,idom,0,i,j)];
				(*par).z.chebd222[idom][i][j] = (*par).z.d222[Index(*par,idom,0,i,j)];
			}
		}
		chebft_Extremes_2D((*par).z.cheb[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).z.chebd1[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).z.chebd2[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).z.chebd11[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).z.chebd12[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).z.chebd22[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).z.chebd111[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).z.chebd112[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).z.chebd122[idom], (*par).ns[idom], (*par).nt);
		chebft_Extremes_2D((*par).z.chebd222[idom], (*par).ns[idom], (*par).nt);
	}

	if(1==(*par).verbose)
	{
		printf("Get_z\n");
	}
}

ftype Get_dphi_drho(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j)
{
	ftype dphidA=phi_A_B.d1[Index(par,idom,ipot,i,j)],
		  dphidB=phi_A_B.d2[Index(par,idom,ipot,i,j)],
		  drhodA = par.rho.d1[Index(par,idom,ipot,i,j)],
		  dzdA = par.z.d1[Index(par,idom,ipot,i,j)],
		  drhodB = par.rho.d2[Index(par,idom,ipot,i,j)],
		  dzdB = par.z.d2[Index(par,idom,ipot,i,j)];

	return -((-dphidB*dzdA + dphidA*dzdB)/(drhodB*dzdA - drhodA*dzdB));
}

ftype Get_dphi_dz(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j)
{
	ftype dphidA=phi_A_B.d1[Index(par,idom,ipot,i,j)],
		  dphidB=phi_A_B.d2[Index(par,idom,ipot,i,j)],
		  drhodA = par.rho.d1[Index(par,idom,ipot,i,j)],
		  dzdA = par.z.d1[Index(par,idom,ipot,i,j)],
		  drhodB = par.rho.d2[Index(par,idom,ipot,i,j)],
		  dzdB = par.z.d2[Index(par,idom,ipot,i,j)];

	return -((dphidB*drhodA - dphidA*drhodB)/(drhodB*dzdA - drhodA*dzdB));
}

ftype Get_d2phi_drho2(parameters par, derivs_2D phi_A_B, ftype dphidrho, ftype dphidz, int idom, int ipot, int i, int j)
{
	ftype dphidAA=phi_A_B.d11[Index(par,idom,ipot,i,j)],
		  dphidAB=phi_A_B.d12[Index(par,idom,ipot,i,j)],
		  dphidBB=phi_A_B.d22[Index(par,idom,ipot,i,j)],

		  drhodA = par.rho.d1[Index(par,idom,ipot,i,j)],
		  drhodB = par.rho.d2[Index(par,idom,ipot,i,j)],
		  drhodAA = par.rho.d11[Index(par,idom,ipot,i,j)],
		  drhodAB = par.rho.d12[Index(par,idom,ipot,i,j)],
		  drhodBB = par.rho.d22[Index(par,idom,ipot,i,j)],

		  dzdA = par.z.d1[Index(par,idom,ipot,i,j)],
		  dzdB = par.z.d2[Index(par,idom,ipot,i,j)],
		  dzdAA = par.z.d11[Index(par,idom,ipot,i,j)],
		  dzdAB = par.z.d12[Index(par,idom,ipot,i,j)],
		  dzdBB = par.z.d22[Index(par,idom,ipot,i,j)];

	return -((2*dphidAB*dzdA*dzdB - 2*dphidz*dzdA*dzdAB*dzdB - dphidBB*pow(dzdA,2) + dphidz*dzdBB*pow(dzdA,2) - dphidAA*pow(dzdB,2) + dphidz*dzdAA*pow(dzdB,2) + dphidrho*(-2*drhodAB*dzdA*dzdB + drhodBB*pow(dzdA,2) + drhodAA*pow(dzdB,2)))*pow(drhodB*dzdA - drhodA*dzdB,-2));
}

ftype Get_d2phi_drho_dz(parameters par, derivs_2D phi_A_B, ftype dphidrho, ftype dphidz, int idom, int ipot, int i, int j)
{
	ftype dphidAA=phi_A_B.d11[Index(par,idom,ipot,i,j)],
		  dphidAB=phi_A_B.d12[Index(par,idom,ipot,i,j)],
		  dphidBB=phi_A_B.d22[Index(par,idom,ipot,i,j)],

		  drhodA = par.rho.d1[Index(par,idom,ipot,i,j)],
		  drhodB = par.rho.d2[Index(par,idom,ipot,i,j)],
		  drhodAA = par.rho.d11[Index(par,idom,ipot,i,j)],
		  drhodAB = par.rho.d12[Index(par,idom,ipot,i,j)],
		  drhodBB = par.rho.d22[Index(par,idom,ipot,i,j)],

		  dzdA = par.z.d1[Index(par,idom,ipot,i,j)],
		  dzdB = par.z.d2[Index(par,idom,ipot,i,j)],
		  dzdAA = par.z.d11[Index(par,idom,ipot,i,j)],
		  dzdAB = par.z.d12[Index(par,idom,ipot,i,j)],
		  dzdBB = par.z.d22[Index(par,idom,ipot,i,j)];

	return (-(dphidBB*drhodA*dzdA) + dphidAB*drhodB*dzdA - dphidrho*drhodAB*drhodB*dzdA + dphidrho*drhodA*drhodBB*dzdA - dphidz*drhodB*dzdA*dzdAB + dphidAB*drhodA*dzdB - dphidrho*drhodA*drhodAB*dzdB - dphidAA*drhodB*dzdB + dphidrho*drhodAA*drhodB*dzdB + dphidz*drhodB*dzdAA*dzdB - dphidz*drhodA*dzdAB*dzdB + dphidz*drhodA*dzdA*dzdBB)*pow(drhodB*dzdA - drhodA*dzdB,-2);
}

ftype Get_d2phi_dz2(parameters par, derivs_2D phi_A_B, ftype dphidrho, ftype dphidz, int idom, int ipot, int i, int j)
{
	ftype dphidAA=phi_A_B.d11[Index(par,idom,ipot,i,j)],
		  dphidAB=phi_A_B.d12[Index(par,idom,ipot,i,j)],
		  dphidBB=phi_A_B.d22[Index(par,idom,ipot,i,j)],

		  drhodA = par.rho.d1[Index(par,idom,ipot,i,j)],
		  drhodB = par.rho.d2[Index(par,idom,ipot,i,j)],
		  drhodAA = par.rho.d11[Index(par,idom,ipot,i,j)],
		  drhodAB = par.rho.d12[Index(par,idom,ipot,i,j)],
		  drhodBB = par.rho.d22[Index(par,idom,ipot,i,j)],

		  dzdA = par.z.d1[Index(par,idom,ipot,i,j)],
		  dzdB = par.z.d2[Index(par,idom,ipot,i,j)],
		  dzdAA = par.z.d11[Index(par,idom,ipot,i,j)],
		  dzdAB = par.z.d12[Index(par,idom,ipot,i,j)],
		  dzdBB = par.z.d22[Index(par,idom,ipot,i,j)];

	return -((2*dphidAB*drhodA*drhodB - 2*dphidrho*drhodA*drhodAB*drhodB - 2*dphidz*drhodA*drhodB*dzdAB - dphidBB*pow(drhodA,2) + dphidrho*drhodBB*pow(drhodA,2) + dphidz*dzdBB*pow(drhodA,2) - dphidAA*pow(drhodB,2) + dphidrho*drhodAA*pow(drhodB,2) + dphidz*dzdAA*pow(drhodB,2))*pow(drhodB*dzdA - drhodA*dzdB,-2));
}

ftype Get_d3phi_drho3(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j)
{
	ftype dphidAAA=phi_A_B.d111[Index(par,idom,ipot,i,j)],
		  dphidAAB=phi_A_B.d112[Index(par,idom,ipot,i,j)],
		  dphidABB=phi_A_B.d122[Index(par,idom,ipot,i,j)],
		  dphidBBB=phi_A_B.d222[Index(par,idom,ipot,i,j)],

		  dphidrho=Get_dphi_drho(par,phi_A_B,idom,ipot,i,j),
		  dphidz=Get_dphi_dz(par,phi_A_B,idom,ipot,i,j),

		  dphidrhorho = Get_d2phi_drho2(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),
		  dphidrhoz = Get_d2phi_drho_dz(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),
		  dphidzz = Get_d2phi_dz2(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),

		  drhodA = par.rho.d1[Index(par,idom,ipot,i,j)],
		  drhodB = par.rho.d2[Index(par,idom,ipot,i,j)],
		  drhodAA = par.rho.d11[Index(par,idom,ipot,i,j)],
		  drhodAB = par.rho.d12[Index(par,idom,ipot,i,j)],
		  drhodBB = par.rho.d22[Index(par,idom,ipot,i,j)],
		  drhodAAA = par.rho.d111[Index(par,idom,ipot,i,j)],
		  drhodAAB = par.rho.d112[Index(par,idom,ipot,i,j)],
		  drhodABB = par.rho.d122[Index(par,idom,ipot,i,j)],
		  drhodBBB = par.rho.d222[Index(par,idom,ipot,i,j)],

		  dzdA = par.z.d1[Index(par,idom,ipot,i,j)],
		  dzdB = par.z.d2[Index(par,idom,ipot,i,j)],
		  dzdAA = par.z.d11[Index(par,idom,ipot,i,j)],
		  dzdAB = par.z.d12[Index(par,idom,ipot,i,j)],
		  dzdBB = par.z.d22[Index(par,idom,ipot,i,j)],
		  dzdAAA = par.z.d111[Index(par,idom,ipot,i,j)],
		  dzdAAB = par.z.d112[Index(par,idom,ipot,i,j)],
		  dzdABB = par.z.d122[Index(par,idom,ipot,i,j)],
		  dzdBBB = par.z.d222[Index(par,idom,ipot,i,j)];

		  return (-3*dphidABB*dzdB*pow(dzdA,2)+3*dphidrho*drhodABB*dzdB*pow(dzdA,2)+6*dphidrhoz*drhodB*dzdAB*dzdB*pow(dzdA,2)+3*dphidz*dzdABB*dzdB*pow(dzdA,2)+3*dphidrhoz*drhodA*dzdB*dzdBB*pow(dzdA,2)+dphidBBB*pow(dzdA,3)-dphidrho*drhodBBB*pow(dzdA,3)-3*dphidrhoz*drhodB*dzdBB*pow(dzdA,3)-dphidz*dzdBBB*pow(dzdA,3)+3*dphidAAB*dzdA*pow(dzdB,2)-3*dphidrho*drhodAAB*dzdA*pow(dzdB,2)-3*dphidrhoz*drhodB*dzdA*dzdAA*pow(dzdB,2)-3*dphidz*dzdA*dzdAAB*pow(dzdB,2)-6*dphidrhoz*drhodA*dzdA*dzdAB*pow(dzdB,2)-3*dphidrhorho*(drhodB*dzdA-drhodA*dzdB)*(-2*drhodAB*dzdA*dzdB+drhodBB*pow(dzdA,2)+drhodAA*pow(dzdB,2))-dphidAAA*pow(dzdB,3)+dphidrho*drhodAAA*pow(dzdB,3)+3*dphidrhoz*drhodA*dzdAA*pow(dzdB,3)+dphidz*dzdAAA*pow(dzdB,3))*pow(drhodB*dzdA-drhodA*dzdB,-3);
}

ftype Get_d3phi_drho2dz(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j)
{
	ftype dphidAAA=phi_A_B.d111[Index(par,idom,ipot,i,j)],
		  dphidAAB=phi_A_B.d112[Index(par,idom,ipot,i,j)],
		  dphidABB=phi_A_B.d122[Index(par,idom,ipot,i,j)],
		  dphidBBB=phi_A_B.d222[Index(par,idom,ipot,i,j)],

		  dphidrho=Get_dphi_drho(par,phi_A_B,idom,ipot,i,j),
		  dphidz=Get_dphi_dz(par,phi_A_B,idom,ipot,i,j),

		  dphidrhorho = Get_d2phi_drho2(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),
		  dphidrhoz = Get_d2phi_drho_dz(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),
		  dphidzz = Get_d2phi_dz2(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),

		  drhodA = par.rho.d1[Index(par,idom,ipot,i,j)],
		  drhodB = par.rho.d2[Index(par,idom,ipot,i,j)],
		  drhodAA = par.rho.d11[Index(par,idom,ipot,i,j)],
		  drhodAB = par.rho.d12[Index(par,idom,ipot,i,j)],
		  drhodBB = par.rho.d22[Index(par,idom,ipot,i,j)],
		  drhodAAA = par.rho.d111[Index(par,idom,ipot,i,j)],
		  drhodAAB = par.rho.d112[Index(par,idom,ipot,i,j)],
		  drhodABB = par.rho.d122[Index(par,idom,ipot,i,j)],
		  drhodBBB = par.rho.d222[Index(par,idom,ipot,i,j)],

		  dzdA = par.z.d1[Index(par,idom,ipot,i,j)],
		  dzdB = par.z.d2[Index(par,idom,ipot,i,j)],
		  dzdAA = par.z.d11[Index(par,idom,ipot,i,j)],
		  dzdAB = par.z.d12[Index(par,idom,ipot,i,j)],
		  dzdBB = par.z.d22[Index(par,idom,ipot,i,j)],
		  dzdAAA = par.z.d111[Index(par,idom,ipot,i,j)],
		  dzdAAB = par.z.d112[Index(par,idom,ipot,i,j)],
		  dzdABB = par.z.d122[Index(par,idom,ipot,i,j)],
		  dzdBBB = par.z.d222[Index(par,idom,ipot,i,j)];

		  return (-2*dphidrho*drhodA*drhodABB*dzdA*dzdB-2*dphidAAB*drhodB*dzdA*dzdB+2*dphidrho*drhodAAB*drhodB*dzdA*dzdB+2*dphidz*drhodB*dzdA*dzdAAB*dzdB-2*dphidz*drhodA*dzdA*dzdABB*dzdB+dphidABB*dzdA*(drhodB*dzdA+2*drhodA*dzdB)-2*dphidrhorho*drhodBB*dzdA*dzdB*pow(drhodA,2)-2*dphidrhoz*dzdA*dzdB*dzdBB*pow(drhodA,2)+2*dphidrhorho*drhodAA*dzdA*dzdB*pow(drhodB,2)+2*dphidrhoz*dzdA*dzdAA*dzdB*pow(drhodB,2)-dphidBBB*drhodA*pow(dzdA,2)-dphidrho*drhodABB*drhodB*pow(dzdA,2)+2*dphidrhorho*drhodA*drhodB*drhodBB*pow(dzdA,2)+dphidrho*drhodA*drhodBBB*pow(dzdA,2)-dphidz*drhodB*dzdABB*pow(dzdA,2)+2*dphidrhoz*drhodAB*drhodB*dzdB*pow(dzdA,2)+dphidrhoz*drhodA*drhodBB*dzdB*pow(dzdA,2)+2*dphidzz*drhodB*dzdAB*dzdB*pow(dzdA,2)+2*dphidrhoz*drhodA*drhodB*dzdBB*pow(dzdA,2)+dphidzz*drhodA*dzdB*dzdBB*pow(dzdA,2)+dphidz*drhodA*dzdBBB*pow(dzdA,2)-2*dphidrhorho*drhodAB*pow(drhodB,2)*pow(dzdA,2)-2*dphidrhoz*dzdAB*pow(drhodB,2)*pow(dzdA,2)-dphidrhoz*drhodB*drhodBB*pow(dzdA,3)-dphidzz*drhodB*dzdBB*pow(dzdA,3)-dphidAAB*drhodA*pow(dzdB,2)+dphidrho*drhodA*drhodAAB*pow(dzdB,2)+dphidAAA*drhodB*pow(dzdB,2)-2*dphidrhorho*drhodA*drhodAA*drhodB*pow(dzdB,2)-dphidrho*drhodAAA*drhodB*pow(dzdB,2)-2*dphidrhoz*drhodA*drhodAB*dzdA*pow(dzdB,2)-dphidrhoz*drhodAA*drhodB*dzdA*pow(dzdB,2)-2*dphidrhoz*drhodA*drhodB*dzdAA*pow(dzdB,2)-dphidzz*drhodB*dzdA*dzdAA*pow(dzdB,2)-dphidz*drhodB*dzdAAA*pow(dzdB,2)+dphidz*drhodA*dzdAAB*pow(dzdB,2)-2*dphidzz*drhodA*dzdA*dzdAB*pow(dzdB,2)+2*dphidrhorho*drhodAB*pow(drhodA,2)*pow(dzdB,2)+2*dphidrhoz*dzdAB*pow(drhodA,2)*pow(dzdB,2)+dphidrhoz*drhodA*drhodAA*pow(dzdB,3)+dphidzz*drhodA*dzdAA*pow(dzdB,3))*pow(drhodB*dzdA-drhodA*dzdB,-3);
}

ftype Get_d3phi_drhodz2(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j)
{
	ftype dphidAAA=phi_A_B.d111[Index(par,idom,ipot,i,j)],
		  dphidAAB=phi_A_B.d112[Index(par,idom,ipot,i,j)],
		  dphidABB=phi_A_B.d122[Index(par,idom,ipot,i,j)],
		  dphidBBB=phi_A_B.d222[Index(par,idom,ipot,i,j)],

		  dphidrho=Get_dphi_drho(par,phi_A_B,idom,ipot,i,j),
		  dphidz=Get_dphi_dz(par,phi_A_B,idom,ipot,i,j),

		  dphidrhorho = Get_d2phi_drho2(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),
		  dphidrhoz = Get_d2phi_drho_dz(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),
		  dphidzz = Get_d2phi_dz2(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),

		  drhodA = par.rho.d1[Index(par,idom,ipot,i,j)],
		  drhodB = par.rho.d2[Index(par,idom,ipot,i,j)],
		  drhodAA = par.rho.d11[Index(par,idom,ipot,i,j)],
		  drhodAB = par.rho.d12[Index(par,idom,ipot,i,j)],
		  drhodBB = par.rho.d22[Index(par,idom,ipot,i,j)],
		  drhodAAA = par.rho.d111[Index(par,idom,ipot,i,j)],
		  drhodAAB = par.rho.d112[Index(par,idom,ipot,i,j)],
		  drhodABB = par.rho.d122[Index(par,idom,ipot,i,j)],
		  drhodBBB = par.rho.d222[Index(par,idom,ipot,i,j)],

		  dzdA = par.z.d1[Index(par,idom,ipot,i,j)],
		  dzdB = par.z.d2[Index(par,idom,ipot,i,j)],
		  dzdAA = par.z.d11[Index(par,idom,ipot,i,j)],
		  dzdAB = par.z.d12[Index(par,idom,ipot,i,j)],
		  dzdBB = par.z.d22[Index(par,idom,ipot,i,j)],
		  dzdAAA = par.z.d111[Index(par,idom,ipot,i,j)],
		  dzdAAB = par.z.d112[Index(par,idom,ipot,i,j)],
		  dzdABB = par.z.d122[Index(par,idom,ipot,i,j)],
		  dzdBBB = par.z.d222[Index(par,idom,ipot,i,j)];

		  return (2*dphidrho*drhodA*drhodABB*drhodB*dzdA+2*dphidz*drhodA*drhodB*dzdA*dzdABB+2*dphidAAB*drhodA*drhodB*dzdB-2*dphidrho*drhodA*drhodAAB*drhodB*dzdB-2*dphidz*drhodA*drhodB*dzdAAB*dzdB-dphidABB*drhodA*(2*drhodB*dzdA+drhodA*dzdB)+dphidBBB*dzdA*pow(drhodA,2)-dphidrhorho*drhodB*drhodBB*dzdA*pow(drhodA,2)-dphidrho*drhodBBB*dzdA*pow(drhodA,2)+dphidrho*drhodABB*dzdB*pow(drhodA,2)-2*dphidrhorho*drhodAB*drhodB*dzdB*pow(drhodA,2)-2*dphidrhoz*drhodBB*dzdA*dzdB*pow(drhodA,2)-2*dphidrhoz*drhodB*dzdAB*dzdB*pow(drhodA,2)+dphidz*dzdABB*dzdB*pow(drhodA,2)-dphidrhoz*drhodB*dzdA*dzdBB*pow(drhodA,2)-2*dphidzz*dzdA*dzdB*dzdBB*pow(drhodA,2)-dphidz*dzdA*dzdBBB*pow(drhodA,2)+dphidrhorho*drhodBB*dzdB*pow(drhodA,3)+dphidrhoz*dzdB*dzdBB*pow(drhodA,3)+dphidAAB*dzdA*pow(drhodB,2)-dphidrho*drhodAAB*dzdA*pow(drhodB,2)+2*dphidrhorho*drhodA*drhodAB*dzdA*pow(drhodB,2)-dphidz*dzdA*dzdAAB*pow(drhodB,2)+2*dphidrhoz*drhodA*dzdA*dzdAB*pow(drhodB,2)-dphidAAA*dzdB*pow(drhodB,2)+dphidrhorho*drhodA*drhodAA*dzdB*pow(drhodB,2)+dphidrho*drhodAAA*dzdB*pow(drhodB,2)+2*dphidrhoz*drhodAA*dzdA*dzdB*pow(drhodB,2)+dphidrhoz*drhodA*dzdAA*dzdB*pow(drhodB,2)+2*dphidzz*dzdA*dzdAA*dzdB*pow(drhodB,2)+dphidz*dzdAAA*dzdB*pow(drhodB,2)-dphidrhorho*drhodAA*dzdA*pow(drhodB,3)-dphidrhoz*dzdA*dzdAA*pow(drhodB,3)+2*dphidrhoz*drhodA*drhodB*drhodBB*pow(dzdA,2)+2*dphidzz*drhodA*drhodB*dzdBB*pow(dzdA,2)-2*dphidrhoz*drhodAB*pow(drhodB,2)*pow(dzdA,2)-2*dphidzz*dzdAB*pow(drhodB,2)*pow(dzdA,2)-2*dphidrhoz*drhodA*drhodAA*drhodB*pow(dzdB,2)-2*dphidzz*drhodA*drhodB*dzdAA*pow(dzdB,2)+2*dphidrhoz*drhodAB*pow(drhodA,2)*pow(dzdB,2)+2*dphidzz*dzdAB*pow(drhodA,2)*pow(dzdB,2))*pow(drhodB*dzdA-drhodA*dzdB,-3);
}

ftype Get_d3phi_dz3(parameters par, derivs_2D phi_A_B, int idom, int ipot, int i, int j)
{
	ftype dphidAAA=phi_A_B.d111[Index(par,idom,ipot,i,j)],
		  dphidAAB=phi_A_B.d112[Index(par,idom,ipot,i,j)],
		  dphidABB=phi_A_B.d122[Index(par,idom,ipot,i,j)],
		  dphidBBB=phi_A_B.d222[Index(par,idom,ipot,i,j)],

		  dphidrho=Get_dphi_drho(par,phi_A_B,idom,ipot,i,j),
		  dphidz=Get_dphi_dz(par,phi_A_B,idom,ipot,i,j),

		  dphidrhorho = Get_d2phi_drho2(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),
		  dphidrhoz = Get_d2phi_drho_dz(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),
		  dphidzz = Get_d2phi_dz2(par,phi_A_B,dphidrho,dphidz,idom,idom,i,j),

		  drhodA = par.rho.d1[Index(par,idom,ipot,i,j)],
		  drhodB = par.rho.d2[Index(par,idom,ipot,i,j)],
		  drhodAA = par.rho.d11[Index(par,idom,ipot,i,j)],
		  drhodAB = par.rho.d12[Index(par,idom,ipot,i,j)],
		  drhodBB = par.rho.d22[Index(par,idom,ipot,i,j)],
		  drhodAAA = par.rho.d111[Index(par,idom,ipot,i,j)],
		  drhodAAB = par.rho.d112[Index(par,idom,ipot,i,j)],
		  drhodABB = par.rho.d122[Index(par,idom,ipot,i,j)],
		  drhodBBB = par.rho.d222[Index(par,idom,ipot,i,j)],

		  dzdA = par.z.d1[Index(par,idom,ipot,i,j)],
		  dzdB = par.z.d2[Index(par,idom,ipot,i,j)],
		  dzdAA = par.z.d11[Index(par,idom,ipot,i,j)],
		  dzdAB = par.z.d12[Index(par,idom,ipot,i,j)],
		  dzdBB = par.z.d22[Index(par,idom,ipot,i,j)],
		  dzdAAA = par.z.d111[Index(par,idom,ipot,i,j)],
		  dzdAAB = par.z.d112[Index(par,idom,ipot,i,j)],
		  dzdABB = par.z.d122[Index(par,idom,ipot,i,j)],
		  dzdBBB = par.z.d222[Index(par,idom,ipot,i,j)];

		  return (3*dphidABB*drhodB*pow(drhodA,2)-3*dphidrho*drhodABB*drhodB*pow(drhodA,2)-3*dphidrhoz*drhodB*drhodBB*dzdA*pow(drhodA,2)-3*dphidz*drhodB*dzdABB*pow(drhodA,2)-6*dphidrhoz*drhodAB*drhodB*dzdB*pow(drhodA,2)-6*dphidzz*drhodB*dzdAB*dzdB*pow(drhodA,2)-3*dphidzz*drhodB*dzdA*dzdBB*pow(drhodA,2)-dphidBBB*pow(drhodA,3)+dphidrho*drhodBBB*pow(drhodA,3)+3*dphidrhoz*drhodBB*dzdB*pow(drhodA,3)+3*dphidzz*dzdB*dzdBB*pow(drhodA,3)+dphidz*dzdBBB*pow(drhodA,3)-3*dphidAAB*drhodA*pow(drhodB,2)+3*dphidrho*drhodA*drhodAAB*pow(drhodB,2)+6*dphidrhoz*drhodA*drhodAB*dzdA*pow(drhodB,2)+3*dphidz*drhodA*dzdAAB*pow(drhodB,2)+6*dphidzz*drhodA*dzdA*dzdAB*pow(drhodB,2)+3*dphidrhoz*drhodA*drhodAA*dzdB*pow(drhodB,2)+3*dphidzz*drhodA*dzdAA*dzdB*pow(drhodB,2)+dphidAAA*pow(drhodB,3)-dphidrho*drhodAAA*pow(drhodB,3)-3*dphidrhoz*drhodAA*dzdA*pow(drhodB,3)-3*dphidzz*dzdA*dzdAA*pow(drhodB,3)-dphidz*dzdAAA*pow(drhodB,3))*pow(drhodB*dzdA-drhodA*dzdB,-3);
}


void CreateFFTWplansLobatto(parameters *par)
{
	//For a Lobatto-Grid the to and back transformation is FFTW_REDFT00
	//For a Gauss-Grid the transformation to Coefficients is FFTW_REDFT10 and to Points is FFTW_REDFT01
	int n;
	double *in, *out;

	for(int idom=0; idom < NDOM; idom++)
	{
		//direction 1
		n=(*par).ns[idom];
	    in =  (double *) fftw_malloc(sizeof(ftype)*n);
	    out = (double *) fftw_malloc(sizeof(ftype)*n);
		(*par).FFTWplanLobatto[idom][0] = fftw_plan_r2r_1d(n, in, out, FFTW_REDFT00, FFTW_PATIENT);
	    fftw_free(in);
	    fftw_free(out);
		//direction 2
		n=(*par).nt;
	    in =  (double *) fftw_malloc(sizeof(ftype)*n);
	    out = (double *) fftw_malloc(sizeof(ftype)*n);
		(*par).FFTWplanLobatto[idom][1] = fftw_plan_r2r_1d(n, in, out, FFTW_REDFT00, FFTW_PATIENT);
	    fftw_free(in);
	    fftw_free(out);
	}
	//TODO die Pläne für die inverse Transformation einführen (obwohl sie identisch sind) um Flüchtigkeits-
	//fehler zu vermeiden
}

void DestroyFFTWplansLobatto(parameters *par)
{
	for(int idom=0; idom < NDOM; idom++)
	{
		//direction 1
		fftw_destroy_plan((*par).FFTWplanLobatto[idom][0]);
		//direction 2
		fftw_destroy_plan((*par).FFTWplanLobatto[idom][1]);
	}
	fftw_cleanup();
}

ftype vector_transition(parameters par, int i, int j, int component)
{
	ftype psi = par.psi.d0[Index(par,0,0,i,j)],
		  eta1= par.eta1, eta2 = par.eta2, a0=par.a0;

	if(1 == component)
	{
		return (pow(a0,-1)*pow(-2*eta2*psi*pow(eta1,0.5) + eta1*(eta2 - pow(psi,2)) + eta2*pow(psi,2),-1)*
			     (eta1*pow(-eta2,0.5) - eta1*cos(2*psi*pow(eta1,-0.5)*(-psi + pow(eta1,0.5))*pow(-eta2,0.5))*
			        cosh(pow(psi,2) + eta2*pow(eta1,-1)*pow(-psi + pow(eta1,0.5),2))*pow(-eta2,0.5) +
			       (eta1*psi - eta2*psi + eta2*pow(eta1,0.5))*sin(2*psi*pow(eta1,-0.5)*(-psi + pow(eta1,0.5))*pow(-eta2,0.5))*
			        sinh(pow(psi,2) + eta2*pow(eta1,-1)*pow(-psi + pow(eta1,0.5),2))))/2.;
	}
	if(2 == component)
	{
		return (pow(a0,-1)*pow(-2*eta2*psi*pow(eta1,0.5) + eta1*(eta2 - pow(psi,2)) + eta2*pow(psi,2),-1)*
			     (cos(2*psi*pow(eta1,-0.5)*(-psi + pow(eta1,0.5))*pow(-eta2,0.5))*cosh(pow(psi,2) + eta2*pow(eta1,-1)*pow(-psi + pow(eta1,0.5),2))*
			        (eta1*psi - eta2*psi + eta2*pow(eta1,0.5)) - (eta1*psi - eta2*psi + eta2*pow(eta1,0.5))*
			        pow(cosh(pow(psi,2) + eta2*pow(eta1,-1)*pow(-psi + pow(eta1,0.5),2)),2) +
			       sinh(pow(psi,2) + eta2*pow(eta1,-1)*pow(-psi + pow(eta1,0.5),2))*
			        (eta1*pow(-eta2,0.5)*sin(2*psi*pow(eta1,-0.5)*(-psi + pow(eta1,0.5))*pow(-eta2,0.5)) +
			          (eta1*psi - eta2*psi + eta2*pow(eta1,0.5))*sinh(pow(psi,2) + eta2*pow(eta1,-1)*pow(-psi + pow(eta1,0.5),2)))))/2.;
	}
	if(3 == component)
	{
		return 0;
	}

	return NAN;
}

void ChangeResolution_X(parameters parI, parameters parG, ftype *XI, ftype *XG)
{
	int Ns, Nt;
	ftype A,B;
	derivs_2D vI, vG;

	allocate_derivs_2D(&vI, parI.ntotal);
	allocate_derivs_2D(&vG, parG.ntotal);

	copy_dvector(vI.d0,XI,0,parI.ntotal);
	Get_Chebyshev_Coefficients_v(parI,&vI,0);

	for(int idom = 0; idom<NDOM; idom++)
	{
		Ns = parG.ns[idom] - 1;
		Nt = parG.nt - 1;
		for(int i=0;i<parG.ns[idom]; i++)
		{
			for(int j=0; j<parG.nt; j++)
			{
				A = sqr(sin(Pih * i / cast<ftype>(Ns)));
				B = -cos(Pi * j / cast<ftype>(Nt));

				XG[Index(parG,idom,0,i,j)] = chebevxy(0,1,-1,1,vI.cheb[idom],parI.ns[idom],parI.nt,A,B);
			}
		}
	}

	free_derivs_2D(&vI);
	free_derivs_2D(&vG);

	for(int idom=0; idom<NDOM; idom++)
	{
		free_dmatrix(vI.cheb[idom],0,parI.ns[idom],0,parI.nt);
	}
}

void Average_at_scri(parameters par, derivs_2D v)
{
	int idom=1;
	ftype Rp=par.Rp, phi_average, Omega_average, phi, omega, rho, z;

	phi_average=0;
	Omega_average = 0;
	for(int i=0; i<par.ns[idom]; i++)
	{
			rho = par.rho.d0[Index(par,idom,0,i,0)];
			z = par.z.d0[Index(par,idom,0,i,0)];
			omega = 1 - (rho*rho + z*z)/(Rp*Rp);
			phi = v.d0[Index(par,idom,0,i,0)];

			phi_average+=phi;
			Omega_average+=omega/(phi*phi);
	}

	phi_average = phi_average/par.ns[idom];
	Omega_average = Omega_average/par.ns[idom];

	cout << "j=0" << " \t <phi>=" << phi_average << " \t <Omega>=" << Omega_average << endl;


	phi_average=0;
	Omega_average = 0;
	for(int i=0; i<par.ns[idom]; i++)
	{
			rho = par.rho.d0[Index(par,idom,0,i,par.nt-1)];
			z = par.z.d0[Index(par,idom,0,i,par.nt-1)];
			omega = 1 - (rho*rho + z*z)/(Rp*Rp);
			phi = v.d0[Index(par,idom,0,i,par.nt-1)];

			phi_average+=phi;
			Omega_average+=omega/(phi*phi);
	}

	phi_average = phi_average/par.ns[idom];
	Omega_average = Omega_average/par.ns[idom];

	cout << "j=0" << " \t <phi>=" << phi_average << " \t <Omega>=" << Omega_average << endl;
}
