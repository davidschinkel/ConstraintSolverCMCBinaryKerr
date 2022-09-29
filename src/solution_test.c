#include "main.h"

void solution_test(parameters par, derivs_2D v, int idom, int i, int j)
{
	ftype phi = v.d0[Index(par,idom,0,i,j)],
		  rho = par.rho.d0[Index(par,idom,0,i,j)],
		  z = par.z.d0[Index(par,idom,0,i,j)],
		  A = sqr(sin(Pih * i / cast<ftype>(par.ns[idom]-1))),
		  B = -cos(Pi * j / cast<ftype>(par.nt-1)),
		  eps=0.001,
		  phi_rho_p, phi_rho_m, phi_z_p, phi_z_m,
		  phi_rho_pp, phi_rho_mm, phi_z_pp, phi_z_mm,
		  phi_rho_ppp, phi_rho_mmm, phi_z_ppp, phi_z_mmm,
		  dphidrho, dphidz, d2phidrho2, d2phidz2,
		  phi5 = sqr(sqr(phi))*phi, phim7 = 1/(phi5 * sqr(phi)),
		  H[Hmax+1], sol;

	cout << scientific;
	//derivatives with respect to rho
	Get_A_B_from_rho_z(par, idom, A, B, rho+eps,z);
	phi_rho_p = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);
	Get_A_B_from_rho_z(par, idom, A, B, rho+2*eps,z);
	phi_rho_pp = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);
	Get_A_B_from_rho_z(par, idom, A, B, rho+3*eps,z);
	phi_rho_ppp = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);
	Get_A_B_from_rho_z(par, idom, A, B, rho-eps,z);
	phi_rho_m = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);
	Get_A_B_from_rho_z(par, idom, A, B, rho-2*eps,z);
	phi_rho_mm = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);
	Get_A_B_from_rho_z(par, idom, A, B, rho-3*eps,z);
	phi_rho_mmm = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);

	dphidrho = (phi_rho_p - phi_rho_m)/(2*eps);
	//d2phidrho2 = (phi_rho_p - 2*phi + phi_rho_m)/(eps*eps);
	d2phidrho2 = (2*phi_rho_mmm-27*phi_rho_mm+270*phi_rho_m-490*phi+270*phi_rho_p-27*phi_rho_pp+2*phi_rho_ppp)/(180*eps*eps);

	//cout << "dphidrho: " << dphidrho <<"\t" << Get_dphi_drho(par,v,idom,0,i,j) << "\t" << dphidrho - Get_dphi_drho(par,v,idom,0,i,j) << endl;
	cout << "\t\tdphidrho_FD - dphidrho_PS= "<< dphidrho - Get_dphi_drho(par,v,idom,0,i,j) << endl;

	//derivatives with respect to z
	Get_A_B_from_rho_z(par, idom, A, B, rho,z+eps);
	phi_z_p = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);
	Get_A_B_from_rho_z(par, idom, A, B, rho,z+2*eps);
	phi_z_pp = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);
	Get_A_B_from_rho_z(par, idom, A, B, rho,z+3*eps);
	phi_z_ppp = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);
	Get_A_B_from_rho_z(par, idom, A, B, rho,z-eps);
	phi_z_m = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);
	Get_A_B_from_rho_z(par, idom, A, B, rho,z-2*eps);
	phi_z_mm = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);
	Get_A_B_from_rho_z(par, idom, A, B, rho,z-3*eps);
	phi_z_mmm = chebevxy(0,1,-1,1,v.cheb[idom],par.ns[idom],par.nt,A,B);

	dphidz = (phi_z_p - phi_z_m)/(2*eps);
	//d2phidz2 = (phi_z_p - 2*phi + phi_z_m)/(eps*eps);
	d2phidz2 = (2*phi_z_mmm-27*phi_z_mm+270*phi_z_m-490*phi+270*phi_z_p-27*phi_z_pp+2*phi_z_ppp)/(180*eps*eps);

	//cout <<"dphidz: " << dphidz <<"\t" << Get_dphi_dz(par,v,idom,0,i,j) << "\t" << dphidz - Get_dphi_dz(par,v,idom,0,i,j) << endl;
	cout <<"\t\tdphidz_FD - dphidz_PS=" << dphidz - Get_dphi_dz(par,v,idom,0,i,j) << endl;

	//cout << "d2phidrho2: " << d2phidrho2 << " \t" << Get_d2phi_drho2(par,v,dphidrho,dphidz,idom,0,i,j) << "\t" << d2phidrho2 - Get_d2phi_drho2(par,v,dphidrho,dphidz,idom,0,i,j) << "\t" << 1- (d2phidrho2 / Get_d2phi_drho2(par,v,dphidrho,dphidz,idom,0,i,j)) << endl;
	//cout << "d2phidz2: " <<d2phidz2 << " \t" << Get_d2phi_dz2(par,v,dphidrho,dphidz,idom,0,i,j) << "\t" << d2phidz2 - Get_d2phi_dz2(par,v,dphidrho,dphidz,idom,0,i,j) << "\t" << 1- (d2phidz2 / Get_d2phi_dz2(par,v,dphidrho,dphidz,idom,0,i,j))<< endl;

	cout << "\t\td2phidrho2_FD - d2phidrho2_PS=" << d2phidrho2 - Get_d2phi_drho2(par,v,dphidrho,dphidz,idom,0,i,j) << endl;
	cout << "\t\td2phidz2_FD - d2phidz2_PS=" << d2phidz2 - Get_d2phi_dz2(par,v,dphidrho,dphidz,idom,0,i,j)<< endl;

	for (int i1 = 1; i1 <= Hmax; i1++)
	{
		H[i1] = par.H[indx_H(par, i1, idom, i, j)];
		//cout << "H["<<i1<<"]="<<H[i1]<<endl;
	}

	sol =  H[1] * d2phidrho2 + H[2] * d2phidz2 + H[4] * dphidrho
					+ H[5] * dphidz	+ H[6] * phi + H[7] * phi5 + H[8]*phim7 + H[9];

	cout << "\t\tLösung der Differentialgleichung ist mit " << sol << " bei rho=" <<rho <<" und z="<< z << " ("<<idom<<", "<< i <<", "<< j << ") erfüllt." << endl;
	cout << fixed;

//	dphidrho = Get_dphi_drho(par,v,idom,0,i,j);
//	dphidz = Get_dphi_dz(par,v,idom,0,i,j);
//	d2phidrho2 = Get_d2phi_drho2(par,v,dphidrho,dphidz,idom,0,i,j);
//	d2phidz2 = Get_d2phi_dz2(par,v,dphidrho,dphidz,idom,0,i,j);
//
//	sol =  H[1] * d2phidrho2 + H[2] * d2phidz2 + H[4] * dphidrho
//						+ H[5] * dphidz	+ H[6] * phi + H[7] * phi5 + H[8]*phim7 + H[9];
//
//	cout << "\t\tLösung der Differentialgleichung ist mit " << sol << " bei rho=" <<rho <<" und z="<< z << " erfüllt." << endl;
}
