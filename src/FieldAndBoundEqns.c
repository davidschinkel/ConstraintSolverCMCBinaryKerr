#include "main.h"

void NonLinFieldEqns(parameters par, int idom, int i, int j, derivs_2D v, ftype *F)
{	// Nonlinear field equations
	int ns = par.ns[idom], Ns = ns - 1, nt = par.nt, Nt = nt - 1,
		indx0 = Index(par, idom, 0, i, j), indx1;
	ftype H[Hmax+1];

	ftype phi = v.d0[indx0],  phi2 = sqr(phi), phi5 = sqr(phi2)*phi,
		  phim7 = 1/(phi5 * phi2),
		  dphidrho = Get_dphi_drho(par,v,idom,0,i,j),
		  dphidz = Get_dphi_dz(par,v,idom,0,i,j),
		  d2phidrho2 = Get_d2phi_drho2(par,v,dphidrho,dphidz,idom,0,i,j),
		  d2phidrhodz = Get_d2phi_drho_dz(par,v,dphidrho,dphidz,idom,0,i,j),
		  d2phidz2 = Get_d2phi_dz2(par,v,dphidrho,dphidz,idom,0,i,j);

	//fordere die Gleichungen
	for (int i1 = 1; i1 <= Hmax; i1++)
	{
		H[i1] = par.H[indx_H(par, i1, idom, i, j)];
	}
	if(0==idom)
	{
		F[indx0]= H[1] * d2phidrho2 + H[2] * d2phidz2 + H[3] * d2phidrhodz + H[4] * dphidrho
				+ H[5] * dphidz	+ H[6] * phi + H[7] * phi5 + H[8]*phim7 + H[9];
	}
	if(1==idom)
	{
		F[indx0]= H[1] * d2phidrho2 + H[2] * d2phidz2 + H[3] * d2phidrhodz + H[4] * dphidrho
				+ H[5] * dphidz	+ H[6] * phi + H[7] * phi5 + H[8]*phim7 + H[9];
	}

	//Scri
	if(0==idom && 0==i)
	{
		indx1 = Index(par, idom, 0, 0, j);
		if(0==par.coefficientswitch)
		{
				F[indx0]=v.d0[indx1] - sqrt(6/(par.Rp*par.Ktr));
		}
		if(1==par.coefficientswitch)
		{
			F[indx0]=v.d0[indx1] - LaplaceSolution(par,idom , i,j);
		}
	}

	//scheinbarer Horizont
	if(1==idom && 0==j)
	{//B=-1
		int B=-1;
		indx1 = Index(par, idom, 0, i, j);
		if(0==par.coefficientswitch)
		{
//			F[indx0]= par.trapping_surface*(par.s[indx_ah1(par,1,i,B)]*dphidrho + par.s[indx_ah1(par,2,i,B)]*dphidz + par.a1h[indx_ah0(par,i,B)]*phi)
//					+ par.a2h[indx_ah0(par,i,B)]*phi2*phi + par.a3h[indx_ah0(par,i,B)]/(phi2*phi);

			ftype ktr = par.Ktr, Rp=par.Rp, rho=par.rho.d0[Index(par,idom,0,i,j)], z=par.z.d0[Index(par,idom,0,i,j)],
				  c1=par.c1, c2=par.c2,
				  d1=par.d1, d2=par.d2,
				  Sz1=par.Sz1, Sz2=par.Sz2,
				  Pz1=par.Pz1, Pz2=par.Pz2,
				  d=par.d1;

			F[indx0] =(2*ktr)/3.-4*pow(phi,-3)*pow(Rp,-2)*(phi*(d*z-pow(rho,2)-pow(z,2))-(dphidrho*rho+dphidz*(-d+z))*(-pow(rho,2)+pow(Rp,2)-pow(z,2)))*pow(-2*d*z+pow(d,2)+pow(rho,2)+pow(z,2),-0.5)+2*pow(phi,-2)*(1-pow(Rp,-2)*(pow(rho,2)+pow(z,2)))*pow(-2*d*z+pow(d,2)+pow(rho,2)+pow(z,2),-0.5)-pow(phi,-6)*(pow(-d+z,2)*pow(pow(rho,2)+pow(d-z,2),-1)*(-(c1*(-2*d1*z+pow(d1,2)+pow(rho,2)-2*pow(z,2))*pow(pow(rho,2)+pow(d1-z,2),-2.5))-c2*(-2*d2*z+pow(d2,2)+pow(rho,2)-2*pow(z,2))*pow(pow(rho,2)+pow(d2-z,2),-2.5)-(3*Pz1*z*(-2*d1*z+pow(d1,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3))/2.-(3*Pz2*z*(-2*d2*z+pow(d2,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3))/2.)+pow(rho,2)*pow(pow(rho,2)+pow(d-z,2),-1)*((pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c1*(-4*z*pow(d1,3)+pow(d1,4)-2*pow(rho,4)-pow(rho,2)*pow(z,2)+pow(d1,2)*(-pow(rho,2)+6*pow(z,2))+d1*(2*z*pow(rho,2)-4*pow(z,3))+pow(z,4))+3*Pz1*z*pow(d1-z,2)*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),0.5)))/2.+(pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c2*(-4*z*pow(d2,3)+pow(d2,4)-2*pow(rho,4)-pow(rho,2)*pow(z,2)+pow(d2,2)*(-pow(rho,2)+6*pow(z,2))+d2*(2*z*pow(rho,2)-4*pow(z,3))+pow(z,4))+3*Pz2*z*pow(d2-z,2)*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),0.5)))/2.)+2*rho*(-d+z)*pow(pow(rho,2)+pow(d-z,2),-1)*((-3*rho*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c1*z*(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2))+Pz1*(-2*d1*z+pow(d1,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),0.5)))/2.-(3*rho*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c2*z*(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2))+Pz2*(-2*d2*z+pow(d2,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),0.5)))/2.))*pow(1-pow(Rp,-2)*(pow(rho,2)+pow(z,2)),3);
		}
		if(1==par.coefficientswitch)
		{
			F[indx0]=v.d0[indx1] - LaplaceSolution(par,idom,i,j);
		}
	}
	if(1==idom && Nt==j)
	{//B=1
		int B=1;
		indx1 = Index(par, idom, 0, i, j);
		if(0==par.coefficientswitch)
		{
//			F[indx0]= par.trapping_surface*(par.s[indx_ah1(par,1,i,B)]*dphidrho + par.s[indx_ah1(par,2,i,B)]*dphidz + par.a1h[indx_ah0(par,i,B)]*phi)
//					+ par.a2h[indx_ah0(par,i,B)]*phi2*phi + par.a3h[indx_ah0(par,i,B)]/(phi2*phi);
			ftype ktr = par.Ktr, Rp=par.Rp, rho=par.rho.d0[Index(par,idom,0,i,j)], z=par.z.d0[Index(par,idom,0,i,j)],
				  c1=par.c1, c2=par.c2,
				  d1=par.d1, d2=par.d2,
				  Sz1=par.Sz1, Sz2=par.Sz2,
				  Pz1=par.Pz1, Pz2=par.Pz2,
				  d=par.d2;

			F[indx0] = (2*ktr)/3.-4*pow(phi,-3)*pow(Rp,-2)*(phi*(d*z-pow(rho,2)-pow(z,2))-(dphidrho*rho+dphidz*(-d+z))*(-pow(rho,2)+pow(Rp,2)-pow(z,2)))*pow(-2*d*z+pow(d,2)+pow(rho,2)+pow(z,2),-0.5)+2*pow(phi,-2)*(1-pow(Rp,-2)*(pow(rho,2)+pow(z,2)))*pow(-2*d*z+pow(d,2)+pow(rho,2)+pow(z,2),-0.5)-pow(phi,-6)*(pow(-d+z,2)*pow(pow(rho,2)+pow(d-z,2),-1)*(-(c1*(-2*d1*z+pow(d1,2)+pow(rho,2)-2*pow(z,2))*pow(pow(rho,2)+pow(d1-z,2),-2.5))-c2*(-2*d2*z+pow(d2,2)+pow(rho,2)-2*pow(z,2))*pow(pow(rho,2)+pow(d2-z,2),-2.5)-(3*Pz1*z*(-2*d1*z+pow(d1,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3))/2.-(3*Pz2*z*(-2*d2*z+pow(d2,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3))/2.)+pow(rho,2)*pow(pow(rho,2)+pow(d-z,2),-1)*((pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c1*(-4*z*pow(d1,3)+pow(d1,4)-2*pow(rho,4)-pow(rho,2)*pow(z,2)+pow(d1,2)*(-pow(rho,2)+6*pow(z,2))+d1*(2*z*pow(rho,2)-4*pow(z,3))+pow(z,4))+3*Pz1*z*pow(d1-z,2)*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),0.5)))/2.+(pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c2*(-4*z*pow(d2,3)+pow(d2,4)-2*pow(rho,4)-pow(rho,2)*pow(z,2)+pow(d2,2)*(-pow(rho,2)+6*pow(z,2))+d2*(2*z*pow(rho,2)-4*pow(z,3))+pow(z,4))+3*Pz2*z*pow(d2-z,2)*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),0.5)))/2.)+2*rho*(-d+z)*pow(pow(rho,2)+pow(d-z,2),-1)*((-3*rho*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c1*z*(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2))+Pz1*(-2*d1*z+pow(d1,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),0.5)))/2.-(3*rho*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c2*z*(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2))+Pz2*(-2*d2*z+pow(d2,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),0.5)))/2.))*pow(1-pow(Rp,-2)*(pow(rho,2)+pow(z,2)),3);
		}
		if(1==par.coefficientswitch)
		{
			F[indx0]=v.d0[indx1] - LaplaceSolution(par,idom,i,j);
		}
	}

	//Übergangsbedingungen
	if(0==idom && Ns==i)
	{
		//Fordere die Identität des Potentials
		F[indx0]= v.d0[Index(par,0,0,i,j)] - v.d0[Index(par,1,0,0,j)];
	}
	if(1==idom && 0==i)
	{
		//Fordere die Identität der Ableitung des Potentials
		int Ns_dom0=par.ns[0]-1;
		ftype dphidrho_dom0 = Get_dphi_drho(par,v,0,0,Ns_dom0,j),
			  dphidz_dom0 = Get_dphi_dz(par,v,0,0,Ns_dom0,j),
			  dphidrho_dom1 = Get_dphi_drho(par,v,1,0,0,j),
			  dphidz_dom1 = Get_dphi_dz(par,v,1,0,0,j),
			  s1 = vector_transition(par,i,j,1),
			  s2 = vector_transition(par,i,j,2);
		F[indx0] = s1*(dphidrho_dom0 - dphidrho_dom1) + s2*(dphidz_dom0 - dphidz_dom1);
	}
}

void LinFieldEqns(parameters par, int idom, int i, int j, derivs_2D v, derivs_2D Dv, ftype *JDX)
{	// Linear version of 'Nonlinear field equations'

	int ns = par.ns[idom], Ns = ns - 1, nt = par.nt, Nt = nt - 1,
			indx0 = Index(par, idom, 0, i, j), indx1;
		ftype H[Hmax+1];

		ftype phi = v.d0[indx0],  phi2 = sqr(phi), phi4 = sqr(phi2), phim8 = 1/(phi4 * phi4),
			  Dphi = Dv.d0[indx0],
			  Ddphidrho = Get_dphi_drho(par,Dv,idom,0,i,j),
			  Ddphidz = Get_dphi_dz(par,Dv,idom,0,i,j),
			  Dd2phidrho2 = Get_d2phi_drho2(par,Dv,Ddphidrho,Ddphidz,idom,0,i,j),
			  Dd2phidrhodz = Get_d2phi_drho_dz(par,Dv,Ddphidrho,Ddphidz,idom,0,i,j),
			  Dd2phidz2 = Get_d2phi_dz2(par,v,Ddphidrho,Ddphidz,idom,0,i,j);

	//fordere die Gleichungen
	for (int i1 = 1; i1 <= Hmax; i1++)
	{
		H[i1] = par.H[indx_H(par, i1, idom, i, j)];
	}
	if(0==idom)
	{
		JDX[indx0]= H[1] * Dd2phidrho2 + H[2] * Dd2phidz2 + H[3] * Dd2phidrhodz + H[4] * Ddphidrho
				    + H[5] * Ddphidz + H[6] * Dphi + H[7] * 5*phi4*Dphi + H[8]*(-7*phim8*Dphi);
	}
	if(1==idom)
	{
		JDX[indx0]= H[1] * Dd2phidrho2 + H[2] * Dd2phidz2 + H[3] * Dd2phidrhodz + H[4] * Ddphidrho
					+ H[5] * Ddphidz + H[6] * Dphi + H[7] * 5*phi4*Dphi + H[8]*(-7*phim8*Dphi);
	}

	//Scri
	if(0==idom && 0==i)
	{
		indx1 = Index(par, idom, 0, 0, j);
		if(0==par.coefficientswitch)
		{
			JDX[indx0]=Dv.d0[indx0];
		}
		if(1==par.coefficientswitch)
		{
			JDX[indx0]=Dv.d0[indx0];
		}
	}

	//scheinbarer Horizont
	if(1==idom && 0==j)
	{//B=-1
		int B=-1;
		indx1 = Index(par, idom, 0, i, j);
		if(0==par.coefficientswitch)
		{
			JDX[indx0]= par.s[indx_ah1(par,1,i,B)]*Ddphidrho + par.s[indx_ah1(par,2,i,B)]*Ddphidz
		        + par.a1h[indx_ah0(par,i,B)]*Dphi + par.a2h[indx_ah0(par,i,B)]*3*phi2*Dphi + (-3*par.a3h[indx_ah0(par,i,B)]/(phi4))*Dphi;
		}
		if(1==par.coefficientswitch)
		{
			JDX[indx0]=Dv.d0[indx1];
		}
	}
	if(1==idom && Nt==j)
	{//B=1
		int B=1;
		indx1 = Index(par, idom, 0, i, j);
		if(0==par.coefficientswitch)
		{
			JDX[indx0]= par.s[indx_ah1(par,1,i,B)]*Ddphidrho + par.s[indx_ah1(par,2,i,B)]*Ddphidz
					+ par.a1h[indx_ah0(par,i,B)]*Dphi + par.a2h[indx_ah0(par,i,B)]*3*phi2*Dphi + (-3*par.a3h[indx_ah0(par,i,B)]/(phi4))*Dphi;
		}
		if(1==par.coefficientswitch)
		{
			JDX[indx0]=Dv.d0[indx0];
		}
	}

	//Übergangsbedingungen
	if(0==idom && Ns==i)
	{
		//Fordere die Identität des Potentials
		JDX[indx0]= Dv.d0[Index(par,0,0,i,j)] - Dv.d0[Index(par,1,0,0,j)];
	}
	if(1==idom && 0==i)
	{
		//Fordere die Identität der Ableitung des Potentials
		int Ns_dom0=par.ns[0]-1;
		ftype Ddphidrho_dom0 = Get_dphi_drho(par,Dv,0,0,Ns_dom0,j),
			  Ddphidz_dom0 = Get_dphi_dz(par,Dv,0,0,Ns_dom0,j),
			  Ddphidrho_dom1 = Get_dphi_drho(par,Dv,1,0,0,j),
			  Ddphidz_dom1 = Get_dphi_dz(par,Dv,1,0,0,j),
			  s1 = vector_transition(par,i,j,1),
			  s2 = vector_transition(par,i,j,2);
		JDX[indx0] = s1*(Ddphidrho_dom0 - Ddphidrho_dom1) + s2*(Ddphidz_dom0 - Ddphidz_dom1);
	}
}

