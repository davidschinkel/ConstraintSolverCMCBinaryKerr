#include "main.h"

ftype Get_phi3_R_mu(parameters par, derivs_2D field, ftype R, ftype mu)
{
	ftype A,B;
	int idom, check;

	A=0;
	B=-1;
	idom=0;

	check = Get_A_B_from_rho_z(par,idom,A,B,R*sqrt(1-mu*mu),R*mu);
	if(0!=check)
	{
		exit(1);
	}

	ftype
		dphidAAA = chebevxy(0,1,-1,1,field.chebd111[idom],par.ns[idom],par.nt,A,B),
		dphidAAB = chebevxy(0,1,-1,1,field.chebd112[idom],par.ns[idom],par.nt,A,B),
		dphidABB = chebevxy(0,1,-1,1,field.chebd122[idom],par.ns[idom],par.nt,A,B),
		dphidBBB = chebevxy(0,1,-1,1,field.chebd222[idom],par.ns[idom],par.nt,A,B),

		rho = chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,A,B),
		drhodA = chebevxy(0,1,-1,1,par.rho.chebd1[idom],par.ns[idom],par.nt,A,B),
		drhodB = chebevxy(0,1,-1,1,par.rho.chebd2[idom],par.ns[idom],par.nt,A,B),

		z = chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,A,B),
		dzdA = chebevxy(0,1,-1,1,par.z.chebd1[idom],par.ns[idom],par.nt,A,B),
		dzdB = chebevxy(0,1,-1,1,par.z.chebd2[idom],par.ns[idom],par.nt,A,B);
	
	return pow(drhodB*dzdA-drhodA*dzdB,-3)*(3*z*(2*(dphidABB*drhodA-dphidAAB*drhodB)*dzdA*dzdB+(-(dphidBBB*drhodA)+dphidABB*drhodB)*pow(dzdA,2)+(-(dphidAAB*drhodA)+dphidAAA*drhodB)*pow(dzdB,2))*pow(rho,2)+(-3*dphidABB*dzdB*pow(dzdA,2)+dphidBBB*pow(dzdA,3)+3*dphidAAB*dzdA*pow(dzdB,2)-dphidAAA*pow(dzdB,3))*pow(rho,3)-3*rho*(dzdB*(-2*dphidAAB*drhodA*drhodB+dphidABB*pow(drhodA,2)+dphidAAA*pow(drhodB,2))-dzdA*(-2*dphidABB*drhodA*drhodB+dphidBBB*pow(drhodA,2)+dphidAAB*pow(drhodB,2)))*pow(z,2)+(3*dphidABB*drhodB*pow(drhodA,2)-dphidBBB*pow(drhodA,3)-3*dphidAAB*drhodA*pow(drhodB,2)+dphidAAA*pow(drhodB,3))*pow(z,3))*pow(pow(rho,2)+pow(z,2),-1.5);
}

void BondiMass(parameters *par, derivs_2D v)
{
	ftype d1=(*par).d1, d2=(*par).d2, c1=(*par).c1, c2=(*par).c2, Pz1=(*par).Pz1, Pz2=(*par).Pz2,
		  ktr=(*par).Ktr, Rp=(*par).Rp,
		  mu, R=Rp, phi0, *phi3, *c, Integral_phi3, Integral_A;

	phi3 = new ftype[(*par).nt];
	c = new ftype[(*par).nt];

	for(int j=0; j<(*par).nt; j++)
	{
		mu = -cos(Pi * j / cast<ftype>((*par).nt-1));

		phi3[j] = Get_phi3_R_mu((*par),v,R,mu);

		//cout << "phi3[" << j << "]=" << phi3[j] << endl;
	}

	Chebyshev_Coefficients_Lobatto(phi3, c, (*par).nt-1);

	Integral_phi3 = Chebyshev_Definite_Integral_Coefficients(c, (*par).nt-1);

	cout << "\tIntegral_phi3=" << Integral_phi3 << endl;

	Integral_A = (3*Pz1*pow(d1,-2)*pow(R,-1)*(-atanh(d1*pow(R,-1))+d1*R*(pow(d1,2)+pow(R,2))*pow(pow(d1,2)-pow(R,2),-2))+3*Pz2*pow(d2,-2)*pow(R,-1)*(-atanh(d2*pow(R,-1))+d2*R*(pow(d2,2)+pow(R,2))*pow(pow(d2,2)-pow(R,2),-2))+4*c1*pow(R*pow(d1,2)-pow(R,3),-1)+4*c2*pow(R*pow(d2,2)-pow(R,3),-1))/2.;

	cout << "\tIntegral_A=" << Integral_A << endl;

	phi0 = sqrt(6/(ktr*Rp));

	//(*par).MBondi = 0.5*Rp*phi0*phi0*(1/9.)*ktr*Rp*Rp*Rp*(ktr*Integral_A - 2*phi0*Integral_phi3);
	//(*par).MBondi = Third*ktr*pow(Rp,3)*Integral_A - 2*sqrt(TwoThird)*sqrt(pow(Rp,7)/ktr)*Integral_phi3;
	(*par).MBondi = (1/6.)*pow(Rp,3)*(Integral_A*ktr - 2*sqrt(6*Rp/ktr)*Integral_phi3);

	cout << "BondiMass=" << (*par).MBondi << endl;

	cout << "\tKrümmungsanteil: " << (1/6.)*pow(Rp,3)*(Integral_A*ktr) << endl;
	cout <<"\t phi-Anteil: " << (1/6.)*pow(Rp,3)*(- 2*sqrt(6*Rp/ktr)*Integral_phi3) << endl;
//cout << fixed << endl;
}

void Horizon_Area(parameters *par, derivs_2D v, derivs_1D h)
{
	int nH = (*par).nH, check, idom;
	ftype *Integrand, *c_Integrand, mu, R, A, B, phi, hd1,
		  omega, Rp=(*par).Rp, sqrtdetq, Area, d=(*par).d, rho, z;

	Integrand = new ftype[nH];
	c_Integrand = new ftype[nH];

	A=0.5;
	B=-1;
	idom=0;

	for(int j=0; j<nH; j++)
	{
		mu = -cos(Pi * j / cast<ftype>(nH-1));
		R = h.d0[j];
		hd1 = h.d1[j];

		check = Get_A_B_from_rho_z(*par,idom,A,B,R*sqrt(1-mu*mu),R*mu + d);
		if(0!=check)
		{
			exit(1);
		}

		phi = chebevxy(0,1,-1,1,v.cheb[idom],(*par).ns[idom],(*par).nt,A,B);
		rho= chebevxy(0,1,-1,1,(*par).rho.cheb[idom],(*par).ns[idom],(*par).nt,A,B);
		z= chebevxy(0,1,-1,1,(*par).z.cheb[idom],(*par).ns[idom],(*par).nt,A,B);

		omega = 1 - (rho*rho + z*z)/(Rp*Rp);
		sqrtdetq = pow(R,3)*pow(pow(-(pow(hd1,2)*(-1 + pow(mu,2))) + pow(R,2),-1),0.5);

		Integrand[j] = (phi*phi*phi*phi/(omega*omega))*sqrtdetq;
	}

	Chebyshev_Coefficients_Lobatto(Integrand, c_Integrand, nH-1);

	Area = Chebyshev_Definite_Integral_Coefficients(c_Integrand, nH-1);


	if(0==(*par).d)
	{
		(*par).Area_common = Area;
		(*par).proper_circumference_common = sqrt(Area*Pi);
	}
	if((*par).d>0)
	{
		(*par).Area_up = Area;
		(*par).proper_circumference_up = sqrt(Area*Pi);
	}
	if((*par).d<0)
	{
		(*par).Area_down = Area;
		(*par).proper_circumference_down = sqrt(Area*Pi);
	}
	cout << "Horizon_Area A=" << Area << endl;
	cout << "proper_circumference u=" << sqrt(Area*Pi) << endl;
	cout << "proper_radius r=" << sqrt(Area/(4*Pi)) << endl;

	delete[] Integrand;
	delete[] c_Integrand;
}

void Local_Mass(parameters *par)
{
	ftype R, J;

	if(0==(*par).d)
	{
		R = sqrt((*par).Area_common/(4*Pi));
		J = (*par).J_common;
		(*par).Mass_C_common = sqrt(R*R/4. + J*J/(R*R));
		//(*par).Mass_C_common = (R + pow(R, TwoThird)*pow(pow(R,4) + 6*J*(9*J + pow(3,0.5)*pow(27*pow(J,2) + pow(R,4),0.5)),
       //-Third) + pow(pow(R,-1)*(pow(R,4) + 6*J*(9*J + pow(3,0.5)*pow(27*pow(J,2) + pow(R,4),0.5))),Third))/6.;
		(*par).Mass_irr_common = sqrt((*par).Area_common/(16*Pi));
	}
	if((*par).d>0)
	{
		R = sqrt((*par).Area_up/(4*Pi));
		J = (*par).J_up;
		(*par).Mass_C_up = sqrt(R*R/4. + J*J/(R*R));
		//(*par).Mass_C_up = (R + pow(R, TwoThird)*pow(pow(R,4) + 6*J*(9*J + pow(3,0.5)*pow(27*pow(J,2) + pow(R,4),0.5)),
       //-Third) + pow(pow(R,-1)*(pow(R,4) + 6*J*(9*J + pow(3,0.5)*pow(27*pow(J,2) + pow(R,4),0.5))),Third))/6.;
		(*par).Mass_irr_up = sqrt((*par).Area_up/(16*Pi));
	}
	if((*par).d<0)
	{
		R = sqrt((*par).Area_down/(4*Pi));
		J = (*par).J_down;
		(*par).Mass_C_down = sqrt(R*R/4. + J*J/(R*R));
		//(*par).Mass_C_down = (R + pow(R, TwoThird)*pow(pow(R,4) + 6*J*(9*J + pow(3,0.5)*pow(27*pow(J,2) + pow(R,4),0.5)),
       //-Third) + pow(pow(R,-1)*(pow(R,4) + 6*J*(9*J + pow(3,0.5)*pow(27*pow(J,2) + pow(R,4),0.5))),Third))/6.;
		(*par).Mass_irr_down = sqrt((*par).Area_down/(16*Pi));

	}
}
