#include "main.h"

void AngularMomentum(parameters *par, derivs_2D v, derivs_1D h)
{
	int nH=(*par).nH, check, idom;
	ftype sqrtqc, Omega, As, mu, R, hd1, *Integrand, *c_Integrand, J,
		  d=(*par).d, d1=(*par).d1, d2=(*par).d2, Sz1=(*par).Sz1, Sz2=(*par).Sz2,
		  A, B, phi, rho, z, omega, Rp=(*par).Rp;

	Integrand = new ftype[nH];
	c_Integrand = new ftype[nH];

	Derivatives_1D(*par,h,-1,nH);

	A=0;
	B=0;
	idom=1;

	for(int j=0; j<nH; j++)
	{
		mu = -cos(Pi * j / cast<ftype>((*par).nH-1));
		R = h.d0[j];
		hd1= h.d1[j];

		check = Get_A_B_from_rho_z(*par,idom,A,B,R*sqrt(1-mu*mu),R*mu + d);
		if(0!=check) exit(1);

		As = -3*R*(R*(d*mu+R)+d*hd1*(-1+pow(mu,2)))*(-1+pow(mu,2))*pow(1-pow(hd1,2)*(-1+pow(mu,2))*pow(R,-2),-0.5)*pow(-2*d*d1+2*d*mu*R-2*d1*mu*R+pow(d,2)+pow(d1,2)+pow(R,2),-2.5)*pow(-2*d*d2+2*d*mu*R-2*d2*mu*R+pow(d,2)+pow(d2,2)+pow(R,2),-2.5)*(Sz2*pow(-2*d*d1+2*d*mu*R-2*d1*mu*R+pow(d,2)+pow(d1,2)+pow(R,2),2.5)+Sz1*pow(-2*d*d2+2*d*mu*R-2*d2*mu*R+pow(d,2)+pow(d2,2)+pow(R,2),2.5));
		sqrtqc = pow(R,3)*pow(pow(-(pow(hd1,2)*(-1 + pow(mu,2))) + pow(R,2),-1),0.5);

		phi = chebevxy(0,1,-1,1,v.cheb[idom],(*par).ns[idom],(*par).nt,A,B);
		rho= chebevxy(0,1,-1,1,(*par).rho.cheb[idom],(*par).ns[idom],(*par).nt,A,B);
		z= chebevxy(0,1,-1,1,(*par).z.cheb[idom],(*par).ns[idom],(*par).nt,A,B);

		omega = 1 - (rho*rho + z*z)/(Rp*Rp);
		Omega = omega/(phi*phi);

		//Integrand[j] = Omega*Omega*Omega*As*sqrtqc;
		Integrand[j] = As*sqrtqc;
	}

	Chebyshev_Coefficients_Lobatto(Integrand, c_Integrand, nH-1);

	J=0.25*Chebyshev_Definite_Integral_Coefficients(c_Integrand, nH-1);

	cout << "J=" << J << endl;

	if(0==d)
	{
		(*par).J_common=J;
	}
	if(d>0)
	{
		(*par).J_up=J;
	}
	if(d<0)
	{
		(*par).J_down=J;
	}

	delete[] Integrand;
	delete[] c_Integrand;
}
