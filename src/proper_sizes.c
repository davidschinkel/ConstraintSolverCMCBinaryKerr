#include "main.h"

//void proper_circumference(parameters *par, derivs_2D v, derivs_1D h)
//{
//	int nH = (*par).nH, check, idom;
//	ftype r, mu, A, B, *Integrand, *c_Integrand, proper_circumference, omega, phi, rho, z,
//		  d=(*par).d, Rp=(*par).Rp;
//
//	Integrand = new ftype[nH];
//	c_Integrand = new ftype[nH];
//
//	for(int j=0; j<nH; j++)
//	{
//		mu = -cos(Pi * j / cast<ftype>(nH-1));
//		r = h.d0[j];
//
//		check = Get_A_B_from_rho_z((*par),idom,A,B,r*sqrt(1-mu*mu),r*mu + d);
//		phi = chebevxy(0,1,-1,1,v.cheb[idom],(*par).ns[idom],(*par).nt,A,B);
//		rho= chebevxy(0,1,-1,1,(*par).rho.cheb[idom],(*par).ns[idom],(*par).nt,A,B);
//		z= chebevxy(0,1,-1,1,(*par).z.cheb[idom],(*par).ns[idom],(*par).nt,A,B);
//		if(0!=check)
//		{
//			exit(1);
//		}
//
//		omega = 1 - (rho*rho + z*z)/(Rp*Rp);
//
//		Integrand[j] = (phi*phi/omega)*r;
//
//	}
//
//		Chebyshev_Coefficients_Lobatto(Integrand, c_Integrand, nH-1);
//
//		proper_circumference = 2*Chebyshev_Definite_Integral_Coefficients(c_Integrand, nH-1);
//
//		delete[] Integrand;
//		delete[] c_Integrand;
//
//		cout <<"proper_circumference=" << proper_circumference<<endl;
//
//		if(0==(*par).d)
//		{
//			(*par).proper_circumference_common = proper_circumference;
//		}
//		if((*par).d>0)
//		{
//			(*par).proper_circumference_up = proper_circumference;
//		}
//		if((*par).d<0)
//		{
//			(*par).proper_circumference_down = proper_circumference;
//		}
//}

void proper_distance(parameters *par, derivs_2D v)
{
	int nt = (*par).nt;
	ftype *Integrand, *c_Integrand, proper_distance, omega, phi, rho, z, zmax, zmin,
		  Rp=(*par).Rp;

	Integrand = new ftype[nt];
	c_Integrand = new ftype[nt];

	for(int j=0; j<nt; j++)
	{
		phi = v.d0[Index((*par),1,0,(*par).ns[1]-1,j)];
		rho = 0;
		z = (*par).z.d0[Index((*par),1,0,(*par).ns[1]-1,j)];

		omega = 1 - (rho*rho + z*z)/(Rp*Rp);

		Integrand[j] = (phi*phi/omega);
	}
	zmax=(*par).z.d0[Index((*par),1,0,1,0)];
	zmin=(*par).z.d0[Index((*par),1,0,1,nt-1)];
//cout << zmax << endl;
//cout << zmin << endl;
	Chebyshev_Coefficients_Lobatto(Integrand, c_Integrand, nt-1);

	proper_distance = 0.5*(zmax-zmin)*Chebyshev_Definite_Integral_Coefficients(c_Integrand, nt-1);

	delete[] Integrand;
	delete[] c_Integrand;

	cout<<"proper_distance="<< proper_distance << "\t proper_distance/Rp="<< proper_distance/Rp << endl;
	(*par).proper_distance = proper_distance;
}
