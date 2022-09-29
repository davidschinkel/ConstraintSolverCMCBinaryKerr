#include "main.h"

void Getr(parameters *par)
{
	int nt=(*par).nt, Nt=nt-1, iter=0, itmax=20;
	ftype B, facB, norm, r, f, df, max_norm=0, *c, *dc, *ddc;
	ofstream fp;

	(*par).r = new ftype[nt];
	(*par).dr = new ftype[nt];
	(*par).ddr = new ftype[nt];
	c = new ftype[nt];
	dc = new ftype[nt];
	ddc = new ftype[nt];

	for(int t=0; t<nt;t++)
	{
		B = -cos(Pi * t / cast<ftype>(Nt));

		r = (0==t)? sqrt(acosh(1/(*par).nu)) : (*par).r[t-1]; // (condition) ? true : false sets initial guess
		norm = implicit_r(*par,r,B);

		iter=0;
		while(iter < itmax)
		{

			f=implicit_r(*par,r,B);
			df=implicit_r_dr(*par,r,B);

			r= r  - f/df;

			norm = fabs(implicit_r(*par,r,B));
			iter++;
		}
		max_norm = (norm>max_norm)? norm : max_norm;
		(*par).r[t]=r;
	}
	cout << "Berechnung von r(B) mit maximalem Fehler von " << max_norm << endl;
	facB = ftype(-1);
	Chebyshev_Coefficients_Lobatto((*par).r,c,nt-1);
	Chebyshev_Coefficients_Derivative(c,dc,nt-1);
	Chebyshev_Coefficients_Derivative(dc,ddc,nt-1);
	Chebyshev_Collocations_Lobatto((*par).dr,dc,nt-1);
	Chebyshev_Collocations_Lobatto((*par).ddr,ddc,nt-1);
	for(int j=0;j<(*par).nt;j++)
	{
		(*par).dr[j]=(*par).dr[j]*facB;
		(*par).ddr[j]=(*par).ddr[j]*facB*facB;
	}
	fp.open("../run/test/r");
	fp.precision(PREC_OUTPUT);
	for(int t=0;t<nt;t++)
	{
		B = -cos(Pi * t / cast<ftype>(Nt));
		fp << B << "\t" << (*par).r[t] << "\t" << t <<"\t"<< fabs(c[t]) << "\t" << (*par).dr[t] << endl;
	}
	fp.close();
	delete[] c;
	delete[] dc;
	delete[] ddc;
}


ftype implicit_r(parameters par, ftype r, ftype B)
{
	ftype Rp=par.Rp, a0=par.a0;

	return pow(a0,2) - pow(Rp,2) - 2*cos(cos((B*Pi)/ftype(2))*pow(r,2))*pow(a0,2)*pow(cos(cos((B*Pi)/ftype(2))*pow(r,2)) - cosh(pow(r,2)*sin((B*Pi)/ftype(2))),-1);
}

ftype implicit_r_dr(parameters par, ftype r, ftype B)
{
	ftype a0=par.a0;

	return -4*r*pow(a0,2)*pow(cos(cos((B*Pi)/ftype(2))*pow(r,2)) - cosh(pow(r,2)*sin((B*Pi)/ftype(2))),-2)*(cos((B*Pi)/ftype(2))*cosh(pow(r,2)*sin((B*Pi)/ftype(2)))*sin(cos((B*Pi)/ftype(2))*pow(r,2)) + cos(cos((B*Pi)/ftype(2))*pow(r,2))*sin((B*Pi)/ftype(2))*sinh(pow(r,2)*sin((B*Pi)/ftype(2))));
}
