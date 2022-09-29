#include "main.h"

int indx_H(parameters par, int i1, int idom, int i, int j)
{
	int ns0 = par.ns[0], nt = par.nt;
	return i1 - 1 + Hmax * (j + nt * (i + ns0 * idom));
}

int indx_ah0(parameters par, int i, int B)
{//Indexfunktion für die Koeffizienten des scheinbaren Horizontes ohne Index
//Um die beiden Horizonte zu unterscheiden wird der Wert der spektralen Koordinate B
//mit übergeben +1 ist der obere und -1 der untere Horizont
	int ns=par.ns[1];
	return i + ns*((B+1)/2);
}

int indx_ah1(parameters par, int i1, int i, int B)
{//Indexfunktion für die Koeffizienten des scheinbaren Horizontes mit einem Index
//Um die beiden Horizonte zu unterscheiden wird der Wert der spektralen Koordinate B
//mit übergeben +1 ist der obere und -1 der untere Horizont
	return i1 -1 + n1_h1*((B+1)/2 + 2*i);
}

void Get_Coefficients(parameters *par)
{
	int ns[NDOM], nt=(*par).nt,
		nH_max = Hmax * (*par).nt * ((*par).ns[0] + (*par).ns[1]),
		na1h_max = 2*(*par).ns[1], //idom=1
		na3h1_max =n1_h1-1 + n1_h1*(1 + 2*(*par).ns[1]); //idom=1

	ftype Rp=(*par).Rp, ktr=(*par).Ktr,
		  c1=(*par).c1, c2=(*par).c2,
		  d1=(*par).d1, d2=(*par).d2,
		  Sz1=(*par).Sz1, Sz2=(*par).Sz2,
		  Pz1=(*par).Pz1, Pz2=(*par).Pz2,
		  rho , z, AA, omega, omega6, d;
	for(int i=0; i<NDOM; i++) ns[i]=(*par).ns[i];

	//Hamiltonzwangsbedingung
	(*par).H = new ftype[nH_max];
	for(int idom=0; idom<NDOM; idom++)
	{

		for(int i=0; i<ns[idom]; i++)
		{
			for(int j=0; j<nt; j++)
			{
				rho = (*par).rho.d0[Index(*par,idom,0,i,j)];
				z = (*par).z.d0[Index(*par,idom,0,i,j)];
				omega = 1 - (rho*rho + z*z)/(Rp*Rp);

				//cout << "\tidom = " << idom << "\ti=" << i <<"\tj=" << j <<endl;
				if(0==(*par).coefficientswitch)
				{
					//H1
					(*par).H[indx_H(*par,1,idom,i,j)] = rho*(8*omega*omega);
					//H2
					(*par).H[indx_H(*par,2,idom,i,j)] = rho*(8*omega*omega);
					//H3
					(*par).H[indx_H(*par,3,idom,i,j)] = 0;
					//H4
					(*par).H[indx_H(*par,4,idom,i,j)] = 8*(Rp*Rp*Rp*Rp - 2*Rp*Rp*z*z + z*z*z*z - rho*rho*rho*rho)/(Rp*Rp*Rp*Rp);
					//H5
					(*par).H[indx_H(*par,5,idom,i,j)] = rho*(-16*z*(-Rp*Rp + z*z + rho*rho))/(Rp*Rp*Rp*Rp);
					//H6
					(*par).H[indx_H(*par,6,idom,i,j)] = rho*(24/(Rp*Rp));
					//H7
					(*par).H[indx_H(*par,7,idom,i,j)] = rho*(-TwoThird*ktr*ktr);
					//H8
					AA = pow(c1*(-2*d1*z+pow(d1,2)+pow(rho,2)-2*pow(z,2))*pow(pow(rho,2)+pow(d1-z,2),-2.5)+c2*(-2*d2*z+pow(d2,2)+pow(rho,2)-2*pow(z,2))*pow(pow(rho,2)+pow(d2-z,2),-2.5)+(3*Pz1*z*(-2*d1*z+pow(d1,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3))/2.+(3*Pz2*z*(-2*d2*z+pow(d2,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3))/2.,2)+18*pow(rho,4)*pow(Sz1*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-2.5)+Sz2*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-2.5),2)+18*pow(rho,2)*pow(z,2)*pow(Sz1*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-2.5)+Sz2*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-2.5),2)+pow(3*Pz1*z*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-2)-2*c1*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-1.5)+3*Pz2*z*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-2)-2*c2*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-1.5),2)/4.+pow(pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c1*(-4*z*pow(d1,3)+pow(d1,4)-2*pow(rho,4)-pow(rho,2)*pow(z,2)+pow(d1,2)*(-pow(rho,2)+6*pow(z,2))+d1*(2*z*pow(rho,2)-4*pow(z,3))+pow(z,4))+3*Pz1*z*pow(d1-z,2)*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),0.5))+pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c2*(-4*z*pow(d2,3)+pow(d2,4)-2*pow(rho,4)-pow(rho,2)*pow(z,2)+pow(d2,2)*(-pow(rho,2)+6*pow(z,2))+d2*(2*z*pow(rho,2)-4*pow(z,3))+pow(z,4))+3*Pz2*z*pow(d2-z,2)*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),0.5)),2)/4.+(9*pow(rho,2)*pow(pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3.5)*(2*c1*z*(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2))-Pz1*(-2*d1*z+pow(d1,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),0.5))+pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3.5)*(2*c2*z*(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2))-Pz2*(-2*d2*z+pow(d2,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),0.5)),2))/2.;
					omega6=omega*omega*omega*omega*omega*omega;
					(*par).H[indx_H(*par,8,idom,i,j)] = rho*(AA*omega6);
					//H9
					(*par).H[indx_H(*par,9,idom,i,j)] = 0;

				}
				if(1==(*par).coefficientswitch)
				{
					//H1
					(*par).H[indx_H(*par,1,idom,i,j)] = rho;
					//H2
					(*par).H[indx_H(*par,2,idom,i,j)] = rho;
					//H3
					(*par).H[indx_H(*par,3,idom,i,j)] = 0;
					//H4
					(*par).H[indx_H(*par,4,idom,i,j)] = 1;
					//H5
					(*par).H[indx_H(*par,5,idom,i,j)] = 0;
					//H6
					(*par).H[indx_H(*par,6,idom,i,j)] = 0;
					//H7
					(*par).H[indx_H(*par,7,idom,i,j)] = 0;
					//H8
					(*par).H[indx_H(*par,8,idom,i,j)] = 0;
					//H9
					(*par).H[indx_H(*par,9,idom,i,j)] = 0;
				}
			}
		}
	}
	//scheinbarer Horizont
	//*a1h, *a2h, *a3h
	(*par).a1h = new ftype[na1h_max];
	(*par).a2h = new ftype[na1h_max];
	(*par).a3h = new ftype[na1h_max];
	(*par).s = new ftype[na3h1_max];

	for(int i=0; i<(*par).ns[1];i++)
	{
		d=d1;
		rho = (*par).rho.d0[Index(*par,1,0,i,0)];
		z = (*par).z.d0[Index(*par,1,0,i,0)];
		omega = 1 - (rho*rho + z*z)/(Rp*Rp);
		(*par).a1h[indx_ah0(*par,i,-1)] = ((2*d*z - pow(rho,2) - pow(Rp,2) - pow(z,2))*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-0.5)*pow(pow(rho,2) - pow(Rp,2) + pow(z,2),-1))/2.;
		(*par).a2h[indx_ah0(*par,i,-1)] = ktr/(6*omega);
		(*par).a3h[indx_ah0(*par,i,-1)] = omega*omega*((-(pow(-d+z,2)*pow(pow(rho,2)+pow(d-z,2),-1)*(-(c1*(-2*d1*z+pow(d1,2)+pow(rho,2)-2*pow(z,2))*pow(pow(rho,2)+pow(d1-z,2),-2.5))-c2*(-2*d2*z+pow(d2,2)+pow(rho,2)-2*pow(z,2))*pow(pow(rho,2)+pow(d2-z,2),-2.5)-(3*Pz1*z*(-2*d1*z+pow(d1,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3))/2.-(3*Pz2*z*(-2*d2*z+pow(d2,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3))/2.))-pow(rho,2)*pow(pow(rho,2)+pow(d-z,2),-1)*((pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c1*(-4*z*pow(d1,3)+pow(d1,4)-2*pow(rho,4)-pow(rho,2)*pow(z,2)+pow(d1,2)*(-pow(rho,2)+6*pow(z,2))+d1*(2*z*pow(rho,2)-4*pow(z,3))+pow(z,4))+3*Pz1*z*pow(d1-z,2)*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),0.5)))/2.+(pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c2*(-4*z*pow(d2,3)+pow(d2,4)-2*pow(rho,4)-pow(rho,2)*pow(z,2)+pow(d2,2)*(-pow(rho,2)+6*pow(z,2))+d2*(2*z*pow(rho,2)-4*pow(z,3))+pow(z,4))+3*Pz2*z*pow(d2-z,2)*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),0.5)))/2.)-2*rho*(-d+z)*pow(pow(rho,2)+pow(d-z,2),-1)*((-3*rho*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c1*z*(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2))+Pz1*(-2*d1*z+pow(d1,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),0.5)))/2.-(3*rho*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c2*z*(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2))+Pz2*(-2*d2*z+pow(d2,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),0.5)))/2.))/4.);
		(*par).s[indx_ah1(*par,1,i,-1)] = rho*pow(pow(rho,2) + pow(d - z,2),-0.5);
		(*par).s[indx_ah1(*par,2,i,-1)] = (-d + z)*pow(pow(rho,2) + pow(d - z,2),-0.5);
		(*par).s[indx_ah1(*par,3,i,-1)] = 0;

		d=d2;
		rho = (*par).rho.d0[Index(*par,1,0,i,(*par).nt-1)];
		z = (*par).z.d0[Index(*par,1,0,i,(*par).nt-1)];
		omega = 1 - (rho*rho + z*z)/(Rp*Rp);
		(*par).a1h[indx_ah0(*par,i,1)] = ((2*d*z - pow(rho,2) - pow(Rp,2) - pow(z,2))*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-0.5)*pow(pow(rho,2) - pow(Rp,2) + pow(z,2),-1))/2.;
		(*par).a2h[indx_ah0(*par,i,1)] = ktr/(6*omega);
		(*par).a3h[indx_ah0(*par,i,1)] = omega*omega*((-(pow(-d+z,2)*pow(pow(rho,2)+pow(d-z,2),-1)*(-(c1*(-2*d1*z+pow(d1,2)+pow(rho,2)-2*pow(z,2))*pow(pow(rho,2)+pow(d1-z,2),-2.5))-c2*(-2*d2*z+pow(d2,2)+pow(rho,2)-2*pow(z,2))*pow(pow(rho,2)+pow(d2-z,2),-2.5)-(3*Pz1*z*(-2*d1*z+pow(d1,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3))/2.-(3*Pz2*z*(-2*d2*z+pow(d2,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3))/2.))-pow(rho,2)*pow(pow(rho,2)+pow(d-z,2),-1)*((pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c1*(-4*z*pow(d1,3)+pow(d1,4)-2*pow(rho,4)-pow(rho,2)*pow(z,2)+pow(d1,2)*(-pow(rho,2)+6*pow(z,2))+d1*(2*z*pow(rho,2)-4*pow(z,3))+pow(z,4))+3*Pz1*z*pow(d1-z,2)*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),0.5)))/2.+(pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c2*(-4*z*pow(d2,3)+pow(d2,4)-2*pow(rho,4)-pow(rho,2)*pow(z,2)+pow(d2,2)*(-pow(rho,2)+6*pow(z,2))+d2*(2*z*pow(rho,2)-4*pow(z,3))+pow(z,4))+3*Pz2*z*pow(d2-z,2)*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),0.5)))/2.)-2*rho*(-d+z)*pow(pow(rho,2)+pow(d-z,2),-1)*((-3*rho*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c1*z*(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2))+Pz1*(-2*d1*z+pow(d1,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d1*z+pow(d1,2)+pow(rho,2)+pow(z,2),0.5)))/2.-(3*rho*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),-3.5)*(-2*c2*z*(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2))+Pz2*(-2*d2*z+pow(d2,2)+pow(rho,2)+2*pow(z,2))*pow(-2*d2*z+pow(d2,2)+pow(rho,2)+pow(z,2),0.5)))/2.))/4.);
		(*par).s[indx_ah1(*par,1,i,1)] = rho*pow(pow(rho,2) + pow(d - z,2),-0.5);
		(*par).s[indx_ah1(*par,2,i,1)] = (-d + z)*pow(pow(rho,2) + pow(d - z,2),-0.5);
		(*par).s[indx_ah1(*par,3,i,1)] = 0;
	}

}

void Free_Coefficients(parameters *par)
{
	delete[] (*par).H;
	delete[] (*par).s;
	delete[] (*par).a1h;
	delete[] (*par).a2h;
	delete[] (*par).a3h;
	delete[] (*par).r;
	delete[] (*par).dr;
	delete[] (*par).ddr;
	delete[] (*par).psi.d0;
	delete[] (*par).kappa.d0;
	delete[] (*par).MOTS_common;
	delete[] (*par).MOTS_up;
	delete[] (*par).MOTS_down;
	delete[] (*par).MITS_common;
	delete[] (*par).MITS_up;
	delete[] (*par).MITS_down;
	free_derivs_2D(&(*par).rho);
	free_derivs_2D_3rd_derivatives(&(*par).rho);
	free_Chebyshev_Coefficients_v(*par, &(*par).rho,3);
	free_derivs_2D(&(*par).z);
	free_derivs_2D_3rd_derivatives(&(*par).z);
	free_Chebyshev_Coefficients_v(*par, &(*par).z,3);
}
