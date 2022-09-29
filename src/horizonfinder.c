#include "main.h"

int Get_A_B_from_rho_z(parameters par, int &dom, ftype &A, ftype &B, ftype rho, ftype z)
{
	ftype norm, normmin = 1e-8, x[2], dx[2], f[2], J[2][2], det;
	int idom,iter=0, itmax=50, restart=0;
	//TODO extra parameter für die Koordinateninversion

	x[0]=A;
	x[1]=B;
	idom = dom;

	f[0] = rho - chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,A,B);
	f[1] = z - chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,A,B);

	norm = sqrt(f[0]*f[0] + f[1]*f[1])/2.;

	//cout << "norm=" << norm << " \titer=" << iter <<"\tx[0]=" << x[0] << "\tx[1]="<< x[1] << "\trho=" <<chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,x[0],x[1]) << "\tz=" <<chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,x[0],x[1])<<endl;

	while(norm > normmin)
	{
		f[0] = rho - chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,x[0],x[1]);
		f[1] = z - chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,x[0],x[1]);

		J[0][0] = chebevxy(0,1,-1,1,par.rho.chebd1[idom],par.ns[idom],par.nt,x[0],x[1]);
		J[1][0] = chebevxy(0,1,-1,1,par.z.chebd1[idom],par.ns[idom],par.nt,x[0],x[1]);
		J[0][1] = chebevxy(0,1,-1,1,par.rho.chebd2[idom],par.ns[idom],par.nt,x[0],x[1]);
		J[1][1] = chebevxy(0,1,-1,1,par.z.chebd2[idom],par.ns[idom],par.nt,x[0],x[1]);

		det = J[0][1]*J[1][0] - J[0][0]*J[1][1];

		dx[0] = -(f[0]*J[1][1] - f[1]*J[0][1])/det;
        dx[1] = -(f[1]*J[0][0] - f[0]*J[1][0])/det;

        x[0]+=dx[0];
        x[1]+=dx[1];

        f[0] = rho - chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,x[0],x[1]);
        f[1] = z - chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,x[0],x[1]);

        norm = sqrt(f[0]*f[0] + f[1]*f[1])/2.;

        //cout << "\t\tnorm=" << norm << " \titer=" << iter <<"\tx[0]=" << x[0] << "\tx[1]="<< x[1] << "\trho=" <<chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,x[0],x[1]) << "\tz=" <<chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,x[0],x[1])<<endl;

        //check if domain has changed
        if(0==idom && x[0]>1)
        {
        	x[0] = x[0]-1;
        	idom = 1;
        }
        if(1==idom && x[0]<0)
        {
        	x[0] = 1+x[0];
           	idom = 0;
        }
        //check if newton raphson went to far
        x[0]=(0==idom && x[0]<0)? x[0]=0 : x[0];
        x[0]=(1==idom && x[0]>1)? x[0]=1 : x[0];
        x[1]=(x[1]>1)? x[1]=1 : x[1];
        x[1]=(x[1]<-1)? x[1]=-1 : x[1];

        iter ++;
        if(iter>itmax)
        {
        	//Wenn die Iteration fehlschlägt gibt es einen neuen Versuch mit neuen Startbedingungen
        	//die in der Mitte des Gebiets liegen
        	//cout << "restart of coordinate-inversion was necessary the setup hang up at x[0]=" << x[0] << "\tx[1]=" << x[1] <<"\tidom=" << idom <<"." << endl;
        	switch(restart)
        	{
				case 0:
				{
					restart++;
					x[0]=1;
					x[1]=0;
					idom=0;
					iter=0;
					break;
				}
				case 1:
				{
					restart++;
					x[0]=0;
					x[1]=1;
					idom=0;
					iter=0;
					break;
				}
				case 2:
				{
					restart++;
					x[0]=0;
					x[1]=-1;
					idom=0;
					iter=0;
					break;
				}
				case 3:
				{
					restart++;
					x[0]=1;
					x[1]=0;
					idom=1;
					iter=0;
					break;
				}
				case 4:
				{
					restart++;
					x[0]=0;
					x[1]=0;
					idom=0;
					iter=0;
					break;
				}
				case 5:
				{
					restart++;
					x[0]=0.5;
					x[1]=1;
					idom=1;
					iter=0;
					break;
				}
				case 6:
				{
					restart++;
					x[0]=0.5;
					x[1]=-1;
					idom=1;
					iter=0;
					break;
				}
				default:
				{
					//cout << "rho=" << rho <<"\tz=" << z << "\tx[0]=" << x[0] << "\tx[1]="<< x[1] << "\tidom=" << idom << "\tnorm=" << norm << " \titer=" << iter <<" \trestart=" << restart << endl;
					//cout << "Koordinateninversion erreicht vorgegebene Norm nicht!\n Breche ab!" << endl;
					return -1;
					break;
				}
        	}
        }
	}
    //cout << "norm=" << norm << " \titer=" << iter <<" \trestart=" << restart << endl;
	A = x[0];
	B = x[1];
	dom = idom;
	return 0;
}

void Get_Chebyshev_Coefficients_v(parameters par, derivs_2D *v, int level)
{
	for(int idom=0; idom<NDOM; idom++)
	{
		if(level >= 0)
		{
			(*v).cheb[idom] = dmatrix(0,par.ns[idom],0,par.nt);
		}
		if(level >= 1)
		{
			(*v).chebd1[idom] = dmatrix(0,par.ns[idom],0,par.nt);
			(*v).chebd2[idom] = dmatrix(0,par.ns[idom],0,par.nt);
		}
		if(level >= 2)
		{
			(*v).chebd11[idom] = dmatrix(0,par.ns[idom],0,par.nt);
			(*v).chebd12[idom] = dmatrix(0,par.ns[idom],0,par.nt);
			(*v).chebd22[idom] = dmatrix(0,par.ns[idom],0,par.nt);
		}
		if(level >= 3)
		{
			(*v).chebd111[idom] = dmatrix(0,par.ns[idom],0,par.nt);
			(*v).chebd112[idom] = dmatrix(0,par.ns[idom],0,par.nt);
			(*v).chebd122[idom] = dmatrix(0,par.ns[idom],0,par.nt);
			(*v).chebd222[idom] = dmatrix(0,par.ns[idom],0,par.nt);
		}

		for(int i=0; i<par.ns[idom];i++)
		{
			for(int j=0;j<par.nt;j++)
			{
				if(level >= 0)
				{
					(*v).cheb[idom][i][j] = (*v).d0[Index(par,idom,0,i,j)];
				}
				if(level >= 1)
				{
					(*v).chebd1[idom][i][j] = (*v).d1[Index(par,idom,0,i,j)];
					(*v).chebd2[idom][i][j] = (*v).d2[Index(par,idom,0,i,j)];
				}
				if(level >= 2)
				{
					(*v).chebd11[idom][i][j] = (*v).d11[Index(par,idom,0,i,j)];
					(*v).chebd12[idom][i][j] = (*v).d12[Index(par,idom,0,i,j)];
					(*v).chebd22[idom][i][j] = (*v).d22[Index(par,idom,0,i,j)];
				}
				if(level >= 3)
				{
					(*v).chebd111[idom][i][j] = (*v).d111[Index(par,idom,0,i,j)];
					(*v).chebd112[idom][i][j] = (*v).d112[Index(par,idom,0,i,j)];
					(*v).chebd122[idom][i][j] = (*v).d122[Index(par,idom,0,i,j)];
					(*v).chebd222[idom][i][j] = (*v).d222[Index(par,idom,0,i,j)];
				}
			}
		}
		if(level >= 0)
		{
			chebft_Extremes_2D((*v).cheb[idom], par.ns[idom], par.nt);
		}
		if(level >= 1)
		{
			chebft_Extremes_2D((*v).chebd1[idom], par.ns[idom], par.nt);
			chebft_Extremes_2D((*v).chebd2[idom], par.ns[idom], par.nt);
		}
		if(level >= 2)
		{
			chebft_Extremes_2D((*v).chebd11[idom], par.ns[idom], par.nt);
			chebft_Extremes_2D((*v).chebd12[idom], par.ns[idom], par.nt);
			chebft_Extremes_2D((*v).chebd22[idom], par.ns[idom], par.nt);
		}
		if(level >= 3)
		{
			chebft_Extremes_2D((*v).chebd111[idom], par.ns[idom], par.nt);
			chebft_Extremes_2D((*v).chebd112[idom], par.ns[idom], par.nt);
			chebft_Extremes_2D((*v).chebd122[idom], par.ns[idom], par.nt);
			chebft_Extremes_2D((*v).chebd222[idom], par.ns[idom], par.nt);
		}
	}
}

void free_Chebyshev_Coefficients_v(parameters par, derivs_2D *v, int level)
{
	for(int idom=0; idom<NDOM; idom++)
	{
		if(level >= 0)
		{
			free_dmatrix((*v).cheb[idom],0,par.ns[idom],0,par.nt);
		}
		if(level >= 1)
		{
			free_dmatrix((*v).chebd1[idom],0,par.ns[idom],0,par.nt);
			free_dmatrix((*v).chebd2[idom],0,par.ns[idom],0,par.nt);
		}
		if(level >= 2)
		{
			free_dmatrix((*v).chebd11[idom],0,par.ns[idom],0,par.nt);
			free_dmatrix((*v).chebd12[idom],0,par.ns[idom],0,par.nt);
			free_dmatrix((*v).chebd22[idom],0,par.ns[idom],0,par.nt);
		}
		if(level >= 3)
		{
			free_dmatrix((*v).chebd111[idom],0,par.ns[idom],0,par.nt);
			free_dmatrix((*v).chebd112[idom],0,par.ns[idom],0,par.nt);
			free_dmatrix((*v).chebd122[idom],0,par.ns[idom],0,par.nt);
			free_dmatrix((*v).chebd222[idom],0,par.ns[idom],0,par.nt);

		}
	}
}

void Get_Add(parameters par, ftype &A11, ftype &A12, ftype &A21, ftype &A22, ftype R, ftype mu)
{
	ftype d1=par.d1, c1=par.c1, Pz1=par.Pz1,
		  d2=par.d2, c2=par.c2, Pz2=par.Pz2,
		  d=par.d;
	//die Sz1 & Sz2 erscheinen nur in A*3-Komponenten
	A11 = (3*Pz1*(d*pow(-d+d1,2)+d*(-5*d+4*d1)*R*pow(mu,3)-mu*R*(pow(d,2)+pow(d1,2)+2*pow(R,2))+pow(mu,2)*(d*(4*d*d1-3*pow(d,2)-2*pow(d1,2))+2*(-3*d+d1)*pow(R,2)))*pow((-d+d1)*(-d+d1-2*mu*R)+pow(R,2),-3))/2.+c1*(2*(2*d+d1)*mu*R-pow(-d+d1,2)+3*pow(d,2)*pow(mu,2)+2*pow(R,2))*pow((-d+d1)*(-d+d1-2*mu*R)+pow(R,2),-2.5)+(3*Pz2*(d*pow(-d+d2,2)+d*(-5*d+4*d2)*R*pow(mu,3)-mu*R*(pow(d,2)+pow(d2,2)+2*pow(R,2))+pow(mu,2)*(d*(4*d*d2-3*pow(d,2)-2*pow(d2,2))+2*(-3*d+d2)*pow(R,2)))*pow((-d+d2)*(-d+d2-2*mu*R)+pow(R,2),-3))/2.+c2*(2*(2*d+d2)*mu*R-pow(-d+d2,2)+3*pow(d,2)*pow(mu,2)+2*pow(R,2))*pow((-d+d2)*(-d+d2-2*mu*R)+pow(R,2),-2.5);
	A12 = (-3*R*pow((-d+d1)*(-d+d1-2*mu*R)+pow(R,2),-3.5)*(-2*c1*d*(d*mu+R)*((-d+d1)*(-d+d1-2*mu*R)+pow(R,2))+Pz1*(-2*d1*(2*d*mu+R)*(d+mu*R)+3*mu*pow(d,3)+(2*d*mu+R)*pow(d1,2)+R*pow(d,2)*(2+5*pow(mu,2))+5*d*mu*pow(R,2)+pow(R,3))*pow((-d+d1)*(-d+d1-2*mu*R)+pow(R,2),0.5)))/2.-(3*R*pow((-d+d2)*(-d+d2-2*mu*R)+pow(R,2),-3.5)*(-2*c2*d*(d*mu+R)*((-d+d2)*(-d+d2-2*mu*R)+pow(R,2))+Pz2*(-2*d2*(2*d*mu+R)*(d+mu*R)+3*mu*pow(d,3)+(2*d*mu+R)*pow(d2,2)+R*pow(d,2)*(2+5*pow(mu,2))+5*d*mu*pow(R,2)+pow(R,3))*pow((-d+d2)*(-d+d2-2*mu*R)+pow(R,2),0.5)))/2.;
	A21 = A12;
	A22 = (pow(R,2)*pow(-1+pow(mu,2),-1)*pow((-d+d1)*(-d+d1-2*mu*R)+pow(R,2),-3.5)*(2*c1*((-d+d1)*(-d+d1-2*mu*R)+pow(R,2))*(2*d*mu*R-2*d1*(d+mu*R)+pow(d1,2)+pow(d,2)*(-2+3*pow(mu,2))+pow(R,2))-3*Pz1*(pow(d,3)*(-2+3*pow(mu,2))+mu*R*pow(d,2)*(-2+5*pow(mu,2))-2*d1*(d+mu*R)*(mu*R+d*(-1+2*pow(mu,2)))+pow(d1,2)*(mu*R+d*(-1+2*pow(mu,2)))+d*(-1+4*pow(mu,2))*pow(R,2)+mu*pow(R,3))*pow((-d+d1)*(-d+d1-2*mu*R)+pow(R,2),0.5)))/2.+(pow(R,2)*pow(-1+pow(mu,2),-1)*pow((-d+d2)*(-d+d2-2*mu*R)+pow(R,2),-3.5)*(2*c2*((-d+d2)*(-d+d2-2*mu*R)+pow(R,2))*(2*d*mu*R-2*d2*(d+mu*R)+pow(d2,2)+pow(d,2)*(-2+3*pow(mu,2))+pow(R,2))-3*Pz2*(pow(d,3)*(-2+3*pow(mu,2))+mu*R*pow(d,2)*(-2+5*pow(mu,2))-2*d2*(d+mu*R)*(mu*R+d*(-1+2*pow(mu,2)))+pow(d2,2)*(mu*R+d*(-1+2*pow(mu,2)))+d*(-1+4*pow(mu,2))*pow(R,2)+mu*pow(R,3))*pow((-d+d2)*(-d+d2-2*mu*R)+pow(R,2),0.5)))/2.;
}

int FieldEquations_HorizonFinder(parameters par, ftype *F, derivs_1D hh, derivs_2D field, ftype *E)
{
	ftype ktr, Rp, mu, R, rho, z, A, B,
	      phi, dphidA, dphidB, drhodA, drhodB, dzdA, dzdB,
		  dRdrho, dRdz, dmudrho, dmudz, dphidR, dphidmu,
		  Add11, Add12, Add21, Add22, hd1, hd11,
		  d=par.d, omega,
		  coeff_h1, coeff_h2, coeff_h3, su1, su2;
	int idom, check;

	ktr=par.Ktr;
	Rp=par.Rp;

	Derivatives_1D(par,hh,-1,par.nH);

	//j=0 -> B=-1
	A=0;
	B=-1;
	idom=0;

	for(int j=0; j<par.nH; j++)
	{
		mu =-cos(Pi * j / cast<ftype>(par.nH-1));
		R = hh.d0[j];

		hd1 = hh.d1[j];
		hd11 = hh.d11[j];

		check = Get_A_B_from_rho_z(par,idom,A,B,R*sqrt(1-mu*mu),R*mu + d);
		if(0!=check) return -1;
		Get_Add(par,Add11, Add12, Add21, Add22,R ,mu);

		phi = chebevxy(0,1,-1,1,field.cheb[idom],par.ns[idom],par.nt,A,B);
		dphidA = chebevxy(0,1,-1,1,field.chebd1[idom],par.ns[idom],par.nt,A,B);
		dphidB = chebevxy(0,1,-1,1,field.chebd2[idom],par.ns[idom],par.nt,A,B);

		rho= chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,A,B);
		drhodA = chebevxy(0,1,-1,1,par.rho.chebd1[idom],par.ns[idom],par.nt,A,B);
		drhodB = chebevxy(0,1,-1,1,par.rho.chebd2[idom],par.ns[idom],par.nt,A,B);
		z= chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,A,B);
		dzdA = chebevxy(0,1,-1,1,par.z.chebd1[idom],par.ns[idom],par.nt,A,B);
		dzdB = chebevxy(0,1,-1,1,par.z.chebd2[idom],par.ns[idom],par.nt,A,B);

		dRdrho =  rho*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-0.5);
		dRdz = (-d + z)*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-0.5);
		dmudrho = rho*(d - z)*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-1.5);
		dmudz = pow(rho,2)*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-1.5);

		dphidR = (dmudrho*dphidB*drhodA - dmudrho*dphidA*drhodB + dmudz*dphidB*dzdA - dmudz*dphidA*dzdB)/((dmudz*dRdrho - dmudrho*dRdz)*(drhodB*dzdA - drhodA*dzdB));
		dphidmu = (-(dphidB*(dRdrho*drhodA + dRdz*dzdA)) + dphidA*(dRdrho*drhodB + dRdz*dzdB))/((dmudz*dRdrho - dmudrho*dRdz)*(drhodB*dzdA - drhodA*dzdB));

		omega = (1 - (2*d*mu*R + pow(d,2) + pow(R,2))*pow(Rp,-2));

		if(1==par.horizon_finder_surface)
		{
			F[j] = (2*ktr)/3.+pow(phi,-2)*pow(R,-1)*(-3*R*pow(hd1,2)*(-1+pow(mu,2))+pow(hd1,3)*(mu-pow(mu,3))+2*hd1*mu*pow(R,2)+(2*R+hd11*(-1+pow(mu,2)))*pow(R,2))*(1-(2*d*mu*R+pow(d,2)+pow(R,2))*pow(Rp,-2))*pow(-(pow(hd1,2)*(-1+pow(mu,2)))+pow(R,2),-1.5)-4*pow(phi,-3)*pow(R,-1)*pow(Rp,-2)*(-(phi*R*(R*(d*mu+R)+d*hd1*(-1+pow(mu,2))))+(dphidmu*hd1*(-1+pow(mu,2))+dphidR*pow(R,2))*(2*d*mu*R+pow(d,2)+pow(R,2)-pow(Rp,2)))*pow(-(pow(hd1,2)*(-1+pow(mu,2)))+pow(R,2),-0.5)+pow(phi,-6)*pow(R,-2)*pow(Rp,-6)*((Add12+Add21)*hd1*(-1+pow(mu,2))*pow(R,2)+Add11*pow(R,4)+Add22*pow(hd1,2)*pow(-1+pow(mu,2),2))*pow(-(pow(hd1,2)*(-1+pow(mu,2)))+pow(R,2),-1)*pow(2*d*mu*R+pow(d,2)+pow(R,2)-pow(Rp,2),3) - E[j];
		}
		if(-1==par.horizon_finder_surface)
		{

			F[j] = (2*ktr)/3.+pow(phi,-2)*pow(R,-1)*(-3*R*pow(hd1,2)*(-1+pow(mu,2))+pow(hd1,3)*(mu-pow(mu,3))+2*hd1*mu*pow(R,2)+(2*R+hd11*(-1+pow(mu,2)))*pow(R,2))*pow(Rp,-2)*(2*d*mu*R+pow(d,2)+pow(R,2)-pow(Rp,2))*pow(-(pow(hd1,2)*(-1+pow(mu,2)))+pow(R,2),-1.5)+4*pow(phi,-3)*pow(R,-1)*pow(Rp,-2)*(-(phi*R*(R*(d*mu+R)+d*hd1*(-1+pow(mu,2))))+(dphidmu*hd1*(-1+pow(mu,2))+dphidR*pow(R,2))*(2*d*mu*R+pow(d,2)+pow(R,2)-pow(Rp,2)))*pow(-(pow(hd1,2)*(-1+pow(mu,2)))+pow(R,2),-0.5)+pow(phi,-6)*pow(R,-2)*pow(Rp,-6)*((Add12+Add21)*hd1*(-1+pow(mu,2))*pow(R,2)+Add11*pow(R,4)+Add22*pow(hd1,2)*pow(-1+pow(mu,2),2))*pow(-(pow(hd1,2)*(-1+pow(mu,2)))+pow(R,2),-1)*pow(2*d*mu*R+pow(d,2)+pow(R,2)-pow(Rp,2),3) - E[j];
		}
		//cout << j<<"\t" <<coeff_h1 <<"\t"<<coeff_h2<<"\t"<<coeff_h3<< "\t"<<su1 <<"\t"<<su2<<endl;
		if(0==j || par.nH-1==j)
		{
			F[0] = hh.d1[0] - E[j];
			F[par.nH-1] = hh.d1[par.nH-1] - E[j];
		}
	}
return 0;
}

int Check_Domain(parameters par, ftype *h)
{//This function checks whether the given h is still in the domain or already inside the Black Holes
// or beyond scri
	int nH = par.nH;
	ftype rho, z, R, mu, Rp=par.Rp, d1=par.d1, d2=par.d2, r1=par.rH1, r2=par.rH2;

	for (int j=0; j < nH; j++)
	{
		mu = -cos(Pi * j / cast<ftype>(nH-1));
		R = h[j];
		rho = R*sqrt(1-mu*mu);
		z = R * mu + par.d;

//		cout << (rho*rho + (z-d1)*(z-d1)) <<">"<<r1*r1<<"\t";
//		cout << (rho*rho + (z-d2)*(z-d2)) <<">"<<r2*r2<<"\t";
//		cout << (rho*rho + z*z) <<"<"<<Rp*Rp<<endl;

		//inside black holes
		if((rho*rho + (z-d1)*(z-d1)) < r1*r1)
		{
			return -1;
		}
		if((rho*rho + (z-d2)*(z-d2)) < r2*r2)
		{
			return -2;
		}
		//beyond scri
		if((rho*rho + z*z) > Rp*Rp)
		{
			return -3;
		}
		//R must be positive
		if(R < 0)
		{
			return -4;
		}
	}
	return 0;
}

int guess_common_horizon(parameters par, ftype &h0, derivs_2D field)
{
	ftype ktr, Rp, mu, R, rho, z, A, B,
	      phi, dphidA, dphidB, drhodA, drhodB, dzdA, dzdB,
		  dRdrho, dRdz, dmudrho, dmudz, dphidR, dphidmu,
		  Add11, Add12, Add21, Add22, h,
		  d=par.d, step, delta,
		  theta,
		  coeff_h1, coeff_h2, coeff_h3, su1, su2;
	int idom, check;

	ktr=par.Ktr;
	Rp=par.Rp;

	h0=-1;

	//j=0 -> B=-1
	A=0;
	B=-1;
	idom=0;

	//black hole 1
	delta = par.Rp - par.d1+par.rH1;
	step = delta/100.;
	for(h=0.99*par.Rp; h>par.d1+par.rH1; h-=step)
	{
		mu = 1;
		R = h;

		check = Get_A_B_from_rho_z(par,idom,A,B,R*sqrt(1-mu*mu),R*mu + d);
		if(0!=check) return -1;
		Get_Add(par,Add11, Add12, Add21, Add22,R ,mu);

		phi = chebevxy(0,1,-1,1,field.cheb[idom],par.ns[idom],par.nt,A,B);
		dphidA = chebevxy(0,1,-1,1,field.chebd1[idom],par.ns[idom],par.nt,A,B);
		dphidB = chebevxy(0,1,-1,1,field.chebd2[idom],par.ns[idom],par.nt,A,B);

		rho= chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,A,B);
		drhodA = chebevxy(0,1,-1,1,par.rho.chebd1[idom],par.ns[idom],par.nt,A,B);
		drhodB = chebevxy(0,1,-1,1,par.rho.chebd2[idom],par.ns[idom],par.nt,A,B);
		z= chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,A,B);
		dzdA = chebevxy(0,1,-1,1,par.z.chebd1[idom],par.ns[idom],par.nt,A,B);
		dzdB = chebevxy(0,1,-1,1,par.z.chebd2[idom],par.ns[idom],par.nt,A,B);

		dRdrho =  rho*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-0.5);
		dRdz = (-d + z)*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-0.5);
		dmudrho = rho*(d - z)*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-1.5);
		dmudz = pow(rho,2)*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-1.5);

		dphidR = (dmudrho*dphidB*drhodA - dmudrho*dphidA*drhodB + dmudz*dphidB*dzdA - dmudz*dphidA*dzdB)/((dmudz*dRdrho - dmudrho*dRdz)*(drhodB*dzdA - drhodA*dzdB));
		dphidmu = (-(dphidB*(dRdrho*drhodA + dRdz*dzdA)) + dphidA*(dRdrho*drhodB + dRdz*dzdB))/((dmudz*dRdrho - dmudrho*dRdz)*(drhodB*dzdA - drhodA*dzdB));

		coeff_h1=(1/R - (2*(d + R))/((d + R - Rp)*(d + R + Rp)))/2.;
		coeff_h2=(ktr*(1 - pow(d + R,2)*pow(Rp,-2)))/6.;
		coeff_h3=-(Add11*pow(-1 + pow(d + R,2)*pow(Rp,-2),2))/2.;

		su1=1;
		su2=0;

		theta = coeff_h1*phi + coeff_h2*phi*phi*phi + coeff_h3/(phi*phi*phi) + su1*dphidR + su2*dphidmu;

		if(theta < 0)
		{
			h0=h;
			//cout << "h0=" << h0 << endl;
			break;
		}
	}

	//black hole 2
	delta = par.Rp - (-par.d2+par.rH2);
	step = delta/100.;
	for(h=0.99*par.Rp; h>(-par.d2+par.rH2); h-=step)
	{
		mu = -1;
		R = h;

		check = Get_A_B_from_rho_z(par,idom,A,B,R*sqrt(1-mu*mu),R*mu + d);
		if(0!=check) return -1;
		Get_Add(par,Add11, Add12, Add21, Add22,R ,mu);

		phi = chebevxy(0,1,-1,1,field.cheb[idom],par.ns[idom],par.nt,A,B);
		dphidA = chebevxy(0,1,-1,1,field.chebd1[idom],par.ns[idom],par.nt,A,B);
		dphidB = chebevxy(0,1,-1,1,field.chebd2[idom],par.ns[idom],par.nt,A,B);

		rho= chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,A,B);
		drhodA = chebevxy(0,1,-1,1,par.rho.chebd1[idom],par.ns[idom],par.nt,A,B);
		drhodB = chebevxy(0,1,-1,1,par.rho.chebd2[idom],par.ns[idom],par.nt,A,B);
		z= chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,A,B);
		dzdA = chebevxy(0,1,-1,1,par.z.chebd1[idom],par.ns[idom],par.nt,A,B);
		dzdB = chebevxy(0,1,-1,1,par.z.chebd2[idom],par.ns[idom],par.nt,A,B);

		dRdrho =  rho*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-0.5);
		dRdz = (-d + z)*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-0.5);
		dmudrho = rho*(d - z)*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-1.5);
		dmudz = pow(rho,2)*pow(-2*d*z + pow(d,2) + pow(rho,2) + pow(z,2),-1.5);

		dphidR = (dmudrho*dphidB*drhodA - dmudrho*dphidA*drhodB + dmudz*dphidB*dzdA - dmudz*dphidA*dzdB)/((dmudz*dRdrho - dmudrho*dRdz)*(drhodB*dzdA - drhodA*dzdB));
		dphidmu = (-(dphidB*(dRdrho*drhodA + dRdz*dzdA)) + dphidA*(dRdrho*drhodB + dRdz*dzdB))/((dmudz*dRdrho - dmudrho*dRdz)*(drhodB*dzdA - drhodA*dzdB));

		coeff_h1=(1/R + 1/(d - R - Rp) + 1/(d - R + Rp))/2.;
		coeff_h2=(ktr*(1 - pow(d - R,2)*pow(Rp,-2)))/6.;
		coeff_h3=-(Add11*pow(-1 + pow(d - R,2)*pow(Rp,-2),2))/2.;

		su1=1;
		su2=0;

		theta = coeff_h1*phi + coeff_h2*phi*phi*phi + coeff_h3/(phi*phi*phi) + su1*dphidR + su2*dphidmu;

		if(theta < 0)
		{
			h0=(h>h0)?h:h0;
			//cout << "h0=" << h0 <<  endl;
			break;
		}
	}
	if(h0>0)
	{
		cout << "h0=" << h0 << endl;
		return 0;
	}
	else
	{
		return -1;
	}
	exit(1);
}

int Jacobian_HorizonFinder(parameters par, derivs_1D h, derivs_2D v, ftype **DF, ftype *E)
{
	int nH=par.nH, check=0;
	ftype epDF = 1.e-8;

	#pragma omp parallel for
	for(int j=0; j < nH; j++)
	{
		int check_m = 0, check_p=0;
		ftype *Fp, *Fm;
		derivs_1D hp, hm;
		hp.d0 = new ftype[par.nH];
		hp.d1 = new ftype[par.nH];
		hp.d11 = new ftype[par.nH];
		hm.d0 = new ftype[par.nH];
		hm.d1 = new ftype[par.nH];
		hm.d11 = new ftype[par.nH];
		Fp = new ftype[par.nH];
		Fm = new ftype[par.nH];

		for(int l=0; l < nH; l++)
		{
			hp.d0[l] = hm.d0[l] = h.d0[l];
		}
		hp.d0[j] += epDF;
		hm.d0[j] -= epDF;

		check_p += FieldEquations_HorizonFinder(par, Fp, hp, v, E);
		if(check_p != 0)
		{//wenn Fp nicht mehr im Gebiet liegt, berechne die Ableitung nur mit Fm und F selbst
			hp.d0[j] = h.d0[j];
			check += FieldEquations_HorizonFinder(par, Fp, hm, v, E);
		}

		check_m += FieldEquations_HorizonFinder(par, Fm, hm, v, E);
		if(check_m != 0)
		{//wenn Fm nicht mehr im Gebiet liegt, berechne die Ableitung nur mit Fp und F selbst
			hm.d0[j] = h.d0[j];
			check += FieldEquations_HorizonFinder(par, Fm, hm, v, E);
		}

		if((check_m!=0) || (check_p!=0))
		{
			for(int l=0; l < nH; l++)
			{
				DF[l][j] = (Fp[l]-Fm[l])/epDF;
			}
		}
		else
		{
			for(int l=0; l < nH; l++)
			{
				DF[l][j] = 0.5*(Fp[l]-Fm[l])/epDF;
			}
		}

		if((check_m!=0) & (check_p!=0)) check +=-1;

		hp.d0[j] = hm.d0[j] = h.d0[j];

		delete[] hp.d0;
		delete[] hp.d1;
		delete[] hp.d11;
		delete[] hm.d0;
		delete[] hm.d1;
	    delete[] hm.d11;
	    delete[] Fp;
	    delete[] Fm;
	}
	if(0!=check) return -1;
	return 0;
}

int newton_HorizonFinder(parameters par, derivs_1D h, derivs_2D v, ftype *E)
{
	int l=1, lmax=par.Newton_itmax, check, nH=par.nH;
	ftype *F, **DF, d, dmin = 1e-5, dmax=1.e5, *h_aux, fac;
	Eigen::MatrixXd A_Eigen(nH,nH);
	Eigen::VectorXd b_Eigen(nH), x_Eigen(nH);
	double dt=0, dtJac=0, dtInv=0, t_start;

	h_aux = new ftype[nH];
	F = new ftype[nH];
	DF = new ftype*[nH];
	for (int i = 0; i < nH; i++)
	{
		DF[i] = new ftype[nH];
	}

	FieldEquations_HorizonFinder(par, F, h, v, E);
	d=normalized_norm(F,nH);

	//cout <<"\t initial norm \t "<< scientific << d << fixed <<endl;

	while(l <= lmax && d > dmin && d < dmax)
	{
		FieldEquations_HorizonFinder(par, F, h, v, E);
		d=normalized_norm(F,nH);

		if (d>0 && d<dmax)
		{
			t_start = omp_get_wtime();
			check = Jacobian_HorizonFinder(par, h, v, DF, E);
			if(0!=check)
			{
				delete[] F;
				delete[] h_aux;
				for (int i = 0; i < nH; i++)
				{
					delete[] DF[i];
				}
				delete[] DF;
				return -1;
			}
			for(int i = 0; i< nH;i++)
			{
				for(int j = 0; j< nH; j++)
				{
					A_Eigen(i,j)=DF[i][j];
				}
				b_Eigen(i)=-F[i];
			}
			dtJac = omp_get_wtime() - t_start;
			x_Eigen=A_Eigen.lu().solve(b_Eigen);

			check = -1;
			fac=1;
			while(check != 0)
			{
				for(int j=0; j<nH; j++)
				{
					h_aux[j] = h.d0[j] + fac*ftype(x_Eigen(j));
				}
				check = Check_Domain(par,h_aux);
				fac = (check != 0)? fac*0.9 : fac;
//				if(check != 0)
//				{
//					cout << "check found: " << check << "\tfac=" << fac << "h="<<h.d0[nH/2] << "h_aux=" << h_aux[nH/2] << "x=" << x_Eigen(nH/2) << endl;
//				}
			}

			for(int j=0; j<nH; j++)
			{
				h.d0[j] +=  fac*ftype(x_Eigen(j));
				//cout << "j="<<j<< "\tX="<<X[j]<<"\tF="<<F[j]<<endl;
			}
			dtInv = omp_get_wtime() - t_start - dtJac;
			dt = omp_get_wtime() - t_start;
		}

		FieldEquations_HorizonFinder(par, F, h, v, E);
		d=normalized_norm(F,nH);

//		cout <<"\r l="<<l<<"\t norm="<<scientific <<d << fixed <<"\t in " << dt << "s\t(Jacobimatrix: "<< dtJac <<"s, Inversion: " << dtInv << "s)" << scientific << flush;

		if(boost::math::isinf(d) || boost::math::isnan(d) || d>dmax)
		{
			delete[] F;
			delete[] h_aux;
			for (int i = 0; i < nH; i++)
			{
				delete[] DF[i];
			}
			delete[] DF;
			return -1;
		}

		if((-1==par.horizon_finder_surface) & (h.d0[par.nH/2]>0.9*par.Rp))
		{
			//Zu nah an Scri macht der MITS-Finder bisweilen Probleme
			cout << "\nMITS zu nah an Scri - setze Horizontfläche auf Scri" << endl;
			for(int j=0; j<nH; j++)
			{
				h.d0[j] = par.Rp;
			}
			FieldEquations_HorizonFinder(par, F, h, v, E);
			d=normalized_norm(F,nH);
			cout << "norm =" << d << endl;
		}

		l++;
	}
	//cout << endl;

	if(l>lmax)
	{
		delete[] F;
		delete[] h_aux;
		for (int i = 0; i < nH; i++)
		{
			delete[] DF[i];
		}
		delete[] DF;
		return -1;
	}

	delete[] F;
	delete[] h_aux;
	for (int i = 0; i < nH; i++)
	{
		delete[] DF[i];
	}
	delete[] DF;
	return 0;
}

int HorizonFinder(parameters *par, derivs_2D v, ftype h0, ftype &h_old)
{
	int check;
	ftype *F, *E, *E0, e;
	derivs_1D h;
	int n_Seq=6, nH=(*par).nH;
	char name[NMAX];

	F = new ftype[(*par).nH];
	E = new ftype[(*par).nH];
	E0 = new ftype[(*par).nH];
	fill0_dvector(E0,0,(*par).nH);
	h.d0 = new ftype[(*par).nH];
	h.d1 = new ftype[(*par).nH];
	h.d11 = new ftype[(*par).nH];
	cout << "---------------------------" << endl;
	cout << "Starting the horizon finder" << endl;
	(*par).horizon_finder_surface = 1;

	//Set initial guess for h
	for(int i=0; i<nH; i++)
	{
		h.d0[i]=h0;
	}

	FieldEquations_HorizonFinder(*par,F, h, v, E0);
	copy_dvector(E0,F,0,(*par).nH);

	for(int i_Seq=0; i_Seq<=n_Seq; i_Seq++)
	{
		e= 1 - i_Seq/cast<ftype>(n_Seq);

		//cout << "e=" << e;

		for(int i=0; i<nH; i++)
		{
			E[i]=e*E0[i];
		}

		check = newton_HorizonFinder(*par, h, v, E);
		//cout << "h[0]=" << h.d0[0] << "\td=" << (*par).d << endl;
		if(0!=check)
		{
			delete[] F;
			delete[] E;
			delete[] E0;
			delete[] h.d0;
			delete[] h.d1;
			delete[] h.d11;
			return -1;
		}
	}
	cout << "print Horizon " << endl;
	if((*par).d>0)
	{
		h_old=-1;
		for(int j=0; j<nH; j++)
		{
			(*par).MOTS_up[j]=h.d0[j];
		}
		sprintf(name,"../run/plots/plot_horizon_up_%ld",(*par).timestamp);
	}
	if((*par).d<0)
	{
		h_old=-1;
		for(int j=0; j<nH; j++)
		{
			(*par).MOTS_down[j]=h.d0[j];
		}
		sprintf(name,"../run/plots/plot_horizon_down_%ld",(*par).timestamp);
	}
	if(0==(*par).d)
	{
		h_old=maximum(h.d0[0],h.d0[nH-1]);
		for(int j=0; j<nH; j++)
		{
			(*par).MOTS_common[j]=h.d0[j];
		}
		sprintf(name,"../run/plots/plot_horizon_common_%ld",(*par).timestamp);
	}
	PrintHorizon(*par, h.d0, name);
	Check_dTheta(*par, h, v);
	AngularMomentum(par,v,h);
	Horizon_Area(par, v, h);
	Local_Mass(par);
	//proper_circumference(par,v,h);

	delete[] F;
	delete[] E;
	delete[] E0;
	delete[] h.d0;
	delete[] h.d1;
	delete[] h.d11;
	return 0;
}

void Check_dTheta(parameters par, derivs_1D hh, derivs_2D field)
{
	int nH=par.nH, check_p, check_m;
	ftype *E, *Fp, *Fm, *dTheta, eps=1e-2;
	derivs_1D h;

	Fp=new ftype[nH];
	Fm=new ftype[nH];
	E=new ftype[nH];
	dTheta=new ftype[nH];
	fill0_dvector(E,0,nH);
	h.d0 = new ftype[par.nH];
	h.d1 = new ftype[par.nH];
	h.d11 = new ftype[par.nH];

	if(-1==par.horizon_finder_surface)
	{
		eps = 1e-1;
	}

	for(int j=0;j<nH; j++)
	{
		h.d0[j]=hh.d0[j] + eps;
	}
	Derivatives_1D(par, h, -1, nH);
	check_p = FieldEquations_HorizonFinder(par,Fp,h,field,E);
	if(0 != check_p)
	{
		for(int j=0;j<nH; j++)
		{
			h.d0[j]=hh.d0[j];
		}
		Derivatives_1D(par, h, -1, nH);
		FieldEquations_HorizonFinder(par,Fp,h,field,E);
	}

	for(int j=0;j<nH; j++)
	{
		h.d0[j]=hh.d0[j] - eps;
	}
	Derivatives_1D(par, h, -1, nH);
	check_m = FieldEquations_HorizonFinder(par,Fm,h,field,E);
	if(0 != check_m)
	{
		for(int j=0;j<nH; j++)
		{
			h.d0[j]=hh.d0[j];
		}
		Derivatives_1D(par, h, -1, nH);
		FieldEquations_HorizonFinder(par,Fm,h,field,E);
	}

	for(int j=0; j<nH; j++)
	{
		//cout << "j=" << j << "\t  Fp=" << Fp[j] << "\t  Fm=" << Fm[j] << "\t dTheta=" << (Fp[j] - Fm[j]) << "\t\th=" << hh.d0[j] << endl;
		dTheta[j]=(Fp[j] - Fm[j]);
	}
	//cout << "check_m=" << check_m << "\tcheck_p="<<check_p << "\teps="<<eps<< endl;
	cout << "Vorzeichen des Normalenableitung des Horizontes = " << Check_horizon_normal_derivative_sign(dTheta, nH) << endl;

	delete[] Fp;
	delete[] Fm;
	delete[] E;
	delete[] dTheta;
	delete[] h.d0;
	delete[] h.d1;
	delete[] h.d11;
}

void MITS_Finder(parameters par, derivs_2D v)
{
	int check;
	ftype *F, *E, *E0, e;
	derivs_1D h;
	int n_Seq=6, nH=par.nH;
	char name[NMAX];

	F = new ftype[par.nH];
	E = new ftype[par.nH];
	E0 = new ftype[par.nH];
	fill0_dvector(E0,0,par.nH);
	h.d0 = new ftype[par.nH];
	h.d1 = new ftype[par.nH];
	h.d11 = new ftype[par.nH];

	cout << "--------------------------------------------" << endl;
	cout << "Search for marginally inner trapped surfaces" << endl;
	par.horizon_finder_surface = -1;

	for(int j=0; j<par.nH; j++)
	{
		if(par.d>0)
		{
			h.d0[j]=par.MOTS_up[j];
		}
		if(par.d<0)
		{
			h.d0[j]=par.MOTS_down[j];
		}
		if(0==par.d)
		{
			if(1==par.MOTS_common_exist)
			{
				h.d0[j]=par.MOTS_common[j];
			}
			else
			{
				h.d0[j]=maximum(par.d1+par.rH1,-par.d2+par.rH2) + 0.1*(par.Rp - maximum(par.d1+par.rH1,-par.d2+par.rH2));
				//h.d0[j] = par.Rp;
			}

		}
		//cout << "j=" << j<< "\th="<<h.d0[j] << endl;
	}

	fill0_dvector(E0,0,nH);
	FieldEquations_HorizonFinder(par,F, h, v, E0);
	copy_dvector(E0,F,0,nH);

	for(int i_Seq=0; i_Seq<=n_Seq; i_Seq++)
	{
		e= 1 - i_Seq/cast<ftype>(n_Seq);

		//cout << "e=" << e;

		for(int i=0; i<nH; i++)
		{
			E[i]=e*E0[i];
		}

		check = newton_HorizonFinder(par, h, v, E);
		//cout << "h[0]=" << h.d0[0] << "\td=" << par.d << endl;
		if(0!=check)
		{
			delete[] F;
			delete[] E;
			delete[] E0;
			delete[] h.d0;
			delete[] h.d1;
			delete[] h.d11;
			cout << "no marginally inner trapped surface found" << endl;
			return;
		}
	}
	cout <<"found marginally inner trapped surface, printing surface" << endl;
	if(par.d>0)
	{
		sprintf(name,"../run/plots/plot_MITS_up_%ld",par.timestamp);
	}
	if(par.d<0)
	{
		sprintf(name,"../run/plots/plot_MITS_down_%ld",par.timestamp);
	}
	if(0==par.d)
	{
		sprintf(name,"../run/plots/plot_MITS_common_%ld",par.timestamp);
	}
	PrintHorizon(par, h.d0, name);

	Check_dTheta(par, h, v);

}

int Check_horizon_normal_derivative_sign(ftype *dTheta, int n)
{
	int signum=0;
	//Die Ränder sind durch Bedingungen an die Ableitung gegeben
	for(int i = 1; i<n-2; i++)
	{
		signum+=sgn(dTheta[i]);
	}
	return signum/(n-3);
}
