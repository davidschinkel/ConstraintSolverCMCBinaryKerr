#include "main.h"

// Scanning and printing data files.

void ScanConfig(parameters *par)
{
	FILE *stream_config;
	char string[NMAX], number[NMAX];
	int i_aux, j_aux, idom, ns, nt, index;  // i_aux: an int auxiliary variable
	ftype d_aux;                    // d_aux: a ftype auxiliary variable

	stream_config = fopen("Config", "r");
	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*par).Newton_itmin            = i_aux;
	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*par).Newton_itmax            = i_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).Newton_tol_factor_steps = d_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).Newton_tol_finish = d_aux;
	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*par).Newton_verb             = i_aux;

	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*par).num_threads             = i_aux;

	for (idom = 0; idom < NDOM; idom++)
	{
		j_aux = fscanf(stream_config, " %d %s ", &ns, string);
		(*par).ns[idom] = ns;
	}
	j_aux = fscanf(stream_config, " %d %s ", &nt, string);
	(*par).nt = nt;

	(*par).nst[0] = 0;
	for (idom = 0; idom < NDOM; idom++)
	{
		ns = (*par).ns[idom];
		(*par).nst[idom + 1] = (*par).nst[idom] + ns * nt;
	}
	//Berechnen der Anzahl an Punkten für alle Potential auf allen Gebieten
	(*par).ntotal=NPOT*(*par).nt*((*par).ns[0] +(*par).ns[1]);
	//Berechnen der Anzahl an Punkten für einen Skalar auf allen Gebieten
	(*par).Scalarlength=((*par).ns[0] + (*par).ns[1]) * nt;

	//Berechnet den Zusammenhang zwischen dem Index "n" in X und
	//"n_v" in v. Dadurch, dass der Funktionswert am Gebiets-
	//übergang nur einmal gespeichert wird verbessert sich die
	//Geschwindigkeit des bicgstab-Verfahrens.
	//Get_Arrays_n_And_n_v(parGoal);
	//In der derzeitigen Implementation werden die Funktionswerte
	//am Gebietsüberlapp getrennt gespeichert, deswegen ist
	//Get_Arrays_n_And_n_v nicht notwendig und n_2D wird auf ntotal gesetzt
	(*par).n_2D=(*par).ntotal;

	j_aux = fscanf(stream_config, " %d  %s ", &i_aux, string);
	(*par).n_Seq = i_aux;

	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).a0 = d_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).Rp = d_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).Ktr = d_aux;

	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).rH1 = d_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).c1 = d_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).Sz1 = d_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).Pz1 = d_aux;

	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).rH2 = d_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).c2 = d_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).Sz2 = d_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*par).Pz2 = d_aux;

	(*par).eta1 = asinh((*par).a0/(*par).rH1);
	(*par).eta2 = -asinh((*par).a0/(*par).rH2);

	(*par).d1 = sqrt((*par).a0*(*par).a0 + (*par).rH1*(*par).rH1);
	(*par).d2 = -sqrt((*par).a0*(*par).a0 + (*par).rH2*(*par).rH2);

	(*par).nu=-((*par).a0*(*par).a0 - (*par).Rp*(*par).Rp)/((*par).a0*(*par).a0 + (*par).Rp*(*par).Rp);

	fclose(stream_config);
}

void PrintToFile(parameters par, ftype *X, const char *filename)
{
	ofstream fp;
	int idom, ipot, i, j, ns, nt, indx;

	fp.open(filename);
	fp << par.a0  << " \ta0"  << endl;
	fp << par.Rp  << " \tRp"  << endl;
	fp << par.Ktr << " \tKtr"  << endl;
	fp << par.rH1 << " \tr1"  << endl;
	fp << par.c1  << " \tc1"  << endl;
	fp << par.Sz1 << " \tSz1" << endl;
	fp << par.Pz1 << " \tPz1" << endl;
	fp << par.rH2 << " \trH2" << endl;
	fp << par.c2  << " \tc2"  << endl;
	fp << par.Sz2 << " \tSz2" << endl;
	fp << par.Pz2 << " \tPz2" << endl;
	fp << par.trapping_surface << " \ttrapping_surface" << endl;

	for (idom = 0; idom < NDOM; idom++)
	{
		fp << par.ns[idom] << "\t ns_domain_" << idom << endl;
	}
	fp << endl << par.nt << "\t nt" << endl;
	fp.precision(PREC_OUTPUT);
	fp << endl;

	for (idom = 0; idom < NDOM; idom++)
	{   // Writing the Potential values
		ns = par.ns[idom];
		nt = par.nt;
		for (ipot = 0; ipot < NPOT; ipot++)
		{
			for (j = 0; j < nt; j++)
			{
				for (i = 0; i < ns; i++)
				{
					indx = Index(par, idom, ipot, i, j);
					fp << X[indx] << "\t (idom,ipot,i,j)=("<<idom<<","<<ipot<<","<<i<<","<<j<<")"<<endl;
				}
			}
		}
	}
	fp.close();
}

void CreateInitialData(parameters par, ftype *X)
{
	int ns, nt;

	for (int idom = 0; idom < NDOM; idom++)
	{
		ns = par.ns[idom];
		nt = par.nt;
		for (int ipot = 0; ipot < NPOT; ipot++)
		{
			for (int j = 0; j < nt; j++)
			{
				for (int i = 0; i < ns; i++)
				{
					ftype aux;
					ftype rho = par.rho.d0[Index(par,idom,0,i,j)],
					z = par.z.d0[Index(par,idom,0,i,j)];
					if(0==par.coefficientswitch)
					{
						aux = sqrt(6/(par.Rp*par.Ktr));
					}
					else if(1==par.coefficientswitch)
					{
						aux = LaplaceSolution(par,idom , i,j) + cos(z);
					}
					else {exit(1);}
					X[Index(par,idom,ipot,i,j)]=aux;
				}
			}
		}
	}
	PrintToFile(par,X,"../run/solution/X_initialguess");
}

void PrintCheb(ftype *c, int length, const char *name)
{
	int i;
	ofstream fp;

	fp.open(name,std::fstream::in | std::fstream::out | std::fstream::trunc);
	fp.precision(PREC_OUTPUT);
	fp << "#" << name << endl;
	for (i = 0; i < length; i++)
	{
		fp << i << "\t" << c[i] << endl;
	}
	fp.close();
}

void Plot_v(parameters par,derivs_2D v, const char *name)
{
	int ns, nt;
	ftype rho,z;
	ofstream fp;

	fp.open(name);
	for(int idom=0;idom<NDOM;idom++)
	{
		ns = par.ns[idom];
		nt = par.nt;

		for(int i=0; i<ns;i++)
		{
			for(int j=0; j<nt; j++)
			{
				rho = par.rho.d0[Index(par,idom,0,i,j)];
				z = par.z.d0[Index(par,idom,0,i,j)];

				fp<< rho <<"\t\t" << z <<"\t\t" << v.d0[Index(par,idom,0,i,j)] << endl;
			}
			fp << endl;
		}
	}
	fp.close();
}

void Plot_data(parameters par,derivs_2D v, const char *name)
{
	int ns, nt;
	ftype rho,z, phi, omega, Rp=par.Rp;
	ofstream fp;

	fp.open(name);

	fp << "# a0=" << par.a0 << ", Rp=" << par.Rp <<", Ktr=" << par.Ktr
		 << ", rH1=" << par.rH1 << ", c1=" << par.c1 << ", Sz1=" << par.Sz1 << ", Pz1=" << par.Pz1 <<", r1=" << par.rH1
		 << ", c2=" << par.c2 << ", Sz2=" << par.Sz2 << ", Pz2=" << par.Pz2 <<", r2=" << par.rH2<< ", timestamp=" << par.timestamp << endl;

	for(int idom=0;idom<NDOM;idom++)
	{
		ns = par.ns[idom];
		nt = par.nt;

		for(int i=0; i<ns;i++)
		{
			for(int j=0; j<nt; j++)
			{
				//Die Division durch R+ ergibt einheitenlose Größen
				rho = par.rho.d0[Index(par,idom,0,i,j)];
				z = par.z.d0[Index(par,idom,0,i,j)];
				omega = 1 - (rho*rho + z*z)/(Rp*Rp);
				phi = v.d0[Index(par,idom,0,i,j)];

				fp<< rho/Rp <<"\t\t" << z/Rp <<"\t\t" << phi << " \t " << omega/sqr(phi) << endl;
			}
			fp << endl;
		}
	}
	fp.close();
}

void ScanFile(parameters *par, const char *name)
{
	int j_aux, ns, nt, indx, i_aux;
	ftype d_aux;
	char string[NMAX], number[NMAX];
	FILE *stream;

	stream = fopen(name, "r");

	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).a0 = d_aux;
	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).Rp = d_aux;
	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).Ktr = d_aux;

	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).rH1 = d_aux;
	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).c1 = d_aux;
	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).Sz1 = d_aux;
	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).Pz1 = d_aux;

	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).rH2 = d_aux;
	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).c2 = d_aux;
	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).Sz2 = d_aux;
	j_aux = fscanf(stream, "%s %s",  number, string);		d_aux=cast<ftype>(number);	(*par).Pz2 = d_aux;

	j_aux = fscanf(stream, " %d %s ", &i_aux, string);
	(*par).trapping_surface = i_aux;

	for (int idom = 0; idom < NDOM; idom++)
	{
		j_aux = fscanf(stream, " %d %s ", &ns, string);
		(*par).ns[idom] = ns;
	}
	j_aux = fscanf(stream, " %d %s ", &nt, string);
	(*par).nt = nt;

	(*par).nst[0] = 0;
	for (int idom = 0; idom < NDOM; idom++)
	{
		ns = (*par).ns[idom];
		(*par).nst[idom + 1] = (*par).nst[idom] + ns * nt;
	}
	//Berechnen der Anzahl an Punkten für alle Potential auf allen Gebieten
	(*par).ntotal=NPOT*(*par).nt*((*par).ns[0] +(*par).ns[1]);
	(*par).Scalarlength=((*par).ns[0] + (*par).ns[1]) * nt;
	(*par).n_2D=(*par).ntotal;

	(*par).eta1 = asinh((*par).a0/(*par).rH1);
	(*par).eta2 = -asinh((*par).a0/(*par).rH2);

	(*par).d1 = sqrt((*par).a0*(*par).a0 + (*par).rH1*(*par).rH1);
	(*par).d2 = -sqrt((*par).a0*(*par).a0 + (*par).rH2*(*par).rH2);

	(*par).nu=-((*par).a0*(*par).a0 - (*par).Rp*(*par).Rp)/((*par).a0*(*par).a0 + (*par).Rp*(*par).Rp);

	(*par).X_aux = new ftype[(*par).ntotal];

	for (int idom = 0; idom < NDOM; idom++)
	{
		ns = (*par).ns[idom];
		nt = (*par).nt;
		for (int ipot = 0; ipot < NPOT; ipot++)
		{
			for (int j = 0; j < nt; j++)
			{
				for (int i = 0; i < ns; i++)
				{
					indx = Index(*par, idom, ipot, i, j);
					j_aux = fscanf(stream, " %s %s ",  number, string);		d_aux=cast<ftype>(number);
					(*par).X_aux[indx] = d_aux;
				}
			}
		}
	}
	fclose(stream);
}

void PrintHorizon(parameters par, ftype *h, const char *name)
{
	ftype rho, z, mu;
	ofstream fp;

	fp.open(name);

	fp << "# a0=" << par.a0 << ", Rp=" << par.Rp <<", Ktr=" << par.Ktr
		 << ", r1=" << par.rH1 << ", c1=" << par.c1 << ", Sz1=" << par.Sz1 << ", Pz1=" << par.Pz1 <<", r1=" << par.rH1
		 << ", c2=" << par.c2 << ", Sz2=" << par.Sz2 << ", Pz2=" << par.Pz2 <<", r2=" << par.rH2<< ", timestamp=" << par.timestamp << endl;


	for(int j=0; j<par.nH;j++)
	{
		mu = -cos(Pi * j / cast<ftype>(par.nH-1));
		z = (h[j]*mu + par.d)/par.Rp;
		rho = (h[j]*sqrt(1-mu*mu))/par.Rp;

		fp << rho << "\t\t" << z << endl;
	}
	fp.close();
}

void Print_physical_data(parameters par, const char *name)
{
	ftype MB = par.MBondi, M_common, A_common, J_common, M_up, A_up, J_up, M_down, A_down, J_down,
		  pd = par.proper_distance, pc_common, pc_up, pc_down;
	ofstream fp;

	if(1==par.MOTS_common_exist)
	{
		M_common = par.Mass_C_common;
		A_common = par.Area_common;
		J_common = par.J_common;
		pc_common = par.proper_circumference_common;
	}
	else
	{
		M_common = NAN;
		A_common = NAN;
		J_common = NAN;
		pc_common = NAN;
	}

	if(1==par.MOTS_up_exist)
	{
		M_up = par.Mass_C_up;
		A_up = par.Area_up;
		J_up = par.J_up;
		pc_up = par.proper_circumference_up;
	}
	else
	{
		M_up = NAN;
		A_up = NAN;
		J_up = NAN;
		pc_up = NAN;
	}

	if(1==par.MOTS_down_exist)
	{
		M_down = par.Mass_C_down;
		A_down = par.Area_down;
		J_down = par.J_down;
		pc_down = par.proper_circumference_down;
	}
	else
	{
		M_down = NAN;
		A_down = NAN;
		J_down = NAN;
		pc_down = NAN;
	}

	fp.open(name);
	cout<<endl;
	cout << "timestamp: " <<par.timestamp<<endl;
	cout << "printing to: " << name<<endl;

	fp << par.timestamp <<" \t " << MB << " \t " << M_common <<" \t " << A_common <<" \t " << J_common << " \t " <<  M_up <<" \t " << A_up <<" \t " << J_up << " \t " << M_down <<" \t " << A_down <<" \t " << J_down << " \t " << pd << " \t " << pc_common <<" \t " << pc_up << " \t " << pc_down << endl;

	fp.close();
}

void PrintCoordinateLines(parameters par)
{
	ofstream fp;
	int ns,nt, idom;
	ftype rho,z, A, B, psi, kappa;
	char name[NMAX];

	//Gebiet 0
	idom = 0;
	//A=const
	for(A=0; A<=1; A+=0.1)
	{
		sprintf(name,"Koordinatenlinien_dom=%d_A=%3.2f",idom, float(A));
		fp.open(name);

		fp << "# z \t rho" << endl;

		for(B = -1; B<=1; B+=0.005)
		{
			// rho = chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,A,B);
			// z = chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,A,B);

			// fp<< z <<" \t " << rho << endl;

			psi = chebevxy(0,1,-1,1,par.psi.cheb[idom],par.ns[idom],par.nt,A,B);
			kappa = chebevxy(0,1,-1,1,par.kappa.cheb[idom],par.ns[idom],par.nt,A,B);

			fp<< psi <<" \t " << kappa << endl;
		}

		fp.close();
	}

	//B=const
	for(B=-1; B<=1; B+=0.1)
	{
		sprintf(name,"Koordinatenlinien_dom=%d_B=%3.2f",idom, float(B));
		fp.open(name);

		fp << "# z \t rho" << endl;

		for(A = 0; A<=1; A+=0.005)
		{
			// rho = chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,A,B);
			// z = chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,A,B);

			// fp<< z <<" \t " << rho << endl;

			psi = chebevxy(0,1,-1,1,par.psi.cheb[idom],par.ns[idom],par.nt,A,B);
			kappa = chebevxy(0,1,-1,1,par.kappa.cheb[idom],par.ns[idom],par.nt,A,B);

			fp<< psi <<" \t " << kappa << endl;
		}

		fp.close();
	}

	//Gebiet 1
	idom = 1;
	//A=const
	for(A=0; A<=1; A+=0.1)
	{
		sprintf(name,"Koordinatenlinien_dom=%d_A=%3.2f",idom, float(A));
		fp.open(name);

		fp << "# z \t rho" << endl;

		for(B = -1; B<=1; B+=0.005)
		{
			// rho = chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,A,B);
			// z = chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,A,B);

			// fp<< z <<" \t " << rho << endl;

			psi = chebevxy(0,1,-1,1,par.psi.cheb[idom],par.ns[idom],par.nt,A,B);
			kappa = chebevxy(0,1,-1,1,par.kappa.cheb[idom],par.ns[idom],par.nt,A,B);

			fp<< psi <<" \t " << kappa << endl;
		}

		fp.close();
	}

	//B=const
	for(B=-1; B<=1; B+=0.1)
	{
		sprintf(name,"Koordinatenlinien_dom=%d_B=%3.2f",idom, float(B));
		fp.open(name);

		fp << "# z \t rho" << endl;

		for(A = 0; A<=1; A+=0.005)
		{
			// rho = chebevxy(0,1,-1,1,par.rho.cheb[idom],par.ns[idom],par.nt,A,B);
			// z = chebevxy(0,1,-1,1,par.z.cheb[idom],par.ns[idom],par.nt,A,B);

			// fp<< z <<" \t " << rho << endl;

			psi = chebevxy(0,1,-1,1,par.psi.cheb[idom],par.ns[idom],par.nt,A,B);
			kappa = chebevxy(0,1,-1,1,par.kappa.cheb[idom],par.ns[idom],par.nt,A,B);

			fp<< psi <<" \t " << kappa << endl;
		}

		fp.close();
	}
}