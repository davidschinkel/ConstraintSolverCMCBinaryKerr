#include "main.h"

int main(int argc, char **argv)
{
	int i_Seq, check, idom, i, j;
	parameters parI, parG, par;
	ftype *X, h, h0, h_old;
	time_t timestamp;
	double t_start;
	char name[NMAX];
	derivs_2D v;

	cout <<" *******************************************************************" << endl;
	cout <<" *                       CMC_Data Buchman2009                      *" << endl;
	cout <<" *******************************************************************" << endl;

	parI.coefficientswitch = 0; //Hamiltonzwangsbedingung
	//parI.coefficientswitch = 1; //Laplacegleichung
	parG.coefficientswitch = par.coefficientswitch = parI.coefficientswitch;

	parG.trapping_surface = +1;  //+ : outgoing, - ingoing
	par.trapping_surface = parG.trapping_surface;

	if(0==par.coefficientswitch)
	{
		cout <<"\tGleichungssystem: Hamiltonzwangsbedingung" << endl;
		if(1==par.trapping_surface)
		{
			cout <<"\tRandbedingung: outgoing"<<endl;
		}
		if(-1==par.trapping_surface)
		{
		cout <<"\tRandbedingung: ingoing" << endl;
		}
	}
	if(1==par.coefficientswitch)
	{
		cout << "\tGleichungssystem: Laplacegleichung" << endl;
	}

	par.verbose = parG.verbose = parI.verbose = 0;

	par.h0_old = -1;

	//lesen der Anfangsdaten und der Config
	ScanConfig(&parG);
	ScanConfig(&par);

	omp_set_num_threads(par.num_threads);

	if(argc > 1)
	{
		ScanFile(&parI,argv[1]);
		cout <<"\tDatei gelesen, gefundene Konfiguration" << endl;
		cout << "\t a0="<<parI.a0 << " \tRp="<<parI.Rp<<" \tK="<<parI.Ktr<<endl;
		cout << "\t r1="<<parI.rH1 << "\t c1="<<parI.c1 << " \tSz1="<<parI.Sz1<<" \tPz1="<<parI.Pz1<<endl;
		cout << "\t r2="<<parI.rH2 << "\t c2="<<parI.c2 << " \tSz2="<<parI.Sz2<<" \tPz2="<<parI.Pz2<<endl;
		cout << "\t ns[0]="<<parI.ns[0] << " \tns[1]="<<parI.ns[1]<<" \tnt="<<parI.nt<<endl;
		cout << "\t trapping_surface " << par.trapping_surface << endl;
		X=new ftype[parG.ntotal];
		if((parI.nt == parG.nt) && (parI.ns[0] == parG.ns[0]) && (parI.ns[1] == parG.ns[1]))
		{
			copy_dvector(X,parI.X_aux,0,parG.ntotal);
			delete[] parI.X_aux;
		}
		else
		{
			ChangeResolution_X(parI, parG, parI.X_aux, X);
			delete[] parI.X_aux;
		}
	}
	else
	{
		ScanConfig(&parI);
		parI.c1=parI.c2=parI.Sz1=parI.Sz2=parI.Pz1=parI.Pz2=0;
		X=new ftype[parI.ntotal];
	}

	//Erzeugen der FFTW-Pläne
	t_start = omp_get_wtime();
	CreateFFTWplansLobatto(&par);
	cout << "Created FFTW-Plan in " << omp_get_wtime() - t_start << "s" << endl;

	//erzeugen der Anfangsdaten
	if(1==argc)
	{
		CreateFFTWplansLobatto(&parI);
		Getr(&parI);
		#pragma omp parallel sections
		{
			#pragma omp section
			{Get_psi(&parI);}
			#pragma omp section
			{Get_kappa(&parI);}
		}
		#pragma omp parallel sections
		{
			#pragma omp section
			{Get_rho(&parI);}
			#pragma omp section
			{Get_z(&parI);}
		}
		Get_Coefficients(&parI);
		ChebyshevTestH(parI);

		// cout << "plotting coordinates" << endl;
		// PrintCoordinateLines(parI);
		// exit(1);

		CreateInitialData(parI,X);

		//Berechnung der Lösung mit c=Sz=Pz=0
		cout << "\nCalculation of initial configuration with c=Sz=Pz=0" << endl;
		cout << " -----------------------------------------------------"<< endl;
		parI.Newton_tol=parI.Newton_tol_finish;
		cout << scientific;
		newton_Eigen(parI,X);
		cout << fixed << "\n";

		DestroyFFTWplansLobatto(&parI);
	}

	//Newton-Raphson-Löser
	for (i_Seq = 0; i_Seq <= par.n_Seq; i_Seq++)
	{
		//Berechnung der Schrittweite h
		if (par.n_Seq > 0)
			h = i_Seq / cast<ftype>(par.n_Seq);
		else
			h = 0;

		cout << " ---------------------------------------------------------------" << endl;
		cout << " Calculation of Sequence element No. \t"<<i_Seq<<"\th="<<h<< endl;
		cout << " ---------------------------------------------------------------" << endl;

		par.a0=(1-h)*parI.a0 + h*parG.a0;
        par.Ktr=(1-h)*parI.Ktr + h*parG.Ktr;
        par.rH1=(1-h)*parI.rH1 + h*parG.rH1;
        par.rH2=(1-h)*parI.rH2 + h*parG.rH2;
		par.c1=(1-h)*parI.c1 + h*parG.c1;
		par.c2=(1-h)*parI.c2 + h*parG.c2;
		par.Sz1=(1-h)*parI.Sz1 + h*parG.Sz1;
		par.Sz2=(1-h)*parI.Sz2 + h*parG.Sz2;
		par.Pz1=(1-h)*parI.Pz1 + h*parG.Pz1;
		par.Pz2=(1-h)*parI.Pz2 + h*parG.Pz2;
		//par.trapping_surface=(1-h)*parI.trapping_surface + h*parG.trapping_surface;

		par.Newton_tol = par.Newton_tol_finish;

		par.eta1 = asinh(par.a0/par.rH1);
		par.eta2 = -asinh(par.a0/par.rH2);
		par.d1 = sqrt(par.a0*par.a0 + par.rH1*par.rH1);
		par.d2 = -sqrt(par.a0*par.a0 + par.rH2*par.rH2);
		par.nu=-(par.a0*par.a0 - par.Rp*par.Rp)/(par.a0*par.a0 + par.Rp*par.Rp);

		cout << "\t a0="<<par.a0 << " \tRp="<<par.Rp<<" \tK="<<par.Ktr<<endl;
		cout << "\t r1="<<par.rH1 << "\t c1="<<par.c1 << " \tSz1="<<par.Sz1<<" \tPz1="<<par.Pz1<<endl;
		cout << "\t r2="<<par.rH2 << "\t c2="<<par.c2 << " \tSz2="<<par.Sz2<<" \tPz2="<<par.Pz2<<endl;
		cout << "\t trapping_surface " << par.trapping_surface << endl;

		Getr(&par);
		#pragma omp parallel sections
		{
			#pragma omp section
			{Get_psi(&par);}
			#pragma omp section
			{Get_kappa(&par);}
		}
		#pragma omp parallel sections
		{
			#pragma omp section
			{Get_rho(&par);}
			#pragma omp section
			{Get_z(&par);}
		}
		Get_Coefficients(&par);
		ChebyshevTestH(par);

		//Newton-Verfahren
		cout << scientific;
		//newton_LU(par, X);
		newton_Eigen(par,X);
		//newton(par, X);
		cout << fixed;

		timestamp = time(NULL);
		par.timestamp = timestamp;
		cout <<"timestamp=" << timestamp << endl;
		allocate_derivs_2D(&v, par.ntotal);
		allocate_derivs_2D_3rd_derivatives(&v,par.ntotal);
		Get_v_From_X(par,X,v);
		Derivatives_2D_3(par,v);
		Get_Chebyshev_Coefficients_v(par,&v,3);

		srand(time(NULL));
		j = (rand() % (par.nt-2)) + 1 ;
		srand(j*time(NULL));
		idom = rand() % 2;
		srand(2*j*time(NULL));
		i = (rand() % (par.ns[idom]-2)) + 1;
		solution_test(par,v,idom,i,j);

		ChebyshevTestv(par,v,timestamp);
		sprintf(name,"./solution/X_%ld",timestamp);
		PrintToFile(par, X, name);
		sprintf(name,"../run/plots/plot_data_%ld",timestamp);
		Plot_data(par,v,name);

		BondiMass(&par,v);
		proper_distance(&par,v);

		par.nH=par.ns[1];
		par.MOTS_common = new ftype[par.nH];
		par.MOTS_up = new ftype[par.nH];
		par.MOTS_down = new ftype[par.nH];
		par.MITS_common = new ftype[par.nH];
		par.MITS_up = new ftype[par.nH];
		par.MITS_down = new ftype[par.nH];
		par.MOTS_down_exist=par.MOTS_up_exist=par.MOTS_common_exist=0;

		//search for common horizon
		par.d=0;
		check = guess_common_horizon(par,h0, v);

		if(-1==check)
		{//es wurde keine gemeinsamer Horizont abgeschätzt
			if(par.h0_old > maximum(fabs(par.d1 + par.rH1),fabs(par.d2-par.rH2)))
			{//wenn es einen alten Horizont gibt und dieser nicht in den Schwarzen Löchern liegt
				h0=par.h0_old;
			}
			else
			{//wenn keinen alten Horizont, oder der alte Horizont unbrauchbar ist
				h0= maximum(fabs(par.d1 + par.rH1),fabs(par.d2-par.rH2)) + 0.1*(par.Rp - maximum(fabs(par.d1 + par.rH1),fabs(par.d2-par.rH2)));
			}
		}

		cout << "par.d=" << par.d << endl;
		check = HorizonFinder(&par,v, h0, h_old);
		if(0==check)
		{
			par.h0_old = h_old;
			cout << "\nCommon horizon found!\n" << endl;
			par.MOTS_common_exist = 1;
			MITS_Finder(par,v);
		}
		else
		{
			cout << "\n No common horizon found!\n" << endl;
			MITS_Finder(par, v);
			cout << "\nSearching for horizon around Black Hole 1" << endl;
			par.d = par.d1;
			cout << "par.d=" << par.d << endl;
			//h0: Das Minimum der folgenden Möglichkeiten: Abstand zu Scri, Abstand zum nächsten Schwarzen Loch, der doppelte Horizontradius
			h0 = minimum(1.5*par.rH1, minimum(par.Rp - fabs(par.d1) - par.rH1, (fabs(par.d1) - par.rH1) + (fabs(par.d2) - par.rH2))/2. + par.rH1 );
			//h0 = par.rH1;
			check = HorizonFinder(&par,v, h0, h_old);
			if(0==check)
			{
				cout << "found horizon" << endl;
				par.MOTS_up_exist = 1;
				MITS_Finder(par, v);
			}
			else
			{
				cout << "problem finding single horizon around Black Hole 1" << endl;
				cout << "the setup was:" << endl;
				cout << "\t a0="<<par.a0 <<  " \t Rp="<<par.Rp<<  " \t K="<<par.Ktr<<endl;
				cout << "\t r1="<<par.rH1 << " \t c1="<<par.c1 << " \t Sz1="<<par.Sz1<<" \t Pz1="<<par.Pz1<<endl;
				cout << "\t r2="<<par.rH2 << " \t c2="<<par.c2 << " \t Sz2="<<par.Sz2<<" \t Pz2="<<par.Pz2<<endl;
				cout << "\t ns[0]="<<par.ns[0] << " \tns[1]="<<par.ns[1]<<" \tnt="<<par.nt<<endl;
				cout << "i_Seq=" << i_Seq << "\ttimestamp=" << timestamp <<  endl;
				free_derivs_2D(&v);
				free_derivs_2D_3rd_derivatives(&v);
				free_Chebyshev_Coefficients_v(par,&v,3);
				Free_Coefficients(&par);
				delete[] X;
				DestroyFFTWplansLobatto(&par);
				cout << "exit to system" << endl;
				exit(1);
			}

			cout << "\nSearching for horizon around Black Hole 2" << endl;
			par.d = par.d2;
			cout << "par.d=" << par.d << endl;
			//h0: Das Minimum der folgenden Möglichkeiten: Abstand zu Scri, Abstand zum nächsten Schwarzen Loch, der doppelte Horizontradius
			h0 = minimum(1.5*par.rH2, minimum(par.Rp - fabs(par.d2) - par.rH2, (fabs(par.d2) - par.rH2) + (fabs(par.d1) - par.rH1))/2. + par.rH2 );
			//h0 = par.rH2;
			check = HorizonFinder(&par,v, h0, h_old);
			if(0==check)
			{
				cout << "found horizon" << endl;
				par.MOTS_down_exist = 1;
				MITS_Finder(par, v);
			}
			else
			{
				cout << "problem finding single horizon around Black Hole 2" << endl;
				cout << "the setup was:" << endl;
				cout << "\t a0="<<par.a0 << " \tRp="<<par.Rp<<" \tK="<<par.Ktr<<endl;
				cout << "\t r1="<<par.rH1 << "\t c1="<<par.c1 << " \tSz1="<<par.Sz1<<" \tPz1="<<par.Pz1<<endl;
				cout << "\t r2="<<par.rH2 << "\t c2="<<par.c2 << " \tSz2="<<par.Sz2<<" \tPz2="<<par.Pz2<<endl;
				cout << "\t ns[0]="<<par.ns[0] << " \tns[1]="<<par.ns[1]<<" \tnt="<<par.nt<<endl;
				cout << "i_Seq=" << i_Seq << "\ttimestamp=" << timestamp << endl;
				free_derivs_2D(&v);
				free_derivs_2D_3rd_derivatives(&v);
				free_Chebyshev_Coefficients_v(par,&v,3);
				Free_Coefficients(&par);
				delete[] X;
				DestroyFFTWplansLobatto(&par);
				cout << "exit to system" << endl;
				exit(1);
			}

		}

		Average_at_scri(par,v);

		free_derivs_2D(&v);
		free_derivs_2D_3rd_derivatives(&v);
		free_Chebyshev_Coefficients_v(par,&v,3);

		//Ausgabe der physikalischen Größen
		if(1==par.MOTS_common_exist)
		{
			cout << "gemeinsamer Horizont: M_Bondi=" << par.MBondi <<"\t M_C=" << par.Mass_C_common <<"\t M_irr=" << par.Mass_irr_common <<"\tA=" <<par.Area_common << "\tJ_common=" << par.J_common << endl;
			cout << "\t\t M_B/R_+= " << par.MBondi/par.Rp <<"\t M_C/R_+=" << par.Mass_C_common/par.Rp <<"\t M_irr/R_+=" << par.Mass_irr_common/par.Rp <<"\tA/R_+^2=" <<par.Area_common/(par.Rp*par.Rp) << "\tj_B=" << par.J_common/pow(par.MBondi,2) << "\tj_C=" << par.J_common/pow(par.Mass_C_common,2)<< "\tj_irr=" << par.J_common/pow(par.Mass_irr_common,2) << endl;		\
			cout << "eps_AB=" << par.Area_common/(8*Pi*(pow(par.MBondi,2) + sqrt(pow(par.MBondi,4) - pow(par.J_common,2))));
		}
		if(1==par.MOTS_up_exist)
		{
			cout << "oberer Horizont: M_Bondi=" << par.MBondi <<"\t M_C=" << par.Mass_C_up <<"\t M_irr=" << par.Mass_irr_up <<"\tA=" <<par.Area_up << "\tJ=" << par.J_up << endl;
			cout << "\t\t M_B/R_+= " << par.MBondi/par.Rp <<"\t M_C/R_+=" << par.Mass_C_up/par.Rp <<"\t M_irr/R_+=" << par.Mass_irr_up/par.Rp <<"\tA/R_+^2=" <<par.Area_up/(par.Rp*par.Rp) << "\tj_B=" << par.J_up/pow(par.MBondi,2) << "\tj_C=" << par.J_up/pow(par.Mass_C_up,2)<< "\tj_irr=" << par.J_up/pow(par.Mass_irr_up,2) << "\tJ=" << par.J_up << endl;
		}
		if(1==par.MOTS_down_exist)
		{
			cout << "unterer Horizont: M_Bondi=" << par.MBondi <<"\t M_C=" << par.Mass_C_down <<"\t M_irr=" << par.Mass_irr_down <<"\tA=" <<par.Area_down <<"\tJ=" << par.J_down << endl;
			cout << "\t\t M_B/R_+= " << par.MBondi/par.Rp <<"\t M_C/R_+=" << par.Mass_C_down/par.Rp <<"\t M_irr/R_+=" << par.Mass_irr_down/par.Rp <<"\tA/R_+^2=" <<par.Area_down/(par.Rp*par.Rp) << "\tj_B=" << par.J_down/pow(par.MBondi,2) << "\tj_C=" << par.J_down/pow(par.Mass_C_down,2)<< "\tj_irr=" << par.J_down/pow(par.Mass_irr_down,2) << "\tJ=" << par.J_up << endl;
		}
		if(2==(par.MOTS_up_exist + par.MOTS_down_exist))
		{
			cout << "eps_AB=" << (par.Area_up + par.Area_down)/(8*Pi*(pow(par.MBondi,2) + sqrt(pow(par.MBondi,4) - pow(par.J_up + par.J_down,2)))) << endl;	
			cout << "MB/(MC1 + MC2) - 1=" << par.MBondi/(par.Mass_C_up + par.Mass_C_down) - 1 << " \t MB/(Mirr1 + Mirr2) - 1=" << par.MBondi/(par.Mass_irr_up + par.Mass_irr_down) - 1 << endl;
		}
		

		sprintf(name,"./physical_data/data_%ld",timestamp);
		Print_physical_data(par, name);

		Free_Coefficients(&par);
	}

	//Vergleichen mit der vorgegebenen Lösung
	if(1==par.coefficientswitch)
	{
		Check_Laplace(par,X);
	}

	delete[] X;
	//DestroyFFTWplansLobatto(&par);

	return 0;
}
