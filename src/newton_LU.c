#include "main.h"

void Jacobian(parameters par, ftype *X, ftype **DF)
{
	int ntotal = par.ntotal;

	#pragma omp parallel for
	for(int j=0; j < ntotal; j++)
	{
		ftype *DX, *JDX;
		DX  = new ftype[ntotal];
		JDX = new ftype[ntotal];
		fill0_dvector(DX,   0, ntotal);
		fill0_dvector(JDX,  0, ntotal);

		DX[j] = 1;
		J_times_DX(par, X, DX, JDX);
		for(int l=0; l < ntotal; l++){
			DF[l][j] = JDX[l];
		}
		DX[j] = 0;
		delete[] DX;
		delete[] JDX;
	}
}

void DF_of_X(parameters par, ftype *X, ftype **DF)
{
	int ntotal = par.ntotal;
	ftype epDF = 1.e-08;

	#pragma omp parallel for
	for(int j=0; j < ntotal; j++)
	{
		ftype *Xp, *Fp, *Xm, *Fm;
		derivs_2D v;
		allocate_derivs_2D(&v, ntotal);
		Xp = new ftype[ntotal];
		Fp = new ftype[ntotal];
		Xm = new ftype[ntotal];
		Fm = new ftype[ntotal];

		for(int l=0; l < ntotal; l++)
		{
			Xp[l] = Xm[l] = X[l];
		}
		Xp[j] += epDF;
		Xm[j] -= epDF;
		F_of_X(par, Xp, v, Fp);
		F_of_X(par, Xm, v, Fm);
		for(int l=0; l < ntotal; l++)
		{
			DF[l][j] = 0.5*(Fp[l]-Fm[l])/epDF;
		}
		Xp[j] = Xm[j] = X[j];

		delete[] Xp;
		delete[] Xm;
		delete[] Fm;
		delete[] Fp;
		free_derivs_2D(&v);
	}
}

void newton_LU(parameters par, ftype *X)
{
	int ntotal=par.ntotal, *indx, j, l=par.Newton_itmin, lmax=par.Newton_itmax;
	ftype *F, d=1.e05, dmin = par.Newton_tol, dmax=1.e15, dd;
	derivs_2D v;
	double dt=0, dtJac=0, dtInv=0, t_start;

	F     = new ftype[ntotal];
	allocate_derivs_2D(&v, ntotal);

	indx  = ivector(0, ntotal-1);

	ftype **DF = new ftype*[ntotal];
	for (int i = 0; i < ntotal; i++)
	{
		DF[i] = new ftype[ntotal];
	}

	F_of_X(par, X, v, F);

	d=normalized_norm(F,ntotal);
	cout <<"\t initial norm"<<"\t d="<<d<<endl;

	while(l <= lmax && d > dmin && d < dmax)
	{

		if (d>0 && d<dmax)
		{
			t_start = omp_get_wtime();
			Jacobian(par, X, DF);
			DF_of_X(par, X, DF);

			dtJac = omp_get_wtime() - t_start;

			ludcmp(DF, ntotal, indx, &dd);
			lubksb(DF, ntotal, indx, F);
			for(j=0; j<ntotal; j++)
			{
				X[j] -=  F[j];
			}
			dtInv = omp_get_wtime() - t_start - dtJac;
			dt = omp_get_wtime() - t_start;
		}

		F_of_X(par, X, v, F);
		d=normalized_norm(F,ntotal);
		cout <<"\t l="<<l<<"\t d="<<d<<"\t in " << fixed << dt << "s\t(Jacobimatrix: "<< dtJac <<"s, Inversion: " << dtInv << "s)" << scientific << endl;
		PrintToFile(par, X,"./solution/X_running");

		if(boost::math::isinf(d) || boost::math::isnan(d) || d>dmax)
		{
//			cout <<"d="<<d << "... printing status quo and exit to system" << endl;
//			PrintToFile(par, X,"./solution/X_abort");
//			for(int i=0;i<ntotal;i++)
//			{
//				if(fabs(F[i])>1e10 || boost::math::isinf(F[i]) || boost::math::isnan(F[i]))
//				{
//					cout << "i="<< i <<"\tF=" << F[i] << endl;
//					InvertIndex(par,i);
//				}
//			}
			exit(1);
		}
		l++;
	}

	if(boost::math::isinf(d) || boost::math::isnan(d) || d>dmax)
	{
		cout <<"d="<<d << "... printing status quo and exit to system" << endl;
		PrintToFile(par, X,"./solution/X_abort");
		exit(1);
	}
	PrintToFile(par, F, "./test/F");
	delete[] F;
	free_derivs_2D(&v);

	free_ivector(indx, 0, ntotal-1);
	for (int i = 0; i < ntotal; i++)
	{
		delete[] DF[i];
	}
	delete[] DF;
}
