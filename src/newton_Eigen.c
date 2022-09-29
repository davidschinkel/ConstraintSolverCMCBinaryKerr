#include "main.h"

void newton_Eigen(parameters par, ftype *X)
{
	int ntotal=par.ntotal, l=1, lmax=par.Newton_itmax;
	ftype *F, *X_aux, d, dmin = par.Newton_tol, dmax=1.e5;
	derivs_2D v;
	Eigen::MatrixXd A_Eigen(ntotal,ntotal);
	Eigen::VectorXd b_Eigen(ntotal), x_Eigen(ntotal);
	double dt=0, dtJac=0, dtInv=0, t_start;

	F = new ftype[ntotal];
	X_aux = new ftype[ntotal];
	allocate_derivs_2D(&v, ntotal);

	ftype **DF = new ftype*[ntotal];
	for (int i = 0; i < ntotal; i++)
	{
		DF[i] = new ftype[ntotal];
	}
	Eigen::setNbThreads(4);

	F_of_X(par, X, v, F);

	d=normalized_norm(F,ntotal);
	cout <<"\t initial norm"<<"\t d="<<d<<endl;

	while(l <= lmax && d > dmin && d < dmax)
	{

		if (d>0 && d<dmax)
		{
			t_start = omp_get_wtime();
			DF_of_X(par, X, DF);
			for(int i = 0; i< ntotal;i++)
			{
				for(int j = 0; j< ntotal; j++)
				{
					A_Eigen(i,j)=DF[i][j];
				}
				b_Eigen(i)=-F[i];
			}
			dtJac = omp_get_wtime() - t_start;
			x_Eigen=A_Eigen.lu().solve(b_Eigen);

			for(int j=0; j<ntotal; j++)
			{
				X[j] +=  ftype(x_Eigen(j));
				//cout << "j="<<j<< "\tX="<<X[j]<<"\tF="<<F[j]<<endl;
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
			cout <<"d="<<d << "... printing status quo and exit to system" << endl;
			PrintToFile(par, X,"./solution/X_abort");
			for(int i=0;i<ntotal;i++)
			{
				if(fabs(F[i])>1e10 || boost::math::isinf(F[i]) || boost::math::isnan(F[i]))
				{
					cout << "i="<< i <<"\tF=" << F[i] << "\tx_Eigen="<<x_Eigen(i)  <<"\tX=" << X[i] << endl;
					InvertIndex(par,i);
				}
			}
			exit(1);
		}

		l++;
	}

	if(l>lmax)
	{
		cout << "Maximale Anzahl an Iterationsschritten erreicht, breche ab!" << endl;
		exit(1);
	}


	if(boost::math::isinf(d) || boost::math::isnan(d) || d>dmax)
	{
		cout <<"d="<<d << "... printing status quo and exit to system" << endl;
		PrintToFile(par, X,"./solution/X_abort");
		exit(1);
	}

	delete[] F;
	delete[] X_aux;
	free_derivs_2D(&v);

	for (int i = 0; i < ntotal; i++)
	{
		delete[] DF[i];
	}
	delete[] DF;
}
