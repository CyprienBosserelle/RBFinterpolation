#include <iostream>
#include <math.h>
#include <armadillo>
#include "RBF.h"
//#define pi 3.14159265

using namespace std;
using namespace arma;



double sign(double x)
{
	double y=0;

	if (abs(x) > datum::eps)
	{
		y = x / abs(x);
	}
	return y;
}



mat RBFgaussian(mat r, double gamma)
{
	//exp(-0.5 * r % r / (gamma * gamma));
	return exp(-0.5 * r % r / (gamma * gamma));
}

double calculateGamma(int ncenters, int ndim, mat x)
{
	
	//Not working
	mat ep=pow(prod(max(x)-min(x))/ncenters,1.0/((double)ndim));

	return (double)ep(0, 0);
}

/*double optimgamma(int ncenters, mat A0, mat RBFcoeff)
{
	//From:
	// Alex Chirokov, alex.chirokov@gmail.com
	// 16 Feb 2006
	mat A, invA;

	A = A0(span(0, ncenters - 1), span(0, ncenters - 1));
	invA = pinv(A);


}*/

double costeps(int ncenters, int ndim, double ep, mat centersnorm, mat data)
{
	//ep is epsilon which means gamma
	mat A(ncenters, ncenters);
	mat r, t, o, P, X, B, b, RBFcoeff;
	mat test, rr, s, f;

	A.zeros();

	for (int i = 0; i < ncenters; i++)
	{
		//
		for (int j = 0; j <= i; j++)
		{
			//
			r = norm(centersnorm.col(i) - centersnorm.col(j));

			t = exp(-0.5 * r % r / (ep * ep));// is 1 by 1 yeah?

			A(i, j) = t(0, 0);
			A(j, i) = t(0, 0);
		}
	}

	o = ones(ncenters, 1);

	P = join_rows(o, centersnorm.t());

	X = join_rows(A, P);

	B = join_rows(P.t(), zeros(ndim + 1, ndim + 1));
	A = join_cols(X, B);

	//A.save("A_mat.txt", raw_ascii);

	b = join_cols(data.t(), zeros(ndim + 1, 1));
	//b.save("b_mat.txt", raw_ascii);

	RBFcoeff = solve(A, b);

	mat iA = A(span(0, ncenters - 1), span(0, ncenters - 1));
	//iA.save("A1.txt", arma_ascii);
	mat invA = pinv(iA,datum::eps,"std");

	mat KK = data - RBFcoeff(ncenters);
	for (int i = 0; i < ndim; i++)
	{
		KK = KK - RBFcoeff(ncenters + i + 1)*centersnorm.row(i);
	}

	mat ceps = (invA*trans(KK)) / diagvec(invA);
	double yy = norm(ceps);

	return yy;


}

mat RBFtrain(int ncenters, int ndim, double gamma, mat centersnorm, mat data)
{
	// Calculate the RBFcoeff 
	mat A(ncenters, ncenters);
	mat r, t, o, P, X, B, b, RBFcoeff;
	mat test, rr, s, f;

	A.zeros();

	for (int i = 0; i < ncenters; i++)
	{
		//
		for (int j = 0; j <= i; j++)
		{
			//
			r = norm(centersnorm.col(i) - centersnorm.col(j));

			t = exp(-0.5 * r % r / (gamma * gamma));// is 1 by 1 yeah?

			A(i, j) = t(0, 0);
			A(j, i) = t(0, 0);
		}
	}

	o = ones(ncenters, 1);

	P = join_rows(o, centersnorm.t());

	X = join_rows(A, P);

	B = join_rows(P.t(), zeros(ndim + 1, ndim + 1));
	A = join_cols(X, B);

	//A.save("A_mat.txt", raw_ascii);

	b = join_cols(data.t(), zeros(ndim + 1, 1));
	//b.save("b_mat.txt", raw_ascii);

	RBFcoeff = solve(A, b);

	return RBFcoeff;

}


double RBFinterp(int ncenters, int ndim, double gamma, mat RBFcoeffw, mat centersnorm, mat pointsnorm)
{
	// perform the RBF interpolation using gamma, RBFcoeff the centers normalised and the new locations to interpolate to;
	

	mat rr, s, x;

	x = centersnorm;

	rr = zeros(1, ncenters);

	rr = (pointsnorm*ones(1, ncenters)) - x;

	rr = sqrt(sum(rr%rr, 0));

	for (int j = 0; j < ncenters; j++)
	{
		rr(j) = norm(pointsnorm - x.col(j));
	}

	


	mat tmp1, tmp2, tmp3;

	tmp1 = exp(-0.5 * rr % rr / (gamma * gamma));
	

	tmp2 = RBFcoeffw.rows(0, ncenters - 1);
	

	tmp3 = tmp2 % tmp1.t();
	

	s = RBFcoeffw(ncenters) + sum(tmp3, 0);
	

	for (int k = 0; k < ndim; k++)
	{
		s = s + RBFcoeffw(ncenters + k + 1)*pointsnorm(k, 0);
	}


	return (double)s(0, 0);
}


int
main(int argc, char** argv)
{

	//
	Param Param;


	vector<double> results;


	mat x, y;
	mat test;
	int twodee = 0;
	int nx = 1;
	int ny = 1;
	int nt; //


	mat RBFcoeff;
	cube RBFcoeffGrid;
	cube yGrid;
	double *xx = NULL;//used to reconstruct the netcdf files
	double *yy = NULL;//used to reconstruct the netcdf files
	double *theta; //used to reconstruct the netcdf files


	//////////////////////////////////////////////////////
	/////             Read Operational file          /////
	//////////////////////////////////////////////////////

	std::ifstream fs("RBF_param.txt");

	if (fs.fail()){
		std::cerr << "RBF_param.txt file could not be opened" << std::endl;
		//write_text_to_log_file("ERROR: XBG_param.txt file could not be opened...use this log file to create a file named XBG_param.txt");
		//SaveParamtolog(XParam);
		exit(1);
	}
	// Read and interpret each line of the XBG_param.txt
	std::string line;
	while (std::getline(fs, line))
	{
		//std::cout << line << std::endl;

		//Get param or skip empty lines
		if (!line.empty() && line.substr(0, 1).compare("#") != 0)
		{
			Param = readparamstr(line, Param);
			//std::cout << line << std::endl;
		}

	}
	fs.close();

	///////////////////////////////////////////////////////////////////////
	// CALCULATE GAMMA ?
	///////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////
	// SANITY CHECK ON THE PARAMETERS
	//////////////////////////////////////////////////////////////////////

	//sanity check
	if (Param.outputfile.empty())
	{
		Param.outputfile = "RBFoutput.txt";
	}

	if (Param.RBFcoefffile.empty() && Param.saveRBFcoeffs == 1)
	{
		Param.RBFcoefffile = "RBFcoeff.txt";
	}


	

	////////////////////////////////////////////////////////////////////
	// LOAD THE CENTERS OF THE RBF
	////////////////////////////////////////////////////////////////////

	
	//load centers
	x = readdatafile(Param.centersfile);
	
	Param.ndim = x.n_rows;

	Param.ncenters = x.n_cols;

	
	// Calculate min and max for normalisation of centers and training data
	std::vector<double> maxdimval, mindimval;
	

	for (int i = 0; i < Param.ndim; i++)
	{
		if (i != Param.isdir)
		{
			maxdimval.push_back(x(i, 0));
			mindimval.push_back(x(i, 0));
		}
		else
		{
			maxdimval.push_back(180.0/datum::pi);
			mindimval.push_back(0.0);
		}

	}
	
	
	
	for (int n = 0; n < Param.ncenters; n++)
	{
		for (int i = 0; i < Param.ndim; i++)
		{
			if (i != Param.isdir)
			{
				maxdimval[i] = max(x(i, n), maxdimval[i]);
				mindimval[i] = min(x(i, n), mindimval[i]);
			}

		}


	}
	// Normalise the centers to min max

	for (int n = 0; n < Param.ncenters; n++)
	{
		for (int i = 0; i < Param.ndim; i++)
		{
			x(i, n) = (x(i, n) - mindimval[i]) / (maxdimval[i] - mindimval[i]);
		}
	
	}

	////////////////////////////////////////////////////////////////////////////
	// TRAIN THE RBF (THIS IS OPTIONAL IF THE TRAINING HAS PRIORILY BEEN DONE)
	////////////////////////////////////////////////////////////////////////////

	
	if (Param.trainRBF == 1 && !Param.trainingfile.empty())
	{
		//first check if the input file is a nc file which impoly a training in 2D
		//if not then it is a qucick 1D stuff

		std::vector<std::string> extvec = split(Param.trainingfile, '.');

		std::string fileext = extvec.back();

		int strcmp = fileext.compare(0, 2, "nc");
		if (strcmp == 0)//IF 2D
		{
			twodee = 1;
		}

		if (twodee == 0)
		{
			//
			// load the training data (should be 1 column with ncenter lines)
			y = readdatafile(Param.trainingfile);

			///////////////////////////////////////////////////////////////////////
			// CALCULATE GAMMA ?
			///////////////////////////////////////////////////////////////////////
			if (Param.gamma<=0.0)
			{
				double ming = 0.001;
				double maxg = 0.900;
				double c = 0.5*(3.0 - sqrt(5.0));
				double a = ming;
				double b = maxg;
				double v = a + c*(b - a);
				double w = v;
				double xf = v;
				double d = 0.0;
				double e = 0.0;
				double xx = xf;
				double starteps = v;
				double fx, fu;
				fx= costeps(Param.ncenters, Param.ndim, starteps, x, y);

				double fv = fx;
				double fw = fx;
				double xm = 0.5*(a + b);
				double tol = 0.0001;
				double tol1 = sqrt(datum::eps)*abs(xf) + tol / 3.0;
				double tol2 = 2.0*tol1;

				double r, p, q, si;

				int gs = 1; //Parabolic fit by default
				int iter = 0;
				while (abs(xf - xm) > (tol2 - 0.5*(b - a)))
				{
					gs = 1;
					//Test if a parabolic fit is possible
					if (abs(e) > tol1)
					{
						//Parabolic fit
						gs = 0;
						r = (xf - w)*(fx - fv);
						q = (xf - v)*(fx - fw);
						p = (xf - v)*q - (xf - w)*r;
						q = 2.0*(q - r);
						if (q > 0.0)
						{
							p = -1.0 * p;
						}
						q = abs(q);
						r = e;  e = d;
						//Is the parabola acceptable
						if ((abs(p) < abs(0.5*q*r)) && (p > q*(a - xf)) && (p < q*(b - xf)))
						{
							d = p / q;
							xx = xf + d;
							// f must not be evaluated too close to ax or bx
							if (((xx - a) < tol2) || ((b - xx) < tol2))
							{
								si = sign(xm) + (((xm - xf) == 0.0) ? 1.0 : 0.0);
								d = tol1*si;

							}
							

						}
						else
						{
							gs = 1;
						}
					}
					if (gs == 1)
					{
						//A golden-section step is required
						if (xf >= xm)
						{
							e = a - xf;
						
						}
						else
						{
							e = b - xf;
						}

						d = c*e;
						
					}
					//The function must not be evaluated too close to xf
					si = sign(d) + ((d == 0) ? 1.0 : 0.0);
					xx = xf + si * max(abs(d), tol1);
					fu = costeps(Param.ncenters, Param.ndim, xx, x, y);

					//Update a, b, v, w, x, xm, tol1, tol2
					if (fu <= fx)
					{
						if (xx >= xf)
						{
							a = xf;
						}
						else
						{
							b = xf;
						}
						v = w;
						fv = fw;
						w = xf;
						fw = fx;
						xf = xx;
						fx = fu;
					}
					else
					{
						if (xx < xf)
						{
							a = xx;
						}
						else
						{
							b = xx;
						}
						if ((fu <= fw) || (w == xf))
						{
							v = w;
							fv = fw;
							w = xx;
							fw = fu;
						}
						else
						{
							if ((fu <= fv) || (v == xf) || (v == w))
							{
								v = xx;
								fv = fu;
							}
						}
						xm = 0.5*(a + b);
						tol1 = sqrt(datum::eps)*abs(xf) + tol / 3.0;
						tol2 = 2.0*tol1;
					}
					iter++;

					//Need a Sanity break here
					if (iter > 200)
					{
						break;
					}
				} //end while
				

				Param.gamma = fx;
			}



			//train the RBF
			RBFcoeff = RBFtrain(Param.ncenters, Param.ndim, Param.gamma, x, y);

			if (Param.saveRBFcoeffs == 1)
			{
				//Convert mat to vector
				std::vector<double> RBFcoeffvec;

				for (int n = 0; n < RBFcoeff.n_rows; n++)
				{
					RBFcoeffvec.push_back(RBFcoeff(n));
				}
				//write data file
				writedatafile(RBFcoeffvec, Param.RBFcoefffile);
			}
		}
		else
		{
			// 2D case
			//read grid size
			readgridncsize(Param.trainingfile, nx, ny, nt);

			xx = (double *)malloc(nx*sizeof(double));
			yy = (double *)malloc(ny*sizeof(double));
			readxync(Param.trainingfile, xx, yy);

			// init RBFCoeffGrid
			//mat RBFcoeffGrid = zeros(nx, ny, nt);
			yGrid = read3Dnc(Param.trainingfile, nx, ny, nt);
			RBFcoeffGrid = RBFcoeffGrid.zeros(Param.ncenters + Param.ndim + 1, ny, nx);
			for (int xi = 0; xi < nx; xi++)
			{
				for (int yi = 0; yi < ny; yi++)
				{
					y = yGrid(span(), span(yi, yi), span(xi, xi));
					RBFcoeff = RBFtrain(Param.ncenters, Param.ndim, Param.gamma, x, y.t());
					RBFcoeffGrid(span(), span(yi, yi), span(xi, xi)) = RBFcoeff;

				}
			}
			
			if (Param.saveRBFcoeffs == 1)
			{
				int ntheta = (Param.ncenters + Param.ndim + 1);
				double *RBFtrained2d;
				
				theta = (double *)malloc(ntheta * sizeof(double));
				RBFtrained2d = (double *)malloc(ntheta*nx*ny*sizeof(double));


				for (int n = 0; n < ntheta; n++)
				{
					theta[n] = n;
				}
				for (int xi = 0; xi < nx; xi++)
				{
					for (int yi = 0; yi < ny; yi++)
					{
						for (int ti = 0; ti < ntheta; ti++)
						{
							//
							RBFtrained2d[xi + yi*nx + ti*ny*nx] = RBFcoeffGrid(ti, yi, xi);
						}
					}
				}

				// Write 3d netcdf
				create3dnc(Param.RBFcoefffile, nx, ny, ntheta,  xx, yy, theta, RBFtrained2d);
				free(RBFtrained2d);
				free(theta);
			}

		}
		
		

		
	}
	else //No training i.e. we already have RBF coeffiscient in a input
	{
		// First look at the size of the grid 
		
		std::vector<std::string> extvec = split(Param.RBFcoefffile, '.');

		std::string fileext = extvec.back();

		int strcmp = fileext.compare(0,2,"nc");// this is safer than comparing the whole string which can have trhe ?varname attached to it
		if (strcmp == 0)//IF 2D
		{
			twodee = 1;
			//read grid size
			readgridncsize(Param.RBFcoefffile, nx, ny, nt);

			xx = (double *)malloc(nx*sizeof(double));
			yy = (double *)malloc(ny*sizeof(double));
			readxync(Param.RBFcoefffile, xx, yy);
			// init RBFCoeffGrid
			//mat RBFcoeffGrid = zeros(nx, ny, nt);
			RBFcoeffGrid = read3Dnc(Param.RBFcoefffile, nx, ny, nt);

			
			
		}
		else //1D
		{
			// Then it must be loaded
			RBFcoeff = readdatafile(Param.RBFcoefffile);
			RBFcoeff = RBFcoeff.t(); // Because the subroutine reads this the wrong way around...
		}
		
	}

	if (Param.interpRBF == 1 && !Param.inputfile.empty())
	{
		
		//Load the test data
		test = readdatafile(Param.inputfile);
		//test = test.t();
		// Normalise the data to the centers max
		for (int n = 0; n < test.n_cols; n++)
		{
			for (int i = 0; i < Param.ndim; i++)
			{
				test(i, n) = (test(i, n) - mindimval[i]) / (maxdimval[i] - mindimval[i]);
			}
			//test(0, n) = (test(0, n) - minHs) / (maxHs - minHs);
			//test(1, n) = (test(1, n) - minT) / (maxT - minT);
			//test(2, n) = test(2, n)*datum::pi / 180.0;
			//test(3, n) = (test(3, n) - minNM) / (maxNM - minNM);
		}

		if (twodee == 0)
		{
			//In this case 
			// Do the interpolation

			for (int n = 0; n < test.n_cols; n++)
			{
				results.push_back(RBFinterp(Param.ncenters, Param.ndim, Param.gamma, RBFcoeff, x, test.col(n)));

			}

			// write data file
			writedatafile(results, Param.outputfile);
		}
		else //2D case
		{
			//
			//cube resultsTD = zeros(RBFcoeffGrid.n_cols, RBFcoeffGrid.n_slices, test.n_cols);
			//
			double * results2d;

			results2d = (double *)malloc(test.n_cols*nx*ny*sizeof(double));
			
			theta = (double *)malloc(test.n_cols*sizeof(double));

			

			for (int n = 0; n < test.n_cols; n++)
			{
				theta[n] = n;
			}

			for (int xi = 0; xi < nx; xi++)
			{
				for (int yi = 0; yi < ny; yi++)
				{
					for (int n = 0; n < test.n_cols; n++)
					{
						results2d[xi+yi*nx+n*ny*nx]=RBFinterp(Param.ncenters, Param.ndim, Param.gamma, RBFcoeffGrid(span(),span(yi,yi),span(xi,xi)), x, test.col(n));

					}
				}
			}
			// Write 3d netcdf
			create3dnc(Param.outputfile, nx, ny, test.n_cols, xx, yy, theta, results2d);

			free(results2d);
			free(xx); free(yy); free(theta);
		}
	}
	


	return 0;
}
