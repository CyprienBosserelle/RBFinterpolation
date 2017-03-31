#include <iostream>
#include <math.h>
#include <armadillo>
#include "RBF.h"
//#define pi 3.14159265

using namespace std;
using namespace arma;







mat RBFgaussian(mat r, double gamma)
{
	//exp(-0.5 * r % r / (gamma * gamma));
	return exp(-0.5 * r % r / (gamma * gamma));
}

double calculateGamma(int ncenters, int dims, mat r)
{
	double variance = -1;
	double gamma;
	for (int i = 0; i < ncenters; ++i)
	{
		for (int j = i + 1; j < ncenters; ++j)
		{
			variance = max(variance, (double)r(i,0) * (double)r(j,0)); // Need to check
		}
	}
		
	variance *= (1.0 / ((double)ncenters));
	gamma = (-1.0 / (2.0 * variance));

	return gamma;
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

mat RBFtrain(int ncenters, int dim, double gamma, mat centersnorm, mat data)
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

	B = join_rows(P.t(), zeros(dim + 1, dim + 1));
	A = join_cols(X, B);

	//A.save("A_mat.txt", raw_ascii);

	b = join_cols(data.t(), zeros(dim + 1, 1));
	//b.save("b_mat.txt", raw_ascii);

	RBFcoeff = solve(A, b);

	return RBFcoeff;

}


double RBFinterp(int ncenters, int dim, double gamma, mat RBFcoeff, mat centersnorm, mat pointsnorm)
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
	

	tmp2 = RBFcoeff.rows(0, ncenters - 1);
	

	tmp3 = tmp2 % tmp1.t();
	

	s = RBFcoeff(ncenters) + sum(tmp3, 0);
	

	for (int k = 0; k < dim; k++)
	{
		s = s + RBFcoeff(ncenters + k + 1)*pointsnorm(k, 0);
	}


	return (double)s(0, 0);
}


int
main(int argc, char** argv)
{

	//

	mat x, y;

	double gamma = 0.2672;

	int ncenters;

	int dim;



	vector<double> results;

	// Set up variables
	//x = zeros(dim, ncenters);
	


	//load centers
	x = readdatafile("test_MDA_300.dat");
	
	dim = x.n_rows;

	ncenters = x.n_cols;

	//y = zeros(1, ncenters);
	//std::ifstream fs;
	//int nline = 0;
	//std::string line;
	//std::vector<std::string> lineelements;
	

	// Calculate min and max for normalisation of centers and training data
	double maxHs = x(0, 1);
	double maxT = x(1, 1);
	double maxNM = x(3, 1);

	double minHs = x(0, 1);
	double minT = x(1, 1);
	double minNM = x(3, 1);
	
	for (int n = 0; n < ncenters; n++)
	{
		maxHs = max(x(0, n), maxHs);
		maxT = max(x(1, n), maxT);
		maxNM = max(x(3, n), maxNM);

		minHs = min(x(0, n), minHs);
		minT = min(x(1, n), minT);
		minNM = min(x(3, n), minNM);
	}
	// Normalise the centers to min max

	for (int n = 0; n < ncenters; n++)
	{
		x(0, n) = (x(0, n) - minHs) / (maxHs - minHs);
		x(1, n) = (x(1, n) - minT) / (maxT - minT);
		x(2, n) = x(2, n)*datum::pi/180.0;
		x(3, n) = (x(3, n) - minNM) / (maxNM - minNM);
	}


	// load the training data (should be 1 column with ncenter lines)
	y = readdatafile("MauiBay_Shore_MDA_zs.txt");
	

	
	//train the RBF
	mat RBFcoeff;
	RBFcoeff = RBFtrain(ncenters, dim, gamma, x, y);
	mat test;
	
	mat readline = zeros(dim, 1);
	
	//Load the test data
	test = readdatafile("Test_data_MB_Shore.txt");

	
	

	// Normalise the data to the centers max
	for (int n = 0; n < test.n_cols; n++)
	{
		test(0, n) = (test(0, n) - minHs) / (maxHs - minHs);
		test(1, n) = (test(1, n) - minT) / (maxT - minT);
		test(2, n) = test(2, n)*datum::pi / 180.0;
		test(3, n) = (test(3, n) - minNM) / (maxNM - minNM);
	}

	// Do the interpolation
	
	for (int n = 0; n < test.n_cols; n++)
	{
		results.push_back(RBFinterp(ncenters, dim, gamma, RBFcoeff, x, test.col(n)));
		
	}



	return 0;
}
