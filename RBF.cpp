#include <iostream>
#include <math.h>
#include <armadillo>

#define pi 3.14159265

using namespace std;
using namespace arma;




void split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss;
	ss.str(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		if (!item.empty())//skip empty tokens
		{
			elems.push_back(item);
		}

	}
}


std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

std::string trim(const std::string& str, const std::string& whitespace)
{
	const auto strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos)
		return ""; // no content

	const auto strEnd = str.find_last_not_of(whitespace);
	const auto strRange = strEnd - strBegin + 1;

	return str.substr(strBegin, strRange);
}




int
main(int argc, char** argv)
{

	//

	mat x, y;

	double cst = 0.2672;

	int ncenters = 300;

	int dim = 4;



	vector<double> results;


	x = zeros(dim, ncenters);
	y = zeros(1, ncenters);

	std::ifstream fs;


	fs.open("test_MDA_300.dat");

	if (fs.fail()){
		std::cerr << "test_MDA_300.dat file could not be openned" << std::endl;
	}

	std::string line;
	std::vector<std::string> lineelements;
	int nline=0;
	while (std::getline(fs, line))
	{
		if (!line.empty())
		{
			lineelements = split(line, '  ');
			if (lineelements.size() < 4)
			{
				std::cerr <<  "ERROR test_MDA_300.dat file format error. only " << lineelements.size() << " where 4 were expected. Exiting." << std::endl;
				exit(1);
			}
			for (int n = 0; n < dim; n++)
			{
				x(n, nline) = std::stod(lineelements[n]);
			}
			
			nline++;
		}
	}
	fs.close();


	// Normalise the centers to min max
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

	for (int n = 0; n < ncenters; n++)
	{
		x(0, n) = (x(0, n) - minHs) / (maxHs - minHs);
		x(1, n) = (x(1, n) - minT) / (maxT - minT);
		x(2, n) = x(2, n)*pi/180.0;
		x(3, n) = (x(3, n) - minNM) / (maxNM - minNM);
	}



	fs.open("MauiBay_Shore_MDA_zs.txt");

	if (fs.fail()){
		std::cerr << "MauiBay_Shore_MDA_zs.txt file could not be openned" << std::endl;
	}

	line.clear();
	lineelements.clear(); // being frugal here
	nline = 0;
	while (std::getline(fs, line))
	{
		if (!line.empty())
		{
			
			y(nline) = std::stod(line);
			

			nline++;
		}
	}
	fs.close();


	
	//train the RBF

	mat A(ncenters, ncenters);
	mat r,t,o,P, X,B,b, RBFcoeff;
	mat test, rr,s,f;

	A.zeros();

	for (int i = 0; i < ncenters; i++)
	{
		//
		for (int j = 0; j <= i; j++)
		{
			//
			r = norm(x.col(i)-x.col(j));

			t= exp(-0.5 * r % r / (cst * cst));// is 1 by 1 yeah?
			
			A(i, j) = t(0,0);
			A(j, i) = t(0,0);
		}
	}

	o = ones(ncenters, 1);

	P = join_rows(o, x.t());

	X = join_rows(A, P);

	B = join_rows(P.t(), zeros(dim + 1, dim + 1));
	A = join_cols(X, B);

	A.save("A_mat.txt", raw_ascii);

	b = join_cols(y.t(), zeros(dim + 1, 1));
	b.save("b_mat.txt", raw_ascii);

	RBFcoeff = solve(A, b);
	RBFcoeff.save("RBFcoeff_mat.txt", raw_ascii);

	mat readline = zeros(dim, 1);
	test = zeros(dim, 49211);/////// !!!!!

	//Load the test data
	fs.open("Test_data_MB_Shore.txt");

	if (fs.fail()){
		std::cerr << "Test_data_MB_Shore.txt file could not be openned" << std::endl;
	}

	line.clear();
	lineelements.clear();
	nline = 0;
	while (std::getline(fs, line))
	{
		if (!line.empty())
		{
			lineelements = split(line, '\t');
			if (lineelements.size() < 4)
			{
				std::cerr << "ERROR Test_data_MB_Shore.txt file format error. only " << lineelements.size() << " where 4 were expected. Exiting." << std::endl;
				exit(1);
			}
			for (int n = 0; n < dim; n++)
			{
				test(n, nline) = std::stod(lineelements[n]);
				
			}

			

			nline++;
		}
	}
	fs.close();

	// Normalise the data to the centers max
	for (int n = 0; n < nline; n++)
	{
		test(0, n) = (test(0, n) - minHs) / (maxHs - minHs);
		test(1, n) = (test(1, n) - minT) / (maxT - minT);
		test(2, n) = test(2, n)*pi / 180.0;
		test(3, n) = (test(3, n) - minNM) / (maxNM - minNM);
	}

	rr = zeros(1, ncenters);

	for (int n = 0; n < nline; n++)
	{
		rr = (test.col(n)*ones(1, ncenters)) - x;

		rr = sqrt(sum(rr%rr, 0));

		for (int j = 0; j < ncenters; j++)
		{
			rr(j) = norm(test.col(n) - x.col(j));
		}

		rr.save("rr_mat.txt", raw_ascii);


		mat tmp1, tmp2, tmp3;

		tmp1 = exp(-0.5 * rr % rr / (cst * cst));
		tmp1.save("tmp1_mat.txt", raw_ascii);

		tmp2 = RBFcoeff.rows(0, ncenters - 1);
		tmp2.save("tmp2_mat.txt", raw_ascii);

		tmp3 = tmp2 % tmp1.t();
		tmp3.save("tmp3_mat.txt", raw_ascii);

		s = RBFcoeff(ncenters) + sum(tmp3 , 0);
		s.save("s_mat.txt", raw_ascii);
		
		for (int k = 0; k < dim; k++)
		{
			s = s + RBFcoeff(ncenters + k + 1)*test(k, n);
		}


		results.push_back( s(0, 0));

	}



	return 0;
}
