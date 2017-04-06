#ifndef RBF_H
#define RBF_H
 

#include <stdio.h>
#include <cmath>
#include <netcdf.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <armadillo>

//using namespace std;


struct Param{
	//Input parameters
	double gamma; // future this should be a vector
	int ncenters;
	int ndim;

	int isdir= -1; //which dimension is directional (index starting at 0!)

	int trainRBF = 1;

	int saveRBFcoeffs = 1;

	int interpRBF = 1;




	//
	std::string centersfile;
	std::string trainingfile;
	std::string RBFcoefffile;
	std::string inputfile;
	std::string outputfile;

};

Param readparamstr(std::string line, Param param);
std::string findparameter(std::string parameterstr, std::string line);
void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
std::string trim(const std::string& str, const std::string& whitespace);
arma::mat readdatafile(std::string filename);
void writedatafile(std::vector<double> outputdata, std::string outputfile);
void readgridncsize(std::string ncfile, int &nx, int &ny, int &nt);
arma::cube read3Dnc(std::string ncfile, int nx, int ny, int nt);
void handle_error(int status);

extern "C" void create3dnc(std::string outfile, int nx, int ny, int nt, double dx, double dy, double dtheta, double totaltime, double *xx, double *yy, double *theta, double * var);


#endif