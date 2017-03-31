#ifndef RBF_H
#define RBF_H
 

#include <stdio.h>
#include <cmath>

#include <iostream>
#include <string.h>
#include <math.h>
#include <armadillo>

//using namespace std;


struct Param{
	//Input parameters
	double gamma; // future this should be a vector

	int ncenters;

	int dim;

	//
	std::string centersfile;
	std::string trainingfile;
	std::string RBFcoefffile;
	std::string inputfile;
};

Param readparamstr(std::string line, Param param);
std::string findparameter(std::string parameterstr, std::string line);
void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
std::string trim(const std::string& str, const std::string& whitespace);
arma::mat readdatafile(std::string filename);


#endif