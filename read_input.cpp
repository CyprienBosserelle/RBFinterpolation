//////////////////////////////////////////////////////////////////////////////////
//XBeach_GPU                                                                    //
//Copyright (C) 2013 Bosserelle                                                 //
//                                                                              //
//This program is free software: you can redistribute it and/or modify          //
//it under the terms of the GNU General Public License as published by          //
//the Free Software Foundation.                                                 //
//                                                                              //
//This program is distributed in the hope that it will be useful,               //
//but WITHOUT ANY WARRANTY; without even the implied warranty of                //    
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 //
//GNU General Public License for more details.                                  //
//                                                                              //
//You should have received a copy of the GNU General Public License             //
//along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//////////////////////////////////////////////////////////////////////////////////


#include <armadillo>
#include <iostream>
#include "RBF.h"

Param readparamstr(std::string line, Param param)
{


	std::string parameterstr, parametervalue;

	///////////////////////////////////////////////////////
	// General parameters
	parameterstr = "centersfile";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.centersfile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}

	parameterstr = "trainingfile";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.trainingfile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}

	parameterstr = "RBRcoefffile";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.RBFcoefffile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}
	
	parameterstr = "inputfile";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.inputfile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}
	
	//
	parameterstr = "ncenters";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ncenters = std::stoi(parametervalue);
	}

	
	parameterstr = "gamma";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.gamma = std::stod(parametervalue);
	}
	
	

	return param;
}


std::string findparameter(std::string parameterstr, std::string line)
{
	std::size_t found, Numberstart, Numberend;
	std::string parameternumber,left,right;
	std::vector<std::string> splittedstr;
	
	// first look fo an equal sign
	// No equal sign mean not a valid line so skip
	splittedstr=split(line, '=' );
	if (splittedstr.size()>1)
	{
		left = trim(splittedstr[0]," ");
		right = splittedstr[1]; // if there are more than one equal sign in the line the second one is ignored
		found = left.compare(parameterstr);// it needs to strictly compare
		if (found == 0) // found the parameter
		{
			//std::cout <<"found LonMin at : "<< found << std::endl;
			//Numberstart = found + parameterstr.length();
			splittedstr = split(right, ';');
			if (splittedstr.size() >= 1)
			{
				parameternumber = splittedstr[0];
			}
			//std::cout << parameternumber << std::endl;

		}
	}
	return trim(parameternumber, " ");
}

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


void writedatafile(std::vector<double> outputdata, std::string outputfile)
{
	FILE * fs;
	errno_t  err= fopen_s(&fs,outputfile.c_str(),"w");
	for (int i = 0; i < outputdata.size(); i++)
	{
		//
		fprintf(fs, "%lf\n", outputdata[i]);
	}
	fclose(fs);
	
}


arma::mat readdatafile(std::string filename)
{
	int ndim = 1;
	int nline = 0;
	

	arma::mat linedata;
	arma::mat filedata;
	

	

	//Open file 
	std::ifstream fs;


	fs.open(filename);

	if (fs.fail()){
		std::cerr << filename << " file could not be openned" << std::endl;
	}

	std::string line;
	std::vector<std::string> lineelements, spacedelements, tabselenments, commaelements, semicolelements;
	//int nline = 0;
	while (std::getline(fs, line))
	{
		if (!line.empty() && line.substr(0, 1).compare("#") != 0)
		{
			nline++;
			// Here space comma and tabs have special meaning and should not be used unless it ois for delimiting fields
			//check which delimiter gives the largest number of elements
			spacedelements = split(line, ' ');
			tabselenments = split(line, '\t');
			commaelements = split(line, ',');
			semicolelements = split(line, ';');

			lineelements = spacedelements;

			if (lineelements.size() < tabselenments.size())
			{
				lineelements.clear();
				lineelements = tabselenments;
			}

			if (lineelements.size() < commaelements.size())
			{
				lineelements.clear();
				lineelements = commaelements;
			}

			if (lineelements.size() < semicolelements.size())
			{
				lineelements.clear();
				lineelements = semicolelements;
			}

			ndim = lineelements.size();

			linedata = arma::zeros(ndim, 1);
			for (int n = 0; n < lineelements.size(); n++)
			{
				linedata(n, 0) = std::stod(trim(lineelements[n]," "));
			}

			filedata = arma::join_horiz(filedata, linedata);
			
		}
		
	}
	fs.close();
	return filedata;
}

