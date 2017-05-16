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

	parameterstr = "RBFcoefffile";
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
	

	//std::string outputfile;
	parameterstr = "outputfile";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.outputfile = parametervalue;
		//std::cerr << "Bathymetry file found!" << std::endl;
	}
	//
	parameterstr = "ncenters";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.ncenters = std::stoi(parametervalue);
	}

	parameterstr = "isdir";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.isdir = std::stoi(parametervalue);
	}

	parameterstr = "trainRBF";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.trainRBF = std::stoi(parametervalue);
	}
	parameterstr = "saveRBFcoeffs";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.saveRBFcoeffs = std::stoi(parametervalue);
	}
	parameterstr = "interpRBF";
	parametervalue = findparameter(parameterstr, line);
	if (!parametervalue.empty())
	{
		param.interpRBF = std::stoi(parametervalue);
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
	std::size_t found;
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
	int ndim = 1; // adjusted later
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


void readgridncsize(std::string ncfile, int &nx, int &ny, int &nt)
{
	//read the dimentions of grid, levels and time 
	int status;
	int ncid, ndimshh;
	
	int varid;




	int dimids[NC_MAX_VAR_DIMS];   /* dimension IDs */
	char varname[NC_MAX_NAME + 1];
	size_t  *ddimhh;
	//char ncfile[]="ocean_ausnwsrstwq2.nc";

	std::vector<std::string> splittedstr;
	std::string mainvarname, filename;
	// first look fo an question mark
	// If present then the user specified the main variable name, if absent we have to figure it out
	splittedstr = split(ncfile, '?');

	if (splittedstr.size() > 1)
	{
		filename = splittedstr[0];
		mainvarname = splittedstr[1];
	}
	else
	{
		// No user specified vaariables
		//Else use a default name
		filename = ncfile;
		status = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
		if (status != NC_NOERR) handle_error(status);
		int nvarinfile;

		status = nc_inq_nvars(ncid, &nvarinfile);
		

		if (nvarinfile == 1)
		{
			status = nc_inq_varname(ncid, nvarinfile-1,varname);
			mainvarname.assign(varname);

		}
		else
		{
			mainvarname = "RBFcoeff";
		}
		status = nc_close(ncid);
		

	}



	//Open NC file
	printf("Open file\n");
	status = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	//printf(" %s...\n", hhvar);
	status = nc_inq_varid(ncid, mainvarname.c_str(), &varid);
	if (status != NC_NOERR)	handle_error(status);



	status = nc_inq_varndims(ncid, varid, &ndimshh);
	if (status != NC_NOERR) handle_error(status);
	//printf("hhVar:%d dims\n", ndimshh);

	status = nc_inq_vardimid(ncid, varid, dimids);
	if (status != NC_NOERR) handle_error(status);

	ddimhh = (size_t *)malloc(ndimshh*sizeof(size_t));

	//Read dimensions nx_u ny_u 
	for (int iddim = 0; iddim < ndimshh; iddim++)
	{
		status = nc_inq_dimlen(ncid, dimids[iddim], &ddimhh[iddim]);
		if (status != NC_NOERR) handle_error(status);

		//printf("dim:%d=%d\n", iddim, ddimhh[iddim]);
	}

	if (ndimshh > 2)
	{
		nt = ddimhh[ndimshh-3];
		ny = ddimhh[ndimshh-2];
		nx = ddimhh[ndimshh-1];
	}
	else
	{
		nt = 1;
		ny = ddimhh[0];
		nx = ddimhh[1];
	}

	


	status = nc_close(ncid);

	free(ddimhh);
	


}
void readxync(std::string ncfile, double *&xx, double *&yy)
{
	//
	//read the dimentions of grid, levels and time 
	int status;
	int ncid, ndimshh;

	int varid;




	int dimids[NC_MAX_VAR_DIMS];   /* dimension IDs */
	char varname[NC_MAX_NAME + 1];
	char xxvarname[NC_MAX_NAME + 1];
	char yyvarname[NC_MAX_NAME + 1];
	size_t  *ddimhh;
	//char ncfile[]="ocean_ausnwsrstwq2.nc";

	std::vector<std::string> splittedstr;
	std::string mainvarname, filename;
	// first look fo an question mark
	// If present then the user specified the main variable name, if absent we have to figure it out
	splittedstr = split(ncfile, '?');

	if (splittedstr.size() > 1)
	{
		filename = splittedstr[0];
		mainvarname = splittedstr[1];
	}
	else
	{
		// No user specified vaariables
		//Else use a default name
		filename = ncfile;
		status = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
		if (status != NC_NOERR) handle_error(status);
		int nvarinfile;

		status = nc_inq_nvars(ncid, &nvarinfile);


		if (nvarinfile == 1)
		{
			status = nc_inq_varname(ncid, nvarinfile - 1, varname);
			mainvarname.assign(varname);

		}
		else
		{
			mainvarname = "RBFcoeff";
		}
		status = nc_close(ncid);


	}



	//Open NC file
	printf("Open file\n");
	status = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	//printf(" %s...\n", hhvar);
	status = nc_inq_varid(ncid, mainvarname.c_str(), &varid);
	if (status != NC_NOERR)	handle_error(status);



	status = nc_inq_varndims(ncid, varid, &ndimshh);
	if (status != NC_NOERR) handle_error(status);
	//printf("hhVar:%d dims\n", ndimshh);

	status = nc_inq_vardimid(ncid, varid, dimids);
	if (status != NC_NOERR) handle_error(status);

	
	// now get the dim name and teh coresponding variable
	int xdimid = dimids[ndimshh - 1];
	int ydimid = dimids[ndimshh - 2];

	int xvarid, yvarid;
	size_t nnx, nny;

	status = nc_inq_dim(ncid, xdimid, xxvarname,&nnx);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_dim(ncid, ydimid, yyvarname, &nny);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, xxvarname, &xvarid);
	if (status != NC_NOERR) handle_error(status);
	
	status = nc_inq_varid(ncid, yyvarname, &yvarid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_var_double(ncid, xvarid, xx);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_var_double(ncid, yvarid, yy);
	if (status != NC_NOERR) handle_error(status);
	
	
}


arma::cube read3Dnc(std::string ncfile, int nx, int ny, int nt)
{
	//read the dimentions of grid, levels and time 
	int status;
	int ncid,varid;

	double * tmpdatastore;
	tmpdatastore = (double *)malloc(nt*nx*ny*sizeof(double));

	int dimids[NC_MAX_VAR_DIMS];   /* dimension IDs */
	char varname[NC_MAX_NAME + 1];
	size_t  *ddimhh;
	//char ncfile[]="ocean_ausnwsrstwq2.nc";

	std::vector<std::string> splittedstr;
	std::string mainvarname, filename;
	// first look fo an question mark
	// If present then the user specified the main variable name, if absent we have to figure it out
	splittedstr = split(ncfile, '?');

	if (splittedstr.size() > 1)
	{
		filename = splittedstr[0];
		mainvarname = splittedstr[1];
	}
	else
	{
		// No user specified vaariables
		//Else use a default name
		filename = ncfile;
		status = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
		if (status != NC_NOERR) handle_error(status);
		int nvarinfile;

		status = nc_inq_nvars(ncid, &nvarinfile);


		if (nvarinfile == 1)
		{
			status = nc_inq_varname(ncid, nvarinfile - 1, varname);
			mainvarname.assign(varname);

		}
		else
		{
			mainvarname = "RBFcoeff";
		}
		status = nc_close(ncid);


	}



	//Open NC file
	printf("Open file\n");
	status = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, mainvarname.c_str(), &varid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_var_double(ncid, varid, tmpdatastore);
	if (status != NC_NOERR) handle_error(status);

	//restructure the data to the mat... So old school need to use the netcdf4 thing
	arma::cube data;
	data = arma::zeros(nt, ny, nx);
	for (int k = 0; k < nt; k++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int i = 0; i < nx; i++)
			{
				data(k, j, i) = tmpdatastore[i+j*nx+k*nx*ny];
			}
		}
	}


	status = nc_close(ncid);
	free(tmpdatastore);
	return data;
}


arma::mat read2Dnc(std::string ncfile, int nx, int ny)
{
	//read the dimentions of grid, levels and time 
	int status;
	int ncid, varid;

	double * tmpdatastore;
	tmpdatastore = (double *)malloc(nx*ny*sizeof(double));

	int dimids[NC_MAX_VAR_DIMS];   /* dimension IDs */
	char varname[NC_MAX_NAME + 1];
	size_t  *ddimhh;
	//char ncfile[]="ocean_ausnwsrstwq2.nc";

	std::vector<std::string> splittedstr;
	std::string mainvarname, filename;
	// first look fo an question mark
	// If present then the user specified the main variable name, if absent we have to figure it out
	splittedstr = split(ncfile, '?');

	if (splittedstr.size() > 1)
	{
		filename = splittedstr[0];
		mainvarname = splittedstr[1];
	}
	else
	{
		// No user specified vaariables
		//Else use a default name
		filename = ncfile;
		status = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
		if (status != NC_NOERR) handle_error(status);
		int nvarinfile;

		status = nc_inq_nvars(ncid, &nvarinfile);


		if (nvarinfile == 1)
		{
			status = nc_inq_varname(ncid, nvarinfile - 1, varname);
			mainvarname.assign(varname);

		}
		else
		{
			mainvarname = "z";
		}
		status = nc_close(ncid);


	}



	//Open NC file
	printf("Open file\n");
	status = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_inq_varid(ncid, mainvarname.c_str(), &varid);
	if (status != NC_NOERR) handle_error(status);

	status = nc_get_var_double(ncid, varid, tmpdatastore);
	if (status != NC_NOERR) handle_error(status);

	//restructure the data to the mat... So old school need to use the netcdf4 thing
	arma::mat data;
	data = arma::zeros( ny, nx);
	for (int j = 0; j < ny; j++)
		{
			for (int i = 0; i < nx; i++)
			{
				data(j, i) = tmpdatastore[i + j*nx ];
			}
		}
	


	status = nc_close(ncid);
	free(tmpdatastore);
	return data;
}



extern "C" void create3dnc(std::string outfile, int nx, int ny, int nt, double *xx, double *yy, double *theta, double * var)
{
	int status;
	int ncid, xx_dim, yy_dim, time_dim, p_dim, tvar_id, gvar_id;

	size_t nxx, nyy, ntt;
	static size_t start[] = { 0, 0, 0 }; // start at first value 
	static size_t count[] = { nt, ny, nx };

	static size_t Gstart[] = { 0, 0 }; // start at first value 
	//static size_t Gcount[] = { ny, nx };

	int time_id, xx_id, yy_id, tt_id;	//
	nxx = nx;
	nyy = ny;
	ntt = nt;

	std::vector<std::string> splittedstr;
	std::string mainvarname, filename;
	// first look fo an question mark
	// If present then the user specified the main variable name, if absent we have to figure it out
	splittedstr = split(outfile, '?');

	if (splittedstr.size() > 1)
	{
		filename = splittedstr[0];
		mainvarname = splittedstr[1];
	}
	else
	{
		// No user specified vaariables
		//Else use a default name
		filename = outfile;
		mainvarname = "3Dvar";



	}


	//create the netcdf dataset
	status = nc_create(filename.c_str(), NC_NOCLOBBER, &ncid);

	//Define dimensions: Name and length

	status = nc_def_dim(ncid, "xx", nxx, &xx_dim);
	status = nc_def_dim(ncid, "yy", nyy, &yy_dim);
	status = nc_def_dim(ncid, "ntheta", ntt, &p_dim);
	//status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	//int tdim[] = { time_dim };
	int xdim[] = { xx_dim };
	int ydim[] = { yy_dim };
	int pdim[] = { p_dim };

	//define variables: Name, Type,...
	int  var_dimids[3];
	int gam_dimids[2];

	//var_dimids[0] = time_dim;
	var_dimids[0] = p_dim;
	var_dimids[1] = yy_dim;
	var_dimids[2] = xx_dim;

	gam_dimids[0] = yy_dim;
	gam_dimids[1] = xx_dim;

	//status = nc_def_var(ncid, "time", NC_DOUBLE, 1, tdim, &time_id);
	status = nc_def_var(ncid, "xx", NC_DOUBLE, 1, xdim, &xx_id);
	status = nc_def_var(ncid, "yy", NC_DOUBLE, 1, ydim, &yy_id);
	status = nc_def_var(ncid, "theta", NC_DOUBLE, 1, pdim, &tt_id);


	status = nc_def_var(ncid, mainvarname.c_str(), NC_DOUBLE, 3, var_dimids, &tvar_id);

	status = nc_enddef(ncid);


	static size_t tst[] = { 0 };
	static size_t xstart[] = { 0 }; // start at first value 
	static size_t xcount[] = { nx };

	static size_t ystart[] = { 0 }; // start at first value 
	static size_t ycount[] = { ny };

	static size_t tstart[] = { 0 }; // start at first value 
	static size_t tcount[] = { nt };


	//Provide values for variables
	//status = nc_put_var1_double(ncid, time_id, tst, &totaltime);
	status = nc_put_vara_double(ncid, xx_id, xstart, xcount, xx);
	status = nc_put_vara_double(ncid, yy_id, ystart, ycount, yy);
	status = nc_put_vara_double(ncid, tt_id, tstart, tcount, theta);

	status = nc_put_vara_double(ncid, tvar_id, start, count, var);
	//status = nc_put_vara_double(ncid, gvar_id, Gstart, Gcount, Gvar);
	status = nc_close(ncid);

}




extern "C" void createTrainingnc(std::string outfile,int nx, int ny, int nt, double *xx, double *yy, double *theta, double * var, double * Gvar)
{
	int status;
	int ncid, xx_dim, yy_dim, time_dim, p_dim, tvar_id, gvar_id;

	size_t nxx, nyy, ntt;
	static size_t start[] = { 0, 0, 0 }; // start at first value 
	static size_t count[] = { nt, ny, nx };

	static size_t Gstart[] = { 0, 0 }; // start at first value 
	static size_t Gcount[] = { ny, nx };

	int time_id, xx_id, yy_id, tt_id;	//
	nxx = nx;
	nyy = ny;
	ntt = nt;

	std::vector<std::string> splittedstr;
	std::string mainvarname, filename;
	// first look fo an question mark
	// If present then the user specified the main variable name, if absent we have to figure it out
	splittedstr = split(outfile, '?');

	if (splittedstr.size() > 1)
	{
		filename = splittedstr[0];
		mainvarname = splittedstr[1];
	}
	else
	{
		// No user specified vaariables
		//Else use a default name
		filename = outfile;
		mainvarname = "3Dvar";
		


	}


	//create the netcdf dataset
	status = nc_create(filename.c_str(), NC_NOCLOBBER, &ncid);

	//Define dimensions: Name and length

	status = nc_def_dim(ncid, "xx", nxx, &xx_dim);
	status = nc_def_dim(ncid, "yy", nyy, &yy_dim);
	status = nc_def_dim(ncid, "ntheta", ntt, &p_dim);
	//status = nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim);
	//int tdim[] = { time_dim };
	int xdim[] = { xx_dim };
	int ydim[] = { yy_dim };
	int pdim[] = { p_dim };

	//define variables: Name, Type,...
	int  var_dimids[3];
	int gam_dimids[2];

	//var_dimids[0] = time_dim;
	var_dimids[0] = p_dim;
	var_dimids[1] = yy_dim;
	var_dimids[2] = xx_dim;

	gam_dimids[0] = yy_dim;
	gam_dimids[1] = xx_dim;

	//status = nc_def_var(ncid, "time", NC_DOUBLE, 1, tdim, &time_id);
	status = nc_def_var(ncid, "xx", NC_DOUBLE, 1, xdim, &xx_id);
	status = nc_def_var(ncid, "yy", NC_DOUBLE, 1, ydim, &yy_id);
	status = nc_def_var(ncid, "theta", NC_DOUBLE, 1, pdim, &tt_id);


	status = nc_def_var(ncid, mainvarname.c_str(), NC_DOUBLE, 3, var_dimids, &tvar_id);
	status = nc_def_var(ncid, "gamma", NC_DOUBLE, 2, gam_dimids, &gvar_id);

	status = nc_enddef(ncid);


	static size_t tst[] = { 0 };
	static size_t xstart[] = { 0 }; // start at first value 
	static size_t xcount[] = { nx };

	static size_t ystart[] = { 0 }; // start at first value 
	static size_t ycount[] = { ny };

	static size_t tstart[] = { 0 }; // start at first value 
	static size_t tcount[] = { nt };


	//Provide values for variables
	//status = nc_put_var1_double(ncid, time_id, tst, &totaltime);
	status = nc_put_vara_double(ncid, xx_id, xstart, xcount, xx);
	status = nc_put_vara_double(ncid, yy_id, ystart, ycount, yy);
	status = nc_put_vara_double(ncid, tt_id, tstart, tcount, theta);

	status = nc_put_vara_double(ncid, tvar_id, start, count, var);
	status = nc_put_vara_double(ncid, gvar_id, Gstart, Gcount, Gvar);
	status = nc_close(ncid);

}
void handle_error(int status) {
	if (status != NC_NOERR) {
		fprintf(stderr, "Netcdf %s\n", nc_strerror(status));
		


		//fprintf(logfile, "Netcdf: %s\n", nc_strerror(status));
		exit(-1);
	}
}
