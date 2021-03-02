/*
 * helper.cpp
 *
 *  Created on: Nov 4, 2011
  *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.

 */

#include "helper.h"
#include<algorithm>


extern ofstream LIN;

void error (const string&  msg)
{
	//first delete 
	cerr<<" ERROR: " <<msg;
	LIN<<" ERROR:"<< msg;
	LIN.close();
	exit(1);
}
// print LIn
void printLIN(const string & msg)
{
	cout<< msg;
	cout.flush();
	LIN<< msg;
	LIN.flush();

}

void printWithSpace(const string& msg, const string& value, const int coutWidth)
{
	cout<<setw(coutWidth)<< msg<< value;
	cout.flush();
	LIN<<setw(coutWidth)<< msg<<value;
	LIN.flush();
	cout << endl;
	LIN<<"\n";
}

string change_int_into_string(const int value)
{
	string count;
	stringstream ss;
	ss<<value;
	count=ss.str();
	return count;
}

string change_double_into_string(const double value)
{
	string count;
	stringstream ss;
	ss<<value;
	count=ss.str();
	return count;
}
// function by return type.
template< class T>
string changeIntoString(const T value)
{
	string count;
	stringstream ss;
	ss<< value;
	count=ss.str();
	return count;
}
// comparing 
bool compfun_with_given_vector(const vector<string> &myVec, const string & s )
	{ 
		 vector<string>nVec=myVec;
		 sort(nVec.begin(),nVec.end()); 
		 bool bv = binary_search(nVec.begin(),nVec.end(),s);
		 return bv;
	}
	
vector<int> pos_in_vec(vector<string>& comVec, const string&  s )
{
	vector<int>myints;
	int ncol =-9; // if not found
	for(vector<string>::iterator it=comVec.begin();it<comVec.end();++it)
	{
		if(*it==s)
			ncol=(int)(it-comVec.begin());
		if(ncol>=0)	
			myints.push_back(ncol);	
		ncol=-9;
	}	
 return myints;

}

vector<int> pos_notin_vec(vector<string>& comVec, const string&  s )
{
	vector<int>myints;
	int ncol =-9; // if not found
	for(vector<string>::iterator it=comVec.begin();it<comVec.end();++it)
	{
		if(*it!=s)
			ncol=(int)(it-comVec.begin());
		if(ncol>=0)	
			myints.push_back(ncol);	
		ncol=-9;
	}	
 return myints;

}
//
void displayHelp(const bool tfValue) //i.e. if hilfe true
{
if(tfValue)
	{
		//string msg= "At the moment this program can only read PLINK and MACH format files and  convert them  into MACH or IMPUTE and PLINK or IMPUTE format respectively.\n";
		//printLIN (msg);
		//
		//printLIN("========================================================================\n");
		printLIN("*->fcGENE's command options can be divided into two parts:\n\t(a)uploading genotype data given in  formats of popular software like PLINK, MACH(MERLIN), IMPUTE and BEAGLE \n\t(b)writing the uploaded data into any format specified by command option \"- -oformat\". \n");
		printLIN ("\t*->For a detailed information on fcGENE, please read its pdf-documentation file,\n\twhich can be downloaded from:  \"http://sourceforge.net/projects/fcgene/files/\" \n" );
		printLIN ("*->Some basic commands of fcGENE:\n");
		printLIN("\t**->To read \"plink\" format  files:\n");
		printLIN ( "\t./fcgene --ped  plink.ped --map plink.map \n");
		//
		printLIN("\t**->To read \"mach\" format  files:\n");
		printLIN ( "\t./fcgene --ped  mach.ped --dat mach.dat \n");

		printLIN( "\t**-> Converting from PLINK format to MACH format:\n");
		printLIN ( "\t./fcgene --map fcgene.map --ped fcgene.ped --oformat mach --out plink_mach \n");
		printLIN ( "\t--out command is optional.If this command does not exist, out file will be saved with the name fcgene_out.\n");
		printLIN( "\t**-> From Plink format to IMPUTE format:\n");
		printLIN ( "\t./fcgene --map fcgene.map --ped fcgene.ped --oformat impute --out plink_impute \n");
		//
		printLIN( "\t**-> Converting from MACH format to PLINK format:\n");
		printLIN ( "\t./fcgene --dat fcgene.dat --ped fcgene.ped --oformat plink --out mach_plink \n");
		printLIN ( "\tTo update the pedinfo and snpinfo again: \n");
		printLIN ( "\t. \"--pedinfo pedinfo.txt\" and \"--mapinfo mapinfo.txt\" commands can be used everywhere while converting format \n");
		printLIN( "\t**-> Converting from MACH format to IMPUTE format:\n");
		printLIN ( "\t./fcgene --dat fcgene.dat --ped fcgene.ped --oformat impute --out mach_impute \n");
		printLIN( "\t**-> Converting from PLINK format to IMPUTE format:\n");
		printLIN ( "\t./fcgene --ped plink.ped --map plink.map  --oformat impute --out plink_impute \n");
		//
		printLIN( "\t**-> Converting from IMPUTE format to PLINK format:\n");
		printLIN ( "\t./fcgene --gens fcgene.gens --oformat plink --out impute_plink \n");
		//
		printLIN( "\t**-> Converting from PLINK format to BEAGLE format:\n");
		printLIN ( "\t./fcgene --ped plink.ped --map plink.map  --oformat beagle --out plink_beagle \n");
		//
		printLIN( "\t**-> Converting from PLINK format to Bimbam format:\n");
		printLIN ( "\t./fcgene --map fcgene.map --ped fcgene.ped --oformat bimbam --out plink_bimbam \n");
		printLIN("|===============================================================================|\n");

	}

return ;}
//---------------------------------------------------------------------//

std::vector<std::string> vector_from_string(const std:: string & ss){
	std::string _str="";
	std::istringstream iss (ss);
	std::vector<std::string> tokens;
	while(iss>>_str){  tokens.push_back(_str);}
	return tokens;
}
std::vector<std::string> vector_from_string(const std:: string&ss, const char sep){
	std::string _str="";
	std::istringstream iss (ss);
	std::vector<std::string> tokens;
	while( getline(iss,_str, sep)){  tokens.push_back(_str);}
	return tokens;
}
int str2int(const std::string &in, const int missing_value)
{
	if ((in.size() == 0) || (in == "."))
		return missing_value;
	else
		return atoi(in.c_str());
}

double str2double(const std::string &in, const double missing_value)
{
	if ((in.size() == 0) || (in == "."))
		return missing_value;
	else
		return atof(in.c_str());
}

std::string int2str(const int in, const int missing_value)
{
	if (in == missing_value)
		return ".";
	else
	{
		static std::ostringstream out;
		out.str(""); out.clear();
		out << in;
		return out.str();
	}
}

std::string double2str(const double in, const double missing_value)
{
	if (in == missing_value)
		return ".";
	else
	{
		static std::ostringstream out;
		out.str(""); out.clear();
		out << in;
		return out.str();
	}
}

double round_nplaces(const double&value, const int& to)
{
	double places = pow(10.0, to);
	//std::cout<<"places: "<< places<<std::endl;
	return round(value * places) / places;
}
