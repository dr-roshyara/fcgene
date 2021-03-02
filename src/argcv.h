/*
 *
 *  Created on: Dec 21, 2011
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _ARGCV_H__
#define _ARGCV_H__
#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include <cstdlib>
#include <algorithm>
#include<exception>
#include "helper.h"
#include "bpar.h"
using namespace std;
class Argcv
{
public:
	Argcv(int , char**);
	~Argcv();
	 static int scsFind(const string &s ); 							// single command string finding
	 static void dcsFind(string & coms, const string &s ); 			// double command  command finding
	 static bool compfun(const string & s ){bool bv = binary_search(stcomVec.begin(),stcomVec.end(),s); return bv;}
	 static void printOptions(const Bpar* pBpar);
	 static void print_new_options(const Bpar* pBpar);
	 static void print_given_options(); //this function just prints the given command options but do not check anything
	 static vector<string> parse_subArg(string& subArg, const string & seperator);
	 //void cuc(const  vector<string>& ); // check unused commands
	 //elements
	static vector<Bpar*> bparam_vec;
	static vector<string> comVec;
private:
	static int argCount;
 	static int bparas_vec_sz;
	static vector< bool > parsed;
	static vector<string> original; 			// unsed command. one can use comVec instead of original.
	static vector<bool> option; 				// unused command till now
	static vector<string> stcomVec; 			// sorted command vector it will be used for binary_search
	static string msg;

};

#endif
