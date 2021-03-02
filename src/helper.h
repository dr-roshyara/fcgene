/*
 * helper.h
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

#ifndef _HELPER_H__
#define _HELPER_H__
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<cstdlib>
#include<cstdio>
#include<bitset>
#include<cmath>
#include<iomanip>
#include<ctime>
#include <sys/timeb.h>
#include<new>
#include<iterator>
#include<cstring>
#include<iostream>
#include<cstdlib>
#include<vector>
#include<string>
#include<cstring>
#include<cstdio>
#include<set>
#include<sstream>
#include <cerrno>


#include "config.h"

//extern vector<CSNP*> genInfo;
using namespace std;
void error(const string & msg);
void printLIN (const string& msg);
void printWithSpace (const string& msg, const string& value, const int coutWidth=0);

void dynDellocate (); // dynDellocate;
string change_int_into_string(const int value);
string change_double_into_string(const double value);
bool compfun_with_given_vector(const vector<string> &myVec, const string & s );
vector<int> pos_in_vec(vector<string>& comVec, const string&  s );
vector<int> pos_notin_vec(vector<string>& comVec, const string&  s );

//
template < class T>
string changeIntoString(const T value);

void displayHelp(const bool tfValue);
// functions added on 2014.01.22
std::vector<std::string> vector_from_string(const std:: string & ss);
std::vector<std::string> vector_from_string(const std:: string&ss, const char sep);
int str2int(const std::string &in, const int missing_value=-1);
double str2double(const std::string &in, const double missing_value=-1.0);
std::string int2str(const int in, const int missing_value=-1.0);
std::string double2str(const double in, const double missing_value=-1.0);

double round_nplaces(const double&value, const int& to);


#endif /* _HELPER_H__ */
/***
 * how to ada a new command:
 *  first crate variables related to command in bpar.h and bpar.cpp
 *  	1-
 */
