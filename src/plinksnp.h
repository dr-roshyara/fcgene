/*
 * locus.h
 *
 *  Created on: 14.11.2011
 *      Author: nroshyar
  *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef _LOCUS_H__
#define _LOCUS_H__
#include "helper.h"
#include "machdat.h"
using namespace std;
extern ofstream LIN;
struct CPSNP:public CBSNP
{

	//functions to read plink map file
	// plink way.
		static void readMapFile(const string&  fileName,  vector<CBSNP*>& genInfo);
		static bool checkSNPinfo(const vector<string> & line, const long int& countSNPs);
		static CBSNP* fillUpLocusInfo(const vector<string>& lines);
		//static void extraMapFileLeser(const string& fileName, const vector<CBSNP*> & genInfo);
	//	void write_snpinfoFile( const vector<CBSNP*>& genInfo,const string& outFileName);
		// static void writeOutMachFiles(const vector<Cpedgree*>& pedInfo,const vector<CBSNP*> &genInfo,const string& outFileName);


	//vector<CBSNP*>genInfo;

	//functions
	CPSNP():CBSNP( ){	}
	~CPSNP();
	static void snpinfoFileLeser(const string& fileName, const vector<CBSNP*> & genInfo);
	static bool _given_three_cols;
};

//vector<CBSNP*> CBSNP::genInfo;

// change



#endif /* _LOCUS_H__ */


