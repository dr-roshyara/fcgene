/*
 * pedgree.h
 *
 *  Created on: Dec 9, 2011
  *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _PEDGREE_H__
#define _PEDGREE_H__
#include "plinksnp.h"
#include "helper.h"
#include "machped.h"
//defining base class
using namespace std;
extern ofstream LIN;
struct CPPED:public CBPED
{
	CPPED():CBPED(){}
	 //vector<CBPED*> vPed;
	 ~CPPED();
	// functions to read plink files
	 	 static CBPED*  pedInfoLeser(const vector<string>& zeilenVec, const long int& snpAnzahl, const int& akt_zeile);
	 	static bool zeilenLeser(IFS pf, vector<string>& pZeile); // reading end of line  returns true;
	 	static void pedFileleser(const string& pedFileName,	const vector<CBSNP*>& genInfo, vector<CBPED*> &vPedInfo);
		//plink-dosage file
	 	static void lese_plink_dosage_file(const vector<CBPED*>&pedInfo, vector<CBSNP*>&genInfo,vector<string>&cZeile,const float& thresh, const string&dosageFileName);
	 	static void lese_plink_raw_file(vector<CBPED*>&pedInfo, vector<CBSNP*>&genInfo,vector<string>&cZeile,const string&rawFileName);
	 	static void lese_plink_fam_file(vector<CBPED*>&pedInfo,const string&famFileName);
	 	//
	 	static void lese_bim_file(vector<CBSNP*>&genVec, const string & fileName);
	 	static void lese_bed_file(vector<CBPED*>&pedVec,vector<CBSNP*>&genVec, const string & fileName);
	 	static bool openBinaryFile(string s, ifstream & BIT);
	 	static void write_BITFILE(const vector<CBSNP*>&genVec, const vector<CBPED*>&pedVec, const string &outFile);
	 	//read extra co-variate file This will be used especially for writing SNPTEST format file
	 	static bool _read_recodeA_type;
	 	static bool	_read_recodeAD_type;

	 	//r
	 //functions
	 void displayIndividualSummary(const vector<CBPED*>& vPedInfo);
	 //
	//converting functions
	 	static void writeOutMachFiles(const vector<CBPED*>& pedInfo ,const vector<CBSNP*> &genInfo, const string&outFileName);
	 	static void writeOutImputeFiles(const vector<CBPED*>& pedInfo,const vector<CBSNP*> &genInfo, const string&outFileName);
	 	static void write_pedinfoFiles(vector<CBPED*> pedInfo, const string& outFileName);
	 	//static void write_plink_beagle_genotypeFile(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName);
	 	//static void plink_type_commands(const string& read_code_type);
	 // convert into plink-raw files
	 // this function needs --oformat plink-recodeA plink-recodeAD
	 	static void write_plink_raw_files(const vector<CBPED*>&pedInfo,const vector<CBSNP*>&genInfo, const string&outFileName);
	 	static void write_plink_raw_dose_files(const vector<CBPED*>&pedInfo,const vector<CBSNP*>&genInfo, const string&outFileName, const string& ref_allele);

		//static void readBinData(const string &s);
private:
	static const int coutWidth =50;

};
#endif /* _PEDGREE_H__ */
