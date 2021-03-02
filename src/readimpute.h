/*
 * readimpute.h
 *
 *  Created on: 03.01.2012
  *  Created on: Dec 29, 2011
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef _READIMPUTE_H__
#define _READIMPUTE_H__
#include "machdat.h"
#include "machped.h"
#include  "helper.h"
extern ofstream LIN;
struct CISNP:public CBSNP
{
	CISNP():CBSNP(){ }
	//vector<CBSNP*>genInfo;
	~CISNP();
	static void gensFileLeser(vector<CBSNP*>& genInfo,const string& gensFileName);
	//static void extraGensFileLeser(const string& fileName, const vector<CBSNP*> & genInfo);
	static void lese_imputegens_file(vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName);
	static void lese_imputed_info_score_file(vector<CBSNP*>&genInfo,vector<string>&cZeile, const double& info_thresh, const double& maf_thresh, const string& fileName);
	//helper function
	static void lese_imputegens_file_genInfo( CBSNP* const pCBSNP, const vector<string>& zeilenVec,const int akt_zeile);
	static void read_impute_ref_legendfile(vector<CBSNP*>& genInfo,const string& legendFileName);
	static void lese_shapeit_haps_file(vector<CBPED*>&vPed,vector<CBSNP*>& genInfo, const string& fileName);
};


struct CIPED:public CBPED
{

	//vector<CBPED*> vPed;  // vector of CBPED;+
		//functions
	CIPED():CBPED(){}
	~CIPED();
	static void imputegen_to_pedgen_converter(const vector<CBSNP*>& genInfo,const float& thresh);
	//static void write_impute_plink_Files(const vector<CBPED*>& pedInfo, const vector<CBSNP*> & genInfo, const string&outFileName,const float& thresh);
	//static void write_impute_mach_Files(const vector<CBPED*>& pedInfo, const vector<CBSNP*> & genInfo, const string&outFileName,const float& thresh);
	static void pedFileLeser(vector<CBPED*>&pedInfo,const vector<CBSNP*>genInfo);
	static void lese_hapmap_genotype_vec(CBSNP* const pCBSNP,const vector<string> &zeilenVec, const int& nLines);
	static void lese_impute_hapmap_file(vector<CBPED*>&vPed,vector<CBSNP*>& genInfo,vector<string>& cZeile,const string& fileName);
	//static void lese_impute_hapmap_legend_file(vector<CBPED>&vPed, vector<CBSNP*>& genInfo, vector<string>&cZeile, const string& fileName);
	static void lese_shapeit_sample_file(vector<CBPED*>&vPed,const string& fileName);
};





#endif /* _READIMPUTE_H__ */
