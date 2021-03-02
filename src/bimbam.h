/*
 * bimbam.h
 *
 *  Created on: May 3, 2012
 *	Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 *      These are the basic parameters
 *
 */

#ifndef _BIMBAM_H__
#define _BIMBAM_H__
#include"helper.h"
#include"machdat.h"
#include"machped.h"
#include<cstdio>
#include<cstdlib>
#include<cmath>

class CBIMSNP:public CBSNP
{
public:
	CBIMSNP(){};
	~CBIMSNP();
	//vector<CBSNP*>genInfo;
	//
	static void read_bimbam_pos_data(vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName);
	static void read_bimbam_snpinfo_data(vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName,const float &maf_value);
	static void bimbam_type_commands(const string & read_code_type);
};
// ped class

class CBIMPED:public CBPED
{
public:
	CBIMPED(){};
	~CBIMPED();
	//vector<CBPED*> vPed;
//functions
	static void read_bimbam_data(vector<CBPED*>& pedInfo, vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName);
	static void read_bimbam_gprobs_data(vector<CBPED*>& pedInfo, vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName, const float & thresh);
	static void read_bimbam_bestguess_data(vector<CBPED*>& pedInfo, vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName);


	//help function
	static void lese_bimbam_bestguess_genInfo(const vector<CBSNP*>&genInfo, const vector<string>& zeilenVec,unsigned const int akt_zeile);

};


#endif /* _BIMBAM_H__ */
