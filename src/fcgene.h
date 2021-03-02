/*
 * fcgene.h
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

#ifndef _FCGENE_H__
#define _FCGENE_H__
#include "helper.h"
#include "argcv.h"
#include "bpar.h"
#include "cparse.h"
#include  "plinksnp.h"
#include "plinkped.h"
//#include "mach.h"
//#include "impute.h"
#include "machdat.h"
#include "machped.h"
#include "readimpute.h"
#include "beagle.h"
#include "bimbam.h"
#include "general.h"
#include "matrix_def.h"
class TimeInfo;
//extern vector<CSNP>genInfo;
class Argcv;
class Bpar;
class CSNP;
//void read_data(CGENERAL* const pGENERAL,  Bpar* pBpar);
	void handel_with_general_commands(CGENERAL* const pGENERAL,Bpar* pBpar);
	int handel_with_merging(vector<CGENERAL*> CGEN_VEC,Bpar*pBpar);
	void merge_second_data_into_first(CGENERAL* first, const CGENERAL* second, const float &thresh);
//static void show_data_summary_if_given_merge(const vector<CGENERAL*>& pGEN_VEC);
// This function will
	//void handel_with_fst_calculation(vector<CGENERAL*> CGEN_VEC,Bpar* pBpar);
  // void calculate_fst_btwn_1st_n_2nd(const CGENERAL* first, const CGENERAL* second);
   void handel_with_splitting(const CGENERAL* const pCGENERAL, const Bpar* const pBpar);

class fcGENE{
public:
	fcGENE();


private:

};
class TimeInfo
{
	public:
	TimeInfo();
	inline void setTime(time_t & tValue){ time(&tValue);}
	void displayTime();
	~TimeInfo();
private:
	time_t startTime;
	time_t endTime;
	double diffTime;
	struct tm * timeInfo;
	int a_hours; // hours taken for analysis
	int a_minutes; // minutes taken for analysis
	double a_seconds; // seconds taken for analysis

};


#endif /* _FCGENE_H__ */
