/*
 * beagle.h
 *
 *  Created on: Jan 6, 2012
  *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _BEAGLE_H__
#define _BEAGLE_H__
#include"helper.h"
#include"machdat.h"
#include"machped.h"

class CBEPED:public CBPED
{
	public:
	CBEPED():CBPED(){}
	//vector<CBPED*> vPed;

	static void read_bgl_data(vector<CBPED*>& pedInfo, vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName);
	static void read_bgl_gprobs_data(vector<CBPED*>& pedInfo, vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName, const float & thresh);
	void static write_beagle_gprobsFile(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName);

	~CBEPED();
};
class CBESNP:public CBSNP
{
public:
	CBESNP():CBSNP(){}
	//vector<CBSNP*> genInfo;
	~CBESNP();
	static void lese_beagle_rsq_score_file(vector<CBSNP*>&genInfo,vector<string>&cZeile, const double& rsq_thresh, const string& fileName);
	//static void beagle_type_commands(const string & read_code_type);
	static bool _given_gprobs_header;
};
#endif /* _BEAGLE_H__ */
