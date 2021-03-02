/*
 * minimac.h
 *
 *  Created on: 22.01.2013
*  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MINIMAC_H__
#define _MINIMAC_H__
#include "helper.h"
#include "machdat.h"
#include "machped.h"

extern ofstream LIN;
struct CMINIMACSNP:public CBSNP{
	CMINIMACSNP():CBSNP(){ }
		//vector<CBSNP*>genInfo;
	~CMINIMACSNP();
	static void lese_minimac_info_file(const string& fileName, vector<CBSNP*>&genInfo, const float& rsq_value, const float& maf_value);
};
struct CMINIMACPED:public CBPED{
	CMINIMACPED():CBPED(){ }
	static  void lese_minimac_probFile( vector<CBPED*>& pedInfo, const vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName, const float & thresh);
		//vector<CBSNP*>genInfo;
	~CMINIMACPED();
};



#endif /* _MINIMAC_H__ */
