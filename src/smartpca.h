/*
 * smartpca.h
 * Created on: 17.01.2013
 *  Author: Nab Raj Roshyara
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _SMARTPCA_H__
#define _SMARTPCA_H__
#include "helper.h"
#include "machdat.h"
#include "machped.h"
class CSMARTPCASNP{
public:
	static void write_smartpca_pedsnp_file(const vector<CBSNP*>&genInfo,const string& outFileName);
	static void write_smartpca_command_file(const string& outFileName);
};

class CSMARTPCAPED{
public:
	static void write_smartpca_ped_file(const vector<CBPED*>&pedInfo,const vector<CBSNP*>&genInfo,const string& outFileName);
	static void write_smartpca_pedind_file(const vector<CBPED*>&pedInfo,const string& outFileName);
	static void read_smartpca_groupLabel_file(vector<CBPED*>& pedInfo, const string& inFileName);

};



#endif /* _SMARTPCA_H__ */
