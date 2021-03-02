/*
 * haploview.cpp
 *
 *  Created on: 17.01.2013
 *  Author: Nab Raj Roshyara
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */
#include"haploview.h"
void CHAPPED::write_ped_file(const vector<CBPED*>&pedInfo, const vector<CBSNP*>&genInfo,const string& outFileName){
	CBPED::write_plinkPedFile(pedInfo,genInfo, outFileName);
}
void CHAPSNP::write_info_file(const vector<CBSNP*>&genInfo,const string& outFileName)
{
	CBSNP::write_rsid_basePairposition(genInfo, outFileName);
	printLIN("\t **->Haploview's SNP-info file has been written out and saved as  \""+\
			  			outFileName+".info \".\n");
}



