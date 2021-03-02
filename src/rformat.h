/*
 * 	rformat.h
 *
 *  Created on: 09.01.2013
 *  Author: Nab Raj Roshyara
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _RFORMAT_H__
#define _RFORMAT_H__
#include "helper.h"
#include "machdat.h"
#include "machped.h"
// for locus information
class CRSNP:public CBSNP{
public:
	static void lese_genotype_file(vector<CBPED*>& vPed, vector<CBSNP*>& genVec, const string & inFileName);
	static void lese_transposed_genotype_file(vector<CBPED*>& vPed, vector<CBSNP*>& genVec, const string & inFileName);
	static void write_ref_and_alternative_allele_info_file(const vector<CBSNP*>& genInfo, const string & outFileName, const string & ref_allele);
};

// for pedigree information
class CRPED:public CBPED
{
public:
	static void write_affection_status_file( const vector<CBPED*>& vPed, const string & outFileName);
	static void write_genotype_file( const vector<CBPED*>& vPed, const vector<CBSNP*>& genInfo, const string & outFileName, const string & ref_allele);
	static void write_transposed_genotype_file( const vector<CBPED*>& vPed, const vector<CBSNP*>& genInfo, const string & outFileName, const string & ref_allele);
	static void write_allele_dose_file( const vector<CBPED*>& vPed, const vector<CBSNP*>& genInfo, const string & outFileName, const string& ref_allele);
	//reading


};

#endif /* _RFORMAT_H__ */
