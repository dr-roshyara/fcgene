/*
 * cparse.h
 *
 *  Created on: Nov 4, 2011
  *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 *      This file is created to parse the command line
 */
#ifndef _CPARSE_H__
#define _CPARSE_H__

#include<iostream>
#include<string>
#include "helper.h"
#include "bpar.h"
#include "argcv.h"

using namespace std;
void setCoptions( Bpar *  const pBpar);

void cfe(const string & file, const bool &  bValue ); // check file exits
void parse_force_filter_option(const string& force_option,  vector<string>&firstPart, vector<string>&secondPart, const bool&bValue);
// This function will paste family id and indiviudal ids.
string paste_snp_or_ped_info(const string& paste_option,vector<string>&pasted_strings,const string&indicator);
void parse_option_iid( Bpar * const pBpar );
void parse_plink_commands( Bpar * const pBpar );
void parse_general_commands(Bpar * const pBpar );
void parse_mach_minimach_commands(	 Bpar * const pBpar );
void parse_impute_snptest_commands( Bpar * const pBpar );
void parse_beagle_commands( Bpar * const pBpar );
void parse_bimbam_commands( Bpar * const pBpar );
void parse_haploview_commands( Bpar * const pBpar );
void parse_eigensoft_commands( Bpar * const pBpar );
void parse_summary_statistics_commands(	Bpar * const pBpar );
void parse_fst_calculation_command(Bpar* const pBpar);
void check_input_output_file_names_n_paths(Bpar* const pBpar);
void parse_rformat_commands(Bpar* const pBpar);
void parse_shapeit_commands(Bpar* const pBpar);
void parse_ssplit_isplit_commands(Bpar* const pBpar);
//void parse_filter_snp_option(const string& filter_option,vector<string>&firstPart,vector<string>&secondPart,const bool&bValue);

#endif
