/*
 * general.h
 *
 *  Created on: 04.10.2012
 *      Author: nroshyar
  *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __GENERAL_H__
#define __GENERAL_H__
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include "helper.h"
#include "bpar.h"
#include "machdat.h"
#include "machped.h"
#include "vcf_header.h"
#include "vcf.h"
vector<int> provide_correct_geno_012(CBSNP* const pCBSNP,const string &ref_allele="");
void calculate_nei_wc_hudson_fst(vector<CBPED*>&pedVec, const vector<CBSNP*>&genVec, const string& groupLabelFileName,const string &outFileName, const string& ref_allele="");
void handel_reich_fst_calculation(vector<CBPED*>&pedVec, const vector<CBSNP*>&genVec);
void handel_fst_calculations( vector<CBPED*>&pedVec, const vector<CBSNP*>&genVec,const string& groupLabelFileName, const string & outFileName );
void splitt_snps_iids( const vector<CBPED*>& pedVec, const vector<CBSNP*>& genVec, const Bpar * const pBpar);

struct CGENERAL
{
public:
		//static variables
		vector<CBSNP*> genVec;
		vector<CBPED*> pedVec;
		 VCF CVCF;

		//

		void general_commands( Bpar * const pBpar);
		 //this function modifies pedigree info and snpinfo if forced. e.g. pid =fid_pid ;
		 void modify_snp_ped_info( const Bpar * const pBpar);
		 void mach_type_commands(Bpar * const pBpar);
		 void minimac_type_commands( Bpar * const pBpar);
		 void impute_type_commands( Bpar * const pBpar);
		 void snptest_type_commands( Bpar * const pBpar);
		 void plink_type_commands(Bpar* const pBpar);
		 void bimbam_type_commands(Bpar* const pBpar);
		 void beagle_type_commands(Bpar* const pBpar);
		 void rformat_type_commands(Bpar* const pBpar);
		 void shapeit_type_commands(Bpar* const pBpar);
		 void handle_with_force_option(const vector<string>&firstPart, const vector<string>&secondPart,Bpar* const pBpar);
		 void handle_with_filter_snp_option(const Bpar * const pBpar);
		 void handle_with_filter_indiv_option(const Bpar * const pBpar);
		 void assign_snp_indiv_summary();
		 void display_indiv_summary();
		 void display_snp_summary();
		 void read_data(Bpar* const pBpar);
		 static bool are_they_same_snp(const CBSNP* const pFirst, const CBSNP* const pSecond);
		 static bool check_make_same_alleleOrder( CBSNP* const pFirst,  CBSNP* const pSecond);
		 static void check_n_update_snpInfo(CBSNP* const pFirst, const CBSNP*const pSecond);
		 static void check_n_make_same_genoFormat(const CBSNP* const pFirst, CBSNP* const pSecond, const float &thresh);
		 static void write_phase_formatted_files(const vector<CBPED*> &pedVec, const vector<CBSNP*>genVec, const string& outFile);
		 static void read_popInfo_file(const vector<CBPED*> &pedVec,const string& inFileName);
		// const vector<CBPED*>&pedVec, const vector<CBSNP*>&genVec,const string& groupLabelFileName, const string & outFileName
		 CGENERAL();
		 ~CGENERAL();
		 // genotype formats
		 bool given_geno; 	// hard call genotypes
		 bool given_pgeno;	// probabilities of genotypes
		 bool given_geno_0123;
		 bool given_geno_dose;

private:
	//summary of snps
	 int nSNPs;
	 int ngoodSNPs;
	 int nbadSNPs;
	 //summary of samples
	 int nindivs;
	 int ngoodIndivs;
	 int nbadIndivs;
	 int nMiss_status;
	 int nUnaffected ;
	 int nAffected;
	 int nMale;
	 int nFemale;
	 int nNosex;


};
#endif //__GENERAL_H__
