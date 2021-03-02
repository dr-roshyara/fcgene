/*
 * machPedgree.h
 *
 *  Created on: Dec 21, 2011
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef _MACHPEDGREE_H__
#define _MACHPEDGREE_H__
#include "machdat.h"
#include "helper.h"
#include "gz_stream.h"

#include<new>
using namespace std;
extern ofstream LIN;
/**	while creating a new PED file  what should be given at least ?
 * 	1. assign sex
 * 	2. assign status
 * 	3.
 *
 */
struct CBPED
{ // CBPED = Class MACH PED
	bool sex;
	bool miss_sex;
	float age ;
	//if male, sex=TRUE;miss_sex=false, if female sex=false miss_sex=false,
	//	if not defined then sex=fasle, miss_sex=true;
	bool pheno;	 // disease pheno.
	bool miss_pheno;
	// if pheno= disease, then pheno=true, miss_state=false, if no diseases, then pheno=false; miss_state=false
	// if undefined, then pheno=false, miss_state=true;
	double ind_CR; // person_call_rate
	bool 	quality ; // to assgin quality =true always quality = qualtiy && true;
	string indId;
	string famId;
	string matId;
	string patId;
	string groupLabel;
	static vector<string> save_line;
	vector<string> covariate;
	static const unsigned  int coutWidth; // width of cout
	static unsigned int nMale;
	static unsigned int nFemale;
	static 	unsigned int nNosex;
	static	unsigned int nAffected; // no of unaffected
	static	unsigned int nUnaffected; // no of affected
	static 	unsigned int nMiss_status; //  undefined status
	static  unsigned int ngood_snps;	// no of good quality snps.
	static  vector<string> covariate_names; // read form file
	static  vector<string> _cov_names_frm_argv; // covaraite names given from argv
	static  vector<string> _cov_types_frm_argv; // covaraite types given from argv
	static  bool	 _given_cov_names;
	static  bool	 _given_cov_types;
	static 	vector<CBPED*> vPed;
	// defining necessary functions
		explicit CBPED(const CBPED& pBPED); // copy constructor
	 	CBPED & operator=(const CBPED& pCBPED); // copy assignment
	// if pheno_status is given , then please don't forget to make miss_pheno=false;
	explicit CBPED():sex(false),
					miss_sex(true),
					age(0.0),
					pheno(false),
					miss_pheno(true),
					ind_CR(1.0),
					quality(true),
					indId(""),famId(""), matId(""),patId(""),
					groupLabel(""),
					given_indiv_crate(false)
					{

					}

	// friend functions

	friend int CBSNP::stringLeser(IFS &pf , string& sValue, bool& rowEnd);
	friend  string CBSNP::stringAssigntoRightPlace(const bool rowEnd, const int ccol,\
			int & nrow, int& ncol, const string& stringValue);
	//
	~CBPED();
	//----------------------------------------------------------------//
	// functions to read
	//----------------------------------------------------------------//
	static void genInfoLeser(const vector<CBSNP*>& genInfo, const string& first , const string& second, const int & ncol , const int& nrow, const int& cSNPN );
	static void pedinfoFileLeser(const vector<CBPED*> & pedInfo,const string& fileName);
	//alternatives this is for mach and ped both
	static void lese_pedfile( vector<CBPED*>& pedInfo, const vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName);
	// new functions to read ped file
	static CBPED* lese_pedfile_pedInfo(const vector<string>& zeilenVec, const long int& snpAnzahl, const int& akt_zeile);
	static bool lese_pedfile_zeile(IFS  &pf, vector<string>& pZeile); // reading end of line  returns true;
	static void lese_pedfile_genInfo(const vector<CBSNP*>&genInfo, const vector<string>& zeilenVec,const int akt_zeile);
	static void lese_pedfile_all(vector<CBPED*> &vPedInfo,	const vector<CBSNP*>&genInfo, vector<string>&zeilenVector,const string& pedFileName);
	static void update_pedInfo_given_gProb(vector<CBPED*>&pedInfo,const vector<CBSNP*>genInfo, const bool update_sex=0, const bool update_pheno=0);
	static void lese_plink_covariate_file(const vector<CBPED*>& pedInfo,const string& fileName);
	static void lese_remove_indivs_file(const string& fileName,  const vector<CBPED*> & pedInfo);
	//functions to calculate
	static void calculate_indiv_callrate(const vector<CBPED*>&pedInfo, const vector<CBSNP*> &genInfo);
	static void calculate_snp_hwe(const vector<CBSNP*> & genVec,const vector<CBPED*>& vPed);
	static void calculate_maf( const vector<CBSNP* >& genVec, const vector<CBPED*> &vped);

	//----------------------------------------------------------------//
	// functions to write
	//----------------------------------------------------------------//
	static void schreibe_ind_callrate(const vector<CBSNP*>&genInfo,const vector<CBPED*>&pedInfo, const string & outFileName);
	static void write_plinkPedFile(const vector<CBPED*>& pedInfo, const vector<CBSNP*> & genInfo, const string&outFileName);
	static void writeOutPlinkFiles(const vector<CBPED*>& pedInfo, const vector<CBSNP*> & genInfo, const string&outFileName);
	static void write_impute_gensFile(const vector<CBSNP*>& genInfo,const vector<CBPED*>&pedInfo,
				const string& outFileName,
				const bool _given_extra_snpinfo
				);
	static void write_snptest_sample_file(const vector<CBPED*>&pedInfo, const string& outFileName);
	static void write_beagleGenotypeFile(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName);
	static void write_bimbam_geno_file(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName);
	static void write_mach_files(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName);
	static void write_pedinfoFiles(const vector<CBPED*>& pedInfo, const string& outFileName);
	static void write_pedlist_file(const vector<CBPED*>& pedInfo, const string& outFileName);
	static void write_plink_dosagesFiles(const vector<CBPED*>&pedInfo,const vector<CBSNP*>genInfo, const string& outFileName);
	///----------------------------------------------------------------//
	//helper functions
	//----------------------------------------------------------------//
	static void genoProb_to_geno_converter(const vector<CBSNP*>& genInfo,const float& thresh);
	static void lese_sex_info( CBPED* const pCBPED, const string sexValue, const bool update_sex=0);
	static void lese_phenotype_info(CBPED*const pCBPED, const string pheno, const bool update_pheno=0);
	static void add_sex_n_phenotype_info(CBPED* const pCBPED, const bool update_sex,const bool update_pheno);

	static bool lese_genotype_info(CBSNP*pCBSNP, const string first, const  int akt_zeile, string&ERR_MSG);
	static bool  lese_genotype_info(CBSNP*pCBSNP, const string first, const string second, const  int akt_zeile, string&ERR_MSG);
	static void display_ped_summary(void);
	static void cerr_with_delete(vector<CBPED*>& vPed, vector<CBSNP*>&genInfo, const string& cerr_message);
	//functions related to computation
	static void calculate_callrate(const vector<CBPED*>& pedInfo,const vector<CBSNP*> & genInfo);
	//this is a function will includes general stuff.
	static int  no_of_good_indivs(const vector<CBPED*> &vPed );
	

private:
 	void copy(const CBPED & pCBPED);
 	bool given_indiv_crate;

};

class CMPED:public CBPED
{
public:

	//vector<CBPED*> vPed;  // vector of CBPED;+
	//functions
	CMPED():CBPED() {}
	static void pedFileleser(const string& pedFileName,	 vector<CBPED*> &vPedInfo, const vector<CBSNP*>&genInfo);
	//static void writeOutPlinkFiles(const vector<CBPED*>pedInfo, const vector<CBSNP*> & genInfo, const string&outFileName);
	//static  void lese_mach_mlGenoFile( vector<CBPED*>& pedInfo, const vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName);
	static  void lese_mach_mlProbFile( vector<CBPED*>& pedInfo, const vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName, const float & thresh);
	static void lese_mach_geno_mlGenoFile( vector<CBPED*>& pedInfo, const vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName);
	static void lese_mach_hapmap_file(vector<CBPED*>&vPed,vector<CBSNP*>& genInfo,vector<string>& cZeile,const string& fileName);
	// to write out mach references files
	static void write_mach_ref_files(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName);

	//static void lese1_mach_info_mlinfo_file(const string& fileName, vector<CBSNP*>&genInfo);
public:
	//static void zeilenLeser(const string & , vector<CBPED*>& vPed );
	~CMPED();
};

void fix_ints_4_impute_command(const vector<int>& bp, vector<int>&lower_range_vec, vector<int>&upper_range_vec,vector<float> & interval_size_MB);
inline int CBPED::no_of_good_indivs( const vector<CBPED*>& vPed )
{
	 unsigned long int _tmp_ngindiv=0;
	 unsigned int i=0;
	 const CBPED * pCBPED=vPed[0];
	 for(i=0;i<vPed.size();++i)
	 {
		 pCBPED= vPed[i];
		 if(pCBPED->quality)
			 ++_tmp_ngindiv;
	 }

return _tmp_ngindiv;}
#endif /* _MACHPEDGREE_H__ */
