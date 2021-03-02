/*
 * machdat.h
  *  Created on: Dec 29, 2011
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MACHDAT_H__
#define _MACHDAT_H__
#include "helper.h"
#include "gz_stream.h"
using namespace std;
const int coutWidth 		=50;
/*	*while creating new snp
 *  * assign the followings
 *  * geno1
 *  * geno2
 *  * aorder
 *  *allele1
 *  *
 *
 */
 struct CBSNP{
	friend ostream& operator <<(ostream& out, const CBSNP& pCBSNP);
	//static variables
	// other variables
	string nchr; //chromosome numbers
	string snpName;
	string rsId;
	//allele 1 and allele 2 according as they are saved in vectors geno1 and geno2
	string allele1;
	string allele2;
	//strand == 0 -> unknown ('u')
	//strand == 1 -> plus ('+')
	//strand == 2 -> minus ('-')
	int coding_strand;
	// allele 1 and allele2 according as the minor and major allele . In this case minor allele is always the second allele and major is the first.
	string min_allele;
	string maj_allele;
	string dose_allele1;
	string dose_allele2;
	int   nchrobs; // no of non-missing allele count // see in plink frq
	bool change_0123;
	bool change_dose;
	//to determine the quality.
	bool   quality; // to assgin quality =true always quality = quality && true;
	static unsigned int no_of_persons;
	static  unsigned int ngood_snps;	// no of good quality snps.
	double freq; 						// freq of allele 1
	double pvalue_hwe_exact;			// snphwe
	double cm_pos; 						// CM map position
	int bp; 							// base pair positin
	int nonMiss;						 // no of non missing allele
	double snp_CR;						// SNP Call rate

	//the follwoing two vectors save hard call genotypes. i.e.e AA, AB and BB
	vector<bool> geno1; 				// genotype 1
	 vector<bool> geno2; 				// genotype 2
	 vector<bool> aOrder ;				// allele order
	 vector<bool > strand;				 // + = true -=false;+
	 //the follwoing three are the vectors to save probabilites of AA, AB and BB
	 vector<double>pgeno1;
	 vector<double>pgeno2;
	 vector<double>pgeno3;
	  bool given_pgeno;
	 // the follwoing vector geno_0123  save the genotype having format 0,1 and 2 according as the number of minor allele.
	 //i.e. 0 mean no minor allele 1 means one minor and 2 means two minor allele
	 //3 means missing
	 vector<int>geno_0123;
	// bool given_geno_0123;
	 vector<float>geno_dose; // Difficult to save missings. I used -1.0 for missing i.e.
	// bool given_geno_dose;
	 	 	 	 	 	  // if (p(AA)+p(AB)+p(BB)=0.0) then geno_dose[i]=-1.0
	 static const int GENOTYPE_MISSING =3;
	 // geno_0123 is required in the follwoing cases
	 // 1. to calculate maf, to convert into plink-raw, to determine quality.
	 // 2. to covnert intor R format
	//vectors for reading snpinfo file (i.e. extra file to update snp info)
	 static vector<string>sifo_snpid; 	// snp info rsid
	 static vector<string>sifo_rsid; 	// snp info rsid
	 static vector<string> sifo_nchr; 	//chr no
	 static vector<int> sifo_bp; 		//snp info position
	 static vector<double>sifo_cm_pos; 	// snp info cm_pos
	 static vector<string> sifo_allele1;
	 static vector<string> sifo_allele2;
	 static	bool _read_extramap;
	 static  bool _read_extraped;
	 static string 	_miss_gvalue;
	 static string 	_out_format;

	 //vector<CBSNP*>genInfo1;
	 double maf;
 	 double rsq; // for mach
	 //vector<float>freq1;
	CBSNP(): nchr("0"), snpName(""), rsId(""), allele1(""),
			allele2(""),
			coding_strand(0),
			min_allele(""),
			maj_allele(""),
			dose_allele1(""),
			dose_allele2(""),
			nchrobs(0),
			change_0123(true),
			change_dose(true),
			quality(true),
			freq(0.0), pvalue_hwe_exact(0.0),
			cm_pos(-1.0), bp(-1),nonMiss(0),snp_CR(1.0),
		 	given_pgeno(false),
			maf(0.0),rsq(0.0),
			given_both_crate(false),
 	 	 	given_snp_crate(false),
			given_indiv_crate(false),
			given_snp_hwe(false)
 	 	 	//given_geno_0123(false),
 	 	 	//given_geno_dose(false)
 	 	 	 { }
	//function to convert genotypes into 0, 1,2 and 3 according as  minor allele
	  static void convert_genotype_into_0123(const vector<CBSNP*>&genInfo);
	  static float calculate_freq(const CBSNP* pCBSNP);
	  static void calculate_snp_callrate(CBSNP* const pCBSNP);
	  //static void calculate_snp_hwe(CBSNP* const pCBSNP);
	  static double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);
	  static void convert_snp_genotype_into_0123(CBSNP* const pCBSNP);
	  static void convert_snp_genotype_into_dose(CBSNP* const pCBSNP,const string & ref_allele="");
	  static void convert_hardcalls_into_genoprob(CBSNP* const pCBSNP);
	  static void convert_genoprob_into_hardcalls(CBSNP* const pCBSNP,const float&thresh);
	  // plink way.
	static int stringLeser(IFS & pf , string& sValue, bool& rowEnd);
	static string stringAssigntoRightPlace(const bool rowEnd, const int ccol,\
			int & nrow, int& ncol,  const string& stringValue);
	//the following function is extra snp info file.
	static void extraMapFileLeser(const string& fileName, const vector<CBSNP*> & genInfo);
	static void lese_exclude_snps_file(const string& fileName, const vector<CBSNP*> & genInfo);
	static void allele_info_leser(const string& a1, const string& a2,  class CBSNP* it_pCBSNP, const int& nrow, const string& fileName); 

	// writing files:
	static void write_plinkMapFile(const vector<CBSNP*>& genInfo, const string& outFileName);
	static void write_rsid_basePairposition(const vector<CBSNP*>genInfo, const string& outFileName);
	static void write_snpinfoFile( const vector<CBSNP*>& genInfo,const string& outFileName);
	static void write_snplist_file( const vector<CBSNP*>& genInfo,const string& outFileName);
	static void write_bimbam_snp_pos_file(const vector<CBSNP*>& genInfo, const string& outFileName);
	static void write_SNP_frequency_file(const vector<CBSNP*>& genInfo, const string& outFileName);
	static void schreibe_snp_callrate(const vector<CBSNP*>&genInfo,const string & outFileName);
	static void schreibe_ind_callrate(const vector<CBSNP*>&genInfo,const string & outFileName);
	static void schreibe_snp_hwe (const vector<CBSNP*>&genInfo,const string & outFileName);
	//helper function
	static int no_of_good_snps(const vector<CBSNP*> & genVec);
	// this functions includes all commands which are general for all types format 
	
	//
	static void display_snp_summary(const vector<CBSNP*>&genInfo);
	//static vector<CBSNP*>genInfo;
	~CBSNP();
	explicit CBSNP(const CBSNP &); // copy constructor
	CBSNP& operator=(const CBSNP& pCBSNP ); // copy assignment
	static int col_width;
	bool given_both_crate;
	bool given_snp_crate;
	bool given_indiv_crate;
	bool given_snp_hwe;
private:
	void copy(const CBSNP& pCBSNP);



};
struct CMSNP:public CBSNP
{

public:
	CMSNP():CBSNP()
	{

	}
	//vector<CBSNP*>genInfo;
	~CMSNP();
	static void mapFileLeser(const string&  fileName,  vector<CBSNP*>& genInfo);
	static void mach_type_commands(const string & read_code_type);
	static void lese_mach_info_mlinfo_file(const string& fileName, vector<CBSNP*>&genInfo, const float& rsq_value, const float& maf_value=0.0);
	static void read_mach_ref_snpsFile(vector<CBSNP*>& genInfo,const string& snpsFileName);
	static bool _is_code_mach;
};

inline int CBSNP::no_of_good_snps(const vector<CBSNP*> & genVec)
{
	 unsigned long int _tmp_ngood_snps=0;
	 unsigned int j=0;
	 const CBSNP* pCBSNP= genVec[0];
	 for(j=0;j<genVec.size();++j){
		 pCBSNP=genVec[j];
		 if(pCBSNP->quality)
			 ++_tmp_ngood_snps;
	 }

return _tmp_ngood_snps;}

#endif /* _MACHDAT_H__ */
/**  *QUESTION: How Genotypes are saved?
 *	  *if 	homozygote1  then:  FALSE,	FALSE
 	   *if  homozygote2  then: 	TRUE ,	TRUE
 	   *if  heterozygote then: 	FALSE,	TRUE
 	   *if  missing 	 then: 	TRUE, 	FALSE
 */
