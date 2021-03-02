/*
 * bpar.h
 *
 *  Created on: Nov 4, 2011
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 *      These are the basic parameters
 *
 */
#ifndef _BPAR_H__
#define _BPAR_H__
#include "helper.h"
#include<string>
#include<fstream>
#include<string>
#include<algorithm>
#include<vector>
using namespace std;
//vector of merge files;
/**
 *  if you add a variable to Bpar class then also add it to 1. copy() 2. Bpar(). in Bpar.cpp
 *
 */
class Bpar{

public:
	explicit Bpar(const Bpar&);// copy constructor
	Bpar& operator=(const Bpar&pBpar); //copy assignment
	void copy(const Bpar& pBpar); // copy
	Bpar();
	~Bpar();
	int find_single_opt(const string ss); 	// single command string finding
	void find_double_opts(string&coms , const string& s); // single command string finding
	void find_multiple_opts(vector<string>&cur_options , const string& s,const int &n_args); // single command string finding

	bool search_string(const string & s );


	// static variables
	static int  	cpars_sz; // size of vector of basic prameter classes.

	 bool 	given_new_command;
	 bool 	merge_data;
	 // one class is basic. the rest depends on how many files are merged.
	string miss_geno;

	//create variables for commands
	unsigned  int nArgs;
	//vector<string> stcomVec; 		// sorted command vector
	vector<string >original; 		// original command vector ;
	vector<string> options; 		// options to be done ??
	vector<bool>   parsed_comms;	// this will save if all commands in original are parsed.
	//----------------------------------//
	//for plink
	//----------------------------------//
	 bool 		read_ped;
	 bool 		read_map;
	 string		 ped_file;
	 string		 map_file;
	 string		bed_file;
	 string		bim_file;
	 string		fam_file;
	 bool		read_bed;
	 bool		read_bim;

	 bool 	 	three_columns;
	 bool 		 read_plink_dosage;
	 string 	 plink_dosage_file;
	 bool		 read_fam;
	 string		plink_fam_file;
	 bool 		read_plink_rawA;
	 bool		read_plink_rawAD;
	 string 	plink_rawA_file;
	 string 	plink_rawAD_file;
	 bool		write_recodeA_dosage;
	// bool 	write_recodeAD_dosage;
	
	//for covariate file
	 bool		read_covariate_file;
	 string 	covariate_file;
	 bool 		given_covar_names;
	 string 	covar_names;
	 bool 		given_covar_types;
	 string 	covar_types;
	 vector<string> vec_covar_names;
	 vector<string> vec_covar_types;
	
	 string 	file_root;
	 string 	iFormat;
	 string 	oFormat;
	 bool   	change_format;
	 string 	output_fileName;
	 bool 		 uncompressed; // parameters regarding vcf format

	 bool 		hilfe;
	 bool 		is_code_plink;
	//  vector<string> bcVec;
	//----------------------------------//
	// for mach
	//----------------------------------//
	 string 	dat_file;
	 bool 		read_dat;
	 string 	 extraMap_file;
	 bool 		read_extraPed;
	 string 	extraPed_file;
	 string 	mlgeno_file;
	 bool   	read_mlgeno;
	 string 	mlinfo_file;
	 bool 		read_mlinfo;
	 string  	 mach_info_file;
	 bool 		read_mach_info;
	 string 	mlprob_file;
	 bool 		read_mlprob;
	 string 	mach_geno_file;
	 bool 		read_mach_geno;
	 bool 		default_rsq;
	 float 		rsq_value;
	 bool 		default_maf;
	 float 		maf_value;
	 string		temp_string;
	 bool 		is_code_mach;
	 bool 		read_mach_hapmap;
	 string		mach_hapmap_file;
	 bool 		read_mach_hap_snps;
	 string		mach_hap_snp_file;
	//----------------------------------//
	//minimac
	//----------------------------------//
	  string minimac_probFile;
	  bool   read_minimac_probFile;
	  string minimac_infoFile;
	  bool   read_minimac_infoFile;
	  bool   is_code_minimac;
	//----------------------------------//
	//for impute
	//----------------------------------//
	 string 	 gens_file;
	 bool		 read_gens;
	 string 	thresh;
	 bool 	default_thresh;
	 float	 threshold; // later thresh value which is as string will be converted to threshhold
	 bool 	read_impute_info_file;
	 string 	impute_info_file;
	 double 	info_thresh;
	 bool 	default_impute_info;
	 bool 	is_code_impute;
	 bool 	read_impute_hapmap;
	 string	impute_hapmap_file;
	 bool 	read_impute_hap_snps;
	 string	impute_hap_snp_file;
	//----------------------------------//
	//beagle
	//----------------------------------//
	 string bgl_file;
	 bool read_bgl;
	 bool read_bgl_gprobs;
	 bool given_gprobs_header;
	 string bgl_gprobs_file;
	 float bgl_rsq_thresh;
	 string bgl_rsq_file;
	 bool  read_bgl_rsq;
	 bool  is_code_beagle;
	//----------------------------------//
	//bimbam
	//----------------------------------//
	 bool 		read_bimbam;
	 string 	bimbam_file;
	 bool 		read_bimbam_pos;
	 string		bimbam_pos_file;
	 bool 		read_bimbam_gprobs;
	 string		bimbam_gprobs_file;
	 bool		read_bimbam_bestguess;
	 string 	bimbam_bestguess_file;
	 bool		is_code_bimbam;
	//----------------------------------//
	// for SNPTEST
	//-------------------------------------//
	  bool 		read_snptest_sample;
	  string 	snptest_sample_file;
	  bool 		read_sex_col;
	  int 		sex_col;
	  bool 		read_pheno_col;
	  int		pheno_col;
	  bool		 is_code_snptest;
	  //----------------------------------//
	  	// for shapeit
	  //-------------------------------------//
	   bool  read_shapeit_haptype;
	   string haps_file;
	   bool   read_haps_file;
	   string shapeit_sample_file;
	   bool   read_shapeit_sample_file;
	   bool  is_code_shapeit;

	 //-------------------------------------//
	 //for smartpca
	 //-------------------------------------//
	  string  		group_label_fileName;
	  bool			read_group_label;
	//----------------------------------//
	 // for  converting data into simple R files.
	 /**
	  	 *  data will be converted from other format into R files
	  */
	//----------------------------------//
	   bool  read_rgenotype;
	   string rgeno_file;
	   bool  is_code_rformat;

	 //----------------------------------//
	// general
	 //summary statistics
	 bool				write_hwe;
	 bool				write_snp_callrate;
	 bool				write_indv_callrate;
	 bool 				write_freq; // this function is to calculate freq
	 bool 				calculate_freq;
	 bool				calculate_snp_callrate; // to make sure that call rate is not calculated unnecessary times.
	 bool				calculate_ind_callrate; // to make sure that call rate is not calculated unnecessary times.
	 bool				calculate_hardy; // to make sure that hwe is not calculated unnecessary times.
	 bool				filter_snp;	//filtering snps according as callrate maf and hwe
	 bool				filter_indiv;	// filtering snps according as callrate
	 bool				exclude_snps; 	// exclude SNPS given as a list of snps in a separate file
	 bool				remove_indivs;  // remove individuals given as a list of pids in a separate file
	 string				exclude_snplist_file; // saves the file name of snp list to be excluded
	 string				remove_indivlist_file; // saves the file name of indivs to be removed
	 string				filter_snp_opts;
	 string				filter_indiv_opts;
	 vector<string>		filter_snp_firstPart;
	 vector<string>		filter_snp_secondPart;
	 vector<string>		filter_indiv_firstPart;
	 vector<string>		filter_indiv_secondPart;
	// general
	 string 			code_readType; // type can be mach , impute  plink etc.
	 bool				define_readType;
	 bool 				is_general_command; 	// At the moment, general commands are then true if any  of the input files are given.
	 string 			code_read_subType; // type can be mach , impute  plink etc.
	 bool 				read_extraMap; // extra map file

	 bool 				write_snpinfo;
	 bool 				write_pedinfo;
	 bool 				given_force;
	 string 			force_option;
	 vector<string> 	force_first_part;
	 vector<string> 	force_2nd_part;
	 bool				write_pedlist;
	 bool 				write_snplist;
	 bool 				paste_option;
	 vector<string>		strings_2bpasted;// vector of string to be pasted
	 static vector<string> 	pedInfo_idicators; // to save pedinfo indicators.
	 bool			 	paste_iid;
	 vector<int> 		pos_pedinfo_2bpasted;
	 string				paste_seperator;
	 string 			ref_allele;
	 bool 				transpose;

	 // calculation of fst
	 bool 				calculate_fst ;
	// string				popInfoFileName;
	 //spliting the data snpwise
	 vector<int> 		v_ssplitt;
	 vector<int>		v_isplitt;
	 vector<string>		v_split_suffix;
	 bool				ssplitt;
	 bool				isplitt;
	 bool				bpsplitt;



};


inline bool Bpar::search_string(const string & s){vector<string> myvec =original; sort(myvec.begin(),myvec.end());  bool bv = binary_search(myvec.begin(),myvec.end(),s); return bv;}

#endif /* _BPAR_H__ */
