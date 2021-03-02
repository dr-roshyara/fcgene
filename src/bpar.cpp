/*
 * bpar.cpp
 *
 *  Created on: Nov 4, 2011
  *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.

 */

#include "bpar.h"
#include<string>
#include<vector>
using namespace std;

//string  temp_pedinfo[7]={"fid","FID","iid","IID","groupLabel","patid","matid"};
const char* temp_pedinfo[7]={"fid","FID","iid","IID","groupLabel","patid","matid"};
	 vector<string> Bpar::pedInfo_idicators(temp_pedinfo,temp_pedinfo+7);


int  Bpar::cpars_sz				=1; // size of vector of basic parameter classes.
vector<string> pedInfo_idicators;
//Bpar::pedInfo_indicators.assign(temp_pedinfo,temp_pedInfo+7);



Bpar::~Bpar()
{

}
Bpar::Bpar(const Bpar&){ }// copy constructor
Bpar& Bpar::operator=(const Bpar&pBpar){
	copy(pBpar);
return *this;} //copy assignment
void Bpar::copy(const Bpar& pBpar){ //copy

	//stcomVec 				=pBpar.stcomVec;  	// sorted command vector
	original				=pBpar.original; 	// original command vector ;
	options					=pBpar.options; 	// options to be done
	parsed_comms			=pBpar.parsed_comms; // if all commands pare parsed
	given_new_command		=pBpar.given_new_command;
	merge_data				=pBpar.merge_data;
	cpars_sz				=pBpar.cpars_sz;
	//static commands
	file_root				=pBpar.file_root;
	nArgs		 			=pBpar.nArgs;
	cpars_sz				=pBpar.cpars_sz;
	miss_geno				=pBpar.miss_geno;
	iFormat 				=pBpar.iFormat;
	oFormat 				=pBpar.oFormat;
	change_format			=pBpar.change_format;
	output_fileName 		=pBpar.output_fileName;
	hilfe					=pBpar.hilfe;
	code_readType			=pBpar.code_readType;
	code_read_subType		=pBpar.code_read_subType;
	define_readType			=pBpar.define_readType;
	//---------------------------------------------------//
	//for plink
	//---------------------------------------------------//
	read_ped 				=pBpar.read_ped;
	read_map 				=pBpar.read_map;
	ped_file				=pBpar.ped_file;
	map_file				=pBpar.map_file;
	bed_file				=pBpar.bed_file;
	bim_file				=pBpar.bim_file;
	fam_file				=pBpar.fam_file;
	read_bed				=pBpar.read_bed;
	read_bim				=pBpar.read_bim;

	read_plink_dosage		=pBpar.read_plink_dosage;
	plink_dosage_file		=pBpar.plink_dosage_file;
	read_fam				=pBpar.read_fam;
	plink_fam_file			=pBpar.plink_fam_file;
	read_plink_rawA			=pBpar.read_plink_rawA;
	read_plink_rawAD		=pBpar.read_plink_rawAD;
	plink_rawA_file			=pBpar.plink_rawA_file;
	plink_rawAD_file		=pBpar.plink_rawAD_file;
	//write_recodeAD_dosage		=false;
	write_recodeA_dosage	=pBpar.write_recodeA_dosage;
	//
	three_columns			=pBpar.three_columns;
	read_covariate_file		=pBpar.read_covariate_file;
	covariate_file			=pBpar.covariate_file;
	given_covar_names		=pBpar.given_covar_names;
	covar_names				=pBpar.covar_names;
	given_covar_types		=pBpar.given_covar_types;
	covar_types				=pBpar.covar_types;
	vec_covar_names 		=pBpar.vec_covar_names;
	vec_covar_types 		=pBpar.vec_covar_types;
	is_code_plink			=pBpar.is_code_plink;
	//---------------------------------------------------//
	// for mach
	//---------------------------------------------------//
	dat_file				=pBpar.dat_file;
	read_dat 				=pBpar.read_dat;
	read_extraMap			=pBpar.read_extraMap;
	extraMap_file			=pBpar.extraMap_file;
	read_extraPed			=pBpar.read_extraPed;
	extraPed_file 			=pBpar.extraPed_file;
	mlgeno_file 			=pBpar.mlgeno_file;
	read_mlgeno 			=pBpar.read_mlgeno;
	mlprob_file 			=pBpar.mlprob_file;
	read_mlprob 			=pBpar.read_mlprob;
	mach_info_file	 		=pBpar.mach_info_file;
	read_mach_info	 		=pBpar.read_mach_info;
	mlinfo_file 			=pBpar.mlinfo_file;
	read_mlinfo 			=pBpar.read_mlinfo;
	mach_geno_file 			=pBpar.mach_geno_file;
	read_mach_geno  		=pBpar.read_mach_geno;
	default_rsq		  		=pBpar.default_rsq;
	rsq_value 				=pBpar.rsq_value;
	default_maf 	 		=pBpar.default_maf;
	maf_value	 			=pBpar.maf_value;
	temp_string  			=pBpar.temp_string;
	is_code_mach  			=pBpar.is_code_mach;
	read_mach_hapmap	  	=pBpar.read_mach_hapmap;
	mach_hapmap_file	  	=pBpar.mach_hapmap_file;
	read_mach_hap_snps      =pBpar.read_mach_hap_snps;
	mach_hap_snp_file 		=pBpar.mach_hap_snp_file;
	//----------------------------------//
	//minimac
	//----------------------------------//
	minimac_probFile 		=pBpar.minimac_probFile;
	read_minimac_probFile  	=pBpar.read_minimac_probFile;
	minimac_infoFile 		=pBpar.minimac_infoFile;
	read_minimac_infoFile	=pBpar.read_minimac_infoFile;
	is_code_minimac 		=pBpar.is_code_minimac;
	//---------------------------------------------------//
	//for impute
	//---------------------------------------------------//
	gens_file 				=pBpar.gens_file;
	read_gens 				=pBpar.read_gens;
	default_thresh 			=pBpar.default_thresh;
	thresh 					=pBpar.thresh;
	threshold 				=pBpar.threshold;
	//
	//
	read_impute_info_file	=pBpar.read_impute_info_file;
	impute_info_file 		=pBpar.impute_info_file;
	info_thresh  			=pBpar.info_thresh;
	default_impute_info  	=pBpar.default_impute_info;
	is_code_impute 			=pBpar.is_code_impute;
	read_impute_hapmap 		=pBpar.read_impute_hapmap;
	impute_hapmap_file 		=pBpar.impute_hapmap_file;
	read_impute_hap_snps	=pBpar.read_impute_hap_snps;
	impute_hap_snp_file		=pBpar.impute_hap_snp_file;
	//---------------------------------------------------//
	// for SHAPTEIT
	//-------------------------------------//
	 read_shapeit_haptype		=pBpar.read_shapeit_haptype;
	 haps_file					=pBpar.haps_file;
	 read_haps_file				=pBpar.read_haps_file;
	 shapeit_sample_file			=pBpar.shapeit_sample_file;
	 read_shapeit_sample_file	=pBpar.read_shapeit_sample_file;
	 is_code_shapeit			=pBpar.is_code_shapeit;
	//---------------------------------------------------//
	// for SNPTEST
	//-------------------------------------//
	read_snptest_sample 	=pBpar.read_snptest_sample;
	snptest_sample_file 	=pBpar.snptest_sample_file;
	read_sex_col	 		=pBpar.read_sex_col;
	read_pheno_col 			=pBpar.read_pheno_col;
	sex_col 				=pBpar.sex_col;
	pheno_col 				=pBpar.pheno_col;
	is_code_snptest	 		=pBpar.is_code_snptest;
	//---------------------------------------------------//
	//for beagle
	//---------------------------------------------------//
	bgl_file		 		=pBpar.bgl_file;
	read_bgl  				=pBpar.read_bgl;
	bgl_gprobs_file			=pBpar.bgl_gprobs_file;
	read_bgl_gprobs	 		=pBpar.read_bgl_gprobs;
	bgl_rsq_thresh	 		=pBpar.bgl_rsq_thresh;
	bgl_rsq_file  			=pBpar.bgl_rsq_file;
	read_bgl_rsq			=pBpar.read_bgl_rsq;
	is_code_beagle	 		=pBpar.is_code_beagle;
	//---------------------------------------------------//
	//for bimbam
	//---------------------------------------------------//
	read_bimbam	 			=pBpar.read_bimbam;
	bimbam_file		 		=pBpar.bimbam_file;
	given_gprobs_header		=pBpar.given_gprobs_header;
	read_bimbam_gprobs		=pBpar.read_bimbam_gprobs;
	bimbam_gprobs_file  	=pBpar.bimbam_gprobs_file;
	read_bimbam_bestguess	=pBpar.read_bimbam_bestguess;
	bimbam_bestguess_file	=pBpar.bimbam_bestguess_file;
	read_bimbam_pos		 	=pBpar.read_bimbam_pos;
	bimbam_pos_file	 		=pBpar.bimbam_pos_file;
	is_code_bimbam		 	=pBpar.is_code_bimbam;
	//---------------------------------------------------//
	//for r format
	//---------------------------------------------------//
	read_rgenotype			=pBpar.read_rgenotype;
	rgeno_file				=pBpar.rgeno_file;
	is_code_rformat			=pBpar.is_code_rformat;
	//---------------------------------------------------//
	// general
	//---------------------------------------------------//
	ref_allele				=pBpar.ref_allele;
	write_hwe		 		=pBpar.write_hwe;
	write_snp_callrate	 	=pBpar.write_snp_callrate;
	write_indv_callrate	 	=pBpar.write_indv_callrate;
	write_freq		 		=pBpar.write_freq;;
	write_snpinfo	 		=pBpar.write_snpinfo;
	write_pedinfo			=pBpar.write_pedinfo;
	write_pedlist			=pBpar.write_pedlist;
	calculate_snp_callrate	=pBpar.calculate_snp_callrate;
	calculate_ind_callrate	=pBpar.calculate_ind_callrate;
	calculate_hardy	 		=pBpar.calculate_hardy;
	filter_snp		 		=pBpar.filter_snp;
	filter_indiv		 	=pBpar.filter_indiv;
	exclude_snps			=pBpar.exclude_snps; 	// exclude SNPS given as a list of snps in a separate file
	remove_indivs			=pBpar.remove_indivs;  // remove individuals given as a list of pids in a separate file
	exclude_snplist_file	=pBpar.exclude_snplist_file; // saves the file name of snp list to be excluded
	remove_indivlist_file	=pBpar.remove_indivlist_file; // saves the file name of indivs to be removed
	filter_snp_opts	 		=pBpar.filter_snp_opts;
	filter_indiv_opts		=pBpar.filter_indiv_opts;
	// At the moment, general commands are then true if any  of the input files are given.
	is_general_command		=pBpar.is_general_command;
	given_force		 		=pBpar.given_force;
	force_first_part		=pBpar.force_first_part;
	force_2nd_part			=pBpar.force_first_part;
	filter_snp_firstPart 	=pBpar.filter_snp_firstPart;
	filter_snp_secondPart	=pBpar.filter_snp_secondPart;
	filter_indiv_firstPart	=pBpar.filter_indiv_firstPart;
	filter_indiv_secondPart	=pBpar.filter_indiv_secondPart;
	//
	paste_option	 		=pBpar.paste_option;
	strings_2bpasted		=pBpar.strings_2bpasted;
	pedInfo_idicators  		=pBpar.pedInfo_idicators;
	paste_iid		 		=pBpar.paste_iid;
	pos_pedinfo_2bpasted  	=pBpar.pos_pedinfo_2bpasted;
	strings_2bpasted  		=pBpar.strings_2bpasted;
	paste_seperator			=pBpar.paste_seperator;
	v_ssplitt				=pBpar.v_ssplitt;
	v_isplitt				=pBpar.v_isplitt;
	v_split_suffix			=pBpar.v_split_suffix;
	ssplitt					=pBpar.ssplitt;
	isplitt					=pBpar.isplitt;
	bpsplitt				=pBpar.bpsplitt;
	transpose				=pBpar.transpose;

	 uncompressed 			=pBpar.uncompressed;
	//vector<string> Bpar.bcVec.push_back("--ped");
	//---------------------------------------------------//
	//for smartpca
	//---------------------------------------------------//
	read_group_label	 	=pBpar.read_group_label;
	group_label_fileName	=pBpar.group_label_fileName;

	//Bpar-cparse-argv


} // copy

Bpar::Bpar()
{

	file_root					="fcgene";
	nArgs						=0;
	//cpars_sz					=0;
	miss_geno					="0";
	iFormat 					= "";
	oFormat 					= "";
	change_format				=false;
	output_fileName 			= "fcgene_out";
	hilfe						=false;
	code_readType				="";
	code_read_subType			="";
	define_readType				=true;
	//---------------------------------------------------//
	//for plink
	//---------------------------------------------------//
	read_ped 					=false;
	read_map 					= false;
	ped_file					= "fcgene.ped";
	map_file					= "fcgene.map";
	bed_file					="";
	bim_file					="";
	fam_file					="";
	read_bed					=false;
	read_bim					=false;
	read_plink_dosage			=false;
	plink_dosage_file			="";
	read_fam					=false;
	plink_fam_file				="";
	read_plink_rawA				=false;
	read_plink_rawAD			=false;
	plink_rawA_file				="";
	plink_rawAD_file			="";
	//write_recodeAD_dosage		=false;
	write_recodeA_dosage		=false;
	//
	three_columns				=false;
	read_covariate_file			=false;
	covariate_file				="";
	given_covar_names			=false;
	covar_names					="";
	given_covar_types			=false;
	covar_types					="";
	//vec_covar_names;
	//vec_covar_types;
	is_code_plink				=false;
	//---------------------------------------------------//
	// for mach
	//---------------------------------------------------//
	dat_file					= "";
	read_dat					=false;
	read_extraMap				=false; // extra map file
	extraMap_file				="";
	read_extraPed				=false;
	extraPed_file				="";
	mlgeno_file					="";
	read_mlgeno					=false;
	mlprob_file					="";
	read_mlprob					=false;
	mach_info_file				="";
	read_mach_info				=false;
	mlinfo_file					="";
	read_mlinfo					=false;
	mach_geno_file				="";
	read_mach_geno				=false;
	default_rsq					=true;
	rsq_value					=0.0;
	default_maf					=true;
	maf_value					=0.0;
	temp_string					="";
	is_code_mach				=false;
	read_mach_hapmap			=false;
	mach_hapmap_file			="";
	read_mach_hap_snps			=false;
	mach_hap_snp_file			="";
	//----------------------------------//
		 	//minimac
		 //----------------------------------//
	minimac_probFile			="";
	read_minimac_probFile		=false;
	minimac_infoFile			="";
	read_minimac_infoFile		=false;
	is_code_minimac				=false;
	//---------------------------------------------------//
	//for impute
	//---------------------------------------------------//
	gens_file					="";
	read_gens					=false;
	default_thresh				=true;
	thresh						="";
	threshold					= 0.0;
	//
	read_impute_info_file		=false;
	impute_info_file			="";
	info_thresh					=0.0;
	default_impute_info			=true;
	is_code_impute				=false;
	read_impute_hapmap			=false;
	impute_hapmap_file			="";
	read_impute_hap_snps		=false;
	impute_hap_snp_file			="";
	//---------------------------------------------------//
	// for SHAPTEIT
	//-------------------------------------//
	 read_shapeit_haptype		=false;
	 haps_file					="";
	 read_haps_file				=false;
	 shapeit_sample_file			="";
	 read_shapeit_sample_file	=false;
	 is_code_shapeit			=false;
	//---------------------------------------------------//
	// for SNPTEST
	//-------------------------------------//
	read_snptest_sample		=false;
	snptest_sample_file		="";
	read_sex_col				=false;
	read_pheno_col			=false;
	sex_col					=0;
	pheno_col					=0;
	is_code_snptest			=false;
	//---------------------------------------------------//
	//for beagle
	//---------------------------------------------------//
	bgl_file					="";
	read_bgl					=false;
	bgl_gprobs_file			="";
	read_bgl_gprobs			=false;
	bgl_rsq_thresh			=0.0;
	bgl_rsq_file			="";
	read_bgl_rsq			=false;
	is_code_beagle			=false;
	//---------------------------------------------------//
	//for bimbam
	//---------------------------------------------------//
	read_bimbam				=false;
	bimbam_file				="";
	given_gprobs_header		=true;
	read_bimbam_gprobs		=false;
	bimbam_gprobs_file		="";
	read_bimbam_bestguess	=false;
	bimbam_bestguess_file	="";
	read_bimbam_pos			=false;
	bimbam_pos_file			="";
	is_code_bimbam			=false;
	//---------------------------------------------------//
	//for r format
	//---------------------------------------------------//
		read_rgenotype		=false;
		rgeno_file			="";
		is_code_rformat		=false;
	//---------------------------------------------------//
	// general
	//---------------------------------------------------//
	write_hwe				=false;
	write_snp_callrate		=false;
	write_indv_callrate		=false;
	write_freq				=false;
	write_snpinfo			=false;
	write_pedinfo			=false;
	write_pedlist			=false;
	write_snplist			=false;
	calculate_freq			=true;//
	calculate_snp_callrate	=true; // to make sure that call rate is not calculated unnecessary times.
	calculate_ind_callrate	=true; // to make sure that call rate is not calculated unnecessary times.
	calculate_hardy			=true; // to make sure that hwe is not calculated unnecessary times.
	filter_snp				=false;
	filter_indiv			=false;
	exclude_snps			=false;
	remove_indivs			=false;
	exclude_snplist_file	=""; // saves the file name of snp list to be excluded
	remove_indivlist_file	=""; // saves the file name of indivs to be removed
	filter_snp_opts			="";
	filter_indiv_opts			="";
	// At the moment, general commands are then true if any  of the input files are given.
	is_general_command		=false;
	given_force				=false;
	force_option				="";
	//force_first_part;
	//force_2nd_part;
	paste_option			=false;
	//strings_2bpasted;			// vector of string to be pasted
	paste_iid				=false;
	//pos_pedinfo_2bpasted;
	paste_seperator			="";
	uncompressed			=false;
	vector<int> 			v_ssplitt;
	vector<int>				v_isplitt;
	vector<string>			v_split_suffix;
	ssplitt					=false;
	isplitt					=false;
	bpsplitt				=false;
	transpose				=false;

	//vector<string> Bpar.bcVec.push_back("--ped");
	//---------------------------------------------------//
	//for smartpca
	//---------------------------------------------------//
	read_group_label			=false;
	group_label_fileName		="";
	ref_allele					="";
	//Bpar-cparse-argv
	// calculation of fst
	calculate_fst 				=false ;
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
}

int Bpar::find_single_opt(const string ss)
{
	int dis=0;
	vector<string> _tmp_svec =original;
	bool tfValue =search_string(ss);// binary_search(stcomVec.begin(),stcomVec.end(),s);
	if(tfValue)
	{
		vector<string>::iterator it;
		it = find (_tmp_svec.begin(), _tmp_svec.end(), ss); // finding the string in the argv
		dis = int(distance(_tmp_svec.begin(),it)); // index of the string
		parsed_comms[dis]=tfValue;
		//also do it for the original command
	}
	else
	{
		string msg = "An argument for "+ss+ " is missing.\n";
		error(msg);
	}

			return dis;
}
// double command string finding
 void Bpar::find_double_opts(string&coms , const string& s)
 {
	 vector<string> _tmp_svec 	=original;
	 bool tfValue 				=search_string(s);// binary_search(stcomVec.begin(),stcomVec.end(),s);
	if(tfValue)
	{
		 vector<string>::iterator it;
		 it = find (_tmp_svec.begin(), _tmp_svec.end(), s); // finding the string in the argv
		// for(unsigned int i=0; i<_tmp_svec.size();++i)
		// cout << _tmp_svec[i]<<" ";
		// cout <<endl;
		 unsigned int dis = int(distance(_tmp_svec.begin(),it)); // index of the string
		 //cout << "first: "<<dis << endl;
		 parsed_comms[dis]=tfValue;//after finding the term make sure that this  term has been found out.

		 ++it;
		 dis = int(distance(_tmp_svec.begin(),it)); // looking for the next term
		 //cout << "second: "<<dis << endl;
		if(dis<nArgs && _tmp_svec[dis].substr(0,1)=="-")
		{
			string msg = "Are you sure about the argument after "+ s + "? If yes, we suggest you to rename the argument: " + original[dis]+".\n";
			error(msg); // is giving warning message: control reaches end of non-void function .

		}
		//cout <<  ->nArgs<< endl;
		//cout<< (dis<nArgs)<< ", "<<_tmp_svec[dis]<<endl;
		if((dis<nArgs) && _tmp_svec[dis].substr(0,2)!="--")
		{
			 parsed_comms[dis]=tfValue;  // next value is also true. e.g. --ped test.ped  both strings should exit there.
			 coms = _tmp_svec[dis];
		//	cout <<"test"<<endl;
		}
		else
		{
			string msg = "An argument after \""+ s + "\" is missing.\n";
			error(msg); // is giving warning message: control reaches end of non-void function .
			//return (msg);
			return;
		}
	}
	else{
		string msg = "An argument after "+s+ " is missing.\n";
		error(msg);
	}
	 //return (msg);
	_tmp_svec.clear();
	_tmp_svec.resize(0);
	return;
 }
 //for multiple options
void Bpar::find_multiple_opts(vector<string>&cur_options , const string& s, const int &n_args)
{
	// single command string finding
	 vector<string> _tmp_svec 	=original;
	 string coms				="";
	 bool tfValue 				=search_string(s);// binary_search(stcomVec.begin(),stcomVec.end(),s);
	 int _cur_args				=1;
	 if(tfValue)
	 	{
	 		 vector<string>::iterator it;
	 		 it = find (_tmp_svec.begin(), _tmp_svec.end(), s); // finding the string in the argv
	 		// for(unsigned int i=0; i<_tmp_svec.size();++i)
	 		// cout << _tmp_svec[i]<<" ";
	 		// cout <<endl;
	 		 unsigned int dis = int(distance(_tmp_svec.begin(),it)); // index of the string
	 		 //cout << "first: "<<dis << endl;
	 		 parsed_comms[dis]=tfValue;//	after finding the term make sure that this  term has been found out.

	 		 ++it; //

	 		 dis = int(distance(_tmp_svec.begin(),it)); // looking for the next term
	 		 while(_cur_args<n_args)
	 		 {
	 			 cout << _cur_args<<"th term: "<<dis << endl;
	 			 if(dis<nArgs){
	 				 string _msg="There are not enough arguments in the command option.\n";
	 				 error(_msg);
	 			 }
	 			 if( _tmp_svec[dis].substr(0,1)=="-")
	 			 {
	 				 string msg = "Are you sure about the argument after "+ s + "? If yes, we suggest you to rename the argument: " + original[dis]+".\n";
	 				 error(msg); // is giving warning message: control reaches end of non-void function .

	 			 }
	 			 //cout <<  ->nArgs<< endl;
	 			 //cout<< (dis<nArgs)<< ", "<<_tmp_svec[dis]<<endl;
	 			 if(_tmp_svec[dis].substr(0,2)!="--")
	 			 {
	 				 parsed_comms[dis]=tfValue;  // next value is also true. e.g. --ped test.ped  both strings should exit there.
	 				 coms = _tmp_svec[dis];
	 				 if(coms!="")
	 					cur_options.push_back(coms);
	 				 //	cout <<"test"<<endl;
	 			}
	 			else
	 			{
	 				string msg = "An argument after \""+ s + "\" is missing.\n";
	 				error(msg); // is giving warning message: control reaches end of non-void function .
	 			}
	 			 ++_cur_args;
	 			 ++dis;
	 		 }//end of while loop

	 	}//end of finding
	 	else{
	 		string msg = "An argument after "+s+ " is missing.\n";
	 		error(msg);
	 	}
	 	 //return (msg);
	 	_tmp_svec.clear();
	 	_tmp_svec.resize(0);
	 	return;

}


/*
//bool Bpar::dummy = false;
	//bool Bpar::myfunction = false;
	// unsigned  int Bpar::nArgs				=0;
	 string Bpar::miss_geno					="0";
	 string Bpar::iFormat 					= "";
	 string Bpar::oFormat 					= "";
	 bool  	Bpar::change_format				=false;
	 string Bpar::output_fileName 			= "fcgene_out";
	 bool 	Bpar::hilfe						=false;
	 string Bpar::code_readType				="";
	 string Bpar::code_read_subType			="";
	 bool	Bpar::define_readType			=true;
	 //---------------------------------------------------//
	 //for plink
	 //---------------------------------------------------//
	 bool 	Bpar::read_ped 					=false;
	 bool 	Bpar::read_map 					= false;
	 string Bpar::ped_file					= "fcgene.ped";
	 string Bpar::map_file					= "fcgene.map";

	 bool 	Bpar::read_plink_dosage			=false;
	 string Bpar::plink_dosage_file			="";
	 bool	Bpar::read_fam					=false;
	 string Bpar::plink_fam_file			="";
	 bool 	Bpar::read_plink_rawA			=false;
	 bool	Bpar::read_plink_rawAD			=false;
	 string Bpar::plink_rawA_file			="";
	 string Bpar::plink_rawAD_file			="";
	 //bool	Bpar::write_recodeAD_dosage		=false;
	 bool	Bpar::write_recodeA_dosage		=false;
	 //
	 bool 	Bpar::three_columns				=false;
	 bool 	Bpar::read_covariate_file		=false;
	string  Bpar::covariate_file			="";
	bool 	Bpar::given_covar_names			=false;
	string 	Bpar::covar_names				="";
	bool 	Bpar::given_covar_types			=false;
	string 	Bpar::covar_types				="";
	vector<string>Bpar::vec_covar_names;
	vector<string>Bpar::vec_covar_types;
	bool 	Bpar::is_code_plink				=false;
	
	 //---------------------------------------------------//
	 // for mach
	 //---------------------------------------------------//
	 string Bpar::dat_file					= "";
	 bool 	Bpar:: read_dat					=false;
	 bool 	Bpar::read_extraMap				=false; // extra map file
	 string Bpar::extraMap_file				="";
	 bool 	Bpar::read_extraPed				=false;
	 string Bpar::extraPed_file				="";
	 string Bpar::mlgeno_file				="";
	 bool 	Bpar::read_mlgeno				=false;
	 string Bpar::mlprob_file				="";
	 bool	Bpar::read_mlprob				=false;
	 string Bpar::mach_info_file			="";
	 bool	Bpar::read_mach_info			=false;
	 string Bpar::mlinfo_file				="";
	 bool 	Bpar::read_mlinfo				=false;
	 string Bpar::mach_geno_file			="";
	 bool   Bpar::read_mach_geno			=false;
	 bool  	Bpar::default_rsq				=true;
	 float 	Bpar::rsq_value					=0.0;
	 bool  	Bpar::default_maf				=true;
	 float 	Bpar::maf_value					=0.0;
	 string Bpar::temp_string				="";
	 bool 	Bpar::is_code_mach				=false;
	 bool   Bpar::read_mach_hapmap			=false;
	 string Bpar::mach_hapmap_file			="";
	 bool 	Bpar::read_mach_hap_snps		=false;
	 string Bpar::mach_hap_snp_file			="";
	 //----------------------------------//
	 	//minimac
	 //----------------------------------//
	 string Bpar::minimac_probFile			="";
	 bool   Bpar::read_minimac_probFile		=false;
	 string Bpar::minimac_infoFile			="";
	 bool   Bpar::read_minimac_infoFile		=false;
	 bool 	Bpar::is_code_minimac			=false;
	 //---------------------------------------------------//
	 //for impute
	 //---------------------------------------------------//
	 string Bpar::gens_file					="";
	 bool 	Bpar::read_gens					=false;
	 bool 	Bpar::default_thresh			=true;
	 string	Bpar::thresh					="";
	 float 	Bpar::threshold					= 0.0;
	 //
	 bool 	Bpar::read_impute_info_file		=false;
	string 	Bpar::impute_info_file			="";
	double 	Bpar::info_thresh				=0.0;
	bool 	Bpar::default_impute_info		=true;
	bool 	Bpar::is_code_impute			=false;
	bool 	Bpar::read_impute_hapmap		=false;
	string	Bpar::impute_hapmap_file		="";
	bool 	Bpar::read_impute_hap_snps		=false;
	string 	Bpar::impute_hap_snp_file		="";
	//---------------------------------------------------//
	// for SNPTEST
	//-------------------------------------//
	 bool 	Bpar:: read_snptest_sample		=false;
	 string Bpar::snptest_sample_file		="";
	 bool 	Bpar::read_sex_col				=false;
	 bool 	Bpar::read_pheno_col			=false;
	 int 	Bpar::sex_col					=0;
	 int  	Bpar::pheno_col					=0;
	 bool 	Bpar::is_code_snptest			=false;
	 //---------------------------------------------------//
	 //for beagle
	 //---------------------------------------------------//
	 string Bpar::bgl_file					="";
	 bool 	Bpar::read_bgl					=false;
	 string Bpar::bgl_gprobs_file			="";
	 bool 	Bpar::read_bgl_gprobs			=false;
	 float	Bpar::bgl_rsq_thresh			=0.0;
	 string Bpar::bgl_rsq_file				="";
	 bool 	Bpar::read_bgl_rsq				=false;
	 bool  Bpar::is_code_beagle				=false;
	 //---------------------------------------------------//
	 //for bimbam
	 //---------------------------------------------------//
	 bool 	Bpar::read_bimbam				=false;
	string 	Bpar::bimbam_file				="";
	bool	Bpar::given_gprobs_header		=true;
	bool	Bpar::read_bimbam_gprobs		=false;
	string 	Bpar::bimbam_gprobs_file		="";
	bool	Bpar::read_bimbam_bestguess		=false;
	string 	Bpar::bimbam_bestguess_file		="";
	bool   	Bpar::read_bimbam_pos			=false;
	string 	Bpar::bimbam_pos_file			="";
	bool 	Bpar::is_code_bimbam			=false;
	//---------------------------------------------------//
	// general
	//---------------------------------------------------//
	bool	Bpar::write_hwe					=false;
	bool	Bpar::write_snp_callrate		=false;
	bool	Bpar::write_indv_callrate		=false;
	bool	Bpar::write_freq				=false;
	bool 	Bpar::write_snpinfo				=false;
	bool 	Bpar::write_pedinfo				=false;
	bool	Bpar::write_pedlist				=false;
	bool	Bpar::write_snplist				=false;
	bool 	Bpar::calculate_freq			=true;//
	bool	Bpar::calculate_snp_callrate	=true; // to make sure that call rate is not calculated unnecessary times.
	bool	Bpar::calculate_ind_callrate	=true; // to make sure that call rate is not calculated unnecessary times.
	bool	Bpar::calculate_hardy			=true; // to make sure that hwe is not calculated unnecessary times.
	bool 	Bpar::filter_snp				=false;
	bool 	Bpar::filter_indiv				=false;
	string	Bpar::filter_snp_opts			="";
	string  Bpar::filter_indiv_opts			="";
	vector<string>	Bpar::filter_snp_firstPart;
	vector<string>	Bpar::filter_snp_secondPart;
	vector<string>	Bpar::filter_indiv_firstPart;
	vector<string>	Bpar::filter_indiv_secondPart;

	// At the moment, general commands are then true if any  of the input files are given.
	bool 	Bpar::is_general_command		=false;
	bool 	Bpar::given_force				=false;
	string  Bpar::force_option				="";
	vector<string> 	Bpar::force_first_part;
	vector<string>	Bpar::force_2nd_part;
	bool 	Bpar::paste_option				=false;
	vector<string>	Bpar::strings_2bpasted;// vector of string to be pasted
	 const char* temp_pedinfo[7]={"fid","FID","iid","IID","groupLabel","patid","matid"};
	 vector<string> Bpar::pedInfo_idicators(temp_pedinfo,temp_pedinfo+7);
	 	bool	Bpar::paste_iid			=false;
	 vector<int> Bpar::pos_pedinfo_2bpasted;
	string	Bpar::paste_seperator			="";
	 //vector<string> Bpar.bcVec.push_back("--ped");
	//---------------------------------------------------//
	//for smartpca
	//---------------------------------------------------//
	bool Bpar::read_group_label				=false;
	string Bpar::group_label_fileName		="";
//Bpar-cparse-argv
*/
