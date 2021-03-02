/*
 * 	Author: Nab Raj Roshyara

 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */
#include"general.h"
#include "minimac.h"
#include "beagle.h"
#include "bimbam.h"
#include "rformat.h"
#include "haploview.h"
#include  "plinksnp.h"
#include "plinkped.h"
#include "readimpute.h"
#include "smartpca.h"
#include "genable.h"
#include "bpar.h"
#include "matrix_def.h"
#include<cmath>
//bool CGENERAL::_read_extramap=false; //later they will come into private
//bool CGENERAL::_read_extraped=false; // later they will come into private
//string CGENERAL::_miss_gvalue="";
//string CGENERAL::_out_format="";
//vector<CBPED*> CGENERAL::pedVec;
//vector<CBSNP*> CGENERAL::genVec;
CGENERAL::CGENERAL()
{
	 	 nindivs		=0;
	 	 ngoodIndivs	=0;
	 	 nbadIndivs		=0;
		 nMiss_status	=0;
		 nUnaffected	=0 ;
		 nAffected		=0;
		 nMale			=0;
		 nFemale		=0;
		 nNosex			=0;
		 //
		 ngoodSNPs		=0;
		 nbadSNPs		=0;
		 nSNPs			=0;

		 given_geno		=false;
		 given_pgeno	=false;

}
CGENERAL::~CGENERAL()
{ 

if(pedVec.size()>0)
	{
		for(vector<CBPED*>::iterator it=pedVec.begin();it!=pedVec.end();++it)
				delete *it;
		
		pedVec.clear();
		pedVec.resize(0);
	}	
	
	if(genVec.size()>0)
	{
			for(vector<CBSNP*>::iterator it=genVec.begin();it!=genVec.end();++it)
				delete *it;
		
		genVec.clear();
		genVec.resize(0);
	}	

	//cout <<"Hu Hu from CGENERAL! I have deleted your genVec and pedVec !! \n";
}

//CBSNP::_miss_gvalue	=pBpar->miss_value;
void CGENERAL::assign_snp_indiv_summary()
{
	nindivs					=pedVec.size();
	nMiss_status			=CBPED::nMiss_status;
	nUnaffected				=CBPED::nUnaffected;
	nAffected				=CBPED::nAffected;
	nMale					=CBPED::nMale;
	nFemale					=CBPED::nFemale;
	nNosex					=CBPED::nNosex;
	//snp
	CBPED::nMiss_status		=0;
	CBPED::nUnaffected		=0;
	CBPED::nAffected		=0;
	CBPED::nMale			=0;
	CBPED::nFemale			=0;
	CBPED::nNosex			=0;

	//Assign SNP  summary
	//cout <<" genVec.size(): "<< genVec.size() <<endl;
	int good_snps	=0;
	int bad_snps	=0;
	for(vector<CBSNP*>::const_iterator it=genVec.begin(); it<genVec.end(); ++it )
	{

		if((*it)->quality)
			++good_snps;
		//cout << boolalpha << (*it)->quality << " "; //debug
	}
	bad_snps	=(genVec.size()-good_snps);
	ngoodSNPs	=good_snps;
	nbadSNPs 	=bad_snps;
	nSNPs	 	=genVec.size();
	//assign good indivs and nbad indivs
	for(vector<CBPED*>::const_iterator it=pedVec.begin(); it<pedVec.end(); ++it )
	{
		//cout << boolalpha <<"badsnp: "<< ((*it)->quality)<<endl;
		if(!((*it)->quality))
			++nbadIndivs;
	}
	//cout << "nbadIndivs: "<<nbadIndivs <<endl;
	ngoodIndivs=(pedVec.size()-nbadIndivs);

}
void CGENERAL::display_indiv_summary()
{
		// count affected and unnaffected
			cout <<left;
			cout <<"*->Individual summary:\n";
			// male and female Geschiste
			cout <<"\t"<< setw(coutWidth)<< 	"**->Total female: " <<nFemale<<"\n";
			cout <<"\t"<<setw(coutWidth)<< 	"**->Total male: " <<nMale<< "\n";
			cout <<"\t"<<setw(coutWidth)<< 	"**->Total individuals with undefined sex: " <<nNosex<< "\n";
			cout.flush();
			//write to fcgene.out.txt file
			LIN << left;
			LIN <<						 	"\n*->Individual summary:\n";
			LIN <<"\t"<< setw(coutWidth)<< 		"**->Total females: " << nFemale<<"\n";
			LIN <<"\t"<<setw(coutWidth)<< 			"**->Total males: " << nMale<< "\n";
			LIN <<"\t"<<setw(coutWidth)<< 			"**->Total individuals with undefined sex: " <<nNosex<< "\n";
			LIN.flush()<<endl;
			// for affected and unaffected.
			cout <<"\t"<< setw(coutWidth)<< 	"**->Total cases: " << nAffected<<"\n";
			cout <<"\t"<<setw(coutWidth)<< 	"**->Total controls: " << nUnaffected<< "\n";
			cout <<"\t"<<setw(coutWidth)<< 	"**->Total individuals with undefined status: " <<nMiss_status<< "\n";
			cout.flush()<<endl;
			//write to fcgene.out.txt file
			LIN <<"\t"<< setw(coutWidth)<< 	"**->Total cases: " << nAffected<<"\n";
			LIN<<"\t"<<setw(coutWidth)<< 	"**->Total controls: " <<nUnaffected<< "\n";
			LIN<<"\t"<<setw(coutWidth)<<	"**->Total individuals with undefined status: " <<nMiss_status<< "\n";
			LIN.flush();
			if(nbadIndivs>0)
			{
				printLIN("\t**->There are some removed (or bad quality)individuals in the data.\n"
						"\t**->Individuals marked as removed or bad quality are not included  in data transformation.\n"
						"\t**->Total included (or good quality)individuals: "+change_int_into_string(ngoodIndivs)+".\n"
						"\t**->Total removed (or bad quality)individuals: "+change_int_into_string(nbadIndivs)+".\n"
						);
			}

}
void CGENERAL::display_snp_summary()
{

	if(nbadSNPs>0)
	{
			cout <<left;
			cout <<"\n*->SNP summary:\n";
			// male and female Geschiste
			cout <<"\t"<< setw(coutWidth)<<"**->Total included (or good Quality) SNPs: " << ngoodSNPs<<"\n";
			cout <<"\t"<< setw(coutWidth)<<"**->Total excluded (or bad Quality) SNPs: " << nbadSNPs<<"\n";
			cout <<"\t"<< setw(coutWidth)<<"**->SNPs marked as excluded (or bad Quality) will not be considered for Format converting process.\n";
			cout.flush()<<endl;
			//
			LIN <<left;
			LIN <<"\n*->SNP summary:\n";
			// male and female Geschiste
			LIN <<"\t"<< setw(coutWidth)<< 	"**->Total included (or good Quality) SNPs: " << ngoodSNPs<<"\n";
			LIN <<"\t"<< setw(coutWidth)<< 	"**->Total excluded (or bad Quality) SNPs: " << nbadSNPs<<"\n";
			LIN <<"\t"<< setw(coutWidth)<< 	"**->SNPs marked as excluded (or bad Quality) will not be considered any more for format converting process.\n";
			LIN.flush()<<endl;
	}
	else{

		cout <<left;
		cout <<"\n*->SNP summary:\n";
		cout <<"\t"<< setw(coutWidth)<<	"**->Total number of successfully read SNPs: " << ngoodSNPs<<"\n";
		cout.flush()<<endl;
		LIN <<"\t"<< setw(coutWidth)<< "**->Total number of successfully read: " << ngoodSNPs<<"\n";
		LIN.flush()<<endl;
	}


}
void CGENERAL::general_commands(Bpar* const pBpar)
{

	/**	*General commands :

	    * reading extra snpinfo file or ped info file
	    * displaying summary
	    * writing snp info and ped info
	    * converting format of the data
	    *
	 */
	const string _output_filename =pBpar->output_fileName;
	/*CBSNP::_read_extramap	=pBpar->read_extraMap;
	CBSNP::_read_extraped	=pBpar->read_extraPed;

	if(CBSNP::_read_extramap)
		CBSNP::extraMapFileLeser(pBpar->extraMap_file, genVec);
	if(CBSNP::_read_extraped)
		CBPED::pedinfoFileLeser(pedVec,pBpar->extraPed_file);
    				//CBPED::displayIndividualSummary(pedVec);
    				//here starts the data management of mach format
	modify_snp_ped_info(pBpar);
	if(pBpar->given_force)
		handle_with_force_option(pBpar->force_first_part,pBpar->force_2nd_part);
	if(pBpar->read_covariate_file)
		CBPED::lese_plink_covariate_file(pedVec,pBpar->covariate_file);
	*/
	 //assign_snp_indiv_summary(); // pedVec
	 //display_snp_summary();
 	 //display_indiv_summary();

	// After here to writing different files: No updates are done. If necessary, do make the
	// updates before this position.
	//cout<<boolalpha <<"test"<< pBpar->write_pedinfo<<endl;
	if(pBpar->write_snpinfo)
		CBSNP::write_snpinfoFile(genVec,_output_filename);
    if(pBpar->write_pedinfo)
    	CBPED::write_pedinfoFiles(pedVec, _output_filename);
     // Here starts the process of changing  format to other software
    if(pBpar->write_pedlist)
    	CBPED::write_pedlist_file(pedVec,_output_filename);
    if(pBpar->write_snplist)
    	CBSNP::write_snplist_file(genVec,_output_filename);
    //writing freq
    if(pBpar->write_freq)
    {
    	printLIN("\n*->Calculating allele frequencies: \n");
    	CBPED::calculate_maf(genVec,pedVec);
    	CBSNP::write_SNP_frequency_file(genVec, _output_filename);
    	pBpar->write_freq =false;
    }
    // handle with the --missing command
    if(pBpar->write_indv_callrate ||pBpar->write_snp_callrate)
    {
    	CBPED::calculate_callrate(pedVec,genVec);
    	if(pBpar->write_indv_callrate)
    	{
        	pBpar->write_indv_callrate=false;
        	CBPED::schreibe_ind_callrate(genVec,pedVec,_output_filename);

    	}
    	// snp call rate
    	if(pBpar->write_snp_callrate)
    	{
    		pBpar->write_snp_callrate=false;
        	CBSNP::schreibe_snp_callrate(genVec,_output_filename);
    	}


    }
    // hadle with the --hardy command
    if(pBpar->write_hwe)
    {
    	pBpar->write_hwe =false;
    	CBPED::calculate_snp_hwe(genVec,pedVec);
    	CBSNP::schreibe_snp_hwe(genVec,_output_filename);
    }

    // filtering will be just befor the oformat change.
    if(pBpar->filter_snp|| pBpar->filter_indiv)
    {
    	/*
    	 *  Filtering SNPs and individuals is a  sensitive task. If you filter first SNPs and then calculate call rate, then it may be that you have bad sample callrate
    	 *  On the other hand If you do not filter or filter samples first, then it affects the quality of SNPs. So there is not a staright forward solution for it.
    	 *  What I have done here:
    	 *  	1. I have calculated all quality measures call rate , hwe and maf by considering all possible snps and individuals included in the data. That means there is no bias
    	 *  	In calculating quality measures in this sense.
    	 *  	2. I have then filtered SNPs  and indviduals then at the end.
    	 */
    	CBPED::calculate_callrate(pedVec,genVec);
        CBPED::calculate_snp_hwe(genVec,pedVec);
        CBPED::calculate_maf(genVec,pedVec);
    	if(pBpar->filter_snp)
    	{
    		printLIN("*->Filtering SNPs according as the cutoff  given with option \"--filter-snp\":\n");
    		handle_with_filter_snp_option(pBpar);
    		pBpar->filter_snp=false;
    	}
    	    // filtering will be just befor the oformat change.
    	if(pBpar->filter_indiv)
    	{
    		printLIN("*->Filtering individuals according as the cutoff  given with option \"--filter-indiv\":\n");
    		handle_with_filter_indiv_option(pBpar);
    		pBpar->filter_indiv=false;
    	}

    }


    //#######################
    if(pBpar->change_format)
    {
    	 //vector<string> outputFormats;
    	 CBSNP::_out_format=pBpar->oFormat; //this temporary
    	 const string _oFormat =CBSNP::_out_format;
		// change format of the value;
    	if(_oFormat=="plink")
		{
			//if(pBpar->code_readType!="plink") //
			{
				CBPED::writeOutPlinkFiles(pedVec,genVec,_output_filename);
				CBSNP::write_snpinfoFile(genVec,_output_filename);
				CBPED::write_pedinfoFiles(pedVec,_output_filename);
			}
			/*else
			{
				string msg= "Changing from \"plink type of\" format to \""+_oFormat+"\" format is not possible.\n";
				error(msg);
			}	
			*/		
		}
    	else if(_oFormat=="plink-bed")
    	{
    		//
    		CPPED::write_BITFILE(genVec,pedVec,_output_filename);
    		CBSNP::write_snpinfoFile(genVec,_output_filename);
    		CBPED::write_pedinfoFiles(pedVec,_output_filename);

    	}

		// --oformat plink-doase should write out three files 
			//1. fam file 
			//  map file 
			// dosage file 
		else if (_oFormat=="plink-dosage")
		{
			//cout<< _oFormat <<endl ;//debug
			CBPED::write_plink_dosagesFiles(pedVec, genVec, _output_filename);
			CBSNP::write_snpinfoFile(genVec,_output_filename);
			CBPED::write_pedinfoFiles(pedVec,_output_filename);
		}

		else if (_oFormat=="plink-recodeA" || _oFormat=="plink-recodeAD" )
		{
			//	cout<< _oFormat <<endl ;//debug
			//CBSNP::convert_genotype_into_0123(genVec);
			CBPED::calculate_maf(genVec,pedVec);
			CPPED::write_plink_raw_files(pedVec, genVec, _output_filename);
			CBSNP::write_snpinfoFile(genVec,_output_filename);
			CBPED::write_pedinfoFiles(pedVec,_output_filename);
			CBSNP::write_SNP_frequency_file(genVec, _output_filename);
		}
		else if(_oFormat=="recodeA-dose") //|| _oFormat=="recodeAD-dose")
		{
			// recodeA-dosage
			// recodeAd-dosage
			CBPED::calculate_maf(genVec,pedVec);
			CPPED::write_plink_raw_dose_files(pedVec, genVec, _output_filename,pBpar->ref_allele);
			CBSNP::write_snpinfoFile(genVec,_output_filename);
			CBPED::write_pedinfoFiles(pedVec,_output_filename);
			CBSNP::write_SNP_frequency_file(genVec, _output_filename);
		}
    	// Hier kommt mach
		else if(_oFormat=="mach")
		{
				CBPED::write_mach_files(pedVec, genVec, pBpar-> output_fileName);
				//CBSNP::write_snplist_file(genVec,_output_filename);
				CBSNP::write_snpinfoFile(genVec,_output_filename);
				CBPED::write_pedinfoFiles(pedVec,_output_filename);
		}
		else if(_oFormat=="minimac")
		{
			 //cout <<" Hi from minimac !\n";
			 CBPED::write_mach_files(pedVec, genVec, pBpar-> output_fileName);
			 CBSNP::write_snplist_file(genVec,_output_filename);
			 CBSNP::write_snpinfoFile(genVec,_output_filename);
			 CBPED::write_pedinfoFiles(pedVec,_output_filename);

		}
		else if(_oFormat=="impute")
    	{
    		//if(pBpar->code_readType!="impute") //
			{
				CBPED::write_impute_gensFile(genVec,pedVec,_output_filename,pBpar->read_extraMap);
    							//CBPED::write_impute_gensFile(genVec,pedVec,_output_filename);
				CBSNP::write_snpinfoFile(genVec,_output_filename);
				CBPED::write_pedinfoFiles(pedVec,_output_filename);
				/*To write snptest *.sample file
								 * You need to have personal call rate so calculate it first and then write *.sample file
								 */
				CBPED::calculate_callrate(pedVec,genVec);
				CBPED::write_snptest_sample_file(pedVec, _output_filename);
			}
			/*else
			{
				string msg= "Changing from \"impute type of\" format to \""+_oFormat+"\" format is not possible.\n";
				error(msg);
			}
			*/
    	}
		// here is format changing for SNPTEST
		else if(_oFormat=="snptest")
		{
			//if(pBpar->code_readType!="snptest") //
			{
				//cout <<"test";
				CBPED::write_impute_gensFile(genVec,pedVec, _output_filename,pBpar->read_extraMap);
				CBSNP::write_snpinfoFile(genVec,_output_filename);
				CBPED::write_pedinfoFiles(pedVec, _output_filename);
				/*To write snptest *.sample file
				 * You need to have personal call rate so calculate it first and then write *.sample file
				 */
				CBPED::calculate_callrate(pedVec,genVec);
				CBPED::write_snptest_sample_file(pedVec, _output_filename);
			}
			/*else
			{
				string msg= "Changing from \"snptest type of\" format to \""+_oFormat+"\" format is not possible.\n";
				error(msg);
			}
			*/				
		}
		else if(_oFormat=="beagle")
    	{
			//if(pBpar->code_readType!="beagle") //
			{
				CBPED::write_beagleGenotypeFile(pedVec, genVec, _output_filename);
				if(given_pgeno)
					CBEPED::write_beagle_gprobsFile(pedVec, genVec, _output_filename);
    			CBPED::write_pedinfoFiles(pedVec, _output_filename);
    			CBSNP::write_snpinfoFile( genVec,_output_filename);
			}
			/*else
			{
				string msg= "Changing from \"beagle type of\" format to \""+_oFormat+"\" format is not possible.\n";
				error(msg);
			}
			*/			
    	}
    	else if(_oFormat=="bimbam")
    	{

    		CBPED::write_bimbam_geno_file(pedVec, genVec, _output_filename);
    		CBSNP::write_bimbam_snp_pos_file(genVec, _output_filename);
    		CBPED::write_pedinfoFiles(pedVec, _output_filename);
    		CBSNP::write_snpinfoFile( genVec,_output_filename);
    	}
    	else if(_oFormat=="r"||_oFormat=="R")
    	{
    		printLIN("*->Writing files in R format: \n");
    		//CBSNP::convert_genotype_into_0123(genVec);
    		CBPED::calculate_maf(genVec,pedVec);
    		//cout <<"test"<<endl;
    		CRPED::write_affection_status_file(pedVec,_output_filename);
    		CRSNP::write_ref_and_alternative_allele_info_file(genVec, _output_filename, pBpar->ref_allele);

    		if(pBpar->transpose)
    			CRPED::write_transposed_genotype_file(pedVec,genVec,_output_filename,pBpar->ref_allele);
    		else
    			CRPED::write_genotype_file(pedVec,genVec,_output_filename,pBpar->ref_allele);
    		CBSNP::write_SNP_frequency_file(genVec, _output_filename);
       		CBSNP::write_snpinfoFile(genVec,_output_filename);
       		CBPED::write_pedinfoFiles(pedVec, _output_filename);

    	}
    	else if(_oFormat=="r-dose"||_oFormat=="R-dose")
    	{
    		printLIN("*->Writing files in "+_oFormat+" format: \n");
    		//CBSNP::convert_genotype_into_0123(genVec);
    		CBPED::calculate_maf(genVec,pedVec);
    		CRPED::write_affection_status_file(pedVec,_output_filename);
    		CRSNP::write_ref_and_alternative_allele_info_file(genVec, _output_filename, pBpar->ref_allele);

    		CRPED::write_allele_dose_file(pedVec,genVec,_output_filename,pBpar->ref_allele);
    		CBSNP::write_SNP_frequency_file(genVec, _output_filename);
    		CBSNP::write_snpinfoFile(genVec,_output_filename);
    		CBPED::write_pedinfoFiles(pedVec, _output_filename);

    	}

    	else if(_oFormat=="haploview"){
    		CHAPPED::write_ped_file(pedVec,genVec,_output_filename);
    		CHAPSNP::write_info_file(genVec,_output_filename);
    		CBPED::write_pedinfoFiles(pedVec, _output_filename);
    		CBSNP::write_snpinfoFile( genVec,_output_filename);
    	}
    	else if(_oFormat=="eigensoft"||_oFormat=="smartpca")
    	{
    		if(pBpar->read_group_label)
    			CSMARTPCAPED::read_smartpca_groupLabel_file(pedVec,pBpar->group_label_fileName);
    		printLIN("*->Writing files in "+_oFormat+" format: \n");
    		CSMARTPCAPED::write_smartpca_ped_file(pedVec,genVec,_output_filename);
    		CSMARTPCASNP::write_smartpca_pedsnp_file(genVec,_output_filename);
    		CSMARTPCAPED::write_smartpca_pedind_file(pedVec,_output_filename);
    		CSMARTPCASNP::write_smartpca_command_file(_output_filename);
    		CBPED::write_pedinfoFiles(pedVec, _output_filename);
    		CBSNP::write_snpinfoFile( genVec,_output_filename);
    	}
    	else if(_oFormat=="genabel")
    	{
    		printLIN("*->Writing files in "+_oFormat+" format: \n");
    		//cout << genVec.size()<<endl;

    		//for(unsigned int j=0; j<genVec.size();++j)
    		//	cout << genVec[j]->geno1.size() <<",  "<< genVec[j]->geno1.size()<<endl ;
    		CBPED::calculate_maf(genVec,pedVec);
    		CGPED::write_genable_raw_file(pedVec,genVec,_output_filename);
    		CGPED::write_genable_pheno_file(pedVec, _output_filename);
    		CGPED::write_geable_commands(_output_filename);
    		CBPED::write_pedinfoFiles(pedVec, _output_filename);
    		CBSNP::write_snpinfoFile( genVec,_output_filename);
    	}
    	else if(_oFormat=="probabel")
    	    	{
    	    		printLIN("*->Writing files in "+_oFormat+" format: \n");
    	    		//cout << genVec.size()<<endl;

    	    		//for(unsigned int j=0; j<genVec.size();++j)
    	    		//	cout << genVec[j]->geno1.size() <<",  "<< genVec[j]->geno1.size()<<endl ;
    	    		CBPED::calculate_maf(genVec,pedVec);
    	    		CGSNP::write_probabel_mlinfo_file(genVec,_output_filename);
    	    		CGPED::write_probabel_mlprob_file(pedVec, genVec,_output_filename);
    	    		//CGPED::write_genable_pheno_file(pedVec, _output_filename);
    	    		//CGPED::write_geable_commands(_output_filename);
    	    		CBPED::write_pedinfoFiles(pedVec, _output_filename);
    	    		CBSNP::write_snpinfoFile( genVec,_output_filename);
    	    }

    	else if((_oFormat=="phase")||(_oFormat=="fastphase"))
    	{
    		printLIN("*->Writing files in "+_oFormat+" format: \n");
    		CGENERAL::write_phase_formatted_files(pedVec,genVec,pBpar->output_fileName);
    		CBPED::write_pedinfoFiles(pedVec, _output_filename);
    		CBSNP::write_snpinfoFile( genVec,_output_filename);
    	}
    	else if(_oFormat=="mach-ref" && (pBpar->code_readType=="impute_hapmap" || pBpar->code_readType=="shapeit_haplotype"))
    	{
    		printLIN("*->Writing files in "+_oFormat+" format: \n");
    		CMPED::write_mach_ref_files(pedVec, genVec, pBpar-> output_fileName);

    	}
    	else if(_oFormat=="vcf")
    	{

    		//written should be the following
    		// before calculate Allele frequency
    		//
    		 CBPED::calculate_maf(genVec,pedVec);
    		CVCF.write_vcf_format(pedVec,genVec,_output_filename,  VCF_HEADER::meta_id,pBpar->uncompressed);
    		CBSNP::write_snpinfoFile(genVec,_output_filename);
    		CBPED::write_pedinfoFiles(pedVec,_output_filename);

    	}

    	else
    	{
			string msg= "Changing into \" "+_oFormat+"\" format is not possible.\n";
    		error(msg);

    	}
	}

   // calculate fst
    // handle with fst command
  if(pBpar->calculate_fst){
	 // handel_fst_calculations(pedVec, genVec,  pBpar->group_label_fileName,pBpar->output_fileName);
	  calculate_nei_wc_hudson_fst(pedVec, genVec,  pBpar->group_label_fileName,pBpar->output_fileName); //const string& ref_allele last argument put it later 
	  printLIN("\n *->NOTE:\t Negative values of Fst should be considered to be zero even though I let  them to be as they are obtained from the calculation.\n");	
	  //string ss ="ss";
	//  handel_fst(ss); //(pedVec, genVec,pBpar->group_label_fileName, pBpar->output_fileName);
	}
    	//const vector<CBPED*> &pedVec, const vector<CBSNP*>genVec
}
//plink type
void CGENERAL::modify_snp_ped_info(const Bpar* const pBpar)
{
	if(pBpar->paste_iid)
	{
		vector<string> temp_ped_info;
		string temp_iid="";
		unsigned int i,j=0;
		const unsigned int temp_sz= pBpar->pos_pedinfo_2bpasted.size();
		//cout <<"size: "<< temp_sz << " \n";
		//temp_ped_info.push_back(p)
		CBPED* ppedInfo =pedVec[0];
		//cout << pBpar->pos_pedinfo_2bpasted.size()<<"\n ";
		//for(unsigned i=0;i<pBpar->pos_pedinfo_2bpasted.size();++i)
		//		cout<< pBpar->pos_pedinfo_2bpasted[i] <<" ";
		//		cout<< "\n";
		//const char* temp_pedinfo[7]={"fid","FID","iid","IID","groupLabel","patid","matid"};
		for(i=0; i<pedVec.size();++i)
		{
			ppedInfo =pedVec[i];
			//
			temp_ped_info.push_back(ppedInfo->famId);
			temp_ped_info.push_back(ppedInfo->indId);
			temp_ped_info.push_back(ppedInfo->groupLabel);
			temp_ped_info.push_back(ppedInfo->patId);
			temp_ped_info.push_back(ppedInfo->matId);
			temp_ped_info.push_back(ppedInfo->groupLabel);
			while(j<temp_sz)
			{
				if(j<(temp_sz-1))
				temp_iid=temp_iid+temp_ped_info[pBpar->pos_pedinfo_2bpasted[j]]+pBpar->paste_seperator;
				else
					temp_iid=temp_iid+temp_ped_info[j];
				j++;
			}
			ppedInfo->indId=temp_iid;
			//cout << temp_iid << " \n"; //debug
			temp_iid="";
			j=0;
			temp_ped_info.clear();

			//ppedInfo->iid =pBpar->
		}
	}
}
void CGENERAL::plink_type_commands(Bpar* const pBpar)
{
	const string  read_code_type 	=pBpar->code_readType;
	CPPED::_read_recodeA_type		=pBpar->read_plink_rawA;
	CPPED::_read_recodeAD_type		=pBpar->read_plink_rawAD;
	CBPED::_given_cov_names			=pBpar->given_covar_names;
	CBPED::_given_cov_types			=pBpar->given_covar_types;
	CBPED::_cov_names_frm_argv	    =pBpar->vec_covar_names;
	CPSNP::_given_three_cols		=pBpar->three_columns;
	if(pBpar->code_readType=="plink")
	{
		CPSNP::readMapFile(pBpar->map_file,  genVec);
		CPPED::pedFileleser(pBpar->ped_file,genVec,pedVec);

	}
	else if(pBpar->code_readType=="bplink"){
		//read bim file
		CPPED::lese_bim_file(genVec, pBpar->bim_file);
		//read fam file
		CPPED::lese_plink_fam_file(pedVec,pBpar->fam_file);
		//read bed file
		CPPED::lese_bed_file(pedVec,genVec, pBpar->bed_file);
	}
	else if(pBpar->code_readType=="plink-dosage")
	{
		// cout <<"Hi hier ist plink-dosage type \n";//debug
		 if(pBpar->read_map)
			 CPSNP::readMapFile(pBpar->map_file,genVec);
		 CPPED::lese_plink_fam_file(pedVec,pBpar->plink_fam_file);
		 CPPED::lese_plink_dosage_file(pedVec,genVec,CBPED::save_line,pBpar->threshold,pBpar->plink_dosage_file);
		 //CBSNP::given_pgeno=true;
		  given_pgeno=true;
	}
	else if(pBpar->code_readType=="plink-rawA"||pBpar->code_readType=="plink-rawAD")
	{
		 //if(pBpar->read_map)
		CPSNP::readMapFile(pBpar->map_file,genVec);
		// for extra snpinfo file, it is compulsuary now to update allele infos.
		// So no conditional
		CBSNP::extraMapFileLeser(pBpar->extraMap_file, genVec);
		if(pBpar->read_plink_rawA)
		 CPPED::lese_plink_raw_file(pedVec,genVec,
				 CBPED::save_line,pBpar->plink_rawA_file);
		 else
			 CPPED::lese_plink_raw_file(pedVec,genVec,CBPED::save_line,pBpar->plink_rawAD_file);
		 CBSNP::_read_extramap =false; // not to re-read the extra file.

	}


	//PED.displayIndividualSummary(pedVec);

return;}

//read data
void CGENERAL::read_data(Bpar* const pBpar)
{
	//CGENERAL GENERAL();
	//cout << "pBpar->code_readType: " << pBpar->code_readType<<endl;
	const string _name_data_type =pBpar->code_readType;
	if(_name_data_type!="")
	{
		// cout << "pBpar->is_code_impute: "<< pBpar->is_code_impute <<endl;
		if(pBpar->is_code_plink)
			plink_type_commands(pBpar);
		else if(pBpar->is_code_mach)
			mach_type_commands(pBpar);
		else if(pBpar->is_code_minimac)
			minimac_type_commands(pBpar);
		else if(pBpar->is_code_impute)
			impute_type_commands(pBpar);
		//----------------------------------------------------//
		//starting for SNPTEST
		// even though I have writen SNPTEST as a seperate group in the front page
		// SNPTEST belongs to IMPUTE category while creating the codes for SNPTEST
		// So same function should be started again.
		//----------------------------------------------------//
		else if(pBpar->is_code_snptest)
			snptest_type_commands(pBpar);
		else if(pBpar->is_code_beagle)
			beagle_type_commands(pBpar);
		else if(pBpar->is_code_bimbam)
			bimbam_type_commands(pBpar);
		else if(pBpar->is_code_rformat)
			rformat_type_commands( pBpar);
		else if(pBpar->is_code_shapeit)
		{
			shapeit_type_commands(pBpar);
		}
		else
		{
			string msg=  "This type of  format files can't be accepted. Please choose another command.\n";
			error(msg);
		}
		//until now only the files has been read according as its form. From now on general things are done.
	/**	*General commands :

	    * reading extra snpinfo file or ped info file
	    * displaying summary
	    * writing snp info and ped info
	    * converting format of the data
	    *
	 */
		CBSNP::_read_extramap	=pBpar->read_extraMap;
		CBSNP::_read_extraped	=pBpar->read_extraPed;
		if(CBSNP::_read_extramap)
			CBSNP::extraMapFileLeser(pBpar->extraMap_file, genVec);
		if(CBSNP::_read_extraped)
			CBPED::pedinfoFileLeser(pedVec,pBpar->extraPed_file);
						//CBPED::displayIndividualSummary(pedVec);
						//here starts the data management of mach format
		modify_snp_ped_info(pBpar);
		if(pBpar->given_force)
			handle_with_force_option(pBpar->force_first_part,pBpar->force_2nd_part,pBpar);
		if(pBpar->read_covariate_file)
			CBPED::lese_plink_covariate_file(pedVec,pBpar->covariate_file);
		//handle with exclude snps and remove indivis commands
		if(pBpar->exclude_snps)
			CBSNP::lese_exclude_snps_file(pBpar->exclude_snplist_file,genVec);
		if(pBpar->remove_indivs)
			CBPED::lese_remove_indivs_file(pBpar->remove_indivlist_file,pedVec);
		assign_snp_indiv_summary(); // pedVec
		display_snp_summary();
		display_indiv_summary();
	}//end of if  for code_readType.
	 //-----------------------------------------------------------------------------//

}

void CGENERAL::rformat_type_commands(Bpar* const pBpar)
{


	CRSNP::lese_genotype_file(pedVec, genVec, pBpar->rgeno_file);
	//	--snpinfo fcgene_out_alleleInfo.txt
}

//mach type 
void CGENERAL::minimac_type_commands( Bpar* const pBpar)
{
	const string read_code_type=pBpar->code_readType;
	//cout <<"hi from minimac_type_commands\n";
	if(read_code_type=="minimac.prob"){
	//	cout <<" i am a step forward \n";
		vector<string> saveLine;
		CMINIMACSNP::lese_minimac_info_file(pBpar->minimac_infoFile,genVec,pBpar->rsq_value,pBpar->maf_value);
		CMINIMACPED::lese_minimac_probFile(pedVec,genVec,saveLine,pBpar->minimac_probFile, pBpar->threshold);
		//CBSNP::given_pgeno=true;
		given_pgeno=true;
	}
}
 void CGENERAL::mach_type_commands( Bpar* const pBpar)
    {
    	const string read_code_type=pBpar->code_readType;
	 	 //CMPED MPED;
    	//CMSNP MSNP;
		 // mach hapmap file
		 if(read_code_type=="mach_hapmap")
		 {
			 CMSNP::read_mach_ref_snpsFile(genVec,pBpar->mach_hap_snp_file);
			 CMPED::lese_mach_hapmap_file(pedVec,genVec,CMPED::save_line,pBpar->mach_hapmap_file);
			// CBSNP::extraMapFileLeser(pBpar->mach_hap_snp_file, genVec);
		}	
		 // mach mlprob 
    	if(read_code_type=="mach.mlprob")
    	{
    		// MSNP.mapFileLeser(pBpar->mlinfo_file ,genVec);
    		 CMSNP::lese_mach_info_mlinfo_file(pBpar->mlinfo_file, genVec,pBpar->rsq_value,pBpar->maf_value);
    		 CMPED::lese_mach_mlProbFile(pedVec,genVec,CMPED::save_line,pBpar->mlprob_file,pBpar->threshold);
    		// CBSNP::given_pgeno=true;
    		 given_pgeno=true;
    	}

    	else if(read_code_type=="mach")
    	{
    		CMSNP::_is_code_mach=true;
    		CMSNP::mapFileLeser(pBpar->dat_file ,genVec);
			CMPED::lese_pedfile(pedVec,genVec,CMPED::save_line,pBpar->ped_file);
			//for(unsigned int j=0; j<genVec.size();++j)
			// cout <<"genVec[j]->bp: " << genVec[j]->bp << " \n";
			
    	}
		else if(read_code_type=="mach.mlgeno")
    	{
    		CMSNP::lese_mach_info_mlinfo_file(pBpar->mlinfo_file, genVec,pBpar->rsq_value,pBpar->maf_value);
    		//	CMPED::lese_mach_geno_mlGenoFile(pedVec,genVec,CMPED::save_line,pBpar->mlgeno_file);
    		CMPED::lese_mach_geno_mlGenoFile(pedVec,genVec,CMPED::save_line,pBpar->mlgeno_file);
    	}
    	else if(read_code_type=="mach.geno")
    	{
    		 //CMSNP::lese_mach_info_mlinfo_file(pBpar->mach_info_file,genVec);
    		 CMSNP::lese_mach_info_mlinfo_file(pBpar->mach_info_file,genVec,pBpar->rsq_value,pBpar->maf_value);
    		 CMPED::lese_mach_geno_mlGenoFile(pedVec,genVec,CMPED::save_line,pBpar->mach_geno_file);
    		//CMPED::lese_mach_geno_mlGenoFile(pedVec,genVec,CMPED::save_line,pBpar->mach_geno_file);
    	}
    	
		//-------------------------------------------------------------------------------------------------//

		
  }
//impute type 
void CGENERAL::impute_type_commands( Bpar* const pBpar)
{
	const string read_code_type=pBpar->code_readType;
	//CIPED IPED;
		//CISNP ISNP;

		if(read_code_type=="impute_hapmap")
		{
			CISNP::read_impute_ref_legendfile(genVec,pBpar->impute_hap_snp_file);
			CIPED::lese_impute_hapmap_file(pedVec,genVec,CIPED::save_line,pBpar->impute_hapmap_file);
			//CBSNP::extraMapFileLeser(pBpar->impute_hap_snp_file, genVec);
		}	 
		//cout << "test1"<<endl;	 
		else if(read_code_type=="impute")
		{
			//cout <<"genvec size(): "<< genVec.size() <<endl;
			//CISNP::gensFileLeser(genVec,pBpar->gens_file);
			CISNP::lese_imputegens_file(genVec,CIPED::save_line, pBpar->gens_file);
			//cout<<pedVec.size()<<" "<<pedVec.size()<<endl;
			CIPED::pedFileLeser(pedVec,genVec);
			//cout <<"genvec size(): "<< genVec.size() <<endl;
		}
		//if(CBSNP::_read_extramap)
		//	ISNP.extraMapFileLeser(pBpar->extraMap_file, genVec);
		//if(pBpar->read_extraPed)
		//	CBPED::pedinfoFileLeser(pedVec,pBpar->extraPed_file);
		if(pBpar->read_impute_info_file)
			CISNP::lese_imputed_info_score_file(genVec,CIPED::save_line, pBpar->info_thresh, pBpar->maf_value,pBpar->impute_info_file);
		//ISNP.display_snp_summary(genVec);
		//IPED.display_ped_summary();
		CIPED::imputegen_to_pedgen_converter(genVec,pBpar->threshold);
		//CBSNP::given_pgeno=true; // this shows that probababilites are given.
		given_pgeno=true;
return;}

void CGENERAL::snptest_type_commands(Bpar* const pBpar)
{
	const string read_code_type=pBpar->code_readType;
	cout << "This command is under construction. Please contact the software developer if you want to use it. \n";
	//cout <<"SNPTEST: "<< read_type_snptest << endl ;  // debug
	//CIPED::impute_type_commands(pBpar->code_readType);

}
//bimbam
void CGENERAL::bimbam_type_commands( Bpar* const pBpar)
{
	const string read_code_type=pBpar->code_readType;
	//CBIMPED BIMPED;
	//CBIMSNP BIMSNP;
	if(pBpar->code_readType=="bimbam")
	{
		CBIMPED:: read_bimbam_data(pedVec, genVec, CBIMPED:: save_line,pBpar->bimbam_file);
		if(pBpar->read_bimbam_pos)
			CBIMSNP::read_bimbam_pos_data(genVec,  CBIMPED:: save_line,pBpar->bimbam_pos_file);

	}
	else  if(pBpar->code_readType=="bimbam.gprob")
	{
		 //cout <<"pBpar->read_bimbam_gprobs"<<endl;//debug
		if(pBpar->read_bimbam_gprobs)
		{
			CBIMPED:: read_bimbam_gprobs_data(pedVec, genVec, CBIMPED:: save_line,pBpar->bimbam_gprobs_file, pBpar->threshold);
			CBIMSNP::read_bimbam_snpinfo_data(genVec,  CBIMPED:: save_line,pBpar->bimbam_pos_file,pBpar->maf_value);
			//CBSNP::given_pgeno=true;
			given_pgeno=true;

		}
	}
			else if(pBpar->read_bimbam_bestguess)
			 {
				CBIMPED:: read_bimbam_bestguess_data(pedVec, genVec,  CBIMPED:: save_line, pBpar->bimbam_bestguess_file);
				CBIMSNP::read_bimbam_snpinfo_data(genVec,  CBIMPED:: save_line,pBpar->bimbam_pos_file,pBpar->maf_value);
			}
			 else
			 {
				//genotype probability distribution file: "+pBpar->bimbam_gprobs_file;
				// error(probleme);
			 }
			//CBIMSNP::display_snp_summary(genVec);
			//CBIMPED:: display_ped_summary();
			//	CBIMSNP::bimbam_type_commands(pBpar->code_readType);


}
//beagle
void CGENERAL::beagle_type_commands( Bpar* const pBpar)
{
	const string read_code_type=pBpar->code_readType;
	//CBEPED BEPED;
	//CBESNP BESNP;
	CBESNP::_given_gprobs_header =pBpar->given_gprobs_header;
	if(pBpar->code_readType=="beagle" && pBpar->read_bgl)
		CBEPED::read_bgl_data(pedVec, genVec,  CBEPED::save_line,pBpar->bgl_file);
	else if(pBpar->code_readType=="beagle_gprobs" && pBpar->read_bgl_gprobs)
	{
		CBEPED::read_bgl_gprobs_data(pedVec, genVec,  CBEPED::save_line,pBpar->bgl_gprobs_file,pBpar->threshold);
		//CBSNP::given_pgeno=true;
		given_pgeno=true;
	}
	else
	{	string probleme="Problem with reading Beagle type file";
				error(probleme);
	}	
	if(pBpar->read_bgl_rsq)
		CBESNP::lese_beagle_rsq_score_file(genVec,CBEPED::save_line, pBpar->bgl_rsq_thresh, pBpar->bgl_rsq_file);
	// writing summary
	

}

//write a function which update phenotype 
void CGENERAL::handle_with_force_option(const vector<string>& firstPart, const vector<string>&secondPart, Bpar* const pBpar)
{
		vector<string> _firstPart=firstPart;
		vector<string> _2ndPart=secondPart;
		//---------------------------------------------------//
		//handelling with phoenotype
		//---------------------------------------------------//
		vector<string>::iterator it=find(_firstPart.begin(),_firstPart.end(), "pheno");
		//vector<string>::iterator vecEnd=_firstPart.begin();
		//
		//vector<string>::iterator its=_2ndPart.begin();
		//vector<string>::iterator s_vend=_2ndPart.begin();
		// cout<< _firstPart.size() << " " ; //debug 
		if(it!=_firstPart.end())
		{
			const int _cint=(int)(it-_firstPart.begin());
			const string nextPart=_2ndPart[_cint];
			const bool update_pheno=true;	
			
				if(nextPart=="0"|| nextPart=="1"||nextPart=="2" ||nextPart=="-9" ||nextPart=="aff" ||nextPart=="unaff")
				{
					
					 CBPED* pCBPED =pedVec[0];
					for(unsigned int i=0;i<pedVec.size();++i)
					{
						pCBPED =pedVec[i];
						CBPED::lese_phenotype_info(pCBPED, nextPart,update_pheno);
						
					}
				
				}
				else
				{
					string pr_msg ="This command is still in construction.\n";
					error(pr_msg);
				}
				_firstPart.erase(it);
				_2ndPart.erase(_2ndPart.begin()+(_cint));
				
				//vecEnd=remove(_firstPart.begin(),_firstPart.end(),*it);
				//s_vend =remove(_2ndPart.begin(),_2ndPart.end(),nextPart);
			
			
		}
		
		//cout << _firstPart.size()<< endl; //debug 
		//cout << _2ndPart.size()<< endl; //debug
		//---------------------------------------------------//
		//handelling with sex
		//---------------------------------------------------//
		it=find(_firstPart.begin(),_firstPart.end(), "sex");
		if(it!=_firstPart.end())
		{
			const int _cint=(int)(it-_firstPart.begin());
			const string nextPart=_2ndPart[_cint];
			const bool update_sex=true;
			if(nextPart=="0"|| nextPart=="1"||nextPart=="2" ||nextPart=="M"||nextPart=="F"||nextPart=="m"||nextPart=="f")
				{

					 CBPED* pCBPED =pedVec[0];
					for(unsigned int i=0;i<pedVec.size();++i)
					{
						pCBPED =pedVec[i];
						CBPED::lese_sex_info( pCBPED, nextPart,update_sex);
						//CBPED::lese_phenotype_info(pCBPED, nextPart,update_pheno);

					}

				}
				else
				{
					string pr_msg ="This command is still in construction.\n";
					error(pr_msg);
				}
			_firstPart.erase(it);
			_2ndPart.erase(_2ndPart.begin()+(_cint));

		}
		 //cout << _firstPart.size() <<endl;

		//cout <<"_rem:"<< _rem <<endl;
		//---------------------------------------------------//
				//handelling with reference allele
			//---------------------------------------------------//
		it=find(_firstPart.begin(),_firstPart.end(), "ref-allele");
		if(it!=_firstPart.end())
		{
			const int _cint=(int)(it-_firstPart.begin());
			const string nextPart=_2ndPart[_cint];
			//const bool update_ref_allele=true;
			if(nextPart=="allele1")
			{
				pBpar->ref_allele ="allele1";
			}
			else if(nextPart=="allele2")
			{
				pBpar->ref_allele ="allele2";
			}
			else if((nextPart=="major-allele")||(nextPart=="major_allele"))
			{
				pBpar->ref_allele ="major-allele";
			}
			else if((nextPart=="minor-allele")||(nextPart=="minor_allele"))
			{
				pBpar->ref_allele ="minor-allele";
			}
			else{
				printLIN( "The following option given with \"--force\" command, can not be accepted:\n ");
								cerr<< _firstPart[_cint]<< "="<<_2ndPart[_cint] <<endl;
								LIN << _firstPart[_cint]<< "="<<_2ndPart[_cint] <<endl;
							const string _msg = "Please correct comand option given after \"--force\".\n";
							error(_msg);
			}
			_firstPart.erase(it);
			_2ndPart.erase(_2ndPart.begin()+(_cint));
		}
		it=find(_firstPart.begin(),_firstPart.end(), "nchr");
		if(it!=_firstPart.end())
			{
				const int _cint=(int)(it-_firstPart.begin());
				const string nextPart=_2ndPart[_cint];
				CBSNP* pCBSNP =genVec[0];
				for(unsigned int i=0;i<genVec.size();++i)
				{
					pCBSNP =genVec[i];
					pCBSNP->nchr=nextPart;
					//CBPED::lese_phenotype_info(pCBPED, nextPart,update_pheno);

				}
				_firstPart.erase(it);
				_2ndPart.erase(_2ndPart.begin()+(_cint));
			}

		//---------------------------------------------------//
			//handelling with----
		//---------------------------------------------------//
		const int _rem=_firstPart.size();




		if(_rem>0)
		{
			printLIN( "The following options given with \"--force\" command, can not be accepted:\n ");
			for(unsigned int i=0;i<_firstPart.size();++i)
			{
				cerr<< _firstPart[i]<< "="<<_2ndPart[i] <<endl;
				LIN << _firstPart[i]<< "="<<_2ndPart[i] <<endl;
			}
			const string _msg = "Please correct the \"--force\" command.\n";
			error(_msg);
		
		}
		

}

void CGENERAL::handle_with_filter_snp_option(const Bpar* const pBpar)
{
		const vector<string> firstPart=pBpar->filter_snp_firstPart;
		const vector<string> secondPart= pBpar->filter_snp_secondPart;
		string badSNPs_fileName =pBpar->output_fileName+"_removed_snps.txt";
		ofstream ofs;
		ofs.clear(); ofs.close();
		ofs.open(badSNPs_fileName.c_str(),ios::out);
		ofs<<"NCHR"<<" "<<"SNP" <<" "<<"BP";
		ofs<<" "<<"MIN_ALLELE"<<" "<<"MAJ_ALLELE"<<" " <<"MAF";
		ofs<<" "<< "HWE_EXACT"<< " "<<"CRATE";
		ofs<<"\n"; //end of first line
		vector<string> _firstPart			=firstPart;
		vector<string> _2ndPart				=secondPart;
		static unsigned  long int _nbadSNPS	=0;
		bool _ok_maf 						=true;
		bool _ok_hwe 						=true;
		bool _ok_crate 						=true;
		double maf_cutoff 					=0.0;
		double hwe_cutoff 					=0.0;
		double crate_cutoff			 		=0.0;
		CBSNP* pCBSNP =genVec[0];
		vector<string>::iterator it=find(_firstPart.begin(),_firstPart.end(), "maf");
		//cout <<"test \n";
		//cout <<"test2 \n";
		//cout << _firstPart.size()<< endl; //debug
		//cout << _2ndPart.size()<< endl; //debug
		if(it!=_firstPart.end())
		{
			//cout <<"test1\n ";
			const int _cint			=(int)(it-_firstPart.begin());
			 const string nextPart	=_2ndPart[_cint];
			 maf_cutoff 			= atof(nextPart.c_str());
			 _firstPart.erase(it);
			 _2ndPart.erase(_2ndPart.begin()+(_cint));
		}

		//cout <<"test2 \n";
		//cout << _firstPart.size()<< endl; //debug
		//cout << _2ndPart.size()<< endl; //debug
		it=find(_firstPart.begin(),_firstPart.end(), "hwe");
		if(it!=_firstPart.end())
		{
			const int _cint				=(int)(it-_firstPart.begin());
			const string nextPart		=_2ndPart[_cint];
			hwe_cutoff 					= atof(nextPart.c_str());
			_firstPart.erase(it);
			_2ndPart.erase(_2ndPart.begin()+(_cint));
		}
		//
		it=find(_firstPart.begin(),_firstPart.end(), "crate");
		if(it!=_firstPart.end())
		{
				const int _cint				=(int)(it-_firstPart.begin());
				const string nextPart		=_2ndPart[_cint];
				crate_cutoff 				= atof(nextPart.c_str());
				_firstPart.erase(it);
				_2ndPart.erase(_2ndPart.begin()+(_cint));
			}


		  for(unsigned int j=0;j<genVec.size();++j)
			  {
				  pCBSNP =genVec[j];
				  // change into 0123
				  //CBSNP::convert_snp_genotype_into_0123(pCBSNP);
				  //CBSNP::calculate_maf(pCBSNP);
				 // CBSNP::calculate_snp_callrate(pCBSNP); // this calculation is done seperately.
				  _ok_maf 	=(pCBSNP->maf>=maf_cutoff);
				  _ok_hwe 	=(pCBSNP->pvalue_hwe_exact>=hwe_cutoff);
				  _ok_crate =(pCBSNP->snp_CR >=crate_cutoff);
				  pCBSNP->quality =(pCBSNP->quality &&_ok_maf&&_ok_hwe && _ok_crate);
				  if(!pCBSNP->quality)
				  {
					  ++_nbadSNPS;
					  ofs<<pCBSNP->nchr <<" "<<pCBSNP->rsId<<" "<<pCBSNP->bp;
					  if(pCBSNP->min_allele!="")
						ofs<<" "<<  pCBSNP->min_allele;
					  else
						 ofs<<" "<<"0";
					  if(pCBSNP->maj_allele!="")
						  ofs<<" " << pCBSNP->maj_allele;
					  else
						  ofs<<" "<<"0";
					  ofs<<" "<<pCBSNP->maf<< " "<<pCBSNP->pvalue_hwe_exact;
					  ofs<< " "<<pCBSNP->snp_CR;
					  ofs<<"\n"; //end of line
				  }
			  }
		  ofs.clear();ofs.close();

		const int _rem=_firstPart.size();
		if(_rem>0)
		{
			printLIN( "The follwoing options given with \"--filter-snp\" command, can not be accepted:\n ");
			for(unsigned int i=0;i<_firstPart.size();++i)
			{
				cerr<< _firstPart[i]<< "="<<_2ndPart[i] <<endl;
				LIN << _firstPart[i]<< "="<<_2ndPart[i] <<endl;
			}
			const string _msg = "Please correct the \"--filter-snp\" command.\n";
			error(_msg);

		}

	printLIN("\t*->Total no of disqualified SNPs: "+change_int_into_string(_nbadSNPS)+". \n");
	printLIN("\t*->A file containing all disqualified SNPs has been saved with the name\""+badSNPs_fileName+"\"  \n");
return; }
//
//next
void CGENERAL::handle_with_filter_indiv_option(const Bpar* const  pBpar)
{
	const vector<string> firstPart 	=pBpar->filter_indiv_firstPart;
	const vector<string> secondPart	=pBpar->filter_indiv_secondPart;
	string badIndiv_fileName 		=pBpar->output_fileName+"_removed_indiv.txt";
		ofstream ofs;
		ofs.clear(); ofs.close();
		ofs.open(badIndiv_fileName.c_str(),ios::out);
		ofs<<"FID"<<" "<<"IID" <<" "<<"PAT";
		ofs<<" "<<"MAT" " "<<"CRATE";
		ofs<<"\n"; //end of first line
		vector<string> _firstPart			=firstPart;
		vector<string> _2ndPart				=secondPart;
		static unsigned  long int _nbadIndiv=0;
		bool _ok_crate 						=true;
		double crate_cutoff			 		=0.0;
		CBPED* pCBPED =pedVec[0];
		vector<string>::iterator it=find(_firstPart.begin(),_firstPart.end(), "crate");
		if(it!=_firstPart.end())
		{
			//cout <<"test1\n ";
			const int _cint			=(int)(it-_firstPart.begin());
			 const string nextPart	=_2ndPart[_cint];
			 crate_cutoff 			= atof(nextPart.c_str());
			 _firstPart.erase(it);
			 _2ndPart.erase(_2ndPart.begin()+(_cint));
		}

		//CBPED::calculate_indiv_callrate(pedVec, genVec);
		// seperately calculated.
		for(unsigned int j=0;j<pedVec.size();++j)
			  {
				  pCBPED =pedVec[j];
				  _ok_crate =(pCBPED->ind_CR >=crate_cutoff);
				  pCBPED->quality =(pCBPED->quality &&_ok_crate);
				  if(!pCBPED->quality)
				  {
					  ++_nbadIndiv;
					  ofs<<pCBPED->famId <<" "<<pCBPED->indId;
					  if(pCBPED->patId!="")
						ofs<<" "<<  pCBPED->patId;
					  else
						 ofs<<" "<<"0";
					  if(pCBPED->matId!="")
						  ofs<<" " << pCBPED->matId;
					  else
						  ofs<<" "<<"0";
					  ofs<< " "<<pCBPED->ind_CR;
					  ofs<<"\n"; //end of line
				  }
			  }
		  ofs.clear();ofs.close();

		const int _rem=_firstPart.size();
		if(_rem>0)
		{
			printLIN( "The follwoing options given with \"--filter-indiv\" command, can not be accepted:\n ");
			for(unsigned int i=0;i<_firstPart.size();++i)
			{
				cerr<< _firstPart[i]<< "="<<_2ndPart[i] <<endl;
				LIN << _firstPart[i]<< "="<<_2ndPart[i] <<endl;
			}
			const string _msg = "Please correct the \"--filter-indiv\" command.\n";
			error(_msg);

		}

	printLIN("\t*->Total no of disqualified individuals: "+change_int_into_string(_nbadIndiv)+". \n");
	printLIN("\t*->A file containing all disqualified individuals has been saved with the name\""+badIndiv_fileName+"\"  \n");
return; }
bool CGENERAL::are_they_same_snp(const CBSNP* const _pFirst, const CBSNP* const _pSecond)
{
	bool tfVec=false;
	string _snpid1 		 =_pFirst->snpName;
	string _rsid1		 =_pFirst->rsId;
	int 	_bp1		= _pFirst->bp;
	//string  _A11		=_pFirst->allele1;
	//string _A12			=_pFirst->allele2;
	double cm_pos1		=_pFirst->cm_pos;
	/*//int coding_strand1	=_pFirst->coding_strand;
	//int qualty1			=_pFirst->quality;
	string dose_allle11	=_pFirst->dose_allele1;
	string dose_allle12	=_pFirst->dose_allele2;
	string min_allle1	=_pFirst->min_allele;
	string maj_allle1	=_pFirst->maj_allele;
	*/
	//2nd
	string _snpid2 		 =_pSecond->snpName;
	string _rsid2		 =_pSecond->rsId;
	int 	_bp2		= _pSecond->bp;
	//string  _A21		=_pSecond->allele1;
	//string _A22			=_pSecond->allele2;
	double cm_pos2		=_pSecond->cm_pos;
	/*
	int coding_strand2	=_pSecond->coding_strand;
	int qualty2			=_pSecond->quality;
	string dose_allle21	=_pSecond->dose_allele1;
	string dose_allle22	=_pSecond->dose_allele2;
	string min_allle2	=_pSecond->min_allele;
	string maj_allle2	=_pSecond->maj_allele;
	bool _is_ok =false;
	*/
	//cout << _snpid1 <<", "<< _rsid1 <<endl;
	//cout << _snpid2 <<", "<< _rsid2 <<endl;
	if(_snpid1!="")
	{
		//cout << boolalpha<<(_snpid1==_snpid2) <<endl;
		if(_snpid1==_snpid2)
		{
			//comapre rsid
			//cout << boolalpha << "test "<< ((_rsid1!="") && (_rsid2!="")&& (_snpid1!=_rsid2)&&(_snpid1!=_rsid1));
			if((_rsid1!="") && (_rsid2!="")&& (_snpid1!=_rsid2)&&(_snpid1!=_rsid1))
			{
					 tfVec= (_rsid1==_rsid2);
			}
			else
				tfVec =true;

		}// if snp not found, I hope to find according as rsid ;
		else if((_rsid1==_rsid2)&& _snpid2=="")
					tfVec= true;
	}
		// Here I check for rsid
	else if(_rsid1!="")
	{
		if(_rsid1==_rsid2)
			tfVec =true;
		else tfVec =false;

	}
	//check bp
	if((_bp1!=-1)&& (_bp2!=-1))
	{
		if(tfVec&&(_bp1!=_bp2) )
		{

			printLIN("First Data->: snpid: "+ _snpid1+", rsid: "+_rsid1 +", basepair position: "+change_int_into_string(_bp1)+". \n"
					"Second Data->: snpid: "+ _snpid2+", rsid: "+_rsid2 +", basepair position: "+change_int_into_string(_bp2)+". \n"
					);
			error("Base pair positions are not same at both data sets.\n");
		}
		else
		tfVec= (tfVec && (_bp1 ==_bp2));
	}
		//check are two allele the same
	if((cm_pos1!=-1.0)&& (cm_pos2!=-1.0))
	{
		if(tfVec&&(cm_pos1!=cm_pos2) )
		{

			printLIN("First Data->: snpid: "+ _snpid1+", rsid: "+_rsid1 +", basepair position: "+change_int_into_string(_bp1)+", centi morgen: "+change_double_into_string(cm_pos1)+". \n"
					"Second Data->: snpid: "+ _snpid2+", rsid: "+_rsid2 +", basepair position: "+change_int_into_string(_bp2)+", centi morgen: "+change_double_into_string(cm_pos2)+". \n"
					);
			error("Base pair positions are not same at both data sets.\n");
		}
		else
			tfVec= (tfVec && (cm_pos1 ==cm_pos2));
	}

return tfVec;}
bool CGENERAL::check_make_same_alleleOrder( CBSNP* const pFirst,  CBSNP* const pSecond)
{
	//pFist is SNP from original data and pSecond is SNP from other data
	bool _tfvec =false;
	bool _is_ok =true;
	string _A11 =pFirst->allele1;
	string _A12 =pFirst->allele2;
	string _A21 =pSecond->allele1;
	string _A22 =pSecond->allele2;
	if(_A11==_A21)
	{

		_tfvec=true;
		if((_A12!=""&& _A22!="")&& _A12!=_A22)
		{
			_is_ok =false;
			_tfvec =false;
		}
	}
	else if ((_A11==_A22) )
	{
		if(((_A12!="") && (_A12!=_A21) ))
			_is_ok=false;

	}
	if(!_is_ok)
	{
	error("first SNP "+pFirst->rsId+" has Allles: "+_A11+", "+_A12+"\n"
		 "and second SNP "+pSecond->rsId+" has alleles: "+_A21+", "+_A22+".\n"
		 "These SNPs must have the same alleles.\n");
	}
	//cout << "TESt1\n";
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	if(!_tfvec)
	{

			//--------------------------------------------------------------------------------------------------------------------------------------------------------------//
			bool cond_1=pSecond->allele1!="" && pFirst->allele1!="";

			// && pSecond->allele1!="0";
			//bool cond_2=pSecond->allele2!="" && pSecond->allele2!="0";
			//cout << "tmpb2: "<< tmpb2<<endl; // test
			// cout << cond_1 << endl; // test
			if(cond_1)  //this condition is always true
			{
				//1
				bool tmpb1=(pSecond->allele1 ==pFirst->allele1) && (pSecond->allele2 ==pFirst->allele2);
				//2
				bool tmpb2=(pSecond->allele1 ==pFirst->allele2) && (pSecond->allele2 ==pFirst->allele1);
				bool tmpb4=(pSecond->allele1 ==pFirst->allele2) && ((pSecond->allele2 ==""));//
				bool tmpb5=(pSecond->allele2 ==pFirst->allele1) && (pFirst->allele2 =="");
				//3
				bool tmpb3=(pSecond->allele1 ==pFirst->allele1)&& (pSecond->allele2 =="");//||pSecond->allele2 =="0");
				bool tmpb6=(pSecond->allele1 ==pFirst->allele1)&& (pFirst->allele2 =="");//||pSecond->allele2 =="0");
				//4

				//5
				//bool tmpb5=(pSecond->allele2 ==pFirst->allele1) && (pFirst->allele2 =="");
				// this condition does not come!!
				//bool tmpb6=(pSecond->allele2 ==pFirst->allele2) && ((pSecond->allele1 =="") || (pSecond->allele1 =="0"));
				// this condition does not also come!!
				//bool tmpb7=(pSecond->allele1 ==pFirst->allele2) && (( pFirst->allele2 =="" ) || (pFirst->allele2 =="0" ) );
				// this condition is ok nothing to change
				//bool tmpb8=(pSecond->allele1 ==pFirst->allele2) && ((pFirst->allele1 =="" ));// || (pFirst->allele1 =="0" ) );
				//bool tmpb9=(pSecond->allele2 ==pFirst->allele1) && ((pFirst->allele2 =="" ));// || (pFirst->allele2 =="0" ) );
				//cout << "tmpb2: "<< tmpb2<<endl; // test
				if(tmpb1)
				{
					//do nothing
				}
				else if(tmpb2||tmpb4 ||tmpb5 )
				{
					//first check for tmpb4 i.e. if pSecond->allele2==""
					if(pSecond->allele2=="")
					{
						for(unsigned int i=0; i<pSecond->geno1.size();++i)
						{

							//if geno1=T and geno2=T or geno1 =F and geno2=T  (i.e. case homo2 and hetero
							if (((pSecond->geno1[i])&& (pSecond->geno2[i] ))|| (!(pSecond->geno1[i])&& (pSecond->geno2[i])))
							{
								cout <<"Allele in other data: " << 	pSecond->allele1 << " " <<pSecond->allele2 <<endl;
								cout <<"Allele in first data: " << 	pFirst->allele1 << " " <<pFirst->allele2 <<endl;
								string msg= " Allele2 \""+pSecond->allele2+\
									"\"   is empty but there are Homozygous 2  types and heterozygoous genotypes exisit in SNP  \""+\
								pSecond->rsId +"\" \n";
								error(msg);

							}
						}
					}

					//now check for tmpb5
					if(pFirst->allele2=="")
					{
						for(unsigned int i=0; i<pFirst->geno1.size();++i)
						{

							//if geno1=T and geno2=T or geno1 =F and geno2=T  (i.e. case homo2 and hetero
							if (((pFirst->geno1[i])&& (pFirst->geno2[i] ))|| (!(pFirst->geno1[i])&& (pFirst->geno2[i])))
							{
								cout <<"Allele in other data: " << 	pSecond->allele1 << " " <<pSecond->allele2 <<endl;
								cout <<"Allele in first data: " << 	pFirst->allele1 << " " <<pFirst->allele2 <<endl;
								string msg= " Allele2 \""+pFirst->allele2+\
									"\"   is empty but there are Homozygous 2  types and heterozygoous genotypes exisit in SNP  \""+\
									pFirst->rsId +"\" \n";
								error(msg);

							}
						}
					}


					// chnage whole allele 1 and allele 2 with each other .
					//cout << "change of allele! allele 1 is allele2 and allele 2 will be allele1 \n"; //ONly for test
					for(unsigned int i=0; i<pSecond->geno1.size();++i)
					{
						//homo1->homo2
						if(!(pSecond->geno1[i])&& !(pSecond->geno2[i]))
						{
							pSecond->geno1[i]=true;
							pSecond->geno2[i]=true;
						}
						//homo2 ->homo1
						else if((pSecond->geno1[i])&& (pSecond->geno2[i]))
						{
							pSecond->geno1[i]=false;
							pSecond->geno2[i]=false;
						}
						//geno1 =F geno2 =T => heterozygote . here only the order should be changed.
						else if(!(pSecond->geno1[i])&& (pSecond->geno2[i]))
						{
							if(pSecond->aOrder[i])
								pSecond->aOrder[i]=false;
							if(!(pSecond->aOrder[i]))
								pSecond->aOrder[i]=true;
						}
						// change probability genotypes also
						if(pSecond->pgeno1.size()>0 && pSecond->pgeno3.size()>0)
						{
							//cout << no_of_persons <<endl;
							for(unsigned int j=0;j<pSecond->pgeno1.size(); ++j)
								swap(pSecond->pgeno1[j],pSecond->pgeno3[j]);
									 //swap will exchange prob(AA) and prob(BB).prob(AB) remains the same
						}

					}
					// Now change the alleles also but be careful of empty allele
					string _tmp_str =pSecond->allele1;
					pSecond->allele1 =pFirst->allele1;
					if(pFirst->allele2!="")
						pSecond->allele2 =pFirst->allele2;
					else
					{
						pFirst->allele2=pSecond->allele1;
					}	pSecond->allele2=pSecond->allele1;

					//only for test

			}
			else if(tmpb3|| tmpb6)
			{
				// it assumes that there is no allele 2 in the data so check if there is not the case and
				// if this is the case throw error . if not just assign allele2 =second allele2
				//first check for tmpb4 i.e. if pSecond->allele2==""
				if(pSecond->allele2=="")
				{
					for(unsigned int i=0; i<pSecond->geno1.size();++i)
					{

						//if geno1=T and geno2=T or geno1 =F and geno2=T  (i.e. case homo2 and hetero
						if (((pSecond->geno1[i])&& (pSecond->geno2[i] ))|| (!(pSecond->geno1[i])&& (pSecond->geno2[i])))
						{
							cout <<"Allele in other data: " << 	pSecond->allele1 << " " <<pSecond->allele2 <<endl;
							cout <<"Allele in first data: " << 	pFirst->allele1 << " " <<pFirst->allele2 <<endl;
							string msg= " Allele2 \""+pSecond->allele2+\
									"\"   is empty but there are Homozygous 2  types and heterozygoous genotypes exisit in SNP  \""+\
									pSecond->rsId +"\" \n";
							error(msg);
						}
					}
				}
				//now check for tmpb6
				if(pFirst->allele2=="")
				{
					for(unsigned int i=0; i<pFirst->geno1.size();++i)
					{
						//if geno1=T and geno2=T or geno1 =F and geno2=T  (i.e. case homo2 and hetero
						if (((pFirst->geno1[i])&& (pFirst->geno2[i] ))|| (!(pFirst->geno1[i])&& (pFirst->geno2[i])))
						{
							cout <<"Allele in other data: " << 	pSecond->allele1 << " " <<pSecond->allele2 <<endl;
							cout <<"Allele in first data: " << 	pFirst->allele1 << " " <<pFirst->allele2 <<endl;
							string msg= " Allele2 \""+pFirst->allele2+\
									"\"   is empty but there are Homozygous 2  types and heterozygoous genotypes exisit in SNP  \""+\
									pFirst->rsId +"\" \n";
							error(msg);
						}
					}
				}
				if(pFirst->allele2=="")
					pFirst->allele2 =pSecond->allele2;
				else if(pSecond->allele2=="")
					pSecond->allele2 =pFirst->allele2;
			}
			else
			{

				cout <<"Allele in second data file: " << 	pSecond->allele1 << " " <<pSecond->allele2 <<endl;
				cout <<"Allele in first data file: " << 	pFirst->allele1 << " " <<pFirst->allele2 <<endl;
				string msg= "SNP: "+pFirst->rsId+" with Allele1 \""+pSecond->allele1+ "\" or  Allele2 \""+\
				pSecond->allele2+ "\" in  second data  does not match with the SNP "+pFirst->rsId+" with allele1 \""+\
				pFirst->allele1+"\" allele2 \""+pFirst->allele2+"\"	given in  first data.\n";
				error(msg);
			}
			//


		 }
			//
		 else
		 {
				cout <<"Allele in second data file: " << 	pSecond->allele1 << " " <<pSecond->allele2 <<endl;
				cout <<"Allele in first data file: " << 	pFirst->allele1 << " " <<pFirst->allele2 <<endl;
				string msg= "Problem in updating allele Information: \n"
						"SNP: "+pFirst->rsId+" with Allele1 \""+pSecond->allele1+ "\" or  Allele2 \""+\
						pSecond->allele2+ "\" in  second data  does not match with the SNP "+pFirst->rsId+" with allele1 \""+\
						pFirst->allele1+"\" allele2 \""+pFirst->allele2+"\"	given in  first data.\n";
				error(msg);

		 }

	}
return _tfvec;}
//xxxxxxxxxxxxxxxxxxxx
void CGENERAL::check_n_update_snpInfo(CBSNP* const _pFirst, const CBSNP*const _pSecond)
{
		string _snpid1 		=_pFirst->snpName;
		string _rsid1		=_pFirst->rsId;
		int 	_bp1		= _pFirst->bp;
		double cm_pos1		=_pFirst->cm_pos;
		//for second
		string _snpid2 		=_pSecond->snpName;
		string _rsid2		=_pSecond->rsId;
		int 	_bp2		= _pSecond->bp;
		double cm_pos2		=_pSecond->cm_pos;
		//cond1
		bool cond1			= (_snpid1 =="") && (_snpid2!="");
		bool cond2			= (_rsid1 =="") && (_rsid2!="");
		bool cond3 			= (_bp1==-1) && (_bp2!=-1);
		bool cond4			= (cm_pos1==-1) && (cm_pos2!=-1);
		if(cond1)
			_pFirst->snpName= _snpid2;
		if(cond2)
			_pFirst->rsId 	= _rsid2;
		if(cond3)
			_pFirst->bp 	= _bp2;
		if(cond4)
			_pFirst->cm_pos = cm_pos2;
return;}
//xxxxxxxxxxxxx
void CGENERAL::check_n_make_same_genoFormat(const CBSNP* const pFirst, CBSNP* const pSecond, const float &thresh)
{
	bool fgeno		=false;
	bool fpgeno		=false;
	bool sgeno		=false;
	bool spgeno		=false;
	if(pFirst->geno1.size()>0&&(pFirst->geno1.size()==pFirst->geno2.size()) )
		fgeno=true;
	if(pSecond->geno1.size()>0&&(pSecond->geno1.size()==pSecond->geno2.size()) )
		sgeno=true;

	//pgeno
	bool cond_1 =(pFirst->pgeno1.size()>0);
	bool cond_2 =(pFirst->pgeno1.size()==pFirst->pgeno2.size());
	bool cond_3 =(pFirst->pgeno1.size()==pFirst->pgeno3.size());
	//cout <<"cond_1&&(cond_2&&cond_3) "<<(cond_1&&(cond_2&&cond_3))<<", fpgeno "<<fpgeno<<endl;
	if(cond_1&&(cond_2&&cond_3))
	{
		fpgeno=true;
	}
	//second
	//pgeno
	cond_1 =(pSecond->pgeno1.size()>0);
	cond_2 =(pSecond->pgeno1.size()==pSecond->pgeno2.size());
	cond_3 =(pSecond->pgeno1.size()==pSecond->pgeno3.size());
	if(cond_1 &&cond_2&&cond_3)
		spgeno=true;
	if(fgeno && !sgeno)
	{
		//cout << "change geno of 2nd \n";
		CBSNP::convert_genoprob_into_hardcalls(pSecond, thresh);
	}
	//cout << fpgeno << " second: "<< spgeno <<endl;
	if(fpgeno &&!spgeno)
	{
		//cout << "change pgeno of 2nd \n";
		//CBSNP::convert_hardcalls_into_genoprob(pSecond);
		CBSNP::convert_hardcalls_into_genoprob(pSecond);

	}

}

void CGENERAL::write_phase_formatted_files(const vector<CBPED*>&pedVec, const vector<CBSNP*> genVec,
		const string &outFile){

	const int NINDIVS =CBPED::no_of_good_indivs(pedVec);
	const int NSNPS	 =CBSNP::no_of_good_snps(genVec);
	unsigned int i=0;
	unsigned int j=0;
	const CBSNP* pCBSNP=genVec[0];
	const CBPED* pCBPED=pedVec[0];
	const string outName =outFile+".inp";
	ofstream ofs;
	ofs.clear(); ofs.close();
	//opening file
	ofs.open(outName.c_str());
	if(!ofs){
		error("file \""+outFile+"\" can not be written. Please check   whether you have given the correct"
				"file name and file path\n");
	}
	// write no of good individuals individuals
	ofs << NINDIVS<<"\n";
	ofs<<NSNPS<<"\n"; // no of good individuals
	ofs<<"P"; // writing positions

	for(j=0;j<genVec.size();++j)
	{
		pCBSNP =genVec[j];
		if(pCBSNP->quality)
			ofs<<" "<< pCBSNP->bp ;
	}
	//
	ofs<<"\n"; //new line
	//type of makers: all are SNPS so S
	for(j=0;j<genVec.size();++j)
	{
		pCBSNP =genVec[j];
		if(pCBSNP->quality)
			ofs<<"S" ;
	}
	ofs<<"\n"; //new line
	//writing id of an individuals
	for(i=0;i<pedVec.size();++i)
	{
		pCBPED=pedVec[i];
		if(pCBPED->quality){
			ofs<<pCBPED->indId<<"\n"; // write sampleid
			//next line will be the first genotypes
			vector<string>_tmp_2nd_coding(genVec.size());
			for(j=0;j<genVec.size();++j)
			{
				pCBSNP =genVec[j];
				if(pCBSNP->quality)
				{
					if(pCBSNP->geno1.size()!=pedVec.size() || pCBSNP->geno2.size()!=pedVec.size() )
					{
						cout <<"geno1.size() : "<<pCBSNP->geno1.size() <<"geno2.size():  "<<pCBSNP->geno1.size();
						//for(unsigned int i=0;i<pCBSNP->geno1.size();++i)
						//	ofs<< pCBSNP->geno1[i] << " " <<pCBSNP->geno2[i] << " ";
						ofs.close();
						error("Problem in writing out "+change_int_into_string(j) +"th SNP. geno1 and geno2 has no equal size\n");
					}
					else if(pCBSNP->allele1!="" )
					{
						if( !pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if homozygote1  00
						{	ofs <<pCBSNP->allele1;
							_tmp_2nd_coding[j]=pCBSNP->allele1;
						}
						else if( !pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if heterozygote
						{
							//ofs << " "<<pCBSNP->allele1 << " "<<pCBSNP->allele2;
							if(pCBSNP->aOrder[i])
							{
								ofs <<pCBSNP->allele1;
								_tmp_2nd_coding[j]=pCBSNP->allele2;
							}
							if(!pCBSNP->aOrder[i])
							{
								ofs <<pCBSNP->allele2;
								_tmp_2nd_coding[j]=pCBSNP->allele1;
							}
						}
						else if( pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if homozygote2 11
							{
								ofs <<pCBSNP->allele2;
								_tmp_2nd_coding[j]=pCBSNP->allele2;
							}
						else if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if missing
							{
								ofs <<"?";
								_tmp_2nd_coding[j]="?";
							}

						else
						{
							error("Problem in writing out "+change_int_into_string(j) +"th SNP.\n");
						}
					}
					else
					{
						error("Problem in writing out "+change_int_into_string(j) +"th SNP.\n");

					}
				}
			}// end of j genVec loop
			ofs<<"\n";//end of line
			// write the second coding
			for(j=0;j<genVec.size();++j){
				pCBSNP=genVec[j];
				if(pCBSNP->quality)
					ofs<<_tmp_2nd_coding[j];
			}
			ofs<<"\n";//end of line
		}//end of if pCBPED->quality loop
	}//end of pedVec i loop
	//closing file
	ofs.clear(); ofs.close();
  	printLIN("\t**->phase-formatted input  file has been written out and saved as \""+outName+\
  			"\".\n");

}

// read popinfo file for calculating fst
void CGENERAL::read_popInfo_file(const vector<CBPED*> &pedVec,const string& inFileName){
	cout << "Hi from popInfo reading file\n";
}

//xxxxxxxxxxxxxxxxxxxx
//cout << "test5"<<endl; //debug
		//allele1
	//cout << col_allele1 << " "<<col_allele2 <<endl;	//test


void calculate_nei_wc_hudson_fst(vector<CBPED*>&pedVec, const vector<CBSNP*>&genVec, const string& groupLabelFileName,const string &outFileName, const string& ref_allele)
{
	bool  _switch =true; //this can be used later for anything so I made loop . 
	if(_switch)
	{
		//read group label file first 
		CSMARTPCAPED::read_smartpca_groupLabel_file(pedVec,groupLabelFileName);
		printLIN("*->Calculating Fst:\n");
		vector<string>popNames;
		//start to separate the indices according as group label. 
		//take first pid 
		const CBPED* pCBPED 	=pedVec[0];
		int cur_val   	=0;
		vector<string>::iterator it;
		Matrix< int > pop_group_idxs(1,0);
		//pop_group_idxs is matrix of size(NPOP, nindiv) in each pop. that means this may not be a matrix but a list  
		//with different rows
		if(pCBPED->groupLabel!="")
		{
			popNames.push_back(pCBPED->groupLabel);
			//pop_group_idxs.data.resize(1);
			pop_group_idxs.data[0].push_back(0); // this will tell you the index of the group-label that comes for the first time. 
			// the first group label comes of course at the 0 index for c++ i.e. in the 1th place 
		}
		else
		{
			delete pCBPED ;
			string msg= "To calculate F_st, each individual should be assigned to a certain population group.\n"
						"But here the first individual is not assigned. \n";
			error(msg);
		}
		// from here the second indiv 
		for(unsigned int popIdx=1;popIdx<pedVec.size();popIdx++)
		{
			 pCBPED =pedVec[popIdx];
			if(pCBPED->groupLabel!="")
			{
				it=find(popNames.begin(),popNames.end(),pCBPED->groupLabel); // search if the grouplabel is already  there. 
				if(it==popNames.end()) //if  no 
				{
					popNames.push_back(pCBPED->groupLabel); //add group label 
					cur_val =(int) pop_group_idxs.data.size(); // find the index of  first indiv in the data set which has different grouplabel
					pop_group_idxs.data.resize((cur_val+1));
					pop_group_idxs.data[cur_val].push_back(popIdx);
				}else
				{
					//this group Label is given already in popNames
					cur_val =(int) (it-popNames.begin());
					//cout <<boolalpha<< (*it==popNames[cur_val])<<endl;
					pop_group_idxs.data[cur_val].push_back(popIdx);

				}
			}
			else
			{
				delete pCBPED ;
				string msg= "To calculate F_st, each individual should be assigned to a certain population group.\n"
								"But here the"+change_int_into_string(popIdx+1)+"th individual is not assigned. \n";
				error(msg);
			}
			// Till now we have found the  indicies of pids which  belongs to whom. 
			//pop_group_indxs (NPOP, n_indiv)  no of pos times no of indiviudals // each element contains the indices of the individual in pedVec 
			//next job 
			//---------------------------------//
			
		} // end of for popIdx loop
		//cout << "popNames.size(): "<< popNames.size()<<endl;
		//cout << pop_group_idxs.data.size()<<endl;
		//now determine the major allele as first allele for each SNP 
		//first determine 		
		//calculate allele frequencies 
		CBPED::calculate_maf(genVec,pedVec);
    	//basic variables necessary 
		CBSNP* pCBSNP =genVec[0];
		unsigned int gsTotal 	=0; // good snps
		unsigned int gpTotal	=0; // good pedigrees
		unsigned int cur_idx	=0; // this will show the current indx of indiv in popVec 
		unsigned int _tmp_ones 	=0; // total heterozygotes 
		unsigned int _tmp_twos 	=0; //total homo1  homo of reference allele
		//create matrices to save fsts and homozygote values , heterozygote values etc 
		//create  matrixes
		// no of columns are not given.
		//This is because every population can have different number of non missing snps
			// 2 =homo1(ref allele) , 1=hetero 0 =homo2 , 3 =missing
		unsigned const int N_SPOP =popNames.size();
		//these matrices will be of size no of subgroups x no of SNPs 
		Matrix<int> N2_mat(N_SPOP, 0); //matrix of homos denoted by 2. 
		Matrix<int> N1_mat(N_SPOP, 0); //n heterozygotes 
		Matrix<int>	N_mat(N_SPOP, 0); 	// no of persons with non missing genotypes in the snp 
		//Matrix<double> Het(N_SPOP, 0, 0.0);
		Matrix<double> Freq(N_SPOP, 0, 0.0);
		//calculate unbiased heterozygosity
		//Matrix<double> Het_reich(N_SPOP,0,0.0);
		//
		//written on 2014.07.21 
		Matrix<double> fst_wc(N_SPOP,N_SPOP,0.0); // weir and cokerman 
		Matrix<double> fst_hudson(N_SPOP,N_SPOP,0.0);//hudson
		Matrix<double> fst_nei(N_SPOP,N_SPOP,0.0);
		Matrix<double> fst_reich(N_SPOP,N_SPOP,0.0);
		Matrix<double> fst_wc_ori(N_SPOP,N_SPOP,0.0);
		
		for(unsigned int j=0;j<genVec.size();j++)
		{
			
			pCBSNP =genVec[j];
			if(pCBSNP->quality)
			{	
				++gsTotal ; // good snp is added.
				vector<int> _new_012;
				CBSNP::convert_snp_genotype_into_0123(pCBSNP);
				_new_012=provide_correct_geno_012(pCBSNP,ref_allele);
				//do i need to know which allele is reference allele?
				// Now not if yes I have to recognize this also.
				for (unsigned int popIdx=0;popIdx<pop_group_idxs.data.size();popIdx++)
				{
					unsigned const int _cur_sz=pop_group_idxs.data[popIdx].size();
					for( unsigned int i=0;i<_cur_sz;i++)
					{

						cur_idx 	=pop_group_idxs.data[popIdx][i];
						pCBPED		=pedVec[cur_idx];
						// cout <<"\npCBPED->quality: "<<pCBPED->quality<<endl;
						if(pCBPED->quality)
						{

						if(_new_012[cur_idx]!=3)
						  ++gpTotal;
						if(_new_012[cur_idx]==1)
						  ++_tmp_ones;
						else if(_new_012[cur_idx]==2)
						  ++_tmp_twos;
					}
					cur_idx=0;

				}
				// here is  counting for each pops is finished now save
				N2_mat.data[popIdx].push_back(_tmp_twos);
				N1_mat.data[popIdx].push_back(_tmp_ones);
				N_mat.data[popIdx].push_back(gpTotal);
				//Het.data[popIdx].push_back((double)(_tmp_ones*1.0/gpTotal));
				//here comes rest for the SNP
				const int TOTAL_ALLELES		=2*gpTotal;
				const int ALLELECOUNT 		= _tmp_ones +_tmp_twos*2;
				//cout << "ALLELECOUNT: "<<ALLELECOUNT<<", ";
				 float ALLELE_FREQ=0.0;
				if(TOTAL_ALLELES!=0)
				{
					ALLELE_FREQ		=(double) (ALLELECOUNT*1.0)/TOTAL_ALLELES;
				}else
				{
					ALLELE_FREQ		=0.0;
				}	
					//cout << " ALLELE_FREQ: "<<ALLELE_FREQ<< endl;
				Freq.data[popIdx].push_back(ALLELE_FREQ);
					// calculate heterozygosity

				gpTotal 		=0;
				_tmp_ones 		=0;
				_tmp_twos		=0;

			}

			//--------------------------------------------------------//
				
			}//end of pCBSNP->quality loop 
		}
		//Till here 	
		//now check if the matrices are saved correctly 
		for(unsigned int j=0; j<N_SPOP; ++j)
		{
			/*cout << popNames[j]<<": \n";
			for(i=0;i<gsTotal;++i){
			//cout <<" "<<N1_mat.data[j][i]<<endl;
			//cout <<" "<<N_mat.data[j][i]<<endl;
			//cout <<" "<<N1_mat.data[j][i]*1.0/N_mat.data[j][i]<<endl;
			cout <<" "<<Freq.data[j][i];
		}
		cout << "\n";
		*/
			bool b1=(gsTotal==N2_mat.data[j].size());
			bool b2=(gsTotal==N1_mat.data[j].size());
			bool b3=(gsTotal==N_mat.data[j].size() );
			//bool b4=(gsTotal==Het.data[j].size() );
			bool b5=(gsTotal==Freq.data[j].size() );
			bool _b= b1|| b2|| b3 || b5 ;
			if(!_b)
			{
				//cout<<gsTotal << ", "<< N0_mat.data[j].size() <<"\n ";
				//cout << N1_mat.data[j].size() <<endl;
				//cout << N_mat.data[j].size() <<endl ;
				//cout << Het.data[j].size() 	<<endl;
				//cout << Freq.data[j].size() <<endl;
				string msg ="Problem by calculating heteroginity and allele freq for\n"
						"different population. Please contact the software developer\n"
						"problem is in function: handel_with_fst_calculation\n";
				error(msg);
			}
	}
		//now calculate each pair of pops the F_st
		//for every population do something 
		//for all 
		double P_i			=0.0;
		double P_j			=0.0;
	 	//for nei 
		double deno_nei		=0.0;
		double num_nei 		=0.0;
		double P_avg		=0.0; 
		double _FST_num_nei	=0.0; //local for each snp 
		double _FST_deno_nei =0.0; //local for each snp 
		//for wc 
		//double num_wc_first		=0.0;
		//double num_wc_2nd		=0.0;
		//double num_wc_3rd		=0.0;
		double deno_wc_first	=0.0;
		double deno_wc_2nd		=0.0;
		double deno_wc_3rd		=0.0;
		double num_wc 			=0.0;
		double deno_wc			=0.0; 
		//double _FST_wc			=0.0; //local for each SNP	
		double _FST_num_wc			=0.0; //local for each SNP	
		double _FST_deno_wc			=0.0; //local for each SNP	
		//for hudson
		unsigned int N_i		=0;
		unsigned int N_j		=0;
		//double 	 num_first_term	=0.0;
		//double 	 num_2nd_term	=0.0;
		//double 	 num_3rd_term	=0.0;
		double num_hudson		=0.0;
		double deno_hudson		=0.0;
		//double _FST_hudson		=0.0; //local for each SNP
		double _FST_num_hudson		=0.0; //local for each SNP
		double _FST_deno_hudson		=0.0; //local for each SNP
		//for Reich 2009 
		double A_i				=0.0;
		double A_j				=0.0;
		double H_i_hat			=0.0; 
		double H_j_hat			=0.0; 
		double N_hat			=0.0; 
		double D_hat			=0.0;
		double N_hat_total		=0.0;
		double D_hat_total		=0.0;
		double _FST_reich		=0.0;
		//for weir and cockerham original  
		const int R =2; // two population taken at a time
		 double N_bar 		=0.0;
		 double S_square	=0.0;
		 double P_bar 		=0.0;
		 double H_i 		=0.0;
		 double H_j			=0.0;
		 double H_bar		=0.0;
		 double num_wc_ori	=0.0;
		 double deno_wc_ori	=0.0;
		 double num_wc_ori_total 	=0.0;
		 double deno_wc_ori_total	=0.0;
			  
		for (unsigned int i=0; i<N_SPOP;++i )
		{	for (unsigned int j=0; j<N_SPOP;++j )
			{
			//if same population then zero.

			if(i==j)
			{
				fst_wc.data[i][j] 		=0.0; 	// with itself 
				fst_hudson.data[i][j]	=0.0;	//with itself 
				fst_nei.data[i][j]		=0.0;	//with itself 
				fst_reich.data[i][j]	=0.0;	//with itself 
			
			}
			else if(i<j) // calculate only the upper half of the matrix.
			{

				for (unsigned int snpIdx=0; snpIdx<gsTotal;++snpIdx) //actually popIdx is snpidx. but we used here because we didn't want to define another 
				{
					/**
						* Fst has been calculated using the formula in the following paper:
						*1. Estimating and interpreting Fst: The impact of rare variants 
						* Access the most recent version at doi:10.1101/gr.154831.113
						*Genome Res. 2013 23: 1514-1521 originally published online July 16, 2013
						*Gaurav Bhatia, Nick Patterson, Sriram Sankararaman, et al.
						*for Nei : equation  s14
						*for hudson : equation s18 
						*for weir and cockerham:  equation s7
						*For Reich Fst : 
						  Reconstructing Indian population history, 24 September 2009| doi:10.1038/nature08365
						 David Reich1,2*, Kumarasamy Thangaraj3*, Nick Patterson2*, Alkes L. Price2,4* & Lalji Singh3
						 Equation 10 11, and 12 in the supplimentary material 
						*for Original weir and cockerham 1984: 
						Estimating F-Statistics for the Analysis of Population Structure
						B. S. Weir; C. Clark Cockerham
						Evolution, Vol. 38, No. 6. (Nov., 1984), pp. 1358-1370.
						Stable URL: http://links.jstor.org/sici?sici=0014-3820%28198411%2938%3A6%3C1358%3AEFFTAO%3E2.0.CO%3B2-0
						equation after equation 7 on page 1361 thita_hat =.. 
					*/
					//allelle frequencies 
					P_i			=Freq.data[i][snpIdx];
					P_j			=Freq.data[j][snpIdx];
					//stat Nei  
					//----------------//	
					if(P_i==0.0&& P_j==0.0){
						num_nei=0.0; 	
						deno_nei=0.0;
					}else
					{
						P_avg		=(P_i+P_j)/2.0;
						deno_nei	=2.0*P_avg*(1-P_avg);
						num_nei		=1.0*pow(P_i-P_j,0.2e1); // 0.2e1=2.0;
					}	
					//if(deno_nei!=0.0)
					//	_FST_nei +=num_nei/deno_nei;
					//else 	
					//	_FST_nei	+=0.0;
					//
					_FST_num_nei +=num_nei; 
					_FST_deno_nei +=deno_nei;					
					
					//cout <<"Test 1" <<endl ;//-----------------------------------------------------------------------------------------------------------test //
					//end of nei 
					//---------------sa-//	
					//start hudson 
					//average of non missing
					//for large samle size 
					num_hudson =num_nei;
					deno_hudson	=P_i*(1-P_j) +P_j*(1-P_i);
					
					//for general case: delete the above two lines and write it like below.  
					/*N_i 		=(N_mat.data[i][snpIdx]);
					N_j 		=(N_mat.data[j][snpIdx]);
					num_first_term 	= num_nei;
					//second term 
					if(N_i<=1)
					 num_2nd_term 	= P_i*(1.0-P_j);
					 else 
						num_2nd_term 	= (P_i*(1.0-P_j))/(N_i-1);
					//third term 
					if(N_j<=1)
					num_3rd_term 	= P_j*(1.0-P_i);
					 else 
						num_3rd_term 	= (P_j*(1.0-P_i))/(N_j-1);
					num_hudson	=num_first_term + num_2nd_term + num_3rd_term;
					deno_hudson	=P_i*(1-P_j) +P_j*(1-P_i);
					*/
					//if(deno_hudson!=0.0)
					//	_FST_hudson	+= num_hudson/deno_hudson;  //local for each SNP
					//else 
					//_FST_hudson	+= 0.0;  //local for each SNP
					_FST_num_hudson	 +=num_hudson; 
					_FST_deno_hudson +=deno_hudson;  //local for each SNP
					
					//cout <<"Test 2" <<endl ; //-----------------------------------------------------------------------------------------------------------test //
					
					//----------------//	
					//end hudson 
					//start wc 
					//----------------//	
					
					//----------------//	
					//end hudson 
					//start wc 
					//----------------//	
					//for large sample saize 
					num_wc 			=num_nei;
					deno_wc_first 	=num_nei;
					N_i 		=2.0*(N_mat.data[i][snpIdx]);
					N_j 		=2.0*(N_mat.data[j][snpIdx]);
					if(N_j!=0)
					deno_wc_2nd	= 1.0/(N_i/N_j+1.0);
					else 
					deno_wc_2nd	= 0.0;
					
					if(N_j!=0)	
						deno_wc_3rd = N_i/N_j*(P_i*(1.0-P_i))+P_j*(1.0-P_j);
					else 
						deno_wc_3rd = P_j*(1.0-P_j);
					deno_wc		=deno_wc_first+deno_wc_2nd*deno_wc_3rd;
					
					//cout <<"Test 3" <<endl ;//-----------------------------------------------------------------------------------------------------------test //
					//for general case 
					/*	
					num_wc_first 	=2.0*(N_i*N_j)/(N_i+N_j);
					if((N_i+N_j)>2)
						num_wc_2nd	=1.0/(N_i+N_j-2);
					else
						num_wc_2nd  =0.0;
					num_wc_3rd 		=N_i*P_i*(1-P_i)+N_j*P_j*(1-P_j);	
					deno_wc_first	=1.0*(N_i*N_j)/(N_i+N_j)*pow(P_i-P_j,0.2e1);
					deno_wc_2nd		=(2.0*((N_i*N_j)/(N_i+N_j))-1.0);
					deno_wc_2nd		*=num_wc_2nd;
					deno_wc_3rd		=num_wc_3rd;
					num_wc 			=num_wc_first*num_wc_2nd*num_wc_3rd;
					deno_wc			=deno_wc_first+deno_wc_2nd* deno_wc_3rd; 
					
					if(deno_wc!=0)
					_FST_wc			+=num_wc/deno_wc; //local for each SNP	
					else 
					_FST_wc			+=0.0; //local for each SNP	
					*/
					_FST_num_wc	+= num_wc;
					_FST_deno_wc += deno_wc; //local for each SNP	
					
					//----------------//	
					// end wc
					//cout <<"Test 4" <<endl ; //-----------------------------------------------------------------------------------------------------------test //
					//start of Reich 2009
					N_i 		=2.0*(N_mat.data[i][snpIdx]);
					N_j 		=2.0*(N_mat.data[j][snpIdx]);
					A_i 		=2.0*(N2_mat.data[i][snpIdx])+(N1_mat.data[i][snpIdx]);
					A_j 		=2.0*(N2_mat.data[j][snpIdx])+(N1_mat.data[j][snpIdx]);
					if(N_i>1)
					 H_i_hat	=(A_i/N_i)*((N_i-A_i)/(N_i-1));
					else 	
					H_i_hat		=0.0;
					//H_j_hat
					if(N_j>1)
					 H_j_hat	=(A_j/N_j)*((N_j-A_j)/(N_j-1));
					else 	
					H_j_hat		=0.0;
					N_hat		=pow(P_i-P_j,0.2e1); 
					if(N_i>0)
						N_hat 	-=H_i_hat/N_i;
					if(N_j>0)
						N_hat 	-=H_j_hat/N_j;
					//for denominator 	
					D_hat		=N_hat+H_i_hat+H_j_hat; 
					N_hat_total +=N_hat; 
					D_hat_total	+=D_hat;
					//end of Reich
					//-start of wier an cockerham 1994 original 
					//average of non missing
					//Now we implement the formula given  in page 1361 of Weir and Cockerham 1984 
					//Estimating F-Statistics for the Analysis of Population Structure
					//B. S. Weir; C. Clark Cockerham
					//Evolution, Vol. 38, No. 6. (Nov., 1984), pp. 1358-1370.
					//Stable URL: http://links.jstor.org/sici?sici=0014-3820%28198411%2938%3A6%3C1358%3AEFFTAO%3E2.0.CO%3B2-0
					N_bar 		=(N_i+N_j)*1.0/(R*1.0); // average no of population 
					P_i			=Freq.data[i][snpIdx];
					P_j			=Freq.data[j][snpIdx];
					//cout<<"N_bar:"<<N_bar <<endl;
					if(N_bar>0)
						P_bar 	=((N_i*P_i)+(N_j*P_j))/(R*N_bar);
					//cout<< P_i<< " P_bar: "<< P_bar <<endl; 
					//S^2 = 	sum(N_i(P_i-Pbar)^2/((r- l ) n) , the sample variance of allele A frequencies over populations
					S_square	 = N_i * (P_i-P_bar) *(P_i-P_bar); 
					//cout<< "S_square: "<<S_square <<", "<< endl;
					S_square	+= N_j * (P_j-P_bar) *(P_j-P_bar) ;
					//cout<< "S_square: "<<S_square <<", "<< endl;
					S_square	*=1.0/( ((R-1.0)*N_bar));
					//cout<< "S_square: "<<S_square <<", "<< endl;
					
					H_i			 =2.0*P_i *(1-P_i);
					H_j			 =2.0*P_j *(1-P_j);
					//cout<< "S_square: "<<S_square <<", "<< endl;
					//cout <<"Test 5" <<endl ; //-----------------------------------------------------------------------------------------------------------test //
					if(N_bar>0)
					H_bar		=(N_i*H_i+N_j*H_j)/(N_bar*R);
					if(N_bar>1 )
					num_wc_ori =S_square-1.0/(N_bar-1)*(P_bar*(1.0-P_bar)-((R-1.0)/R)*S_square -1.0/4*(H_bar)); 
					else 
					num_wc_ori =S_square; 
					deno_wc_ori	=P_bar*(1-P_bar)+S_square/(R*1.0);
					//cout << " num: "<< num_wc_ori<< ", deno: "<<deno_wc_ori<< endl; 
					num_wc_ori_total +=num_wc_ori; 
					deno_wc_ori_total+=deno_wc_ori;
					
					//cout <<"Test 6" <<endl ;//-----------------------------------------------------------------------------------------------------------test //
					//assigning again to zero
					num_wc_ori		=0.0;
					deno_wc_ori		=0.0; 
					//num_wc_first 	=0.0;
					//num_wc_2nd  	=0.0;
					//num_wc_3rd 		=0.0; 
					deno_wc_first	=0.0;
					deno_wc_2nd		=0.0;
					deno_wc_2nd		=0.0;
					deno_wc_3rd		=0.0;
					num_wc 			=0.0;
					deno_wc			=0.0; 
					P_i=P_j			=0.0;
					N_i =N_j 		=0;
					P_avg			=0.0;
					deno_nei		=0.0;
					num_nei			=0.0;
					//num_first_term 	= 0.0;
					//num_2nd_term 	= 0.0;
					//num_3rd_term 	= 0.0;
					num_hudson		=0.0;
					deno_hudson		=0.0;
					H_i_hat			=0.0; 
					H_j_hat			=0.0; 
					N_hat			=0.0; 
					D_hat			=0.0; 

					//--------------------------------------------//
				}
				
				
				// After calculating a_sum , b_sum and c_sum
				//calculate f_st
				//fst_nei.data[i][j]	 	=_FST_nei*1.0/gsTotal;
				if(_FST_deno_nei!=0)
				fst_nei.data[i][j]	 	=_FST_num_nei/_FST_deno_nei;
				else 
				fst_nei.data[i][j]	 	=0.0;
				// hudson 
				if(_FST_deno_hudson!=0)
				fst_hudson.data[i][j]	=_FST_num_hudson/_FST_deno_hudson;
				else 
				fst_hudson.data[i][j]	=0.0;
				//weir and cockerham 
				if(_FST_deno_wc!=0.0)
				 fst_wc.data[i][j] 		=_FST_num_wc/_FST_deno_wc;
				 else 
				  fst_wc.data[i][j] 		=0.0;
				
				//fst reich
				if(D_hat_total==0.0)
					_FST_reich 				=0.0;
				else 
				_FST_reich 				=N_hat_total/D_hat_total; 
				fst_reich.data[i][j] 	=_FST_reich;
				//fst wc_ori
				if(deno_wc_ori_total==0.0)
					fst_wc_ori.data[i][j] 	=0.0;
				else 
					fst_wc_ori.data[i][j] 	=num_wc_ori_total/deno_wc_ori_total;
				// set to zero again
				num_wc_ori_total	=0.0; 
				deno_wc_ori_total	=0.0;
				//_FST_nei 			=0.0;
				_FST_num_nei		=0.0;
				_FST_deno_nei		=0.0;
				//_FST_hudson 		=0.0;
				_FST_num_hudson 		=0.0;
				_FST_deno_hudson 		=0.0;
				//_FST_wc				=0.0;
				_FST_num_wc				=0.0;
				_FST_deno_wc			=0.0;
				_FST_reich 			=0.0;
			}
			else
			{
				fst_wc.data[i][j] 		=fst_wc.data[j][i];
				fst_hudson.data[i][j]	=fst_hudson.data[j][i];
				fst_nei.data[i][j]		=fst_nei.data[j][i];
				fst_reich.data[i][j]	=fst_reich.data[j][i];
				fst_wc_ori.data[i][j]	=fst_wc_ori.data[j][i];	
			}
			//cout <<"F_st["<<i+1<<","<<j+1<<"]= "<<fst_wc_ori.data[i][j]<< "  ";
			//Population pairwise FST values. A negative symbol (-) indicates that pairwise FST values were not significant.
			//A plus (+) sign indicates a significant FST value adjusted for the number of contrasts using a sequential Bonferroni test
		}
		}
	
	//------------------------------------------------------//
		//write fst now. 
		// open a file to write the output of F_st
		string msg_out = "output file" + outFileName+" does not exits or can't be written out. Please check the file path and file name if they are given correctly.\n";
		ofstream ofs_nei ,ofs_wc,ofs_hudson, ofs_reich,ofs_wc_ori;
		const string wc_fileName= outFileName+"_Fst_modified_Weir_Cockerham_1984.txt";
		const string nei_fileName= outFileName+"_Fst_modified_Nei_1986.txt";
		const string hudson_fileName= outFileName+"_Fst_modified_Hudson_1992.txt";
		const string reich_fileName= outFileName+"_Fst_Reich_2009.txt";
		const string wc_ori_fileName= outFileName+"_Fst_Weir_Cockerham_1984.txt";
		//clear streams ;
		ofs_nei.clear();	 ofs_nei.close();
		ofs_wc.clear(); 	 ofs_wc.close();
		ofs_wc_ori.clear(); 	 ofs_wc_ori.close();
		ofs_hudson.clear();  ofs_hudson.close();
		ofs_reich.clear(); 	 ofs_reich.close();
	// open now
		ofs_nei.open(nei_fileName.c_str(),ios::out);
		ofs_wc.open(wc_fileName.c_str(),ios::out);
		ofs_wc_ori.open(wc_ori_fileName.c_str(),ios::out);
		ofs_hudson.open(hudson_fileName.c_str(),ios::out);
		ofs_reich.open(reich_fileName.c_str(),ios::out);
		if(!ofs_nei)
			error(msg_out);
		if(!ofs_wc)
			error(msg_out);
		if(!ofs_hudson)
			error(msg_out);
		if(!ofs_reich) error(msg_out);	
		if(!ofs_wc_ori)
			error(msg_out);
		
		//
		//write population names
		ofs_nei<<"popName";
		ofs_wc<<"popName";
		ofs_wc_ori<<"popName";
		ofs_hudson<<"popName";
		ofs_reich<<"popName";
		//close clear
		for (unsigned int i=0; i<N_SPOP; ++i)
		{
			ofs_nei << "\t"<<popNames[i] ;
			ofs_wc << "\t"<<popNames[i] ;
			ofs_wc_ori << "\t"<<popNames[i] ;
			ofs_hudson << "\t"<<popNames[i] ;
			ofs_reich <<"\t"<<popNames[i];
		}
		ofs_nei<< "\n";
		ofs_wc<< "\n";
		ofs_wc_ori<< "\n";
		ofs_hudson<< "\n";
		ofs_reich<<"\n";

		// end of the start of file of wright's fst.
		for (unsigned int i=0; i<N_SPOP;++i )
		{
			ofs_nei<<popNames[i] << "\t";
			ofs_wc<<popNames[i] << "\t";
			ofs_wc_ori<<popNames[i] << "\t";
			ofs_hudson<<popNames[i] << "\t";
			ofs_reich<<popNames[i]<<"\t";
			for (unsigned int j=0; j<N_SPOP;++j )
			{
				ofs_nei <<fst_nei.data[i][j]<< "\t";
				ofs_wc <<fst_wc.data[i][j]<< "\t";
				ofs_wc_ori<<fst_wc_ori.data[i][j]<< "\t";
				ofs_hudson <<fst_hudson.data[i][j]<< "\t";
				ofs_reich<<fst_reich.data[i][j]<<"\t";
				//cout <<"F_st["<<i+1<<","<<j+1<<"]= "<<F_st.data[i][j]<< "  ";
				//Population pairwise FST values. A negative symbol (-) indicates that pairwise FST values were not significant.
				//A plus (+) sign indicates a significant FST value adjusted for the number of contrasts using a sequential Bonferroni test
			}
			ofs_nei<<"\n";
			ofs_wc<<"\n";
			ofs_wc_ori<<"\n";
			ofs_hudson<<"\n";
			ofs_reich<<"\n";

			//cout << "\n";
		}
		//###############################################################################
		ofs_nei.clear(); 	ofs_nei.close();
		ofs_wc.clear(); 	ofs_wc.close();
		ofs_wc_ori.clear(); 	ofs_wc_ori.close();
		ofs_hudson.clear(); ofs_hudson.close();
		ofs_reich.clear(); ofs_reich.close();
		string msg = "\t**->Fst using modified Weir and Cockerham Fst nhas been calculated and saved in the file: "+wc_fileName+". \n";
		printLIN(msg);
		msg = "\t**->Fst (Theta_hat)  using  the method (EQUATION 1) of \"Weir & Cockerham 1984, Evolution 38(6):1358-1370\" \n"
				"has been calculated and saved in the file: "+wc_ori_fileName+". \n";
		printLIN(msg);
		//Reich 2009 file
		msg = "\t**->Fst using Nei 1986 definition is written and saved as : "+nei_fileName+". \n";
		printLIN(msg);
		msg = "\t**->Fst defined by Hudson 1992 is written and saved as "+hudson_fileName+". \n";
		printLIN(msg);
		msg = "\t**->Fst defined by Reich 2009 is written and saved as "+reich_fileName+". \n";
		printLIN(msg);

			
		//-----------------------------------------------------------------------------//
			
		
	}//end of _switch loop 

}	
vector<int> provide_correct_geno_012(CBSNP* const pCBSNP,const string &ref_allele)
{
		vector<int> _new_012=pCBSNP->geno_0123;
		//unsigned const long int _nperson=_new_012.size();
	//here no if(pCBSNP->quality)
	{
		// if  reference allele means which is being counted
		// i.e. whose homzygote is saved as 2. 
		if((ref_allele=="")||(ref_allele=="major-allele"))
		{
			
			
			if ((pCBSNP->freq>=0.5) ) 
			// i.e. if the first allele is major
			{ //it is saved as 2. so  change it  
				//change 0 into 2 and 2 into zero
				std::replace (_new_012.begin(),_new_012.end(),0,99);
				std::replace (_new_012.begin(),_new_012.end(),2,0);
				std::replace (_new_012.begin(),_new_012.end(),99,0);
			}
		}
		else if((ref_allele=="minor-allele"))
		{
			if ((pCBSNP->freq<0.5) ) 
				// i.e. if the first allele is major
			{ //it is saved as 2 . so  change it  
					//change 0 into 2 and 2 into zero
				std::replace (_new_012.begin(),_new_012.end(),0,99);
				std::replace (_new_012.begin(),_new_012.end(),2,0);
				std::replace (_new_012.begin(),_new_012.end(),99,0);
			}
		}
		else if((ref_allele=="allele1"))
		{
			// first allele shoudl be saved as 2 
			std::replace (_new_012.begin(),_new_012.end(),0,99);
			std::replace (_new_012.begin(),_new_012.end(),2,0);
			std::replace (_new_012.begin(),_new_012.end(),99,0);
			
		}
		// if ref_allele is allele2 then no change 
		
	}
	return _new_012;	

}
	
void handel_fst_calculations( vector<CBPED*>&pedVec, const vector<CBSNP*>&genVec,const string& groupLabelFileName, const string & outFileName )
{
	//vector<CBPED*>&pedVec, const vector<CBSNP*>&genVec,const string& groupLabelFileName, const string & outFileName

if(true)
	{
		//read popinfo file name
		//CGENERAL::read_popInfo_file(pedVec, pBpar->popInfoFileName);
		CSMARTPCAPED::read_smartpca_groupLabel_file(pedVec,groupLabelFileName);
		printLIN("*->Calculating Fst:\n");

		vector<string>popNames;
		const CBPED* pCBPED 	=pedVec[0];
		//int nSNPs	 	=genVec.size();
		//int nInds		=pedVec.size();
		int cur_val   	=0;
		vector<string>::iterator it;
		Matrix< int > pop_group_idxs(1,0);
		//reading group label and find the corresponding pupulation indices. 
		if(pCBPED->groupLabel!="")
		{
			popNames.push_back(pCBPED->groupLabel);
			//pop_group_idxs.data.resize(1);
			pop_group_idxs.data[0].push_back(0);
		}
		else
		{
			delete pCBPED ;
			string msg= "To calculate F_st, each individual should be assigned to a certain population group.\n"
						"But here the first individual is not assigned. \n";
			error(msg);
		}

		for(unsigned int popIdx=1;popIdx<pedVec.size();popIdx++)
		{
			 pCBPED =pedVec[popIdx];
			if(pCBPED->groupLabel!="")
			{
				it=find(popNames.begin(),popNames.end(),pCBPED->groupLabel);
				if(it==popNames.end())
				{
					popNames.push_back(pCBPED->groupLabel);
					cur_val =(int) pop_group_idxs.data.size();
					pop_group_idxs.data.resize((cur_val+1));
					pop_group_idxs.data[cur_val].push_back(popIdx);
				}else{
					//this group Label is given already in popNames
					cur_val =(int) (it-popNames.begin());
					//cout <<boolalpha<< (*it==popNames[cur_val])<<endl;
					pop_group_idxs.data[cur_val].push_back(popIdx);

				}
			}
			else
			{
				delete pCBPED ;
				string msg= "To calculate F_st, each individual should be assigned to a certain population group.\n"
								"But here the"+change_int_into_string(popIdx+1)+"th individual is not assigned. \n";
				error(msg);
			}

		} // end of for popIdx loop

		unsigned const int N_SPOP =popNames.size();
		 //cout<< count( _pCBSNP->geno_0123.begin(),_pCBSNP->geno_0123.end(),0)<<endl;
		//determinie total alleles etc
		unsigned int gpTotal	 =0; // good pedigrees
		unsigned int gsTotal 	=0; // good snps
		unsigned int i,j,popIdx	=0;
		unsigned int cur_idx	=0; // this will show the current indx of pops
		unsigned int _tmp_ones =0;
		unsigned int _tmp_zeros =0;
		//unsigned int _tmp_threes =0;
		// coput << vPed.size() <<"," <<genVec.size() <<endl; //debug
		//const CBPED* pCBPED=pGENERAL->pedVec[0];
		CBSNP* pCBSNP =genVec[0];
		// cout << pCBSNP->geno_0123.size()<<endl;
		//create  matrixes
		Matrix<int> N0_mat(N_SPOP, 0);
		// no of columns are not given.
		//This is because every population can have different number of non missing snps
		Matrix<int> N1_mat(N_SPOP, 0);
		Matrix<int>	N_mat(N_SPOP, 0);
		//Matrix<double> Het(N_SPOP, 0, 0.0);
		Matrix<double> Freq(N_SPOP, 0, 0.0);
		//calculate unbiased heterozygosity
		//Matrix<double> Het_reich(N_SPOP,0,0.0);
		//
		Matrix<double> F_st(N_SPOP, N_SPOP, 0.0);
		Matrix<double> Fst_Wright(N_SPOP,N_SPOP,0.0);
		Matrix<double> Fst_Reich(N_SPOP,N_SPOP,0.0);
		//written on 2014.07.21 
		Matrix<double> fst_wc(N_SPOP,N_SPOP,0.0); // weir and cokerman 
		Matrix<double> fst_hudson(N_SPOP,N_SPOP,0.0);//hudson
		Matrix<double> fst_nei(N_SPOP,N_SPOP,0.0);
		//written on 2014.07.21 
		//-----------------------------------//

  for(j=0;j<genVec.size();j++)
  {
	  pCBSNP =genVec[j];
	  CBSNP::convert_snp_genotype_into_0123(pCBSNP);
	  //cout << "pCBSNP->quality: "<<pCBSNP->quality<<endl;
	  //cout << j+1<< " ";
	  //for each SNP  and for each population group:
	  // count the correct numbers of 0 1 and 2 if individs are ok
	//cout << "SNP "<<(j+1)<<"\n-----------------------\n"<<endl;
	 //cout<<"show 01023 matrix :\n";
	 //cout <<".....................................#\n";
	 //for(i=0;i<pCBSNP->geno_0123.size();++i)
	//	cout << " "<< pCBSNP->geno_0123[i];
	//	cout << "\n";

	if(pCBSNP->quality)
		{
			++gsTotal;
			for (popIdx=0;popIdx<pop_group_idxs.data.size();popIdx++)
			{
				unsigned const int _cur_sz=pop_group_idxs.data[popIdx].size();
				for( i=0;i<_cur_sz;i++)
				{

					cur_idx =pop_group_idxs.data[popIdx][i];
					pCBPED=pedVec[cur_idx];
					// cout <<"\npCBPED->quality: "<<pCBPED->quality<<endl;
				  if(pCBPED->quality)
				  {

					 if(pCBSNP->geno_0123[cur_idx]!=3)
						  ++gpTotal;
					  if(pCBSNP->geno_0123[cur_idx]==1)
						  ++_tmp_ones;
					  else if(pCBSNP->geno_0123[cur_idx]==0)
						  ++_tmp_zeros;
				  }
				  cur_idx=0;

				}
				// here is  counting for each pops is finished now save
				N0_mat.data[popIdx].push_back(_tmp_zeros);
				N1_mat.data[popIdx].push_back(_tmp_ones);
				N_mat.data[popIdx].push_back(gpTotal);
				//Het.data[popIdx].push_back((double)(_tmp_ones*1.0/gpTotal));
				//here comes rest for the SNP
				const int TOTAL_ALLELES		=2*gpTotal;
				const int ALLELECOUNT 		= _tmp_ones +_tmp_zeros*2;
				//cout << "ALLELECOUNT: "<<ALLELECOUNT<<", ";
				const float ALLELE_FREQ		=(double) (ALLELECOUNT*1.0)/TOTAL_ALLELES;
				//cout << " ALLELE_FREQ: "<<ALLELE_FREQ<< endl;
				Freq.data[popIdx].push_back(ALLELE_FREQ);
				// calculate heterozygosity

				gpTotal 		=0;
				_tmp_ones 		=0;
				_tmp_zeros		=0;

			}

				//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),2,99); // value 99 is used just as a temp coding to swap later .
				//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),0,2);
				//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),99,0);
		}//end of snp quality

	} //End of for genvec loop
	// Check the size of all matrices must be the same
	//cout <<"Het_1: \n...................\n";
	for(j=0; j<N_SPOP; ++j)
	{
		/*cout << popNames[j]<<": \n";
		for(i=0;i<gsTotal;++i){
			//cout <<" "<<N1_mat.data[j][i]<<endl;
			//cout <<" "<<N_mat.data[j][i]<<endl;
			//cout <<" "<<N1_mat.data[j][i]*1.0/N_mat.data[j][i]<<endl;
			cout <<" "<<Freq.data[j][i];
		}
		cout << "\n";
		*/
		bool b1=(gsTotal==N0_mat.data[j].size());
		bool b2=(gsTotal==N1_mat.data[j].size());
		bool b3=(gsTotal==N_mat.data[j].size() );
		//bool b4=(gsTotal==Het.data[j].size() );
		bool b5=(gsTotal==Freq.data[j].size() );
		bool _b= b1|| b2|| b3 || b5 ;
		if(!_b)
		{
			//cout<<gsTotal << ", "<< N0_mat.data[j].size() <<"\n ";
			//cout << N1_mat.data[j].size() <<endl;
			//cout << N_mat.data[j].size() <<endl ;
			//cout << Het.data[j].size() 	<<endl;
			//cout << Freq.data[j].size() <<endl;
			string msg ="Problem by calculating heteroginity and allele freq for\n"
						"different population. Please contact the software developer\n"
						"problem is in function: handel_with_fst_calculation\n";
			error(msg);
		}
	}
	//now calculate each pair of pops the F_st
	 vector<double> N_av(gsTotal,0.0); // average of Non missing.
	 //vector<double>
	  const int R =2; // two population taken at a time
	  double N_bar 		=0.0;
	  int 	 N_i		=0;
	  int 	 N_j		=0;
	  double S_sum		=0.0;
	  double S_square	=0.0;
	  double N_c		=0.0;
	  double P_bar 		=0.0;
	  double P_i		=0.0;
	  double P_j		=0.0;
	  double H_i 		=0.0;
	  double H_j		=0.0;
	  double H_bar		=0.0;
	  double a_val		=0.0;
	  double b_val		=0.0;
	  double c_val		=0.0;
	  double a_sum		=0.0;
	  double b_sum		=0.0;
	  double c_sum		=0.0;
	  double fst_w		=0.0;
	  double fst_w_tmp	=0.0;
	  //for reich fst
	  int a_i			=0; 	// allele count of first allele of fist pop
	  int a_j			=0; 	// allele count of first allele of second pop
	  double hr_i		=0.0;	//heterozygosity according as Reich for first pop
	  double hr_j		=0.0;	//heterozygosity
	  double Nst_r		=0.0;	// NUmerator of Fst  of Reich
	  double Dst_r		=0.0; // Denomerator of Fst of Reich
	  double fst_r_tmp	=0.0; //tmp fst_reich
	  //cout <<"gsTotal: "<< gsTotal <<endl ; //debug
	// open a file to write the output of F_st
	 string msg_out = "output file" + outFileName+" does not exits or can't be written out. Please check the file path and file name if they are given correctly.\n";
	 ofstream ofs ,ofs_wright,ofs_reich;

	//weir fst
	const string fstFileName= outFileName+"_Fst_Weir_Cockerham_1984.txt";
	// write Reich Fst
	const string fst_reich_fileName= outFileName+"_Fst_Reich_2009.txt";
	const string fst_wright_fileName= outFileName+"_Fst_Wright_1951.txt";
	//clear streams ;
	ofs.clear();	  ofs.close();
	ofs_reich.clear(); ofs_reich.close();
	// open now
	ofs.open(fstFileName.c_str(),ios::out);
	ofs_reich.open(fst_reich_fileName.c_str(),ios::out);
	ofs_wright.open(fst_wright_fileName.c_str(),ios::out);

	if(!ofs)
		error(msg_out);
	if(!ofs_reich)
		error(msg_out);
	if(!ofs_wright)
		error(msg_out);
	//
	//write population names
	ofs<<"popName";
	ofs_reich<<"popName";
	ofs_wright<<"popName";
	//close clear
	for (i=0; i<N_SPOP; ++i)
	{
		ofs << "\t"<<popNames[i] ;
		ofs_reich << "\t"<<popNames[i] ;
		ofs_wright << "\t"<<popNames[i] ;

	}
	ofs<< "\n";
	ofs_reich<< "\n";
	ofs_wright<< "\n";

	// end of the start of file of wright's fst.
	for (i=0; i<N_SPOP;++i )
	{
		ofs<<popNames[i] << "\t";
		ofs_reich<<popNames[i] << "\t";
		ofs_wright<<popNames[i] << "\t";
		for (j=0; j<N_SPOP;++j )
		{
			//if same population then zero.

			if(i==j)
			{
				F_st.data[i][j] =0.0;
				Fst_Wright.data[i][j]=0.0;
			}
			else if(i<j) // calculate only the upper half of the matrix.
			{

				for (popIdx=0; popIdx<gsTotal;++popIdx)
				{
					//average of non missing
					N_i 		=(N_mat.data[i][popIdx]);
					N_j 		=(N_mat.data[j][popIdx]);
					N_bar 		=(N_i+N_j)*1.0/R;
					S_sum		= ( (N_i*N_i) + (N_j*N_j))*1.0/(R*N_bar); // sum of square
					N_c	  		=((R*N_bar) -S_sum)*1.0/(R-1) ;
					P_i			=Freq.data[i][popIdx];
					P_j			=Freq.data[j][popIdx];
					P_bar 		 =((N_i*P_i)+(N_j*P_j))/(R*N_bar);
					S_square	 = ((N_i * ((P_i-P_bar)*(P_i-P_bar)) )+ (N_j * ((P_j-P_bar)*(P_j-P_bar))))/( ((R-1)*N_bar));
					H_i			=((N1_mat.data[i][popIdx]*1.0)/(N_mat.data[i][popIdx]));
					H_j			=((N1_mat.data[j][popIdx]*1.0)/(N_mat.data[j][popIdx]));
					H_bar		=(N_i*H_i+N_j*H_j)/(N_bar*R);
					a_val 		=(N_bar/N_c) * ( S_square - ( (1/(N_bar-1)) * ( (P_bar*(1-P_bar)) - (( 1.0*(R-1)/R)*S_square) - ((1.0/4)*H_bar))) );
					b_val 		=(N_bar/(N_bar-1)) * ((P_bar*(1.0-P_bar)) - ((1.0*(R-1)/R)*S_square) - ((((2.0*N_bar)-1)/(4.0*N_bar))*H_bar));
					c_val 		= H_bar/2.0;

					if(!isfinite(a_val))
						a_val=0.0;
					if(!isfinite(b_val))
						b_val=0.0;
					if(!isfinite(c_val))
						c_val=0.0;
					a_sum =a_sum+a_val;
					b_sum =b_sum+b_val;
					c_sum =c_sum+c_val;
					// calculating Reich 2009 F_st.
					/* The follwoing definintion has been taken from David Reich's article suppliment
					*Reconstructing Indian population history  NATURE| Vol 461|24 September 2009 doi: 10.1038/nature08365
					*Choose the variant allele, and suppose that the allele has population frequency p1, p2
					*in populations 1 and 2 respectively.
					*Set q_i=1-pi.
					*Then we can define Wright’s Fst
					*as 	Fst	=N/D
					*where
					*		N =p1(q2-q1)+p2(q1-q2)
					*		D =p1q2+q1p2 = N+p1q1+p2q2
					*This is a definition of Fst, a parameter measuring divergence at a given locus,
					*not a sample statistic.
					*Then Reich et.al further define the F_st estimator as
					* the following definition
					*
					*/
					a_i 		=2*(N0_mat.data[i][popIdx])+(N1_mat.data[i][popIdx]);
					a_j 		=2*(N0_mat.data[j][popIdx])+(N1_mat.data[j][popIdx]);
					hr_i		= a_i*1.0/(2.0*N_i) * ((2.0*N_i -a_i)/((2.0*N_i)-1.0));
					hr_j		= a_j*1.0/(2.0*N_j) * ((2.0*N_j -a_j)/((2.0*N_j)-1.0));
					Nst_r 		=  (((a_i/(2*N_i)) -(a_j/(2*N_j)))* ((a_i/(2*N_i)) -(a_j/(2*N_j)))) - ((hr_i/(2*N_i)) +(hr_j/(2*N_j)) );
					Dst_r		= Nst_r + hr_i + hr_j ;
					fst_r_tmp  +=isfinite((Nst_r/Dst_r)) ? (Nst_r/Dst_r) : 0.0;
					//calculation of Fst according as Wright 1951 :
					/**FST was introduced by Wright (1951) as a measure of
					 *correlation of gene frequencies und suggested the first and simplest
					 *estimator, For one allele at locus k:
					 *F_ST=s^2/(p_bar (1-p_bar))
					 *where s2 = 1/r* sum((p_i-p_bar)^2, r=1..r)
					 *is the observed variance of allele frequencies p_i among the sampled
					 *populations i (i=1, …, r) and p is the mean allele frequency over
					 *all populations. The estimate of FST for multiple loci is calculated
					 *by taking the mean across k loci.
					 *This estimator has a theoretical range between zero and one and is
					 *known to overestimate the level of genetic differentiation especially
					 *at low values.
					 */
					// This is defined in David reich theotirical way.
					//fst_w_tmp	=(P_i*((1-P_j) - (1.0-P_i)) +P_j*((1.0-P_i) - (1.0-P_j)) ) / (P_i*(1.0-P_j) + P_j*(1.0-P_i));
					/**The following formula is derived form solving above definition of Fst.
					 *This definition is given in
					 *Citation: Willing E-M, Dreyer C, van Oosterhout C (2012) Estimates of Genetic Differentiation Measured by FST Do Not Necessarily Require Large Sample Sizes
					 *When Using Many SNP Markers. PLoS ONE 7(8): e42649. doi:10.1371/journal.pone.0042649
					 *Eva-Maria Willing1*, Christine Dreyer1, Cock van Oosterhout2
					 *1 Department of Molecular Biology, Max Planck Institute for Developmental Biology, Tu¨ bingen, Germany, 2 School of Environmental Sciences, University of East Anglia,
					 *Norwich, United Kingdom
					*/
					// 0.2e1=2.0;
					fst_w_tmp	= -1.0*(pow(P_i - P_j, 0.2e1))/ (P_i + P_j)* (1.0/ (-0.2e1 + P_i + P_j));
					if(!isfinite(fst_w_tmp))
						fst_w_tmp=0.0;
					fst_w 		+=fst_w_tmp;
					// assigning all values to zero
					N_i =N_j =0;
					N_bar = S_sum =0.0;
					N_c = P_i = P_j = P_bar =S_square = H_i = H_j = H_bar = 0.0;
					a_val	= b_val =c_val=0.0;
					fst_w_tmp=0.0;


				}
				// After calculating a_sum , b_sum and c_sum
				//calculate f_st
				F_st.data[i][j] =a_sum /(a_sum+b_sum+c_sum);
				a_sum=b_sum=c_sum= 0.0;
				//fst reich
				Fst_Reich.data[i][j]=fst_r_tmp/(gsTotal);
				fst_r_tmp =0.0;
				//wright fst
				Fst_Wright.data[i][j]=fst_w/(gsTotal); //average weir fst.
				fst_w=0.0;
			}
			else
			{
				F_st.data[i][j] =F_st.data[j][i];
				Fst_Reich.data[i][j]=Fst_Reich.data[j][i];
				Fst_Wright.data[i][j]=Fst_Wright.data[j][i];

			}
			ofs <<F_st.data[i][j]<< "\t";
			ofs_reich <<Fst_Reich.data[i][j]<< "\t";
			ofs_wright <<Fst_Wright.data[i][j]<< "\t";

			//cout <<"F_st["<<i+1<<","<<j+1<<"]= "<<F_st.data[i][j]<< "  ";
			//Population pairwise FST values. A negative symbol (-) indicates that pairwise FST values were not significant.
			//A plus (+) sign indicates a significant FST value adjusted for the number of contrasts using a sequential Bonferroni test
		}
		ofs<<"\n";
		ofs_reich<<"\n";
		ofs_wright<<"\n";

		//cout << "\n";
	}
	//###############################################################################
	ofs.clear(); ofs.close();
	ofs_reich.clear(); ofs_reich.close();
	ofs_wright.clear(); ofs_wright.close();

	string msg = "\t**->Fst (Theta_hat)  using  the method (EQUATION 1) of \"Weir & Cockerham 1984, Evolution 38(6):1358-1370\" \n"
				"has been calculated and saved in the file: "+fstFileName+". \n";
	printLIN(msg);
	//Reich 2009 file
	msg = "\t**->Fst  using the definition mentioned in the supplement (Equation 12) of : \"Reconstructing Indian population history  NATURE| Vol 461|24 September 2009 doi: 10.1038/nature08365\"\n"
						 "has been calculated and saved in the file: "+fst_reich_fileName+". \n";

			printLIN(msg);

	msg = "\t**->Fst  using the Wright's definition  mentioned in:\n"
			"1. \"Nei, M. Analysis of gene diversity in subdivided populations. Proc. Nat. Acad. Sci. USA 70, 3321-3321 (1973)\"\n"
			"and in \n"
			"2. \"Kent E. Holsinger, Bruce S. Weir. Genetics in geographically structured populations: defining, estimating, and interpreting FST, 2009\" (equation 2 on  page 4).\n"
			"has been calculated and saved in the file: "+fst_wright_fileName+". \n";

	printLIN(msg);

}// end of calculate_fst

handel_reich_fst_calculation(pedVec, genVec);

}

// splitting snpwise and individual wise
void splitt_snps_iids( const vector<CBPED*>& pedVec, const vector<CBSNP*>&genVec, const Bpar* const pBpar){
	printLIN("->* Splitting the data: \n");
	if(pBpar->ssplitt){

	}


}

//----------------------------------------------//

int numblocks(vector<CBSNP*>snpm, int numsnps, double blocklen)
{
  int n = 0, i ;
  int chrom, xsize, lchrom, olds ;
  double fpos, dis, gpos ;
  	  CBSNP* cupt ;
  //SNP *cupt ;

  lchrom = -1 ; xsize = 0; olds = -1 ;

  	  fpos = 0 ;
  for (i=0; i<numsnps; i++)
  {
	// printf("zz %d %d %d %d\n", i, n, xsize, olds) ;
	  cupt = snpm[i] ;
	 // cupt -> tagnumber = -1 ;
	 // if (cupt -> ignore) continue ;
	  //if (cupt -> isfake) continue ;
   chrom = atoi(cupt -> nchr.c_str()) ;
   gpos = cupt ->cm_pos ; // must be in morgen
   dis = gpos - fpos ;
   if ((chrom != lchrom) || (dis>blocklen))
   {
		if (xsize>0) {
		cupt = snpm[olds] ;
		//   printf("n: %d  %d %d\n", n, cupt->chrom, (int) cupt -> physpos) ;
			++n ;
		}
     lchrom = chrom ;
     fpos = gpos ;
     olds = i ;
     xsize = 1 ;
     continue ;
   }
   ++xsize ;
  }
  return n+1 ;
}

void handel_reich_fst_calculation(vector<CBPED*>&pedVec, const vector<CBSNP*>&genVec)
{
	/*
	cout <<"\n---------------------\n calculting david reich\n";

	 int t1, t2 ;
	 int nblocks, xnblocks;
	 double y, sd;
	 int *blstart, *blsize;
	 double *xfst ;
	 int nrows 	 	=pedVec.size();
	 int ncols 	 	=genVec.size();
	 double blgsize	=0.5; //block size in Morgen
	 int numsnps =ncols; // check later no of good snps
	 // nblocks = numblocks(snpmarkers, numsnps, blgsize) ;
	 */
}

//
void CGENERAL::shapeit_type_commands(Bpar* const pBpar)
{
	//cout << "shapeit command";
	CIPED::lese_shapeit_sample_file(pedVec,pBpar->shapeit_sample_file);
	CISNP::lese_shapeit_haps_file(pedVec,genVec,pBpar->haps_file);


}










































