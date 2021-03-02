/*
 * cparse.cpp
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

#include "cparse.h"
#include "helper.h"
#include<cmath>
#include<stdio.h>
 bool temp_bool =false;
/*	* ALREADY USED command options: 
 	* #########  plink:   ##############
	* 	--ped ,	--map  			// to read plink ped and mapfile
	* 	--file					//combined form
	*	--dosage 				//to read impute dose file
	*	--fam 					//
	*	--reocdeA 				//raw file
	*	--recodeA-dose			// oformat
	*	--recodeAD 				//raw file with AD
	**
	**   #########  mach/minimac:   ##############
	*	--ped 	--dat			// to read mach ped and dat file
	*	--mfile					//ped and dat combined
	*	--mach-info 			// to read mach info file
	* 	--mach-hap	--mach-snp	// to read mach hapmap file
	* 	--mach-geno				// to read mach geno file 
	* 	--mach-mlinfo			// to read mach' mlinfo file
	*	--mach-mlgeno 			// to read mach's mlgeno file
	* 	--mach-mlprob			//to read mach's mlprob file
	* 	--minimac-prob			//minimac's prob file
	* 	--minimac-info 			//minimac's info file
	*	--rsq					//this command gives a thresh value for mach's rsq value
	**
	**   #########  impute:   ##############
	* 	--impute-hap    and
	* 	--impute-legend	 		// to read impute's hapmap format file
	*	--gens					// to read impute gens file
	* 	--thresh				// with this command one can give the thresh value of genotype probabilities
	* 	--info-thresh			// to give info threshold value for mach and impute
	*
	* 	 #########  SNPTEST:   ##############
	*	--covar 					//reads plink covariate file
	*	--covar-name				// with this command you can specify covariate names
	* 	--covar-type				// to specify the type of covariates. // this can be used for impute-snptest  covar file
	*	--sample	 --sex-col 	--pheno-col	// reading sample file of snptest
	*	--group-label   // to read groups or population label for individual while transforming into smartpca
	*
	**   #########  beagle:   ##############
	* 	--bgl 						// to read beagle's bgl file 
	* 	--bgl-gprobs					// to read beagle's gprobs file 
	*	--bgl-rsq  --rsq-thresh			//bg rsq values 
	*
	**   #########  bimbam:   ##############
	* 	--geno  --pos 				//bimbam 
	*	--wgd 					//bimbam 
	*	--wbg --pos  
	*	--maf-thresh			// to give the maf thresh for quality of SNPS
	* 	--force 					// to give a command for forcing  to implement
	* 	--force ref-allele = allele1, allele2 major-allele , minor-allele // this will determine the reference allele for dose
	* 	 e.g. dose =0*AA +1.AB+2.BB
	* 	 #########  rformat:   ##############
	* 	--rgeno 	rgenotype file
	*
	*
	**   #########  general:   ##############
	**   --exclude 			// to exclude SNP list
	**   --remove 			// to remove individual list
	*	--oformat			// to convert oformat files
	*						// current oformats :  plink ,  beagle , mach , impute, snptest, bimbam  plink -dosage
	*	--oformat 			//  probabel
	*	--oformat 			//genabel
	*	--freq 				//  to calcualte allele frequency and to write it down.
	*	--hardy				//hwe
	*	--missing			//callrate
	*	--out 				//  to give output file name
	* 	--help 				//  to use help command
	*	--snpinfo 			//  to read snpinfo file
	*	--pedinfo 			//  to read ped info file
	*	--write-snpinfo		// for writing snp list
	*	--write-pedinfo		//writing ped info
	*	--write-snplist 	// to write snp-list
	*	--write-pedlist 	//to write pedlist
	*	--force 			// to give a command for forcing  to implement
	*	--filter-snp		// to give  filtering options for maf, hwe and callrate
	*	--filter --iid		// to filter individual according as different situtaiton .
	*	--new-start  --new-end // to give new command on fcgene
	*	--merge 			//to merge two files
	*	--uncompress	=false		// this command is used to  write uncompress data
	*	--ssplit	1-10,11-20; 	//this command will split the data snpwise
	*	--isplit 	1-10,11-20 		//this command will split the data indiviudalwise
	*	--bpsplit	1-10,1-20
	*	--transpose				// this is now implimented only for writing genotype data in r or r-dose
	*
	*
	*


	
	
*/
 

//this set option is the main parsing command funciton
void setCoptions(Bpar * const pBpar)
{
	// set command options.
	//for no option
	//pBpar->find_single_opt("--map")
	/*
		// for --file option
	//----------------------------------------------------------------------//
		if(pBpar->search_string("--file"))
	{

		// for mach
		if(pBpar->file_root=="mach")
				{
				pBpar->ped_file 		=pBpar->file_root+".ped";
				pBpar->read_ped 		=true;
				pBpar->read_dat			=true;
				pBpar->dat_file 		=pBpar->file_root+".dat";
				pBpar->code_readType	=pBpar->file_root;
				//check if file exits
				cfe(pBpar->ped_file,pBpar->read_ped);
				cfe(pBpar->dat_file,pBpar->read_dat);

				}
		//for impute
		else if(pBpar->file_root=="impute")
		{
			pBpar->code_readType		=pBpar->file_root;
			pBpar->read_gens			=true;
			pBpar->gens_file			=pBpar->file_root+".gens";
			cfe(pBpar->gens_file,pBpar->read_gens);
			if(pBpar->search_string("--thresh") )
			{
				pBpar->find_double_opts(pBpar->thresh, "--thresh");
				pBpar->default_thresh=false;
				if(!pBpar->default_thresh)
					pBpar->threshold=atof(pBpar->thresh.c_str());
			}

		}
		else if(pBpar->file_root=="beagle")
		{
			pBpar->code_readType		=pBpar->file_root;
			pBpar->read_bgl			=true;
			pBpar->bgl_file			=pBpar->file_root+".bgl";
			cfe(pBpar->bgl_file,pBpar->read_bgl);

		}
		else if(pBpar->file_root=="bimbam")
		{
			cout << "\nThis function is under construction. If you want to use this option please inform us so that  we can  update it for you.\n";
			exit(1);
		}

	}
	// end for --file option
	*/
	//-----------------------------------------------------------------------------//
	// The following is a mix of Plink and mach because they both of ped file
	//-----------------------------------------------------------------------------//
	parse_plink_commands(pBpar);
	// for condition : --map test.map#
	/**------------------------------------------------------------------------
		*mach options. There may be also plink options in this section
		*------------------------------------------------------------------------
	*/
	//--dat option
	parse_mach_minimach_commands(pBpar);
		//
	if(pBpar->search_string("--ped"))
	{
		pBpar->read_ped = true;
		pBpar->find_double_opts(pBpar->ped_file, "--ped");
		cfe(pBpar->ped_file,pBpar->read_ped);
		if(pBpar->read_map)
		{
			pBpar->find_single_opt("--map");
			pBpar->find_double_opts(pBpar->map_file, "--map");
			cfe(pBpar->map_file,pBpar->read_map);
			pBpar->code_readType="plink";


		}
		else if(pBpar->read_dat)
		{
			pBpar->find_double_opts(pBpar->dat_file, "--dat");
			cfe(pBpar->dat_file,pBpar->read_dat);
			if(pBpar->define_readType)
				pBpar->code_readType="mach";
			 pBpar->define_readType=false;

		}
		else
		{
			string msg= "Please give map file with option \"--map\", to read plink format data or dat file with option \"--dat\", to read mach format data.\n";
			error(msg);
		}
	}
	//Impute
	parse_impute_snptest_commands(pBpar);
	parse_shapeit_commands(pBpar);
	parse_beagle_commands(pBpar);
	parse_bimbam_commands(pBpar);
	parse_rformat_commands(pBpar);

	// for the command --out:
// pasting iid and fid
 parse_option_iid(pBpar);

//########################################
 //general commands
 //########################################
 parse_general_commands(pBpar);
 parse_eigensoft_commands(pBpar);
 parse_summary_statistics_commands(pBpar);
// for reading plink fomratted covariate files:
// for making snptest *.sample file

	 // plink covariate file
	 if(pBpar->search_string("--covar"))
	 {
		 if(pBpar->is_general_command)
		 {
			pBpar->find_single_opt("--covar");
			if(pBpar->search_string("--oformat") && pBpar->search_string("snptest"))
				pBpar->find_single_opt("--covar-name"); // first check if th
			pBpar->find_double_opts(pBpar->covariate_file,"--covar");
			 pBpar->read_covariate_file	=true;
			 cfe(pBpar->covariate_file,pBpar->read_covariate_file);
			 if(pBpar->read_covariate_file)
			 {
				 // check if  only certain covariates are supposed to use
				 if(pBpar->search_string("--covar-name"))
				 {
					pBpar->find_single_opt("--covar-name"); // first check if there any value after a command  //1
					pBpar->find_double_opts(pBpar->covar_names,"--covar-name"); //2
					 pBpar->given_covar_names	=true;			//3
					 pBpar->vec_covar_names.resize(0);
					 pBpar->vec_covar_names  =Argcv::parse_subArg(pBpar->covar_names, ",");
					 //for(unsigned int i =0;i<pBpar->vec_covar_names.size();++i)
					//	cout << pBpar->vec_covar_names[i]<< " ";
					//	cout<<endl;	
					// if  --covar -name command is given, then there must be  --covar-type command
					//here must be all things 
					//#########################################start
					ifstream file; file.close();file.clear();
					string fileName	=pBpar->covariate_file;
					string line		="";
					string buffer		="";
					size_t	found;
					string fCovar		="";
					string sCovar		="";
					string tempString	=""	;
					vector<string>temp_header;
					//out << fileName <<endl;	//debug 
					file.open(fileName.c_str(),ios::in);
					if(!file)
					{
								file.setstate(ios_base::failbit);
								file.clear();
								error("Your map file "+fileName +" either does not exit or could not be opend!\n");
		
					}
							 //read the fist rwo of the file which should contain
							 // covariate names
					do
					{
								getline(file,line,'\n');
								if(line[0]=='#')
									continue;
								if(file.eofbit) break;
					 }while(line.empty());
							 file.clear();file.close();
							 stringstream line_parser(line);
					if(line_parser.good())
					{
						while(line_parser>>buffer)
							temp_header.push_back(buffer);
						//for(unsigned int i=0;i<temp_header.size();++i) //debug 
						//	cout<< temp_header[i]<< " "; //debug 
						//	cout <<endl; //debug 
						// first column must be labelled with "FID"
						if(temp_header[0]!="FID")
						{
							string pr_msg="First Column of covariate file "+fileName+ " should be labelled with \"FID\". Please correct the header of first column.\n";
							error(pr_msg);
						}
						// second column should contain indiviudal id 
						if(temp_header[1]!="IID")
							{	
								string pr_msg="Second Column of covariate file "+fileName+ " should be labelled with \"IID\". Please correct the header of second column.\n";
								error(pr_msg);
							}
						//cout <<pBpar->vec_covar_names.size()<<endl;
					}
					else
					{
								string pr_msg = "Problem in parsing the first line of the file. Please contact the software developer.\n";
								error(pr_msg);
					}
					
				//##########################################
					vector<string> temp_vec1;
					vector<string> temp_vec2;
					vector<int>idx_B;
					vector<int>idx_notB;
					if(pBpar->search_string("--oformat") && pBpar->search_string("snptest"))
					{
						pBpar->find_single_opt("--covar-type"); // first check if there any value after a command
						pBpar->find_double_opts(pBpar->covar_types,"--covar-type"); // first check if there any value after a command
						pBpar->given_covar_types=true;
						if(pBpar->given_covar_types)
						{
						//--
							pBpar->vec_covar_types.resize(0);
							pBpar->vec_covar_types  =Argcv::parse_subArg(pBpar->covar_types, ",");
							if(pBpar->vec_covar_names.size()!=pBpar->vec_covar_types.size())
							{
								string msg = "The number of covariates names given with option \"--covar-name\" and the number of covariate types given with option \"--covar-type\" must be same. Please correct or remove unnecessary commands.\n";	
								error(msg);
							}
							// check if vec_covar_type contaions only C, D , P AND B or not 
							 vector<string>total_covar_types;
							 total_covar_types.push_back("C");total_covar_types.push_back("D");
							 total_covar_types.push_back("B");total_covar_types.push_back("P");
							 int ncol		=0;
							 bool tfValue	=0;
							for(unsigned int i=0;i<pBpar->vec_covar_types.size();++i)
							 {
								// coviraite types should be:  C,D,P or B.
								 //  checking if covariate types are given correctly. 
								 // writing 
								tfValue=compfun_with_given_vector(total_covar_types, pBpar->vec_covar_types[i] );
									if(!tfValue)
								{
									printLIN("Covariate type \" "+pBpar->vec_covar_types[i] +"\" is not acceptable.\n");
									string msg = "Covariate types given with option \"--covar-type\"  can be only 'C', 'D', 'B' or 'P', seperated with commas.\n"	;
									error(msg);
								}
							
							}
					// here starts changing  the order of covariate names . IN SNPTEST format files. phenotypes should come after covariates
					// That means P and B types must be at the end  columns So let's put like:   CCDCPPBB 
							//first putting p types phenotypes at the end of the file. 
							//##
							vector<int>idx_P		=pos_in_vec(pBpar->vec_covar_types,"P");
							//vector<string>idx_notP	=pos_notin_vec(pBpar->vec_covar_types,"P");
							//##
							idx_B	=pos_in_vec(pBpar->vec_covar_types,"B");
							
							//idx_notB=pos_notin_vec(pBpar->vec_covar_types,"B");
							for(unsigned int i=0;i<idx_B.size();++i)
							{
								//cout << i << " "; //debug 
								idx_P.push_back(idx_B[i]);
							}
							vector<int>tmp = idx_P;
							sort(tmp.begin(),tmp.end());
							bool _tmp_b_bool=true;
							for(unsigned int i=0;i<pBpar->vec_covar_types.size();++i)
							{
								_tmp_b_bool=binary_search(tmp.begin(),tmp.end(),(int)i);
								if(!_tmp_b_bool)
								{
									idx_notB.push_back(i);
									//cout << "no pheno: "<<i << " "; // debug 
								}
								//else cout<< "pheno: "<<i << " " ;	
							}
							 //cout <<endl; //debug 
							 //Now first put the indices which are in idx_notB and then which are in idx_P (!!not idx_B)
							for(unsigned int i=0;i<idx_notB.size();++i)
							{
								ncol=idx_notB[i];
								temp_vec1.push_back(pBpar->vec_covar_names[ncol]);
								temp_vec2.push_back(pBpar->vec_covar_types[ncol]);
							}
							for(unsigned int i=0;i<idx_P.size();++i)
							{
								ncol=idx_P[i];
								temp_vec1.push_back(pBpar->vec_covar_names[ncol]);
								temp_vec2.push_back(pBpar->vec_covar_types[ncol]);
							}
							if(temp_vec2.size()!=pBpar->vec_covar_types.size()|| temp_vec1.size()!=pBpar->vec_covar_names.size() )
							{
								string msg ="Problem in parsing  --coar-name and covar-type. Please contact the software developer!\n";
								error(msg);
							
							}
							pBpar->vec_covar_names=temp_vec1;
							pBpar->vec_covar_types=temp_vec2;
							temp_vec1.clear();
							temp_vec2.clear();
							idx_B.clear();
							idx_P.clear();
							tmp.clear();
							idx_notB.clear();			
						//--
						}
						
						
						
					}
					//---
					for(unsigned int i=0;i<pBpar->vec_covar_names.size();++i)
					{
						//cout << " index: " << i<<endl; //debug 
						tempString=pBpar->vec_covar_names[i];
						//cout << tempString<<" ";
						found=tempString.find('-');
						if(found!=string::npos)
						{
							fCovar=tempString.substr(0,found);
							sCovar=tempString.substr(found+1);
							//sCovar and fCovar must not be as the other elements vec_covar_names.
							//cout << fCovar <<" " <<sCovar<<endl; //debug 
							idx_B		=pos_in_vec(temp_header, fCovar);
							idx_notB	=pos_in_vec(temp_header, sCovar);
							//cout << idx_B.size()<< idx_notB.size() <<endl; 
							if(idx_B.size()!=1 ||idx_notB.size()!=1)
							{
								printLIN("Covariate name \" "+tempString +"\" is not acceptable.\n");
								string msg = "Covariate "+ sCovar+ " or "+fCovar +" or both  does not match with the covariate names given in "+fileName+". \n"	;
								error(msg);									
							}
							int lowerLimit=idx_B[0];
							int upperLimit=idx_notB[0];
							idx_B.clear();
							idx_notB.clear();
							for(int j=lowerLimit;j<(upperLimit+1);++j)
							{
								//idx_B		=pos_in_vec(pBpar->vec_covar_names, temp_header[j]);
								idx_B		=pos_in_vec(temp_vec1, temp_header[j]);
								if(idx_B.size()>0 )
								{
									printLIN("The range of Covariate name \""+tempString +"\" is not acceptable.\n");
									string msg = "One ore more of the covariates which lie within the range are already given with command option \" --covar-name\".  Please mention a covriate only once.\n";
									error(msg);									
								}
								idx_B.clear();
								temp_vec1.push_back(temp_header[j]);
								if(pBpar->given_covar_types)
									temp_vec2.push_back(pBpar->vec_covar_types[i]);
													
							}	
												
											
						}
						else
						{
							bool tfValue=compfun_with_given_vector(temp_header, tempString);
							//cout << boolalpha << tfValue<<endl;
							if(!tfValue)
							{
								
								printLIN("Covariate name \" "+tempString +"\" is not acceptable.\n");
								string msg = "Covariate names given with option \"--covar-name\"  must match with the covariate names given in "+fileName+". \n"	;
								error(msg);
							}
							idx_B		=pos_in_vec(temp_vec1,tempString);
							if(idx_B.size()>0 )
							{
								printLIN("Covariate name \""+tempString +"\" is not acceptable.\n");
								string msg = "Covariate "+ tempString+ " has  mentioned  more than one times using command option \" --covar-name\".  Please mention a covriate only once.\n";
								error(msg);									
							}
							temp_vec1.push_back(tempString);
							if(pBpar->given_covar_types)
								temp_vec2.push_back(pBpar->vec_covar_types[i]);
							idx_B.clear();
													
						}
										
					}
										
					pBpar->vec_covar_names=temp_vec1;
					if(pBpar->given_covar_types)
						pBpar->vec_covar_types=temp_vec2;
					temp_vec1.clear();
					temp_vec2.clear();
					idx_B.clear();
					idx_notB.clear();
					//--
					//---
					//##########################################
						
				}
				// here must be for the covariate names. 
			 }
		 
		 }else
		 {
			 string msg= " Plink formatted covariate file  can be read only with the plink ped and map file. Please feed these two files with options \"--ped\" and \"--map\".\n";
			 error(msg);

		 }

	 }
	
	
//---------------------------------------------------------------------------------------------//
					// END of  parsed functions
//---------------------------------------------------------------------------------------------//
// Now the only thing is to check whether all input and output files have valid file name and file path.
	 //first check if input file name and output file Names are readable or if they exsits or not.
	 // whether the filepath are given correctly.
	 check_input_output_file_names_n_paths(pBpar);




	// here ends the function setCoptions:
}
// parse general commands
void parse_plink_commands(Bpar * const pBpar ){

	if(pBpar->search_string("--file"))
	{
		if(pBpar->search_string("--ped")|| pBpar->search_string("--map"))
		{
			error( "Use either the --file option or --ped and --map optios.\n ");
		}
		pBpar->find_double_opts(pBpar->file_root,"--file");
		//cout <<"Test\n";
		if(pBpar->file_root!="")
		{
			pBpar->ped_file 		=  pBpar->file_root+".ped";
			pBpar->read_ped 		= true;
			pBpar->read_map		 	=true;
			pBpar->map_file 		= pBpar->file_root+".map";
			pBpar->code_readType	="plink";
			// check if file exits
			cfe(pBpar->ped_file,pBpar->read_ped);
			cfe(pBpar->map_file, pBpar->read_map);
		}

	}
	// bfile
	if(pBpar->search_string("--bfile"))
	{
		if(pBpar->search_string("--bed")|| pBpar->search_string("--bim")|| pBpar->search_string("--fam"))
		{
			error( "Use either the --bfile option or --bed and --bim and --fam optios.\n ");
		}
		pBpar->find_double_opts(pBpar->file_root,"--bfile");
		//cout <<"Test\n";
		if(pBpar->file_root!="")
		{
			pBpar->bed_file 		=  pBpar->file_root+".bed";
			pBpar->bim_file 		= pBpar->file_root+".bim";
			pBpar->fam_file 		= pBpar->file_root+".fam";
			pBpar->read_bed 		= true;
			pBpar->read_bim		 	=true;
			pBpar->read_fam		 	=true;
			pBpar->code_readType	="bplink";
			// check if file exits
			cfe(pBpar->bed_file,pBpar->read_bed);
			cfe(pBpar->bim_file, pBpar->read_bim);
			cfe(pBpar->fam_file, pBpar->read_fam);
		}

	}
	//bfile sepearte
	if(pBpar->search_string("--bed"))
	{
		 pBpar->find_double_opts(pBpar->bed_file,"--bed");
		 pBpar->find_double_opts(pBpar->bim_file,"--bim");
		 pBpar->find_double_opts(pBpar->fam_file,"--fam");
		 pBpar->read_bed	=true;
		 pBpar->read_bim	=true;
		 pBpar->read_fam	=true;
		 pBpar->code_readType	="bplink";
		 // check if file exits
		 cfe(pBpar->bed_file,pBpar->read_bed);
		 cfe(pBpar->bim_file, pBpar->read_bim);
		 cfe(pBpar->fam_file, pBpar->read_fam);

	}

	//map
	if(pBpar->search_string("--map"))
	{
		pBpar->read_map = true;
		//cout <<"test32\n";
		//pBpar->map_file = pBpar->find_double_opts("--map");
		//pBpar->find_double_opts(pBpar->map_file, "--map");
		//pBpar->code_readType="plink";
	}
	if(pBpar->search_string("--fam"))
		pBpar->read_fam=true;

	// plink dosage file
	if(pBpar->search_string("--dosage"))
	{
		if(pBpar->read_fam)
		{
			pBpar->read_plink_dosage =true;
			pBpar->find_double_opts(pBpar->plink_dosage_file,"--dosage");
			if(pBpar->define_readType)
			{
				pBpar->code_readType="plink-dosage";
				pBpar->define_readType =false;
			}
			else{
				 error("Type of input data is already defined. Therefore command option \"--dosage\" is not possible here\n");
			}

			if(pBpar->read_map)
				pBpar->find_double_opts(pBpar->map_file,"--map");
			// may be read_map
			//may be read_fam
		}
		else{
			error("A necessary argument given with \"--fam\" is missing. \n");
		}

	}
	// plink rawA file
	if(pBpar->search_string("--recodeA"))
	{
		// At the moment this is like that you should have snpinfo file to
		// read recode files
		if(pBpar->read_map&& pBpar->search_string("--snpinfo"))
		{
			pBpar->read_plink_rawA=true;
			pBpar->find_double_opts(pBpar->plink_rawA_file,"--recodeA");
			if(pBpar->define_readType)
				pBpar->code_readType="plink-rawA";
			pBpar->define_readType=false;
			pBpar->find_double_opts(pBpar->map_file,"--map");
		}
		else
		{
			error("argument with --map and/or --snpinfo is missing \n");
		}
	}
	// plink rawA file
	if(pBpar->search_string("--recodeAD"))
	{
		if(pBpar->read_map&&  pBpar->search_string("--snpinfo"))
		{
			pBpar->read_plink_rawAD=true;
			pBpar->find_double_opts(pBpar->plink_rawAD_file,"--recodeAD");
			if(pBpar->define_readType)
				pBpar->code_readType="plink-rawAD";
			pBpar->define_readType =false;
			pBpar->find_double_opts(pBpar->map_file,"--map");
		}
		else
		{
			error("argument with --map and/or --snpinfo is missing \n");
		}
	}

	// fam file
	if(pBpar->read_fam)
	{
		if(pBpar->read_plink_dosage)
		{
			pBpar->read_fam=true;
			pBpar->find_double_opts(pBpar->plink_fam_file,"--fam");
		}
	}



} //end of plink commands
//end of parse plink-commands
void parse_mach_minimach_commands(Bpar * const pBpar )
{
	// for mach
	if(pBpar->search_string("--mfile"))
	{
		if(pBpar->search_string("--ped")|| pBpar->search_string("--dat"))
		{
			error( "Use either the --mfile option or --ped and --dat optios.\n ");
		}
		pBpar->find_double_opts(pBpar->file_root,"--mfile");
			//cout <<"Test\n";
		if(pBpar->file_root!="")
		{
			pBpar->ped_file 		=pBpar->file_root+".ped";
			pBpar->read_ped 		=true;
			pBpar->read_dat			=true;
			pBpar->dat_file 		=pBpar->file_root+".dat";
			pBpar->code_readType	="mach";
			//check if file exits
			cfe(pBpar->ped_file,pBpar->read_ped);
			cfe(pBpar->dat_file,pBpar->read_dat);

		}
			//for impute
	}
	if(pBpar->search_string("--dat"))
		{
			pBpar->read_dat = true;
			//pBpar->find_double_opts(pBpar->dat_file, "--dat");
			//pBpar->code_readType="mach";
		}
		// for mach.info or mach.mlinfo option
		if(pBpar->search_string("--mach-info"))
		{
			pBpar->read_mach_info = true;
		}
		//for mach mlinfo
		if(pBpar->search_string("--mach-mlinfo"))
		{
			pBpar->read_mlinfo = true;
		}

		//mach hapmap file
		if(pBpar->search_string("--mach-hap"))
		{
			pBpar->read_mach_hapmap =true;
			pBpar->find_single_opt("--mach-snp"); // this will ask to give --snpinfo command
			pBpar->find_double_opts(pBpar->mach_hapmap_file,"--mach-hap");
			cfe(pBpar->mach_hapmap_file,pBpar->read_mach_hapmap);
			pBpar->code_readType="mach_hapmap";
			// snps
			if(pBpar->search_string("--mach-snp"))
			{
				pBpar->read_mach_hap_snps=true;
				pBpar->find_double_opts(pBpar->mach_hap_snp_file,"--mach-snp");
				cfe(pBpar->mach_hap_snp_file,pBpar->read_mach_hap_snps);

			}
		}
//
		// mach hapmap snps file
		// mach.geno file
		if(pBpar->search_string("--mach-geno"))
		{
			pBpar->read_mach_geno = true;
			pBpar->find_double_opts(pBpar->mach_geno_file, "--mach-geno");
			cfe(pBpar->mach_geno_file,pBpar->read_mach_geno);
			if(pBpar->read_mach_info)
			 {
				pBpar->find_double_opts(pBpar->mach_info_file, "--mach-info");
				cfe(pBpar->mach_info_file,pBpar->read_mach_info);
				pBpar->code_readType="mach.geno";

			 }
			else
			{
				string msg= "Please give  mach.info file with option \"--mach-info\", to read mach.geno format data.\n";
				error(msg);
			}


		}
			//pBpar->find_double_opts(pBpar->dat_file, "--dat");
			//pBpar->code_readType="mach";

	 // for mach mlgeno
		if(pBpar->search_string("--mach-mlgeno"))
		{
			pBpar->read_mlgeno = true;
			pBpar->find_double_opts(pBpar->mlgeno_file, "--mach-mlgeno");
			cfe(pBpar->mlgeno_file,pBpar->read_mlgeno);
			if(pBpar->read_mlinfo)
			 {
				pBpar->find_double_opts(pBpar->mlinfo_file, "--mach-mlinfo");
				cfe(pBpar->mlinfo_file,pBpar->read_mlinfo);
				pBpar->code_readType="mach.mlgeno";
			 }
			else
			{
				string msg= "Please give  mach.mlinfo file with option \"--mach-mlinfo\", to read mach.mlgeno format data.\n";
				error(msg);
			}
			//pBpar->find_double_opts(pBpar->dat_file, "--dat");
			//pBpar->code_readType="mach";
		}
		//--------------------------------------------//
		//for mach mlprob.
		if(pBpar->search_string("--mach-mlprob"))
			{
				pBpar->read_mlprob = true;
				pBpar->find_double_opts(pBpar->mlprob_file, "--mach-mlprob");
				cfe(pBpar->mlprob_file,pBpar->read_mlprob);
				if(pBpar->read_mlinfo)
				{
					pBpar->find_double_opts(pBpar->mlinfo_file, "--mach-mlinfo");
					cfe(pBpar->mlinfo_file,pBpar->read_mlinfo);
					pBpar->code_readType="mach.mlprob";
				}
				else
				{
					string msg= "Please give  mach.mlinfo file with option \"--mach-mlinfo\", to read mach.mlprob format data.\n";
					error(msg);
				}
				//pBpar->find_double_opts(pBpar->dat_file, "--dat");
				//pBpar->code_readType="mach";
			}
		// writing for mach.mlgeno , mach.geno and mach mlprob rsq and maf command
	//-----------------------------------------------------//
		temp_bool = pBpar->search_string("--mach-mlprob") || pBpar->search_string("--mach-mlgeno") ||pBpar->search_string("--mach-geno")||pBpar->search_string("--minimac-prob");
		if(temp_bool)
		{
			//for rsq
				if(pBpar->search_string("--rsq") )
				{
					pBpar->find_double_opts(pBpar->temp_string, "--rsq");
					pBpar->default_rsq=false;
					if(!pBpar->default_rsq)
					{
						pBpar->rsq_value		=atof(pBpar->temp_string.c_str());
						pBpar->temp_string 	="";
					}



				}
			// for maf
			if(pBpar->search_string("--maf-thresh") )
			{
				pBpar->find_double_opts(pBpar->temp_string, "--maf-thresh");
				pBpar->default_maf=false;
				if(!pBpar->default_maf)
				{
					pBpar->maf_value		=atof(pBpar->temp_string.c_str());
					pBpar->temp_string 	="";
				}
			}
			temp_bool=false;

		}
	//-----------------------------------------------------//
		//---------------------------------------------------------------------------------------------//
			//minimac
		//---------------------------------------------------------------------------------------------//
			//for mach mlprob.
			if(pBpar->search_string("--minimac-info"))
				{
					pBpar->read_minimac_infoFile = true;
				}
			if(pBpar->search_string("--minimac-prob"))
			{
				pBpar->read_minimac_probFile = true;
				pBpar->find_double_opts(pBpar->minimac_probFile, "--minimac-prob");
				cfe(pBpar->minimac_probFile,pBpar->read_minimac_probFile);
				if(pBpar->read_minimac_infoFile)
				{
					pBpar->find_double_opts(pBpar->minimac_infoFile, "--minimac-info");
					cfe(pBpar->minimac_infoFile,pBpar->read_minimac_infoFile);
					pBpar->code_readType="minimac.prob";
				}
				else
				{
					string msg= "Please give  minimac.info file with option \"--minimac-info\", to read minimac.prob format data.\n";
					error(msg);
				}
				//pBpar->find_double_opts(pBpar->dat_file, "--dat");
				//pBpar->code_readType="mach";
			}



}//end of function mach_minimach
//
void parse_impute_snptest_commands(Bpar * const pBpar ){
	//---------------------------------------------------------------------------------------------//
		//IMPUTE
	//---------------------------------------------------------------------------------------------//
		//impute hapmap reading
		if(pBpar->search_string("--impute-hap"))
		{
			pBpar->read_impute_hapmap =true;
			pBpar->find_single_opt("--impute-legend"); // this will ask to give --snpinfo command
			pBpar->find_double_opts(pBpar->impute_hapmap_file,"--impute-hap");
			cfe(pBpar->impute_hapmap_file,pBpar->read_impute_hapmap);
			if(pBpar->search_string("--impute-legend"))
			{
				pBpar->read_impute_hap_snps=true;
				pBpar->find_double_opts(pBpar->impute_hap_snp_file,"--impute-legend");
				cfe(pBpar->impute_hap_snp_file,pBpar->read_impute_hap_snps);
			}

			pBpar->code_readType="impute_hapmap";

		}

		// --gens option
		if(pBpar->search_string("--gens"))
		{
			pBpar->find_double_opts(pBpar->gens_file,"--gens");
			pBpar->read_gens =true;
			cfe(pBpar->gens_file,pBpar->read_gens);
			if(pBpar->read_gens)
			{
				if(pBpar->define_readType)
					pBpar->code_readType 	="impute";
				pBpar->define_readType=false;
				if(pBpar->search_string("--thresh") )
				{
					pBpar->find_double_opts(pBpar->thresh, "--thresh");
					pBpar->default_thresh=false;
					if(!pBpar->default_thresh)
					pBpar->threshold=atof(pBpar->thresh.c_str());
					pBpar->thresh="";
				}
				// for impute_info file
				if(pBpar->search_string("--info"))
				{
					pBpar->find_double_opts(pBpar->impute_info_file,"--info");
					pBpar->read_impute_info_file =true;
					cfe(pBpar->impute_info_file,pBpar->read_impute_info_file);
					if(pBpar->read_impute_info_file)
					{//
						if(pBpar->search_string("--info-thresh") )
						{
							pBpar->find_double_opts(pBpar->thresh, "--info-thresh");
							pBpar->default_impute_info=false;
							if(!pBpar->default_impute_info)
							pBpar->info_thresh=atof(pBpar->thresh.c_str());
							pBpar->thresh="";
						}
						else{
							string msg= "Please  give the cutoffs value of impute_info score  with option \"--info-thresh\", so that "+pBpar->impute_info_file+ " file can be read  and implemented for further analysis.\n";
							error(msg);

						}
						//for maf thresh
						// for maf
						if(pBpar->search_string("--maf-thresh") )
						{
							pBpar->find_double_opts(pBpar->temp_string, "--maf-thresh");
							pBpar->default_maf=false;
							if(!pBpar->default_maf)
							{
								pBpar->maf_value		=atof(pBpar->temp_string.c_str());
								pBpar->temp_string 	="";
							}
						}
					}


				}
				// for snp-test sample file
				if(pBpar->search_string("--sample"))
				{
					pBpar->find_double_opts(pBpar->snptest_sample_file,"--sample");
					pBpar->read_snptest_sample=true;
					cfe(pBpar->snptest_sample_file,pBpar->read_snptest_sample);
					if(pBpar->read_snptest_sample)
					{
							// for sex info column
						if(pBpar->search_string("--sex-col") ) // if sex column is specified.
						{
							pBpar->find_double_opts(pBpar->temp_string, "--sex-col");
							pBpar->	read_sex_col		=true;
							if(pBpar->read_sex_col)
								pBpar->  sex_col			=atoi(pBpar->temp_string.c_str());

						}
						// for phenotype (case control) column
						if(pBpar->search_string("--pheno-col") )
						{
							pBpar->find_double_opts(pBpar->temp_string, "--pheno-col");
							pBpar->read_pheno_col		=true;
							if(pBpar->read_pheno_col)
								pBpar-> pheno_col			=atoi(pBpar->temp_string.c_str());

						}
						// cout << "col: "<< pBpar->sex_col << " " << pBpar->pheno_col << endl; //test
						 //cout << (pBpar->sex_col==pBpar->pheno_col)<< endl; //test
						 const int no_allow_cols[] ={1, 2,3};
						 const vector<int>	now_allow_ints(no_allow_cols,no_allow_cols+3);
						 const bool tfValue = binary_search(now_allow_ints.begin(),now_allow_ints.end(),pBpar->sex_col)|| binary_search(now_allow_ints.begin(),now_allow_ints.end(),pBpar->pheno_col);
						if(tfValue)
						{
							string msg= " As described in SNPTEST file format information, first three columns of"+\
							pBpar->snptest_sample_file+ " file should represent pedigree id, sample id and portion of missing. How ever you gave  --sex-col "+\
							change_int_into_string(pBpar->sex_col)+" and --pheno-col "+change_int_into_string(pBpar->pheno_col)+\
							". Please specify these two information correctly.\n";
							error(msg);
						}
						if((pBpar->sex_col==pBpar->pheno_col))
						{
							string msg= "Sex information and phenotype information can not be represented by the same column: "+\
							change_int_into_string(pBpar->sex_col)+". Please give the correct number of columns of sex information and phenotype(case control) information with option \"--sex-col\" and  \"--pheno-col\"  so that "+\
							pBpar->snptest_sample_file+ " file can be read  and implemented for further analysis.\n";
							error(msg);
						}

					}

				}// end of if loop
			}

		}
	//---------------------------------------------------------------------------------------------//
				//SNPTEST
	//---------------------------------------------------------------------------------------------//

}
//

void parse_beagle_commands(Bpar * const pBpar ){


	//---------------------------------------------------------------------------------------------//
			//BEAGLE
	//---------------------------------------------------------------------------------------------//
		// for beagle --bgl option
			if(pBpar->search_string("--bgl"))
			{
				pBpar->find_double_opts(pBpar->bgl_file,"--bgl");
				pBpar->read_bgl			=true;
				pBpar->code_readType 	="beagle";
				cfe(pBpar->bgl_file,pBpar->read_bgl);

			}

		//  for beagle --bgl-grpobs options
				if(pBpar->search_string("--bgl-gprobs"))
				{

					pBpar->find_double_opts(pBpar->bgl_gprobs_file,"--bgl-gprobs");
					pBpar->read_bgl_gprobs			=true;
					pBpar->code_readType 	="beagle_gprobs";
					cfe(pBpar->bgl_gprobs_file,pBpar->read_bgl_gprobs);
					 if(pBpar->search_string("--no-header"))
						 {
							 pBpar->given_gprobs_header=false;
							 pBpar->find_single_opt("--no-header");
						 }


				}//end of if(pBpar->search_string(--bgl-gropus))

				//the following part is for rsq file.

				if(pBpar->read_bgl || pBpar->read_bgl_gprobs)
				{
					//----------------------------------------//
					//the following part is for rsq file.
					if(pBpar->search_string("--bgl-rsq"))
					{
						pBpar->find_double_opts(pBpar->bgl_rsq_file,"--bgl-rsq");
						pBpar->read_bgl_rsq=true;
						cfe(pBpar->bgl_rsq_file, pBpar->read_bgl_rsq);
						if(pBpar->search_string("--rsq-thresh"))
						{
							pBpar->find_double_opts(pBpar->temp_string,"--rsq-thresh");
							pBpar->bgl_rsq_thresh=atof(pBpar->temp_string.c_str());
							if(pBpar->bgl_rsq_thresh>1.0 &&pBpar->bgl_rsq_thresh<0.0)
							{
								string msg= "The value of beagle RSQ cutoff must lie between 0 and 1. please correct the value after --rsq-thresh .\n";
								error(msg);
							}
							pBpar->temp_string="";

						}
						else{
							string msg= "Please  give the cutoffs value of beagle type Rsq  score  with option \"--rsq-thresh\", so that "+pBpar->bgl_rsq_file+ " file can be read  and implemented for further analysis.\n";
							error(msg);

						}

					}
					//--------------------------------------//

				}

}
//
void parse_bimbam_commands(Bpar * const pBpar ){
	//---------------------------------------------------------------------------------------------//
		//bimbam
	//---------------------------------------------------------------------------------------------//
		//for bimbam --bimbam option
		if(pBpar->search_string("--geno"))
		{
			pBpar->find_double_opts(pBpar->bimbam_file,"--geno");
			pBpar->read_bimbam=true;
			pBpar->code_readType		="bimbam";
			cfe(pBpar->bimbam_file,pBpar->read_bimbam);
			if(pBpar->search_string("--pos"))
			{
				// --pos command  reads a different file with  --wbg or --wgd option for bimbam is given.
				// see below
				pBpar->find_double_opts(pBpar->bimbam_pos_file,"--pos");
				pBpar->read_bimbam_pos	=true;
				cfe(pBpar->bimbam_pos_file,pBpar->read_bimbam_pos);
			}else
			{
				string msg= pBpar->bimbam_file +"  can only be read together with \"--pos\" command. Please give SNP information  file (used in BIMBAM ) with \"--pos\" option.\n";
				error(msg);

			}

		}
		//for bimbam --bimbam-gprob option
		if (pBpar->search_string("--testt")){
			cout <<"testting"<<endl;
		}
		if(pBpar->search_string("--wgd"))
		{
			//cout <<"test";
			pBpar->find_double_opts(pBpar->bimbam_gprobs_file,"--wgd");
			pBpar->read_bimbam_gprobs=true;
			cfe(pBpar->bimbam_gprobs_file,pBpar->read_bimbam_gprobs);
			pBpar->code_readType		="bimbam.gprob";
			//cout <<"test1";
			pBpar->find_double_opts(pBpar->bimbam_pos_file,"--pos");
			//cout <<"test2" ;
			if(pBpar->search_string("--pos"))
			{
				pBpar->read_bimbam_pos	=true;
				cfe(pBpar->bimbam_pos_file,pBpar->read_bimbam_pos);
			}else
			{
				string msg= "file "+pBpar->bimbam_gprobs_file +"  can only be read together with \"--pos\" command. Please give snpinfo  file  with \"--pos\" option.\n";
				error(msg);

			}

		}
		//for bimbam --bimbam-best guess reading option
		if(pBpar->search_string("--wbg"))
		{
			pBpar->find_double_opts(pBpar->bimbam_bestguess_file,"--wbg");
			pBpar->read_bimbam_bestguess=true;
			pBpar->code_readType		="bimbam_bestguess";
			cfe(pBpar->bimbam_bestguess_file,pBpar->read_bimbam_bestguess);
			 // --pos option with bimbam --geno reads different file as here. see above for --geno option .
			pBpar->find_double_opts(pBpar->bimbam_pos_file,"--pos");
			if(pBpar->search_string("--pos"))
			{
				pBpar->read_bimbam_pos	=true;
				cfe(pBpar->bimbam_pos_file,pBpar->read_bimbam_pos);
			}else
			{
				string msg= "file "+pBpar->bimbam_file +"  can only be read together with \"--pos\" command. Please give SNP information  file (used in BIMBAM ) with \"--pos\" option.\n";
				error(msg);

			}

		}

	// The following part is for maf-thresh value in bimbam best guess genotype data and bimbam grob data
		// the following part is for maf thresh.
		// for maf
		if(pBpar->read_bimbam_bestguess || pBpar->read_bimbam_gprobs)
		{
			//----------------------------------------//
			if(pBpar->search_string("--maf-thresh") )
			{
			pBpar->find_double_opts(pBpar->temp_string, "--maf-thresh");
			pBpar->default_maf=false;
			if(!pBpar->default_maf)
			{
				pBpar->maf_value		=atof(pBpar->temp_string.c_str());
				pBpar->temp_string 	="";
			}
			}
		}
		//end of maf thresh


}
//
void parse_haploview_commands(Bpar * const pBpar ){

}
//
void parse_eigensoft_commands(Bpar * const pBpar )
{

	//--------------------------------------------//
		//for smartpca
	//--------------------------------------------//
		//cout<<boolalpha<<(pBpar->oFormat=="")<< endl;
		//cout << "----------------------"<<endl;
		 if(pBpar->oFormat=="smartpca"|| pBpar->oFormat=="eigensoft")
		 {
			 if(pBpar->search_string("--group-label"))
			 {
				 pBpar->find_double_opts(pBpar->group_label_fileName,"--group-label");
				  pBpar->read_group_label=true;
				  cfe(pBpar->group_label_fileName,pBpar->read_group_label);
			 }
		 }
	//---------------------------------------------------------------------------------------------//
}
void parse_fst_calculation_command(Bpar* const pBpar){
	// parse fst command
		if((pBpar->search_string("--fst"))){
			 bool  given_popInfo = pBpar->search_string("--group-label");
			 if(given_popInfo)
			 {
					pBpar->find_single_opt("--fst");
					pBpar->calculate_fst =true;
					pBpar->find_double_opts(pBpar->group_label_fileName ,"--group-label");
			 }
			 else{
				 string msg= " You have asked to calculate Fst with \"--fst\" command option.\n"
						 "To calculate Fst among different populations,  individuals should be"
						 "assigned to different population groups  with option \"--group-label\".\n"
						 "So please assign the individuals to the correct population groups.\n"
				 	 	 "Note the assignment of indiviudals to different popolation can be done by using \n"
				 	 	 "two columns, namely the first: individual ids; and second: population assignments.\n"
						 "e.g.\n"
						 "indiv1\t pop1\n"
						 "indiv2\t pop1\n"
						 "indiv3\t pop2\n";

				 error(msg);
			 }


		}

}
//end of fst_calculation commands

void parse_general_commands(Bpar * const pBpar )
{
	//snpinfo and pedinfo is correct for all type of read
		// for --snpinfo
		if(pBpar->search_string("--snpinfo"))
		{
			string code="--snpinfo";
			pBpar->read_extraMap			=true; // extra map file
			bool temp_bool=(pBpar->code_readType!="");
			//bool temp_bool = pBpar->read_map || pBpar->read_dat ||pBpar->read_gens || (pBpar->read_mlgeno) ||(pBpar->read_mlprob);
			if(temp_bool)
			{
				pBpar->find_double_opts(pBpar->extraMap_file, code);
				cfe(pBpar->extraMap_file,pBpar->read_extraMap);
			}
			else
			{
				string msg= " Command \""+code+"\" can only be used to update SNP information after a genotype file has been already read by fcGENE. Please give genotype file first\n";
				error(msg);
			}


		}
		//for --pedinfo
		if(pBpar->search_string("--pedinfo"))
		{
			string code	="--pedinfo";
			pBpar->read_extraPed=true;
			bool temp_bool=(pBpar->code_readType!="");
			if(temp_bool)
			{
				pBpar->find_double_opts(pBpar->extraPed_file, code);
				cfe(pBpar->extraPed_file,pBpar->read_extraPed);
			}
			else
			{
			 string msg ="command \""+code+"\"  can only be used to update pedgree info. \n";
			 error(msg);

			}
		}
		// for condition --iformat ;
		if(pBpar->search_string("--iformat"))
		{
			//pBpar->iFormat = pBpar->find_double_opts("--iformat");
			pBpar->find_double_opts(pBpar->iFormat,"--iformat");
		}
		if(pBpar->search_string("--exclude"))
		{
			 string code="--exclude";
			 pBpar->exclude_snps =true;
			 bool __tmp_bool		 =(pBpar->code_readType!="");
			 if(__tmp_bool)
			 {
				 pBpar->find_double_opts(pBpar->exclude_snplist_file,code);
				 cfe(pBpar->exclude_snplist_file,pBpar->exclude_snps);
			 }
			 else
			 {
				 string msg= " Command \""+code+"\" can only be used to exclude SNPs  after a genotype file has been already read by fcGENE. Please give genotype files first\n";
				 error(msg);
			 }
		}
		//remove pids
		if(pBpar->search_string("--remove"))
		{
			 string code="--remove";
			 pBpar->remove_indivs =true;
			 bool __tmp_bool		 =(pBpar->code_readType!="");
			 if(__tmp_bool)
			 {
				 pBpar->find_double_opts(pBpar->remove_indivlist_file,code);
				 cfe(pBpar->remove_indivlist_file,pBpar->remove_indivs);
			 }
			 else
			 {
				 string msg= " Command \""+code+"\" can only be used to remove individuals  after a genotype file has been already read by fcGENE. Please give genotype files first\n";
				 error(msg);
			 }
		}

		// for general extra commands
		//---------------------------------------------------------------------------------------------//
				pBpar->is_code_plink		= (pBpar->code_readType=="plink")
											||(pBpar->code_readType=="plink-dosage")
											||(pBpar->code_readType=="plink-rawA")
											||(pBpar->code_readType=="bplink")
											||(pBpar->code_readType=="plink-rawAD");
				pBpar->is_code_mach 		=(pBpar->code_readType=="mach")
											||(pBpar->code_readType=="mach.geno")
											||(pBpar->code_readType=="mach.mlgeno")
											||(pBpar->code_readType=="mach.mlprob")
											||(pBpar->code_readType=="mach_hapmap");
				pBpar->is_code_minimac	=pBpar->code_readType=="minimac.prob";
				pBpar->is_code_impute	=(pBpar->code_readType=="impute") ||(pBpar->code_readType=="impute_hapmap");
				pBpar->is_code_snptest	=(pBpar->code_readType=="snptest");
				pBpar->is_code_shapeit	=(pBpar->code_readType=="shapeit_haplotype");
				pBpar->is_code_beagle	=(pBpar->code_readType=="beagle") ||(pBpar->code_readType=="beagle_gprobs");
				pBpar->is_code_bimbam	=(pBpar->code_readType=="bimbam") ||(pBpar->code_readType=="bimbam.gprob")||(pBpar->code_readType=="bimbam_bestguess");
				pBpar->is_code_rformat	=(pBpar->code_readType=="rformat");
				// At the moment, general commands are then true if any  of the input files are given.
				pBpar->is_general_command=pBpar->is_code_plink||pBpar->is_code_mach ||pBpar->is_code_minimac
						||pBpar->is_code_impute ||pBpar->is_code_snptest||pBpar->is_code_beagle
						||pBpar->is_code_rformat
						||pBpar->is_code_bimbam
						||pBpar->is_code_shapeit;


			//if commands are given with --new-start --new-end
				if(pBpar->given_new_command&& pBpar->is_general_command)
				{
					if(pBpar->search_string("--merge"))
					{
						pBpar->find_single_opt("--merge");
							pBpar->merge_data=true;
					}

				}

	//for uncompressed format
				if(pBpar->search_string("--uncompress"))
				{
					 pBpar->find_single_opt("--uncompress");
					 pBpar->uncompressed =true;
				}
	// for splitting snpwise and individual wise
		parse_ssplit_isplit_commands(pBpar);
				// for outout format --oformat:
			//1. first ofrmat should be recognized.
			//2.   the word after --oformat should be  fed in
			// the rest is in general.cpp file
			if(pBpar->search_string("--oformat"))
			{
				//pBpar->oFormat = pBpar->find_double_opts("--oformat");
				pBpar->find_double_opts(pBpar->oFormat,"--oformat");
				//supported formats
				bool _of_bool =((pBpar->oFormat=="plink")
								||(pBpar->oFormat=="plink-bed")
								||(pBpar->oFormat=="plink-recodeA")
								||(pBpar->oFormat=="plink-recodeAD")
								||(pBpar->oFormat=="plink-dosage")||
								(pBpar->oFormat=="recodeA-dose")||
								//(pBpar->oFormat=="recodeAD-dosage")||
								(pBpar->oFormat=="mach")||
								(pBpar->oFormat=="minimac")||
								(pBpar->oFormat=="impute")||
								(pBpar->oFormat=="snptest")||
								(pBpar->oFormat=="beagle")||
								 (pBpar->oFormat=="bimbam")||
								(pBpar->oFormat=="haploview")||
								(pBpar->oFormat=="eigensoft")||
								(pBpar->oFormat=="smartpca")||
								(pBpar->oFormat=="genabel")|| 		//genabel
								(pBpar->oFormat=="probabel")|| 		// probabel
								(pBpar->oFormat=="r")|| (pBpar->oFormat=="R")||
								(pBpar->oFormat=="r-dose")||  (pBpar->oFormat=="R-dose")||
								(pBpar->oFormat=="phase")||  (pBpar->oFormat=="fastphase")||
								(pBpar->oFormat=="mach-ref")||
								(pBpar->oFormat=="vcf")
								);
				if(!_of_bool)
					error("To convert genotype data into \""+pBpar->oFormat+"\" format is not possible.\n");
				if(pBpar->code_readType!="")
				{

					pBpar->change_format=true;
				}
				else
				{
					string msg="Without giving any input files, the command  \"--oformat "+pBpar->oFormat+"\" is not be accepted.\n";
					error(msg);

				}
				// If pBpar->oFormat is "r" or something like this
				if((pBpar->oFormat=="r")|| (pBpar->oFormat=="R")){
					if(	pBpar->search_string("--transpose")){
						pBpar->find_single_opt("--transpose");
						pBpar->transpose=true;
					}
				}

			}

	 //for help
	 if(pBpar->search_string("--help"))
	 {
		 pBpar->hilfe=true;
		 pBpar->find_single_opt("--help");
	 }
	 // for writing snp list
	 if( (pBpar->search_string("--write-snpinfo"))  )
	 {
		 if( pBpar->code_readType!="")
		 {
			 pBpar->write_snpinfo=true;
			 pBpar->find_single_opt("--write-snpinfo");
		 }
		 else
		 {
			 string msg="Without giving any input files, the command  \"--write-snpinfo \" is not be accepted.\n";
			 error(msg);

		 }

	 }
	// force option
	 if(pBpar->search_string("--force"))
	 {
		 if(pBpar->code_readType!="")
		 {
			 pBpar->given_force=true;
			 pBpar->find_double_opts(pBpar->force_option, "--force");
			 parse_force_filter_option(pBpar->force_option,  pBpar->force_first_part, pBpar->force_2nd_part, pBpar->given_force);

		 }
		 else
		 {
			 string msg="Without giving any input files, the command  \"--force \" is not be accepted.\n";
			 error(msg);

		}

	 }
	//next
		// for writing ped ifnormation
		if((pBpar->search_string("--write-pedinfo")))
		{

			 if( pBpar->code_readType!="")
			 {
				pBpar->write_pedinfo=true;
				pBpar->find_single_opt("--write-pedinfo");
			 }
			 else
			 {
					string msg="Without giving any input files, the command  \"--write-pedinfo \" is not be accepted.\n";
					error(msg);

			 }

		}
	// write ped-list
		if((pBpar->search_string("--write-pedlist")))
		{
			 if( pBpar->code_readType!="")
			 {
				pBpar->write_pedlist=true;
				pBpar->find_single_opt("--write-pedlist");
			 }
			 else
			 {
					string msg="Without giving any input files, the command  \"--write-pedlist \" is not be accepted.\n";
					error(msg);
			 }
		}
	// write snp-list
		if((pBpar->search_string("--write-snplist")))
		{
			 if( pBpar->code_readType!="")
			 {
				pBpar->write_snplist=true;
				pBpar->find_single_opt("--write-snplist");
			 }
			 else
			 {
					string msg="Without giving any input files, the command  \"--write-snplist \" is not be accepted.\n";
					error(msg);
			 }
		}
		parse_fst_calculation_command(pBpar);



 }//end of function
//
void parse_summary_statistics_commands(Bpar * const pBpar )
{

	//for freq
		 if(pBpar->search_string("--freq") )
		 {
			 if( pBpar->code_readType!="")
			 {
				 pBpar->write_freq=true;
				pBpar->find_single_opt("--freq");
			 }
			 else
			 {
				 string msg="Without giving any input files, the command  \"--freq \" is not be accepted.\n";
				 error(msg);
			 }
		 }
		 // next
		 //hwe
		if(pBpar->search_string("--hardy")){
			if(pBpar->code_readType!="")
			{
				pBpar->write_hwe=true;
				pBpar->find_single_opt("--hardy");
			}
			else
			{
				string msg="Without giving any input files, the command  \"--hardy\" is not be accepted.\n";
				error(msg);
			}

		}
		//missing
		if(pBpar->search_string("--crate"))
		{
			if(pBpar->code_readType!="")
			{
				pBpar->write_indv_callrate=true;
				pBpar->write_snp_callrate=true;
				pBpar->find_single_opt("--crate");
			}
			else
			{
				string msg="Without giving any input files, the command  \"--crate\" is not be accepted.\n";
				error(msg);
			}
		}
	// for --filter-snp command option
	if(pBpar->search_string("--filter-snp"))
	{
		if(pBpar->code_readType!="")
		{
			pBpar->filter_snp =true;
			pBpar->find_double_opts(pBpar->filter_snp_opts,"--filter-snp");
			//cout<< pBpar->filter_snp_opts <<"  ";
			parse_force_filter_option(pBpar->filter_snp_opts,pBpar->filter_snp_firstPart,pBpar->filter_snp_secondPart,pBpar->filter_snp);
		}
		else
		{
			string msg="Without giving any kind of input files, command option  \"--filter-snp \" can not be e accepted.\n";
			error(msg);

		}
	}
	// for --filter-ind command option
	if(pBpar->search_string("--filter-indiv"))
	{
		if(pBpar->code_readType!="")
		{
			pBpar->filter_indiv =true;
			pBpar->find_double_opts(pBpar->filter_indiv_opts,"--filter-indiv");
			parse_force_filter_option(pBpar->filter_indiv_opts,pBpar->filter_indiv_firstPart,pBpar->filter_indiv_secondPart,pBpar->filter_indiv);
		}
		else
		{
			string msg="Without giving any kind of input files, command option  \"--filter-indiv \" can not be e accepted.\n";
			error(msg);

		}
	}

}//end of function

void cfe(const string & file, const bool &  bValue ){ // check if file exists.

	//bValue will be used for pBpar->map:
	//cout << "bValue"<<bValue<< endl; // only for test
	if(!bValue)
	{
		size_t pos=file.find_last_of(".");
					string msg="Your command  does not include a "+file.substr(pos+1)+\
						" file which has  default name ["+file +"].\n";
					error(msg);
	}
	else
	{
		ifstream inp;
		inp.open(file.c_str(),ifstream::in);
		if(inp.fail())
		{
			inp.close();
			inp.clear(ios::failbit);
			string msg="There does not exist any file named "+file +".\n";
			error(msg);
		}
		inp.close();
	}
	return;
}
// this function is written especially for creating SNPTEST sample file from plink formatted covariate data
void parse_force_filter_option(const string& force_option,  vector<string>&firstPart, vector<string>&secondPart, const bool&bValue)

{
	if(bValue)
	{
		vector<string>myoptions; // this will be all sub arguments
		myoptions.clear();
		string _option=force_option;
		myoptions= Argcv::parse_subArg(_option, ","); // this will separate all sub arguments
		//cout << "option size: "<< myoptions.size() <<endl; //debug 
		_option="";
		firstPart.clear();
		secondPart.clear();	
		size_t found;
		for(unsigned int i=0;i<myoptions.size();++i)
		{
			found=myoptions[i].find_first_of('=');
			if(found!=string::npos)
			{
				firstPart.push_back(myoptions[i].substr(0, found));
				secondPart.push_back(myoptions[i].substr(found+1));
				//std::cout << myoptions[i].substr(0, found) << " " << myoptions[i].substr(found+1) <<endl; //debug 
			}
			else
			{
			  string pr_msg ="Your command option \""+myoptions[i]+"\" given after command option \"--force\" should contain a \'=\' symbol.\n";
			  printLIN(pr_msg);
			  pr_msg="Please correct the command option.\n  Argument before '=' should represent the name of variable you want to change and the argument after '=' should represent a value or a fileName where the values you want to include, are mentioned.\n";
			  error(pr_msg);
			}
		}
		
		if(firstPart.size()!=secondPart.size())
		{
			string pr_msg= "problem in parsing --force option. Please contact the software developer.\n";
			error(pr_msg);
		}
		vector<string>::iterator it1;
		vector<string>::iterator it;
		for( it=firstPart.begin();it!=firstPart.end();++it)
		{
			_option =*it;
			it1=find(it+1,firstPart.end(), _option);
			if((it1!=firstPart.end()))
			{
				string pr_msg= _option+ " is given more than 1 times  after command option  \"--force\". Please use  each option only once.\n";
				//cout << "found at " << (int)(it1+1-firstPart.begin())<<endl; //debug 
				error(pr_msg);
				
			}
		}
	}

return;}
//

string paste_snp_or_ped_info(const string& paste_option,vector<string>&pasted_strings,const string&indicator){

	vector<string>myoptions; // this will be all sub aruguments
	myoptions.clear();
	string _option=paste_option;
	string seperator ="";
	myoptions= Argcv::parse_subArg(_option, ","); // this will separate all sub arguments
	if(myoptions.size()==2)
	{
		pasted_strings.push_back(myoptions[0]);
		pasted_strings.push_back(myoptions[1]);
		//cout <<pasted_strings[0]<<":-> "<<pasted_strings[1];
	}
	else if(myoptions.size()>2)
	{
			// add all sub strings as to be pasted except the last one.
			//if last one has "sep= ", then take string after sep= as seperator
			//otherwise take sep="" and last sub argument also as string to be pasted.
			size_t found;
			for(unsigned int i=0;i<(myoptions.size()-1);++i)
				pasted_strings.push_back(myoptions[i]);
			seperator=myoptions[(myoptions.size()-1)];
			found= seperator.find_first_of('=');
			if(found!=string::npos && (seperator.substr(0, found)=="sep"))
				seperator=seperator.substr(found+1);
			else
			{
				pasted_strings.push_back(seperator);
				seperator="";
			}
	}
	else{
			string pr_msg ="Your command option indicated by  \""+indicator+ "\" should contain  two or more than two sub arguments each of which are seperated by \",\".\n";
			printLIN(pr_msg);
			error(pr_msg);

	}

return seperator;}
/*
void parse_paste_option(const string& paste_option,vector<string>&pasted_strings,string&seperator,const bool to_paste){
	if(to_paste)
	{
		vector<string>myoptions; // this will be all sub aruguments
		myoptions.clear();
		string _option=paste_option;
		myoptions= pBpar->parse_subArg(_option, ","); // this will seperate all sub arguments
		//cout << "option size: "<< myoptions.size() <<endl; //debug
		_option="";
		paste_strings.clear();
		size_t found;
		string sep_firstPart;
		string sep_secondPart;
		bool found_once=true;
		for(unsigned int i=0;i<myoptions.size();++i)
		{
			found=myoptions[i].find_first_of('=');
			if(found!=string::npos)
			{
				if(found_once)
				{
					firstPart.myoptions[i].substr(0, found);
					secondPart.myoptions[i].substr(found+1);
				}
				else
				{
					string pr_msg ="Your command option \""+myoptions[i]+"\" given after command option \"--paste\" should contain symbol \'=\' only once.\n";
					printLIN(pr_msg);
					pr_msg="Please correct the command option.\n  Argument after '=' should represent seperator in between the paste arguments.\n";
					error(pr_msg);
				}
			}
			else{

			}
		}

	}
}
*/
void parse_option_iid(Bpar * const pBpar ){
	// pasting iid and fid
	if(pBpar->search_string("--iid"))
{
		string string_to_paste="";
		if(pBpar->code_readType!="")
		{
			pBpar->paste_iid=true;
			pBpar->find_double_opts(string_to_paste, "--iid");
			//cout <<string_to_paste << " \n";//debug
			//parse and find substrings
			pBpar->paste_seperator=paste_snp_or_ped_info(string_to_paste,pBpar->strings_2bpasted,"--iid");
			//cout << "seperator: "<< pBpar->paste_seperator <<" \n";//debug
			vector<string> temp_strings_2bpasted =pBpar->pedInfo_idicators;
			sort(temp_strings_2bpasted.begin(),temp_strings_2bpasted.end());
			vector<string>::iterator it =pBpar->strings_2bpasted.begin();
			bool temp_tfvalue=false;
			int unsigned info_pos=0;
			for(unsigned int i=0; i<pBpar->strings_2bpasted.size();++i)
			{
				temp_tfvalue= binary_search(temp_strings_2bpasted.begin(),temp_strings_2bpasted.end(),pBpar->strings_2bpasted[i]);
				if(temp_tfvalue)
				{
					it =find(pBpar->pedInfo_idicators.begin(),pBpar->pedInfo_idicators.end(),pBpar->strings_2bpasted[i]);
					info_pos =(int)(it-pBpar->pedInfo_idicators.begin());
					//cout << info_pos <<" ";
						pBpar->pos_pedinfo_2bpasted.push_back(info_pos);
				}else
				{

					string msg="Command option \""+ pBpar->strings_2bpasted[i]+"\" given after option \"--iid\" ca not be  be accepted.\n";
					error(msg);

			 }
		}
	}else
	{

		string msg="Without giving any input files, the command option \"--iid\" is not be accepted.\n";
		error(msg);

	}

}
}
// the following function checks if different input and output files eixits
void check_input_output_file_names_n_paths(Bpar* const pBpar)
{
	ifstream ifs ;
	ofstream ofs;
	//check output files first
	const string msg_out ="Above file can't be opened to write. Please check whether the file path and/or file name are correct.\n";
	const string msg_in ="Above file can't be opened to read. Please check whether the file path and/or file name are correct.\n";
	//for output files :
	/*
	ofs.clear(); ofs.close();
	ofs.open(pBpar->output_fileName.c_str(),ios::out);
	if(!ofs){
		printLIN("Invalid file path or file name: "+pBpar->output_fileName+":\n");
		error(msg_out);
	}else{
		ofs.clear();
		ofs.close();
		char buffer [101];
		puts(buffer,pBpar->output_fileName.c_str());
		remove(buffer);
	}
	*/

	// now different input file names:
	//plink
		const string  read_code_type 	=pBpar->code_readType;
		if(read_code_type =="plink")
		{
			ifs.open(pBpar->map_file.c_str(),ios::in);
			if(!ifs){
				printLIN("Invalid file path or file name: "+pBpar->map_file+":\n");
				error(msg_in);
			}else{
				ifs.clear();
				ifs.close();
			}
			//ped file
			ifs.open(pBpar->ped_file.c_str(),ios::in);
			if(!ifs){
				printLIN("Invalid file path or file name: "+pBpar->ped_file+":\n");
				error(msg_in);
			}else{
				ifs.clear();
				ifs.close();
			}


		}
		if(read_code_type =="plink-dosage")
		{
			//map file
			ifs.open(pBpar->map_file.c_str(),ios::in);
			if(!ifs){
				printLIN("Invalid file path or file name: "+pBpar->map_file+":\n");
				error(msg_in);
			}else{
				ifs.clear();
				ifs.close();
			}
			//plink_fam_file file
			ifs.open(pBpar->plink_fam_file.c_str(),ios::in);
			if(!ifs){
				printLIN("Invalid file path or file name: "+pBpar->plink_fam_file+":\n");
				error(msg_in);
			}else{
				ifs.clear();
				ifs.close();
			}
			//plink_dosage_file
			ifs.open(pBpar->plink_dosage_file.c_str(),ios::in);
			if(!ifs){
				printLIN("Invalid file path or file name: "+pBpar->plink_dosage_file+":\n");
				error(msg_in);
			}else{
				ifs.clear();
				ifs.close();
			}

		}
		if(read_code_type =="plink-rawA"||read_code_type=="plink-rawAD")
		{
			//map file
			ifs.open(pBpar->map_file.c_str(),ios::in);
			if(!ifs){
				printLIN("Invalid file path or file name: "+pBpar->map_file+":\n");
				error(msg_in);
			}else{
				ifs.clear();
				ifs.close();
			}
			//extra map file
			ifs.open(pBpar->extraMap_file.c_str(),ios::in);
			if(!ifs){
				printLIN("Invalid file path or file name: "+pBpar->extraMap_file+":\n");
				error(msg_in);
			}else{
				ifs.clear();
				ifs.close();
			}
			//extra map file
			if(pBpar->read_plink_rawA)
			{
				ifs.open(pBpar->plink_rawA_file.c_str(),ios::in);
				if(!ifs){
					printLIN("Invalid file path or file name: "+pBpar->plink_rawA_file+":\n");
					error(msg_in);
					}else{
						ifs.clear();
						ifs.close();
				}
			}
			else
			{
				ifs.open(pBpar->plink_rawAD_file.c_str(),ios::in);
				if(!ifs){
					printLIN("Invalid file path or file name: "+pBpar->plink_rawAD_file+":\n");
					error(msg_in);
				}else{
					ifs.clear();
					ifs.close();
				}
			}
		}
		if(read_code_type =="rformat")
		{
			ifs.open(pBpar->rgeno_file.c_str(),ios::in);
			if(!ifs){
				printLIN("Invalid file path or file name: "+pBpar->rgeno_file+":\n");
				error(msg_in);
			}else{
				ifs.clear();
				ifs.close();
			}
		}
	//mach

	//impute

	//beagle

	//bimbam

	//shapeit
	if(read_code_type =="shapeit_haplotype")
	{
		ifs.open(pBpar->haps_file.c_str(),ios::in);
		if(!ifs)
		{
			printLIN("Invalid file path or file name: "+pBpar->haps_file+":\n");
			error(msg_in);
		}else
		{
			ifs.clear();
			ifs.close();
		}
		//check for sample file
		ifs.open(pBpar->shapeit_sample_file.c_str(),ios::in);
		if(!ifs)
		{
			ifs.clear();
			ifs.close();
			printLIN("Invalid file path or file name: "+pBpar->shapeit_sample_file+":\n");
			error(msg_in);
		}

		ifs.clear();
		ifs.close();


	}

	// different general files snp info pedinfo etc
	//grouLabel file name for smartpca and fst calculation
	if(pBpar->group_label_fileName!="")
	{
		ifs.clear(); ifs.close();
		ifs.open(pBpar->group_label_fileName.c_str(),ios::in);
		if(!ifs){
			printLIN("Invalid file path or file name: "+pBpar->group_label_fileName+":\n");
			error(msg_in);
		}else{
			ifs.clear();
			ifs.close();
		}

	}
	//	fst popinfo file


}//end of funciton

void parse_rformat_commands(Bpar* const pBpar){
	//---------------------------------------------------//
		//for r format
		//---------------------------------------------------//
		if(pBpar->search_string("--rgeno"))
		{
			pBpar->find_double_opts(pBpar->rgeno_file,"--rgeno");
			pBpar->read_rgenotype		=true;
			pBpar->code_readType 		="rformat";

			//
			if(pBpar->search_string("--snpinfo"))
			{
				string code="--snpinfo";
				pBpar->read_extraMap			=true; // extra map file
				bool temp_bool=(pBpar->code_readType!="");
				//bool temp_bool = pBpar->read_map || pBpar->read_dat ||pBpar->read_gens || (pBpar->read_mlgeno) ||(pBpar->read_mlprob);
				if(temp_bool)
				{
					pBpar->find_double_opts(pBpar->extraMap_file, code);
						cfe(pBpar->extraMap_file,pBpar->read_extraMap);
				}
				else
				{
					string msg= " Command \""+code+"\" can only be used to update SNP information after a genotype file has been already read by fcGENE. Please give genotype file first\n";
						error(msg);
				}


			}else{
				 string _msg ="R format file can only be read  if the allele info file is given with --snpinfo option.\n"
						 "Allele mentioned  in the column with header \"allele1\" is the one whose homozygote is mentioned as 2, and \"allele2\" is the one, whose homozygote is\n"
						 "mentioned as \"0\" in the genotype file. In other words, R-format genotype data is the counts of first alleles at each SNP. \n";
				 error(_msg);
			}
		}

}
//-----------------------------------//
void parse_ssplit_isplit_commands(Bpar* const pBpar)
{
	string _split_option="";
	vector<string> snp_low_up;
	vector<string> ind_low_up;
	int _tmp_first 	=0;
	int _tmp_2nd	=0;

	if(pBpar->search_string("--ssplit") || pBpar->search_string("--bpsplit") )
	{
		bool _tmp_b1 =pBpar->search_string("--ssplit") && pBpar->search_string("--bpsplit");
		if(!_tmp_b1) //eitehr of them should be true. both cant come at same time
		{
			if(pBpar->search_string("--ssplit") )
			{
				pBpar->find_double_opts(_split_option,"--ssplit");
				//cout << _split_option<<endl;
				pBpar->ssplitt =true;
			}else{

				pBpar->find_double_opts(_split_option,"--bpsplit");
				//cout << _split_option<<endl; //debug
				pBpar->bpsplitt =true;
			}
			//pBpar->v_ssplitt= Argcv::parse_subArg(_split_option, ","); // this will separate all sub arguments
			snp_low_up= Argcv::parse_subArg(_split_option, ","); // this will separate all sub arguments
			_split_option="";
		}else{
			string _msg ="Options  \"--ssplit\" and \"--bpsplit\" can't be given at the same time. Please use either \"--ssplit\" or \"--bpsplit\" .\n";
			error(_msg);
		}
	}
	if(pBpar->search_string("--isplit"))
	{
			pBpar->find_double_opts(_split_option,"--isplit");
			pBpar->isplitt =true;
		//	pBpar->v_isplitt= Argcv::parse_subArg(_split_option, ","); // this will separate all sub arguments
			ind_low_up	= Argcv::parse_subArg(_split_option, ","); // this will separate all sub arguments

	}
	//check
	 unsigned int _sSize =snp_low_up.size();
	 unsigned int _iSize =ind_low_up.size();
	string _suffix ="";
	string _cur_snp_opt ="";
	string _cur_ind_opt="";
	vector<string> tmp_low_up;
	//vector<int>int_snp_low_up;
	//vector<int>ind_low_up;
	if(_sSize>1&& _iSize>1)
	{
		if(_sSize!=_iSize){
			string _msg ="error\n";
			error(_msg);
		}
	}
	//find the larger size
	unsigned int _nSize = max(_sSize,_iSize);
	// const unsigned int SZ_bparam		=Argcv::bparam_vec.size(); // size of basic parameters group
	//this loop creates the new parameter class and new General class
	//	cout << _nSize <<endl;
	for(unsigned int j=0;j<_nSize; ++j)
	{
		if(_nSize==_sSize)
		{
			 //snpwise split
				_cur_snp_opt	=snp_low_up[j];
			 _suffix 			="_snp"+_cur_snp_opt;
			 if(_iSize==1)
				 _cur_ind_opt	=ind_low_up[0];
			 else if(_iSize>1)
				 _cur_ind_opt 	=ind_low_up[j];
			 if(_iSize>0)
			 _suffix			+="_indiv"+_cur_ind_opt;
			 else
				 _suffix+="_indiv_all";
		}
		 else
		 {
			 //individual wise split
				_cur_ind_opt=ind_low_up[j];
				_suffix ="_indiv"+_cur_ind_opt;
				if(_sSize==1)
					_cur_snp_opt=snp_low_up[0];
				else if(_sSize>1)
					_cur_snp_opt =snp_low_up[j];
				 if(_iSize>0)
					 _suffix+="_snp"+_cur_snp_opt;
				 else
					 _suffix+="_snp_all";
		 }
		pBpar->v_split_suffix.push_back(_suffix);
		//find the lower and upper index of _snpidx
		if(_cur_snp_opt!="")
		{
			//cout<<_cur_snp_opt<<endl;
			tmp_low_up =vector_from_string(_cur_snp_opt,'-');
			if(tmp_low_up.size()!=2)
			{
				string _msg ="Command option given after --ssplit --bpsplit is incorrect.\n"
						 	 "Command-option --ssplit must contain a lower and a upper index that lie in between any two indices of snplist.\n"
							"Command-option --bpslit must contain a lower and a upper index base pair position.\n"
							"Then a new data can be generated with the SNPs that lie between the two indices.\n"
							"Lower and upper indices must be separated by \"-\" sign. e.g. 1-10.\n"
							"Lower must be greater than zero and smaller than the upper.\n"
							"Multiple separation tasks can also be given to fcGENE.\n"
							"An example of \"--ssplit\" command option is : \"--ssplit 1-10,11-20,20-30\".\n";
					error(_msg);

			}
			_tmp_first =atoi(tmp_low_up[0].c_str());
			_tmp_2nd	=atoi(tmp_low_up[1].c_str());

			if((_tmp_first>=_tmp_2nd) )
			{
				string _msg ="Command option given after --ssplit --bpsplit is incorrect.\n"
					 	 "Command-option --ssplit must contain a lower and a upper index that lie in between any two indices of snplist.\n"
						"Command-option --bpslit must contain a lower and a upper index base pair position.\n"
						"Then a new data can be generated with the SNPs that lie between the two indices.\n"
						"Lower and upper indices must be separated by \"-\" sign. e.g. 1-10.\n"
						"Lower must be greater than zero and smaller than the upper.\n"
						"Multiple separation tasks can also be given to fcGENE.\n"
						"An example of \"--ssplit\" command option is : \"--ssplit 1-10,11-20,20-30\".\n";
				error(_msg);

			}
			if(pBpar->ssplitt&& _tmp_first<=0)
			{
				string _msg ="Command option given after --ssplit --bpsplit is incorrect.\n"
					 	 "Command-option --ssplit must contain a lower and a upper index that lie in between any two indices of snplist.\n"
						"Command-option --bpslit must contain a lower and a upper index base pair position.\n"
						"Then a new data can be generated with the SNPs that lie between the two indices.\n"
						"Lower and upper indices must be separated by \"-\" sign. e.g. 1-10.\n"
						"Lower must be greater than zero and smaller than the upper.\n"
						"Multiple separation tasks can also be given to fcGENE.\n"
						"An example of \"--ssplit\" command option is : \"--ssplit 1-10,11-20,20-30\".\n";
				error(_msg);

			}

			pBpar->v_ssplitt.push_back( _tmp_first);
			pBpar->v_ssplitt.push_back(_tmp_2nd);
			tmp_low_up.clear();
			_tmp_first	=0;
			_tmp_2nd 	=0;

		}//end of checking if _cur_snp_opt is not empty
		if(_cur_ind_opt!="")
		{
			//cout <<_cur_snp_opt<<endl;
			tmp_low_up =vector_from_string(_cur_ind_opt,'-');
			//	cout <<tmp_low_up.size()<<endl;
			//cout << (tmp_low_up[0]>=tmp_low_up[1])<<endl;

			if(tmp_low_up.size()!=2)
			{
				string _msg =" command option given after --isplitt is incorrect.\n"
							"The command-option must contain a lower and a upper index that lie in between any two indices of "
							"individual-list.\n"
							"Then a new data can be generated with the SNPs that lie between the two indices.\n"
							"Lower and upper indices must be separated by \"-\" sign. e.g. 1-10.\n"
							"Lower index must be greater than 0 and  smaller than the upper.\n"
							"Multiple separation tasks can also be given to fcGENE.\n"
							"An example of \"--isplitt\" command option is : \"--isplitt 1-10,11-20,20-30\".\n";
					error(_msg);

			}
			_tmp_first  =atoi(tmp_low_up[0].c_str());
			_tmp_2nd	=atoi(tmp_low_up[1].c_str());

			if((_tmp_first<=0) ||(_tmp_first>=_tmp_2nd))
			{
				// the given index must be greater than zero
				string _msg =" command option given after --isplitt is incorrect.\n"
						"The command-option must contain a lower and a upper index that lie in between any two indices of "
						"individual-list.\n"
						"Then a new data can be generated with the SNPs that lie between the two indices.\n"
						"Lower and upper indices must be separated by \"-\" sign. e.g. 1-10.\n"
						"Lower index must be greater than 0 and smaller than the upper.\n"
						"Multiple separation tasks can also be given to fcGENE.\n"
						"An example of \"--isplitt\" command option is : \"--isplitt 1-10,11-20,20-30\".\n";
				error(_msg);
			}
			pBpar->v_isplitt.push_back(_tmp_first);
			pBpar->v_isplitt.push_back(_tmp_2nd);
			tmp_low_up.clear();
			_tmp_first	=0;
			_tmp_2nd 	=0;

		}//end of checking if _cur_snp_opt is not empty


	}//end of for loop
	//checking for v_splitt and v_isplitt
		_sSize 	=pBpar->v_ssplitt.size();
		_iSize	=pBpar->v_isplitt.size();
		_nSize	=max(_sSize,_iSize);
		if(_sSize>0||_iSize>0){
			//first character
			 if(fmod((double)_sSize,2)!=0.0|| fmod((double)_iSize,2)!=0.0)
			 {
				 string _msg =" problem in Saving  --splitt and --splitt option.\n";
				 error(_msg);
			 }
			//
			 if(_nSize!=(2*pBpar->v_split_suffix.size())){
				 string _msg =" problem in Saving  --splitt and --splitt option.\n"
				 	 	 	 	"_nSize =(2*pBpar->v_split_suffix.size())"	;
				 error(_msg);
			 }
		}

}

//parse shapeit commands
void parse_shapeit_commands(Bpar* const pBpar){
	//cout <<"hi I am from  shapeit\n";
	if(pBpar->search_string("--haps"))
	{
		pBpar->find_double_opts(pBpar->haps_file,"--haps");
		pBpar->read_haps_file			=true;
		pBpar->find_double_opts(pBpar->shapeit_sample_file,"--sample");
		pBpar->read_shapeit_sample_file 	=true;
		if(pBpar->read_shapeit_sample_file){
			pBpar->read_shapeit_haptype		=true;
			pBpar->code_readType			="shapeit_haplotype";
		}
	}

 }
