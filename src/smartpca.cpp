/*
 * smartpca.cpp
 *  Created on: 17.01.2013
 *  Author: Nab Raj Roshyara
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 *   */
#include"smartpca.h"
#include "helper.h"
#include <sstream>
using namespace std;
//Description of smartpca file format
/*
	*Before  I start writing codes for pca format conversion, the following facts must be
	*considered and read carefully first.
	*fact1:
	*******
	*PED format:
  	*genotype file: see example.ped    *** file name MUST end in .ped ***
  	*snp file:      see example.pedsnp *** file name MUST end in .pedsnp ***
    *indiv file:    see example.pedind *** file name MUST end in .pedind ***
	*fact2:
	*******
	*Question: When I run I get an error message about "idnames too long". What should I do?
	*Answer: The software supports sample ID names up to a max of 39 characters.
	*Longer sample ID names must be shortened. In addition, if your data is in PED format, the default is to concatenate the family ID and sample ID names so that their total length must meet this limit;
	*however, you can set "familynames:	NO" so that only the sample ID name will be used and must meet the 39 character limit.
	fact 3:
	*******
	*Question: Should regions of long-range LD in the genome be removed prior to PCA?
	*Answer: Yes, to avoid principal components that are artifacts of long-range LD
	*it is ideal to remove such regions. See Table 1 of Price et al. 2008 AJHG.
	*However, EIGENSTRAT can subsequently be run to compute disease association statistics
	*using the full set of SNPs.
	fact 4:
 	*******
	Question: Can I use EIGENSTRAT in studies involving imputed SNPs?
	Answer: At the moment, our code does not support probabilistic genotypes
	*that may be produced by imputation programs. This is algorithmically straightforward
	*but due to our limited resources, it may be awhile before we can provide this upgrade.
	*In the meantime, a possible solution is to first run PCA on non-imputed SNPs
	*(this will indicate whether there are ancestry differences between cases and controls)
	*and then run EIGENSTRAT to compute disease association statistics for all SNPs by sampling
	*integer-valued genotypes in the case of imputed SNPs.
	fact 5:
	********
	Question: I'm running on an extremely large number of samples and the software runs
	out of memory. Why?
	Answer: The software uses memory proportional to the square of the number of samples.
	In the case of an extremely large number of samples (e.g. >10,000), the software may
	run out of memory. The fast eigenvector approximation described above would actually
	solve this problem, but is not yet implemented.
	fact 6:
	************
	PED format:
  	 genotype file: 	see example.ped    *** file name MUST end in .ped ***
  	  snp file:     	 see example.pedsnp *** file name MUST end in .pedsnp ***
    convertf also supports .map suffix for this input file name
     indiv file:    see example.pedind *** file name MUST end in .pedind ***
                 convertf also supports the full .ped file (example.ped)
		 for this input file
	Note that
	Mandatory suffix names enable our software to recognize this file format.
	The indiv file contains the first 6 or 7 columns of the genotype file.
	The genotype file is 1 line per individual.  Each line contains 6 or 7 columns
  of information about the individual, plus two genotype columns for
  each SNP in the order the SNPs are specified in the snp file.
  Genotype format MUST be either 0ACGT or 01234, where 0 means missing data.
  	 The first 6 or 7 columns of the genotype file are:
    1st column is family ID.
    2nd column is sample ID.
    3rd and 4th column are sample IDs of parents.
    5th column is gender (male is 1, female is 2)
    6th column is case/control status (1 is control, 2 is case) OR
      quantitative trait value OR population group label.
    7th column (this column is optional) is always set to 1.
    [Note: this release *changed* to output .ped files in 6-column format,
     not in 7-column format.  Also see sevencolumnped parameter below.]
  convertf does not support pedigree information, so 1st, 3rd, 4th columns are
    ignored in convertf input and set to arbitrary values in convertf output.
  In the two genotype columns for each SNP, missing data is represented by 0.
The snp file contains 1 line per SNP.  There are 6 columns (last 2 optional):
  1st column is chromosome.  Use X for X chromosome.
    Note: SNPs with illegal chromosome values, such as 0, will be removed
  2nd column is SNP name
  3rd column is genetic position (in Morgans)
  4th column is physical position (in bases)
  Optional 5th and 6th columns are reference and variant alleles.
    For monomorphic SNPs, the variant allele can be encoded as X.
The indiv file contains the first 6 or 7 columns of the genotype file.
The PED format is used by the PLINK package of Shaun Purcell.
  See http://pngu.mgh.harvard.edu/~purcell/plink/.

*/


void CSMARTPCASNP::write_smartpca_pedsnp_file(const vector<CBSNP*>&genInfo,const string& outFileName){

	 ofstream ofs;
	 ofs.clear(); ofs.close();
	 string pedsnp_fileName =outFileName+".pedsnp";
	 ofs.open(pedsnp_fileName.c_str(),ios::out);
	 if(!ofs)
		 error("out file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");
	 CBSNP* pgenInfo= genInfo[0];
	 for(unsigned int i=0; i<genInfo.size();++i)
	   	{
		 	 pgenInfo= genInfo[i];
	   		if(pgenInfo->quality)
	   		{
	 			ofs<< pgenInfo->nchr << " ";
	 			if(pgenInfo->rsId!="") ofs <<pgenInfo->rsId<< " ";
	 			else if(pgenInfo->snpName!="") ofs <<pgenInfo->snpName<< " ";
	 			else if(pgenInfo->rsId==""&& pgenInfo->snpName=="")
	 			{
	 				string snpNewName="snp"+change_int_into_string(i)+"";
	 				ofs<< snpNewName << " ";
	 			}
	 			if((pgenInfo->cm_pos!=(-1.0))&& pgenInfo->cm_pos!=0.0)
	 				ofs<<pgenInfo->cm_pos << " ";
	 			// this is default value;
	 			else
	 			{
	 	 			if(pgenInfo->bp!=-1  ) {
	 	 				ofs<< (double)(pgenInfo->bp)/100.0/10000.0<<" ";
	 	 			}
	 	 			else
	 				ofs<< "0.0"<< " ";

	 			}
	 			if(pgenInfo->bp!=-1  ) ofs<< (long int)pgenInfo->bp<<" ";
	 			else
	 			{
	 				ofs << i << " ";

	 			}
	 			ofs << pgenInfo->allele1 << " ";
	 			if(pgenInfo->allele2!="")
	 				ofs<< pgenInfo->allele2;
	 			else ofs<<"X"; // for monomorph SNP smartpca accepts X
				ofs <<"\n"; // end of line
	 		}
	   	}
	   	ofs.clear();ofs.close();
	   	printLIN("\t**->pedsnp  file has been written out and saved as  \""+\
	 	  			pedsnp_fileName+"\".\n");

}

//
void CSMARTPCAPED::write_smartpca_ped_file(const vector<CBPED*>&pedInfo,const vector<CBSNP*>&genInfo,const string& outFileName)
{
	ofstream ofs ;
		ofs.clear(); ofs.close();
		const string ped_fileName =outFileName+".ped";
		/*
		 * The genotype file contains 1 line per SNP.
  	  	  Each line contains 1 character per individual:
  	  	  0 means zero copies of reference allele.
  	  	  1 means one copy of reference allele.
  	  	  2 means two copies of reference allele.
  	  	  9 means missing data.
		 */
		ofs.open(ped_fileName.c_str(),ios::out);
		 if(!ofs)
				 error("out file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");
		 CBPED* ppedInfo= pedInfo[0];
		 for(unsigned int i=0; i<pedInfo.size(); ++i)
		 {
			 ppedInfo =pedInfo[i];
			 if(ppedInfo->quality)
			 {
				 if((ppedInfo->famId.size() +ppedInfo->indId.size())>=39)
					 ofs<< i+1<< " "<<ppedInfo->indId <<" ";
				 else
					 ofs<< ppedInfo->famId<< " "<<ppedInfo->indId<< " ";
				 if(ppedInfo->patId!="")ofs<<ppedInfo->patId << " ";
					else ofs<< 0<<" ";
				 if(ppedInfo->matId!="") ofs <<ppedInfo->matId <<" ";
				 else ofs<< 0 << " ";
				 if(ppedInfo->sex) ofs<<"1"<< " "; // 1 is male
				 else if(!ppedInfo->sex && !ppedInfo->miss_sex) ofs<<"2"<<" "; // 2 is female
				 else if(ppedInfo->miss_sex) ofs<<"0"<< " "; // 0 is missing
				 // in plink format here comes phenotype but in smartpca
				 // here comes the grouplabel or traits.
				 if(ppedInfo->groupLabel!="")
					 ofs<<ppedInfo->groupLabel <<" ";
				 else
				 {
					 if(ppedInfo->pheno)
						 ofs<<"2"<< " ";
					 else if(!ppedInfo->pheno && !ppedInfo->miss_pheno)
						 ofs<<"1"<< " ";
					 else if(!ppedInfo->pheno && ppedInfo->miss_pheno)
						 ofs<<"-99.0"<< " ";
				 }
				 CBSNP* pCBSNP=genInfo[0];
				 for(unsigned int j=0;j<genInfo.size(); ++j)
				 {
					 //-----------------------------------------------------------------//
					 pCBSNP=genInfo[j];
					 if(pCBSNP->quality) // default value of quality is true;
					 {
						 if(pCBSNP->geno1.size()!=pedInfo.size() || pCBSNP->geno2.size()!=pedInfo.size() )
						 {
							 cout <<"geno1.size() : "<<pCBSNP->geno1.size() <<"geno2.size():  "<<pCBSNP->geno1.size();
							 //for(unsigned int i=0;i<pCBSNP->geno1.size();++i)
							 //	ofs<< pCBSNP->geno1[i] << " " <<pCBSNP->geno2[i] << " ";
							 ofs.close();
							 error("Problem in writing out "+change_int_into_string(j) +"th SNP. geno1 and geno2 has no equal size\n");
						 }
						 else if(pCBSNP->allele1!="" || pCBSNP->allele2!="")
						 {
							 if( !pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if homozygote1  00
								 ofs <<" "<< pCBSNP->allele1 << " "<<pCBSNP->allele1;
							 else if( pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if homozygote2 11
								 ofs <<" " << pCBSNP->allele2 << " "<<pCBSNP->allele2;
							 else if( !pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if heterozygote
							 {
								 //ofs << " "<<pCBSNP->allele1 << " "<<pCBSNP->allele2;
								 if(pCBSNP->aOrder[i])
									 ofs <<" "<< pCBSNP->allele1 << " " <<pCBSNP->allele2;
								 if(!pCBSNP->aOrder[i])
									 ofs<<" " <<pCBSNP->allele2 << " " <<pCBSNP->allele1;
							 }
							 else if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if missing
								 ofs <<" "<<CBSNP::_miss_gvalue<< " "<<CBSNP::_miss_gvalue;
							 else
							 {
								 error("Problem in writing out "+change_int_into_string(j) +"th SNP.\n");
							 }
						 }
						 else if (pCBSNP->allele1!=""&& pCBSNP->allele2=="")
						 {
							 if(!pCBSNP->geno1[i] && !pCBSNP->geno2[i])
								 ofs<< pCBSNP->allele1  <<pCBSNP->allele1; //<<",  ";
							 else if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if missing
								 ofs <<" "<<CBSNP::_miss_gvalue<< " "<<CBSNP::_miss_gvalue;
							 else
							 {
								 string prmsg= pCBSNP->rsId +" can have only homozygoze majors because no other allele is detected for it.\n";
								 error(prmsg)	;
							 }
						 }
						 else
						 {
							 // cout << pCBSNP->geno1[i] << " "<< pCBSNP->geno2[i] <<endl; //debug
							 if(!pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if homozygote1  00
								 ofs <<" "<< "A" << " "<<"A";
							 else if( pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if homozygote2 11
								 ofs <<" " << "B" << " "<<"B";
							 else if( !pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if heterozygote
							 {
								 //ofs << " "<<pCBSNP->allele1 << " "<<pCBSNP->allele2;
								 if(pCBSNP->aOrder[i])
									 ofs <<" "<< "A" << " " <<"B";
								 if(!pCBSNP->aOrder[i])
									 ofs<<" " <<"B"<< " " <<"A";
							 }
							 else if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if missing
								 ofs <<" "<<CBSNP::_miss_gvalue<< " "<<CBSNP::_miss_gvalue;
							 else
							 {
								 error("Problem in writing out "+change_int_into_string(j) +"th SNP.\n");
							 }
						 }
					 }
				 }
							//-----------------------------------------------------------------//
				 ofs << "\n"; //end of line
				 //++ppedInfo; // make an increament to
			 }	//if quality pedinfo
		 } //end of  pedinfo.size for loop

		ofs.clear(); ofs.close();
	   	printLIN("\t**->ped  file has been written out and saved as  \""+\
		 	  			ped_fileName+"\".\n");


}
//

void CSMARTPCAPED::write_smartpca_pedind_file(const vector<CBPED*>&pedInfo,const string& outFileName)
{
	ofstream ofs ;
	ofs.clear(); ofs.close();
	const string pedind_fileName =outFileName+".pedind";
	ofs.open(pedind_fileName.c_str(),ios::out);
	 if(!ofs)
			 error("out file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");
	 CBPED* ppedInfo= pedInfo[0];
	 for(unsigned int i=0; i<pedInfo.size(); ++i)
	 {
		 ppedInfo =pedInfo[i];
		 if(ppedInfo->quality)
		 {
			 if((ppedInfo->famId.size() +ppedInfo->indId.size())>=39)

				 ofs<< i+1<< " "<<ppedInfo->indId <<" ";
			 else
				 ofs<< ppedInfo->famId<< " "<<ppedInfo->indId<< " ";
			 if(ppedInfo->patId!="")ofs<<ppedInfo->patId << " ";
				else ofs<< 0<<" ";
			 if(ppedInfo->matId!="") ofs <<ppedInfo->matId <<" ";
			 else ofs<< 0 << " ";
			 if(ppedInfo->sex) ofs<<"1"<< " "; // 1 is male
			 else if(!ppedInfo->sex && !ppedInfo->miss_sex) ofs<<"2"<<" "; // 2 is female
			 else if(ppedInfo->miss_sex) ofs<<"0"<< " "; // 0 is missing
			 if(ppedInfo->groupLabel!="")
				 ofs<<ppedInfo->groupLabel <<" ";
			 else
			 {
				 if(ppedInfo->pheno)
					 ofs<<"2"<< " ";
				 else if(!ppedInfo->pheno && !ppedInfo->miss_pheno)
					 ofs<<"1"<< " ";
				 else if(!ppedInfo->pheno && ppedInfo->miss_pheno)
					 ofs<<"99"<< " ";
			 }
			 ofs<<"\n"; //end of line
		 }
		}
	ofs.clear(); ofs.close();
   	printLIN("\t**->pedind  file has been written out and saved as  \""+\
	 	  			pedind_fileName+"\".\n");

}
//writing command file and necessary instruction.
void CSMARTPCASNP::write_smartpca_command_file(const string& outFileName )
{
	ofstream ofs, ofls, ofrs;
	ofs.clear(); ofs.close();
	string command_fileName = outFileName+"_smartpca.par";
	//string linux_command_fileName = outFileName+"_LinuxCommands.txt";
	//string r_command_fileName = outFileName+"_Rscript.r";
	//string info_command_fileName = outFileName+"_commandInfo.txt";
	ofs.open(command_fileName.c_str(),ios::out);
	 if(!ofs)
		 error("out file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");
	 ofs<<"#Lines(part of lines) starting with symbol \"#\" are comments. smartpca accepts these lines as comments and ignores when reading the script.\n";
	ofs<<"#Command to run PCA in linux command line:\t ./smartpca -p "+command_fileName+"\n";
	ofs<<"#A detailed information is given in the website: http://www.hsph.harvard.edu/alkes-price/software/\n";
	ofs <<"#Comments in the scripts are also copied form http://helix.nih.gov/Applications/README.popgen.\n";
	ofs<<"genotypename:\t\t "+ outFileName+".ped\n";
	ofs<<"snpname:\t\t "+outFileName+ ".pedsnp\n";
	ofs<<"indivname:\t\t "+ outFileName+".pedind\n";
	ofs<<"evecoutname:\t\t "+outFileName+".evec  #output file of eigenvectors.\n";
	ofs<<"evaloutname:\t\t "+outFileName+".eval  #output file of all eigenvalues\n";
	ofs<<"#optional commands:\n";
	ofs<<"numoutevec:\t\t 10  #number of eigenvectors to output.  Default is 10.\n";
	ofs <<"numoutlieriter:\t\t 5 #maximum number of outlier removal iterations.  Default is 5.  To turn off outlier removal, set this parameter to 0.\n";
	ofs<< "numoutlierevec:\t\t 10 #number of principal components along which to remove outliers during each outlier removal iteration.  Default is 10.\n";
	ofs<<"nsnpldregress:\t\t  0 #If set to a positive integer, then LD correction is turned on,  and input to PCA will be the residual of a regression involving that many  previous SNPs, according to physical location.  See Patterson et al. 2006.   Default is 0 (no LD correction).  If desiring LD correction,  recommend value is: 2.\n";
	ofs<<"outliersigmathresh:\t\t 6.0 #number of standard deviations which an individual must exceed, along one of the top (numoutlierevec) principal components, in order for that individual to be removed as an outlier.  Default is 6.0.\n";
    ofs<< "outlieroutname:\t\t "+outFileName+"_smartpca.log #output logfile of outlier individuals removed. If not specified,  smartpca will print this information to stdout, which is the default.\n";
	ofs<<"usenorm:\t\t YES #Whether to normalize each SNP by a quantity related to allele freq.  Default is YES.  (When analyzing microsatellite data, should be set to NO.\n";
	ofs<<"#outputformat:\t\t EIGENSTRAT #other format types are: ANCESTRYMAP,   PED, PACKEDPED or PACKEDANCESTRYMAP \n";
	ofs.clear(); ofs.close();
	printLIN("\t**->smartpca command parameter  file has been written out and saved as  \""+\
			command_fileName+"\".\n");
	printLIN("\t**->For necessary commands, see \""+outFileName+"_commandInfo.txt\". \n");

	//------------------------------------------------------------------------------------------------------------//
	//next writing  par file for eigenstgrat
	//------------------------------------------------------------------------------------------------------------//
	command_fileName = outFileName+"_eigenstrat.par";
	ofs.open(command_fileName.c_str(),ios::out);
	if(!ofs)
		error("out file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");
	ofs<<"#Lines(part of lines) starting with symbol \"#\" are comments. smartpca and eigenstrat accept these lines as comments and ignores when reading the script.\n";
	ofs<<"#Command to run PCA in linux command line:\t ./smarteigenstrat -p "+command_fileName+"\n";
	ofs<<"#A detailed information is given in the website: http://www.hsph.harvard.edu/alkes-price/software/\n";
	ofs<<"#Comments in the scripts are according as in the readme file of eigenstrat. \n";
	ofs<<"genotypename:		"+outFileName+".ped\n";
	ofs<<"snpname:		"+outFileName+".pedsnp\n";
	ofs<<"indivname:		"+outFileName+".pedind\n";
	ofs<<"pcaname:		"+outFileName+".pca  #modified file of \""+outFileName+".evec\". This file contains eigenvectors. In additional the first line of the file should have the number of eigenvectors  contained in the file.\n";
	ofs<<"outputname:   "+outFileName+"_eigenstrat.txt  #name of output file of chisq association statistics.\n";
	ofs<<"qtmode:			YES #This implies whether group label are given in quantative way or in expression.If quantative , then YES. \n";
	ofs<<"numpc:			10 #Number of eigenvectors to be analysed.\n";
	ofs.clear();ofs.close();
	printLIN("\t**->eigenstrat parameter  file has been written out and saved as  \""+\
			command_fileName+"\".\n");
	//------------------------------------------------------------------------------------------------------------//
	//next writing Linux commands
	//------------------------------------------------------------------------------------------------------------//
	ofls.clear(); ofls.close();
	command_fileName = outFileName+"_LinuxCommands.txt";
	ofls.open(command_fileName.c_str(),ios::out);
		 if(!ofls)
			 error("out file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");
		 ofls<<"#!/bin/bash\n";
		 ofls<<"#Command to run SMARTPCA:\n";
		 ofls <<"#To use this script, make sure that program SMARTPCA,ploteig,smarteigenstrat and twstats  are executable and saved respectively with these names in the same folder from where you run fcGENE.\n";
		 ofls<<"#If this is not the case, you should give its  path and name while excuting above commands.\n";
		 ofls<<"#If necessary, you can also edit (add/delete) the necessary command options in file: "+outFileName+"_smartpca.par\n";
		 ofls<<"./smartpca -p "+outFileName+"_smartpca.par>"+outFileName+"_smartpca_screenOutput.txt \n";
		 ofls<<"#construct a pca plot using ploteig\n";
		 ofls<<"./ploteig -i "+outFileName+".evec -c 1:2 -p Case:Control -x .pdf -o "+outFileName+".xtxt \n";
		 ofls<<"#create an  *.pca file necessary for eigenstrat \n";
		 ofls<<"#Before running the following command, please replace k with  the number of principal components in example.evec file (e.g. 10)\n";
		 ofls<<"cols=$(sed -n 2p smartpca/example.evec |awk '{print  NF  }') \n";
		 ofls<<"k=$(echo $cols-2 | bc) \n";
		 ofls<<"./evec2pca.perl $k "+outFileName+".evec "+outFileName+".pedind "+outFileName+".pca\n";
		 ofls<<"#run eigenstrat using the following command \n";
		 ofls<<"./smarteigenstart -p "+outFileName+"_eigenstrat.par>"+outFileName+"_eigenstrat_screenOutput.txt \n";
		 ofls<<"#To calculate  Tracy-Widom statistics to evaluate the statistical significance of each principal component identified by pca (Patterson et al. 2006): \n";
		 ofls<<"#use twstats program.\n";
		 ofls<<"#Make sure that twtable is contained in current folder. \n";
		 ofls<<"./twstats -t twtable -i "+outFileName+".eval -o "+outFileName+"_twstats.txt \n";
		 ofls.clear();ofls.close();
		 printLIN("\t**->A Linux script file has also been saved as "+command_fileName+"\n");
	//------------------------------------------------------------------------------------------------------------//
	//next writing R script
	//------------------------------------------------------------------------------------------------------------//
		 ofrs.clear(); ofrs.close();
		 command_fileName = outFileName+"_Rscript.r";
		 ofrs.open(command_fileName.c_str(),ios::out);
		 if(!ofrs)
			 error("out file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");
		 ofrs<<"#Command to run SMARTPCA:\n";
		 ofrs <<"#To use this script, make sure that program SMARTPCA  is executable and saved with the name \"smartpca\" in the same folder from where you run fcGENE.\n";
		 ofrs<<"#If this is not the case, you should give its  path and name while excuting above commands.\n";
		 ofrs<<"#If necessary, you can also edit (add/delete) the necessary command options in file: "+outFileName+"_smartpca.par\n";
		 ofrs<<"system(\"./smartpca -p "+outFileName+"_smartpca.par>"+outFileName+"_smartpca_screenOutput.txt\")\n";
		 ofrs <<"#R function to make pca plot: \n";
		 //ofrs<<"#=====================\n";
		 ofrs<<"pca_plot<-function(\n";
		 ofrs<<"		evecFileName, 	#name of *.evec file obtained from eigenstrat\n";
		 ofrs<<"		outFileName	=\"\", #name of plot file\n";
		 ofrs<<"		used_evecs=1:2, #which eigenvectors to plot.  Use 1:2 for top two eigenvectors. \n";
		 ofrs<<"		file_type=\"pdf\" #type of  file to be plotted.You can also write  \"png\".\n";
		 ofrs<<"	 	){\n";
		 ofrs<<" 		if(outFileName==\"\"){\n";
		 ofrs<<" 			outFileName=substring(evecFileName,1,(nchar(evecFileName)-nchar(\".evec\")))\n";
		 ofrs<<"	 		plotName<-paste(outFileName, \".\",file_type,sep=\"\")\n";
		 ofrs<<"		}\n";
		 ofrs<<"	 	evec<-read.table(evecFileName)\n";
		 ofrs<<"	 	ncols<-dim(evec)[2]\n";
		 ofrs<<"		evec<-evec[c(used_evecs+1,ncols)]\n";
		 ofrs<<"	 	colnames(evec)<-c(\"pca1\",\"pca2\",\"groupName\")\n";
		 //ofrs<<"	 	pids<-as.character(evec$names)\n";
		 ofrs<<"	 	groupNames<-evec$groupName\n";
		 ofrs<<"	 	color<-1:length(unique(groupNames))\n";
		 ofrs<<"	 	evec$color<-color[groupNames]\n";
		 ofrs<<"	 	evec$pch<-substring(groupNames,1,5)\n";
		 ofrs<<"		if(file_type==\"pdf\")\n";
		 ofrs<<"	 		pdf(plotName)\n";
		 ofrs<<"		else if(file_type==\"png\")\n";
		 ofrs<<"	 	 	 png(plotName)\n";
		 ofrs<<"		xvalue<-as.numeric(evec$pca1)\n";
		 ofrs<<"	 	yvalue<- as.numeric(evec$pca2)\n";
		 ofrs<<"	 	plot(xvalue,yvalue , pch=\"\",xlab=\"first_component\",ylab=\"second_component\")\n";
		 ofrs<<"	 	text(xvalue,yvalue, labels=evec$pch,col=evec$color)\n";
		 ofrs<<"		dev.off()\n";
		 ofrs<<"	 	cat(\"pca plot can be found in file:\\t\",plotName,\"\\n\")\n";
		 ofrs<<"}\n";
		 ofrs<<"\n";
		 ofrs<<"evecFileName<-\""+outFileName+".evec\"\n";
		 ofrs<<"pca_plot(evecFileName)\n";
		 ofrs<<"#To write file\""+outFileName+".pca\" necessary for eigenstrat: \n";
		 ofrs<<"evec<-as.matrix(read.table(evecFileName,comment.char=\"\",fill=T,header=F))\n";
		 ofrs<<"evec<-evec[,c(-1,-dim(evec)[2])]\n";
		 ofrs<<"outFileName<-\""+outFileName+".pca\"\n";
		 ofrs<<"sink(outFileName) #this command is to open file for writing!\n";
		 ofrs<<"cat(dim(evec)[2],\"\n\")#writes number of eigenvectors in the file \n";
		 ofrs<<"writeLines(evec[1,])# write eigenvalues given in *.evec file \n";
		 ofrs<<"write.table(evec[-1,],file=outFileName,append=T,quot=F,col.names=F,row.names=F,sep=\" \" )\n";
		 ofrs<<"sink() #file closing part \n";
		 ofrs<<"unlist(outFileName)#file closing part\n";
		 ofrs.clear(); ofrs.close();
		 printLIN("\t**->A R-script file has also been saved as "+command_fileName+"\n");
	//------------------------------------------------------------------------------------------------------------//
	//next writing command info
	//------------------------------------------------------------------------------------------------------------//
	command_fileName=outFileName+"_commandInfo.txt";
	ofs.open(command_fileName.c_str(),ios::out);
	 if(!ofs)
		 error("out file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");
	 ofs<<"To make PCA analysis, the following command steps can be recommended. \n";
	 ofs<<"=====================\n";
	 ofs <<"step1:  	run program smartpca using file \""+outFileName+"_smartpca.par\" \n";
	 ofs <<"step2: 		create a plot of the first two (or any specified pair) principal components \n";
	 ofs <<"step3: 		create *.pca file from file \""+outFileName+".evec\"\n";
	 ofs <<"step4: 		run program eigenstrat using file "+outFileName+"_eigenstrat.par\n";
	 ofs <<"step5: 		run program twstats \n";
	 ofs<<"=====================\n";
	 printLIN("\t**->A commandInfo file has also been saved as "+command_fileName+"\n");
	 ofs.clear();ofs.close();


	/*command_fileName=outFileName+"_eigenstrat.par";
	ofs.open(command_fileName.c_str(),ios::out);
	ofs<<"#Command to run convertf in Linux command line:\t ./convertf -p "+command_fileName+"\n";
	ofs <<"#This command will create files for the program eigenstrat.\n";
	ofs <<"#For more information, please visit the link: http://helix.nih.gov/Applications/README.convertf";
	ofs<<"genotypename:\t\t "+ outFileName+".ped\n";
	ofs<<"snpname:\t\t "+outFileName+ ".pedsnp\n";
	ofs<<"indivname:\t\t "+ outFileName+".pedind\n";
	ofs<<"outputformat:\t\t EIGENSTRAT\n";
	ofs<<"genotypeoutname:\t\t "+outFileName+".eigenstratgeno\n";
	ofs<<"snpoutname:\t\t "+outFileName+".snp\n";
	ofs<<"indivoutname:\t\t "+outFileName+".ind\n";
	ofs<<"familynames:\t\t NO\n";
	ofs.clear(); ofs.close();
	printLIN("\t**->command  file for convertf to create files necessary for eigenstrat, has been written out and saved as  \""+\
			command_fileName+"\".\n");
	printLIN("\t**->Command to run program convertf in linux command line:\t ./convertf -p "+command_fileName+"\n");
	*/
}
// read group label file
void CSMARTPCAPED::read_smartpca_groupLabel_file(vector<CBPED*>& pedInfo, const string& inFileName)
{
	 string r_msg	="\n*->Reading file \""+inFileName+ "\":\n";
	 printLIN(r_msg);
	 string line	="";
	 string buffer	="";
	 rFile file;	  file.clear();file.close(); 	// check if it has been already opend.
	 static int extraMapFileLineZaehler	=0;			 // extra map file line zaehler
	 //ios_base::iostate i= file.rdstate();
	 file.open(inFileName.c_str());
	 // this will check if each row has same number of columns.
	 if(!file)
	 {
			file.setstate(ios_base::failbit);
			file.clear();
			error("File "+inFileName+" either does not exit or could not be opend.\n");
	 }
	 // check first header and  the number of columns
	 file.clear();
	 do{
		 getline(file,line, '\n');
		 if(line[0]=='#') // leave the comments;
			 continue;
		 if(file.eofbit) break;
	 }while(line.empty());
	 vector<string> header;
	 header.clear();
	 stringstream line_parser(line);
	 if(line_parser.good())
	 {
		 while(line_parser>>buffer)
			header.push_back(buffer);
		 ++extraMapFileLineZaehler;
	 }
	 else
	 {
		 string msg= "Problem in parsing the first line!\n";
		 error(msg);
	 }
	 //end of parsing first line of the file
	 const unsigned int nCols 	=2;
	 unsigned int curCols 		=0;
	 static int no_of_updates	=0;
	 bool temp_tfvec			=(header.size()==nCols); //&& (header[0]=="indid" || header[0]=="indId" ||header[0]!="indID") &&( header[1]=="group-label"||header[1]=="pop-label");
	 bool found0				=false;
	 unsigned int pos_indId		=0;
	  vector<CBPED*>::const_iterator it_pCBPED=pedInfo.begin();
	 vector<string> temp_vec_indId;
	 vector<string>::iterator it =temp_vec_indId.begin();
	vector<string> sort_temp_vec_indId; //=temp_vec_indId;
	if(temp_tfvec)
	{

 		for(it_pCBPED=pedInfo.begin();it_pCBPED!=pedInfo.end();++it_pCBPED)
			temp_vec_indId.push_back((*it_pCBPED)->indId);
		it_pCBPED=pedInfo.begin();
		sort_temp_vec_indId=temp_vec_indId;
		sort(sort_temp_vec_indId.begin(),sort_temp_vec_indId.end());
		r_msg="\t**->Updating group/population label: \n";
		printLIN(r_msg);
		found0	=binary_search(sort_temp_vec_indId.begin(),sort_temp_vec_indId.end(),header[0]);
		if(found0)
		{
			it=find(temp_vec_indId.begin(),temp_vec_indId.end(),header[0]);
			pos_indId  =(int)(distance(temp_vec_indId.begin(),it));
			it_pCBPED  =pedInfo.begin()+pos_indId;
			(*it_pCBPED)->groupLabel=header[1];
			it_pCBPED	=pedInfo.end();
			++no_of_updates;
		}
	}
	else
	{
		string msg="Problem in reading file \""+inFileName+"\"!\n  File used to update group label or population label should contain two columns having individual ids in the first column and  the grouping-label in the second column.\n";
		error(msg);
	}
	// till now we have just worked for the first line of group label specifying file.
	// now for the second line :
	while(getline(file,line, '\n')) // until eof hit
	{
		//---------------------//
		++extraMapFileLineZaehler;
		if(line.length()==0)
			continue;
			if(line[0]=='#')
				continue;
		//if(!((*it_pCBSNP)->quality))
		//	continue;
		header.resize(0);
		header.clear();
		stringstream line_parser1(line);
		if(line_parser1.good())
		{
			while(line_parser1>>buffer)
			{
				header.push_back(buffer);
				//cout << buffer << " "; // only for test
			}
			curCols=header.size();
			if(header.size()!=nCols)
			{
				string msg= "file \""+inFileName+"\"  should contain  "+change_int_into_string(nCols)+\
						"columns  but in the "+change_int_into_string(extraMapFileLineZaehler)+\
						"th line of the file, there are "+change_int_into_string(curCols)+" columns.\n";
				curCols=0;
				error(msg);
			}
			found0 	=binary_search(sort_temp_vec_indId.begin(),sort_temp_vec_indId.end(),header[0]);
			if(found0)
			{
				it=find(temp_vec_indId.begin(),temp_vec_indId.end(),header[0]);
				pos_indId  =(int)(distance(temp_vec_indId.begin(),it));
				it_pCBPED  =pedInfo.begin()+pos_indId;
				(*it_pCBPED)->groupLabel=header[1];
				it_pCBPED	=pedInfo.end();
				++no_of_updates;
			}

		} //end of if line parser is good();
		else
		{
			r_msg= "Parsing"+ change_int_into_string(extraMapFileLineZaehler)+"th line was not possible.\nPlease contact the software developer.\n";
			error(r_msg);
		}
	}//end of while getline
		//check if first
	r_msg ="\t**->Total no of individuals whose grouping/population-label is updated: "+change_int_into_string(no_of_updates)+".\n";
	printLIN(r_msg);

}//end of function


/* 	*we can also create for
    ind
    snp
    eigenstrat.geno file
    EIGENSTRAT format: used by eigenstrat program
  genotype file: see example.eigenstratgeno
  snp file:      see example.snp (same as above)
  indiv file:    see example.ind (same as above)
	Note that
	The genotype file contains 1 line per SNP.
  Each line contains 1 character per individual:
  0 means zero copies of reference allele.
  1 means one copy of reference allele.
  2 means two copies of reference allele.
  9 means missing data.
  ind:
  pid   grouplabel #no header
  iid1  case
  iid2 control
  snp : #no header
   rs7288834    22        0.020000        16212142 		G	 T
   rs16978746    22        0.130000        20278224	 	T	 X # x means missing
   rs5754387    22        0.130000        20304703 		G	 C
  eigenstrat:
  12121100122111111121
  22222222222292222222
  22221211211222221222
  10112222212222222211
 */
