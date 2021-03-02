/*
 * rformat.cpp
 *
 *  Created on: 09.01.2013
 *      Author: nroshyar
 */
#include "rformat.h"

void CRPED::write_affection_status_file(const vector<CBPED*>& vPed, const string & outFileName){
	//cout << "\n------------------------------------------------\n";
	//cout<< "Hi from write_affection_status_file \n";
	//printLIN("*->Writing files in "+Bpar::oFormat+" format: \n");
	string affstatFileName = outFileName+"_affstat.txt";
	bool pheno_missing=true;
	ofstream ofs ;
	ofs.clear(); ofs.close();
	ofs.open(affstatFileName.c_str());
	if(!ofs)
	{
		error("out file "+ affstatFileName+ " can not be written out.  Please  check whether you have given the correct file Path.\n");
	}
	//ofs <<"test\n";
	ofs << "SMAPLE_ID" << " " <<   "AFFSTAT" <<"\n"; // first line as header

	for(unsigned int i=0; i<vPed.size(); ++i)
	{
		if(vPed[i]->quality)
		{
			if(vPed[i]->indId!="") ofs << vPed[i]->indId << " ";
			else ofs << "indv"<<i+1<< " ";
			if(vPed[i]->pheno)
				ofs<< "1";
			else if(!vPed[i]->pheno&& !vPed[i]->miss_pheno)
				ofs<< "0";
			else if(!vPed[i]->pheno && vPed[i]->miss_pheno)
			{
				ofs<< "NA";
			}
			if(pheno_missing)
			{
				string sout ="\t**->WARNING: Your genotyped data contains individuals with missing phenotype information.\n";
				printLIN(sout);
				sout ="\t It  may affect your analysis.\n";
				printLIN(sout);
			}
			pheno_missing=false ;

			ofs <<"\n";	//end of the line
		}
	}
	ofs.clear();ofs.close();
//	cout << "\n------------------------------------------------\n";
	printLIN("\t**-> Affection status file has been written out and saved as \""+affstatFileName+"\".\n");

}
//
//xxxxxxxxxxxxxxxxxxxxxx
void CRPED::write_allele_dose_file( const vector<CBPED*>& vPed, const vector<CBSNP*>& genInfo, const string & outFileName, const string& ref_allele)
{

	//cout << "\n------------------------------------------------\n";
	//cout << "\n Hi ! from write_genotype_file!\n";
	//printLIN("*->Writing files in "+Bpar::oFormat+" format: \n");
	string genotypeFileName = outFileName+"_dose.txt";
	//string snpinfoFileName  = outFileName+"_genotype.txt";
	//---------------------------------------------//
	ofstream ofs,ofs1 ;
	ofs.clear(); ofs.close();
	ofs.open(genotypeFileName.c_str());
	if(!ofs)
	{
			error("out file "+ genotypeFileName+ " can not be written out.  Please  check whether you have given the correct file Path.\n");
	}
	//ofs<< "test \n";
	// write first the header
	ofs << "SMAPLE_ID" <<" ";
	CBSNP* pCBSNP = genInfo[0];
	for(unsigned int j=0; j<genInfo.size(); ++j)
	{
		//first change genotype into 0123
		pCBSNP=genInfo[j];
		//if(pCBSNP->change_0123)
		//CBSNP::convert_snp_genotype_into_0123(pCBSNP); //convert genotype into 0123
		 //pCBSNP->change_0123=false;
		CBSNP::convert_hardcalls_into_genoprob(pCBSNP);
		CBSNP::convert_snp_genotype_into_dose( pCBSNP,ref_allele);

		 if(pCBSNP->quality)
		 {
			 if(pCBSNP->rsId!="")
				ofs<<pCBSNP->rsId<<" ";
			 else
				ofs<<pCBSNP->snpName<<" ";
		 }
	}
	ofs << "\n";
	//writing genotypes
	for(unsigned int i=0; i<vPed.size(); ++i)
	{
		if(vPed[i]->quality)
		{
			if(vPed[i]->indId!="") ofs << vPed[i]->indId;
					else ofs << "indv"<<i+1;
			//ofs <<"\n"; //debug
			for(unsigned int j=0;j<genInfo.size(); ++j)
			{
				pCBSNP=genInfo[j];
				if(pCBSNP->quality) // default value of quality is true;
				{

					  pCBSNP =genInfo[j];
					  if(!(pCBSNP->geno_dose[i]+1.0))
						  ofs<< " "<<"NA";
					  else
						  ofs<<" "<< pCBSNP->geno_dose[i];
				} // end of pCBSNP->quality
				//cout << (int)(pCBSNP-genInfo.begin())<<endl;

			}
			ofs << "\n"; // end of the line
			pCBSNP =genInfo[0];
		}
	}
	//---------------------------------------------//
	ofs.clear(); ofs.close();
	printLIN("\t**->Genotype dose file has been written out and saved as \""+genotypeFileName+"\".\n"
				"\tGenotype dose \"0.0\" means 0 count (distribution: 0 0 1) of \"ref_allele\" in alleleInfo file, i.e. homozygote alternative allele \n"
				"\tGenotype \"1.0\" 1 count (distribution: 0 1 0) of \"ref_allele\" in alleleInfo file, i.e. heterozygote\n"
				"\tGenotype \"2.0\" means 2 count (distribution: 1 0 0) of \"ref_allele\" in alleleInfo file, i.e. homozygote reference allele \n"
			);
	//cout << "\n------------------------------------------------\n";

 }

//xxxxxxxxxxxxxxxxx

void CRSNP::write_ref_and_alternative_allele_info_file(const vector<CBSNP*>& genInfo, const string & outFileName, const string & ref_allele)
{
	string alleleFileName	= outFileName+"_alleleInfo.txt";
	printLIN("\t**->Writing file:  \""+alleleFileName+"\".\n");
	ofstream  aInfo;
	aInfo.clear(); aInfo.close();
	aInfo.open(alleleFileName.c_str());
	if(!aInfo)
	{
				error("out file "+ alleleFileName+ " can not be written out.  Please  check whether you have given the correct file Path.\n");
	}
	aInfo<<"nchr\t"<<"rsid\t"<<"ref_allele\t"<<"alternative_allele";
	aInfo<<"\n";
	CBSNP* pCBSNP = genInfo[0];
	for(unsigned int j=0; j<genInfo.size(); ++j)
	{
		//first change genotype into 0123
			pCBSNP=genInfo[j];
		 if(pCBSNP->quality)
		 {
			 //write nchr
			 aInfo<<pCBSNP->nchr;
			 //write rsid
			 if(pCBSNP->rsId!="")
				 aInfo<<pCBSNP->rsId<<" ";
			 else
				 aInfo<<pCBSNP->snpName<<" ";
			 	 //write for allele info file
			 	if((ref_allele=="")||(ref_allele=="minor-allele")) //default is minor so;
				{
					if(pCBSNP->min_allele!="")
						aInfo<<"\t"<<pCBSNP->min_allele;
					else
						aInfo<<"\t"<<".";
					if(pCBSNP->maj_allele!="")
						aInfo<<"\t"<<pCBSNP->maj_allele;
					else
						aInfo<<"\t"<<".";
				}else // if reference allelel is given as major
				if(ref_allele=="major-allele")
				{
					if(pCBSNP->maj_allele!="")
						aInfo<<"\t"<<pCBSNP->maj_allele;
					else
						aInfo<<"\t"<<".";
					if(pCBSNP->min_allele!="")
						aInfo<<"\t"<<pCBSNP->min_allele;
					else
						aInfo<<"\t"<<".";
				}else
				if(ref_allele=="allele1")
				{
					if(pCBSNP->allele1!="")
						aInfo<<"\t"<<pCBSNP->allele1;
					else
						aInfo<<"\t"<<".";
					if(pCBSNP->allele2!="")
						aInfo<<"\t"<<pCBSNP->allele2;
					else
						aInfo<<"\t"<<".";

				}else
				{
					if(pCBSNP->allele2!="")
						aInfo<<"\t"<<pCBSNP->allele2;
					else
						aInfo<<"\t"<<".";
					if(pCBSNP->allele1!="")
						aInfo<<"\t"<<pCBSNP->allele1;
					else
						aInfo<<"\t"<<".";

				}
				//write end of the line
				aInfo<<"\n";
		 }

	}
	//close ainfo file
	aInfo.clear();
	aInfo.close();
	printLIN("\t**->Allele Info file has been written out and saved as \""+alleleFileName+"\".\n"
			"\t The allele specified as \"ref_allele\" in this file, is used as the counts of alleles, or as the allele-dose at any particular SNP.\n"
			"\t For example if we convert the genotypes as allele counts (0,1,2), then homozygote of ref_allele is represented  by 2\n"
			"\t and homozygote of alternative allele is represented by 0. This case is similar with the allele-dose file\n"
			);
}
void CRPED::write_genotype_file( const vector<CBPED*>& vPed, const vector<CBSNP*>& genInfo, const string & outFileName, const string & ref_allele)
{
	//cout << "\n------------------------------------------------\n";
	//cout << "\n Hi ! from write_genotype_file!\n";
	//printLIN("*->Writing files in "+Bpar::oFormat+" format: \n");
	string genotypeFileName = outFileName+"_genotype.txt";
	//string snpinfoFileName  = outFileName+"_genotype.txt";
	//---------------------------------------------//
	ofstream ofs;
	ofs.clear(); ofs.close();
	ofs.open(genotypeFileName.c_str());

	if(!ofs)
	{
			error("out file "+ genotypeFileName+ " can not be written out.  Please  check whether you have given the correct file Path.\n");
	}
	//ofs<< "test \n";
	// write first the header
	ofs << "SAMPLE_ID" <<" ";
	CBSNP* pCBSNP = genInfo[0];
	for(unsigned int j=0; j<genInfo.size(); ++j)
	{
		//first change genotype into 0123
		pCBSNP=genInfo[j];
		 if(pCBSNP->quality)
		 {
			 if(pCBSNP->rsId!="")
				ofs<<pCBSNP->rsId<<" ";
			 else
				ofs<<pCBSNP->snpName<<" ";
		 }

	}
	//close ainfo file
	ofs << "\n";
	//writing genotypes
	for(unsigned int i=0; i<vPed.size(); ++i)
	{
		if(vPed[i]->quality)
		{
			if(vPed[i]->indId!="") ofs << vPed[i]->indId;
					else ofs << "indv"<<i+1;
			//ofs <<"\n"; //debug
			for(unsigned int j=0;j<genInfo.size(); ++j)
			{
				pCBSNP=genInfo[j];
				if(pCBSNP->quality) // default value of quality is true;
				{
					//here starts
					// cout<< pCBSNP->geno_0123[i]<<", ";
					if(pCBSNP->geno_0123[i]==3)
					{
						ofs << " "<< "NA";
					}
					else if (pCBSNP->geno_0123[i]==1)
					{
						ofs <<" "<<"1";
					}
					else
					{
						// for default
						 //cout<< pCBSNP->geno_0123[i]<<", ";
						 //cout <<"ref_allele: "<<ref_allele <<",  "<< boolalpha <<"value: "<< ((ref_allele=="")||(ref_allele=="minor_allele"))<<endl;
						 //cout <<"ref_allele: "<<ref_allele <<",  "<< boolalpha <<"value: "<< ((ref_allele=="allele1"))<<endl;
						 //cout <<"ref_allele: "<<ref_allele <<",  "<< boolalpha <<"value: "<< ((ref_allele=="allele2"))<<endl;
						//cout <<"ref_allele: "<<ref_allele <<",  "<< boolalpha <<"value: "<< (((ref_allele=="major-allele")||(ref_allele=="major_allele")))<<endl;
						if((ref_allele=="")||(ref_allele=="minor-allele"))
						{
							if ((pCBSNP->freq>0.5) )
								ofs <<" "<< (pCBSNP->geno_0123[i]);
							else if(pCBSNP->freq<=0.5 )
							{
								if (pCBSNP->geno_0123[i]==0)
									ofs << " "<< "2";
								else if(pCBSNP->geno_0123[i]==2)
									ofs << " "<< "0";
							}
						}
						else if((ref_allele=="allele2"))
						{
							ofs <<" "<< (pCBSNP->geno_0123[i]);
						}
						else if((ref_allele=="allele1"))
						{
							if (pCBSNP->geno_0123[i]==0)
								ofs << " "<< "2";
							else if(pCBSNP->geno_0123[i]==2)
								ofs << " "<< "0";
						}
						else if((ref_allele=="major-allele"))
						{
						 //	cout << "(pCBSNP->freq<=0.5)"<< (pCBSNP->freq<=0.5) <<endl;
							if ((pCBSNP->freq<=0.5) )
								ofs <<" "<< (pCBSNP->geno_0123[i]);
							else if(pCBSNP->freq>0.5 )
							{
								if (pCBSNP->geno_0123[i]==0)
									ofs << " "<< "2";
								else if(pCBSNP->geno_0123[i]==2)
									ofs << " "<< "0";
							}
						}

					}
						//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),2,99); // value 99 is used just as a temp coding to swap later .
						//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),0,2);
						//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),99,0);

				} // end of pCBSNP->quality
				//cout << (int)(pCBSNP-genInfo.begin())<<endl;

			}
			ofs << "\n"; // end of the line
			pCBSNP =genInfo[0];
		}
		//cout <<endl;
	}
	//---------------------------------------------//
	ofs.clear(); ofs.close();
	printLIN("\t**->Genotype file has been written out and saved as \""+genotypeFileName+"\".\n"
			"\tGenotype \"0\" means 0 count of \"ref_allele\" in alleleInfo file, i.e. homozygote alternative allele \n"
			"\tGenotype \"1\" means 1 count of \"ref_allele\" in alleleInfo file, i.e. heterozygote\n"
			"\tGenotype \"2\" means 2 count of \"ref_allele\" in alleleInfo file, i.e. homozygote reference allele \n"
		);
	//cout << "\n------------------------------------------------\n";

 }

 // write transposed genotype file 
void CRPED::write_transposed_genotype_file( const vector<CBPED*>& vPed, const vector<CBSNP*>& genInfo, const string & outFileName, const string & ref_allele)
{
	//cout << "\n------------------------------------------------\n";
	//cout << "\n Hi ! from write_genotype_file!\n";
	//printLIN("*->Writing files in "+Bpar::oFormat+" format: \n");
	string genotypeFileName = outFileName+"_genotype.txt";
	//string snpinfoFileName  = outFileName+"_genotype.txt";
	//---------------------------------------------//
	ofstream ofs;
	ofs.clear(); ofs.close();
	CBSNP* pCBSNP = genInfo[0];
	ofs.open(genotypeFileName.c_str());
	if(!ofs)
	{
			error("out file "+ genotypeFileName+ " can not be written out.  Please  check whether you have given the correct file Path.\n");
	}
	ofs << "rsid";
	//write header 
	for(unsigned int i=0;i<vPed.size();++i)
	{	
		if(vPed[i]->quality)
		{
			if(vPed[i]->indId!="")
			ofs <<" "<< vPed[i]->indId;
					else ofs <<" "<< "indv"<<i+1;
			
		}	
	}
	ofs << "\n";
	
	
	//writing genotypes
	for(unsigned int j=0;j<genInfo.size(); ++j)
	{
		pCBSNP=genInfo[j];
				
						
					
		if(pCBSNP->quality) // default value of quality is true;
		{
					//here starts
			if(pCBSNP->rsId!="")
				ofs<<pCBSNP->rsId<<" ";
			else
				ofs<<pCBSNP->snpName<<" ";
			for(unsigned int i=0; i<vPed.size(); ++i)
			{
				if(vPed[i]->quality)
				{
 	
					// cout<< pCBSNP->geno_0123[i]<<", ";
					if(pCBSNP->geno_0123[i]==3)
					{
						ofs << " "<< "NA";
					}
					else if (pCBSNP->geno_0123[i]==1)
					{
						ofs <<" "<<"1";
					}
					else
					{
						// for default
						 //cout<< pCBSNP->geno_0123[i]<<", ";
						 //cout <<"ref_allele: "<<ref_allele <<",  "<< boolalpha <<"value: "<< ((ref_allele=="")||(ref_allele=="minor_allele"))<<endl;
						 //cout <<"ref_allele: "<<ref_allele <<",  "<< boolalpha <<"value: "<< ((ref_allele=="allele1"))<<endl;
						 //cout <<"ref_allele: "<<ref_allele <<",  "<< boolalpha <<"value: "<< ((ref_allele=="allele2"))<<endl;
						//cout <<"ref_allele: "<<ref_allele <<",  "<< boolalpha <<"value: "<< (((ref_allele=="major-allele")||(ref_allele=="major_allele")))<<endl;
						if((ref_allele=="")||(ref_allele=="minor-allele"))
						{
							if ((pCBSNP->freq>0.5) )
								ofs <<" "<< (pCBSNP->geno_0123[i]);
							else if(pCBSNP->freq<=0.5 )
							{
								if (pCBSNP->geno_0123[i]==0)
									ofs << " "<< "2";
								else if(pCBSNP->geno_0123[i]==2)
									ofs << " "<< "0";
							}
						}
						else if((ref_allele=="allele2"))
						{
							ofs <<" "<< (pCBSNP->geno_0123[i]);
						}
						else if((ref_allele=="allele1"))
						{
							if (pCBSNP->geno_0123[i]==0)
								ofs << " "<< "2";
							else if(pCBSNP->geno_0123[i]==2)
								ofs << " "<< "0";
						}
						else if((ref_allele=="major-allele"))
						{
						 //	cout << "(pCBSNP->freq<=0.5)"<< (pCBSNP->freq<=0.5) <<endl;
							if ((pCBSNP->freq<=0.5) )
								ofs <<" "<< (pCBSNP->geno_0123[i]);
							else if(pCBSNP->freq>0.5 )
							{
								if (pCBSNP->geno_0123[i]==0)
									ofs << " "<< "2";
								else if(pCBSNP->geno_0123[i]==2)
									ofs << " "<< "0";
							}
						}

					}
						//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),2,99); // value 99 is used just as a temp coding to swap later .
						//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),0,2);
						//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),99,0);

				} // end of pCBSNP->quality
				//cout << (int)(pCBSNP-genInfo.begin())<<endl;

			}
			ofs << "\n"; // end of the line
		}
		//cout <<endl;
	}
		//pCBSNP =genInfo[0];
	
	//---------------------------------------------//
	ofs.clear(); ofs.close();
	printLIN("\t**->Genotype file has been written out and saved as \""+genotypeFileName+"\".\n"
			"\tGenotype \"0\" means 0 count of \"ref_allele\" in alleleInfo file, i.e. homozygote alternative allele \n"
			"\tGenotype \"1\" means 1 count of \"ref_allele\" in alleleInfo file, i.e. heterozygote\n"
			"\tGenotype \"2\" means 2 count of \"ref_allele\" in alleleInfo file, i.e. homozygote reference allele \n"
		);
	//cout << "\n------------------------------------------------\n";

 }
 
//---------------------------------------------------------
void CRSNP::lese_genotype_file(vector<CBPED*>& vPed, vector<CBSNP*>& genVec, const string & inFileName)
{
	printLIN("*->Reading file \""+ inFileName + "\": \n");

	IFS file(inFileName.c_str(),"r");
	if(!file)
	{
		string msg= "File \""+ inFileName+ "\" either does not exits or could not be opened.\n"	;
			error(msg);
	}
		// if file is okay then define more varaibles
		string elem					=""; // gelesen von ped File;
		string elem1				="";
		int ncol					=0;
		int ccol					=0;
		static int nrow				=0;
		bool zeilenEnde				=false;
		static unsigned int nInd	=0;
		bool	assign_snpid		=true;
		bool	assign_iid			=true;
		string 	rgeno_type			="iid_snp"; // indivxsnp type of matrix

		//unsigned int nPerson	=0;
		while(!_ifs_eof(file))
		{
				//address the end of the file
			if(_ifs_eof(file))
				{
					zeilenEnde=true;
					if(nrow==0)
					{
						string msg=" There are no elements in "+inFileName +" file.\n";
						error(msg);
					}
					else
					{
						cout << "Total no of Individual read: "<< (nrow-1)<< endl; // only for test
						nrow=0;
					}
					break;
				}
			zeilenEnde=false;
			// define  pointer of CBSNP
	//-----------------------------------------------------------------------------------------------------------//
			//reading first row
	//-----------------------------------------------------------------------------------------------------------//



		if(assign_snpid && assign_iid)
		{

			//read first column
			ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
			elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
			if(_ifs_eof(file)) break; // break is okay because we do not want to read only one column
			if(zeilenEnde) break; // break is also okay because we do not want to read only one  column
			if(elem1.empty() )
				do
				{
					ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
					elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol,elem);
					if(zeilenEnde|| _ifs_eof(file))
						break;
				}while(elem1.empty());
			// If we read the first column and then there is the end of file then that is not correct
			//the follwoing line handels this mistake
			// there should be at least two columns with first the header that is why
			if( zeilenEnde || _ifs_eof(file) )
			{
				if(!elem1.empty())
				{
					string msg= " There is only one column in the first row of "+inFileName+" file. \n" ;
					error(msg);
				}
				else
					break; // only in the first column
			}
			//cout << elem1 <<endl;
			bool _b_tmp = (elem1=="SAMPLE_ID")||((elem1=="id")|| (elem1=="iid")||(elem1=="sample_id"));
			//cout <<boolalpha<<_b_tmp<<endl;
			// This case is true except last column
			if(zeilenEnde || _ifs_eof(file))
						{
							string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+inFileName+"\".\n";
							error(msg);
						}
			if(!_b_tmp)
			{
				//If the first column first has not no element having above then it can be that the type of matrix is
				// given as snp x indiv type of matrix. If this is the condition, it should give as the following
				_b_tmp = (elem1=="rsid")||((elem1=="snpid")|| (elem1=="rsID")||(elem1=="SNPID")||(elem1=="snp_id")||(elem1=="snpID")||(elem1=="rs_id")||(elem1=="rs") ||(elem1=="SNP"));
				//cout <<boolalpha<<_b_tmp<<endl;
				if(!_b_tmp)
				{
					string _msg = "The first element of first column must be specified as either words like:\n"
							"SAMPLE_ID, id, iid, sample_id, if every row contains an individual. or it should contains words like:\n"
							"\"rsid\", \"snpid\",\"rsID\",\"snpID\", if each row contains a SNP. \n";
					error(_msg);
				}else{
					//else reach row contains a snp so
					rgeno_type			="snp_iid";
				}

			}

			//read the rest of first row
			if(rgeno_type=="iid_snp")
			{
				//starting reading second column of first row
					//each column gives the rsid or snpid
				while(!zeilenEnde && !_ifs_eof(file)) //reading start from second column to last column
						{

							//define pointer
							elem1="";
							elem="";
							ccol=0;
							ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
							elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);

							if(!elem1.empty())
							{
															//create SNP
								CBSNP* pCBSNP=new CBSNP;
								if(pCBSNP->rsId=="")
								pCBSNP->rsId	=elem1;
								genVec.push_back(pCBSNP);

							}
							if(zeilenEnde)
								break;

						}//end of 2nd column to last column
						// inform the total number of SNPs
					printLIN("\t**->Total number of SNPs found: "+change_int_into_string( (int)genVec.size() )+".\n ");
			}
			else if(rgeno_type=="snp_iid") // if snp x indiv
			{
				//starting reading second column of first row
					//each column gives the rsid or snpid
				while(!zeilenEnde && !_ifs_eof(file)) //reading start from second column to last column
					{
					//define pointer
					elem1="";
					elem="";
					ccol=0;
					ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
					elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
					if(!elem1.empty())
					{
														//create SNP
							CBPED* pCBPED=new CBPED;
							if(pCBPED->indId=="")
							pCBPED->indId	=elem1;
							//set sex to undefined
							pCBPED->sex				=false;
							pCBPED->miss_sex 		=true;
							++pCBPED->nNosex;
							//set status to undefined
							pCBPED->pheno			=false;
							pCBPED->miss_pheno		=true;
							++pCBPED->nMiss_status; //
							vPed.push_back(pCBPED);

					}

					}//end of 2nd column to last column
					// inform the total number of individuals
			printLIN("\t**->Total number of individuals found: "+change_int_into_string( (int)vPed.size() )+".\n ");

			}
			assign_snpid=false;
			assign_iid	=false;
		}
	if(rgeno_type=="iid_snp")
	{
		//-----------------------------------------------------------------------------------------------------------//
			//reading from second row
		//-----------------------------------------------------------------------------------------------------------//
		zeilenEnde =false;

		while(!_ifs_eof(file))
			{  // do all works till the end of line and then break

				//if(zeilenEnde || _ifs_eof(file))
				//read first column
				elem1="";
				elem="";
				ccol=0;
				ncol=0;
				ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
				elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
				if(_ifs_eof(file)) break; // break is okay because we do not want to read only one column
				if(zeilenEnde) break; // break is also okay because we do not want to read only one  column
				if(elem1.empty() )
					do
					{
						ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
						elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol,elem);
						if(zeilenEnde|| _ifs_eof(file))
							break;
					}while(elem1.empty());

				//this is true except for the last element of the line .
				if(zeilenEnde || _ifs_eof(file))
				{
								string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+inFileName+"\".\n";
								error(msg);
				}

				// read family id in the first column
				//cout << ccol <<", " <<ncol<<",elem: "<<elem1<<endl;
				if(ccol!=0 && ncol==1 ) // this condition must always be true ;
				{
					// The first column contains sample ids so save them
					CBPED* pCBPED 			=new CBPED;
					pCBPED->indId			=elem1;
					//set sex to undefined
					pCBPED->sex				=false;
					pCBPED->miss_sex 		=true;
					++pCBPED->nNosex;
					//set status to undefined
					pCBPED->pheno			=false;
					pCBPED->miss_pheno		=true;
					++pCBPED->nMiss_status; //
					vPed.push_back(pCBPED);
					++nInd;
				}
				else
				{
					string msg= "problem in reading indiviudal ids "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+inFileName+"\". \n";
					error(msg);
				}

				// start of reading from the second column till end
				CBSNP * pCBSNP =genVec[0];
				unsigned int i =0;
				while(true)
				{
					elem1="";
					elem="";
					ccol=0;
					if(elem1.empty() )
						do
						{
							ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
							elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol,elem);
							if(zeilenEnde|| _ifs_eof(file))
								break;
						}while(elem1.empty());
					// start saving the genotypes
					 pCBSNP =genVec[i];
					 if(elem1=="2")
					 {
						 // homozygote allele1: false false
						 pCBSNP->geno1.push_back(false);
						 pCBSNP->geno2.push_back(false);
						 pCBSNP->aOrder.push_back(true);
					  }else if(elem1=="1"){
						  // heterozygote major: false false
						  pCBSNP->geno1.push_back(false);
 						  pCBSNP->geno2.push_back(true);
 						  pCBSNP->aOrder.push_back(true);
					  }else if(elem1=="0"){
						  // homozygote allele2: false false
						  pCBSNP->geno1.push_back(true);
						  pCBSNP->geno2.push_back(true);
						  pCBSNP->aOrder.push_back(true);
					  }else
					  {
						  // homozygote major: false false
						  pCBSNP->geno1.push_back(true);
						  pCBSNP->geno2.push_back(false);
						  pCBSNP->aOrder.push_back(true);
					  }
					 if((vPed.size()!=  pCBSNP->geno1.size() )|| (  vPed.size()!=pCBSNP->geno2.size())){
						 string msg ="Internal error, problem in Saving: \n"
								  	  "\tSize of vPed and geno1 and geno2 must be same. Please contact the software provider\n";
						 error(msg);
					 }
					 ++i;
					 elem1="";
					// cout << ncol<<endl;
					if(zeilenEnde|| _ifs_eof(file))
					{
						// check if i has the equal to snp size
						if(genVec.size()!=i)
						{
							string msg ="Internal error, problem in Saving: \n"
									"There are only "+change_int_into_string(i)+" columns in the "+change_int_into_string(nrow)+"th row of the file.\n"
									"Please check if this line has information for each SNP. Otherwise:\n"
									"Hint: Size of genVec and index i  must have be equal. Please contact the software provider\n";
							error(msg);
						}
							ncol		=0;
							ccol		=0;
							break;
					}


				}
				if(_ifs_eof(file))
					break;

			} //end of the line while(!zeilenEnde && !_ifs_eof(file))

		//-----------------------------------------------------------------------------------------------------------//
			//end of reading
		//-----------------------------------------------------------------------------------------------------------//

		if(_ifs_eof(file))
		{
			zeilenEnde=true;
			if(nrow==0)
			{
				string msg=" There are no elements in "+inFileName +" file.\n";
				error(msg);
			}
			else
			{
				printLIN("\t**->Total no of Individuals read: "+ change_int_into_string(vPed.size())+".\n");
				nrow=0;
			}
			nrow=0;
			break;
		}

	}else if(rgeno_type=="snp_iid")
	{
		//-----------------------------------------------------------------------------------------------------------//
				//reading from second row
		//-----------------------------------------------------------------------------------------------------------//

			while(!_ifs_eof(file))
				{  // do all works till the end of line and then break

					//if(zeilenEnde || _ifs_eof(file))
					//read first column
					elem1="";
					elem="";
					ccol=0;
					ncol=0;
					ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
					elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
					if(_ifs_eof(file)) break; // break is okay because we do not want to read only one column
					if(zeilenEnde) break; // break is also okay because we do not want to read only one  column
					if(elem1.empty() )
						do
						{
							ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
							elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol,elem);
							if(zeilenEnde|| _ifs_eof(file))
								break;
						}while(elem1.empty());

					//this is true except for the last element of the line .
					if(zeilenEnde || _ifs_eof(file))
					{
									string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+inFileName+"\".\n";
									error(msg);
					}

					// read family id in the first column
					//cout << ccol <<", " <<ncol<<",elem: "<<elem1<<endl;
					CBSNP* pCBSNP 			=new CBSNP;

					if(ccol!=0 && ncol==1 ) // this condition must always be true ;
					{
						// The first column contains sample ids so save them
						pCBSNP->rsId			=elem1;
						pCBSNP->snpName			=elem1;

					}
					else
					{
						string msg= "problem in reading individual ids "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+inFileName+"\". \n";
						error(msg);
					}

					// start of reading from the second column till end
					//CBSNP * pCBSNP =genVec[0];
					unsigned int i =0;
					while(true)
					{


						elem1="";
						elem="";
						ccol=0;
						if(elem1.empty() )
							do
							{
								ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
								elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol,elem);
								if(zeilenEnde|| _ifs_eof(file))
									break;
							}while(elem1.empty());
						// start saving the genotypes
						 //pCBSNP =genVec[i];
						 if(elem1=="2")
						 {
							 // homozygote allele1: false false
							 pCBSNP->geno1.push_back(false);
							 pCBSNP->geno2.push_back(false);
							 pCBSNP->aOrder.push_back(true);
						  }else if(elem1=="1"){
							  // heterozygote major: false false
							  pCBSNP->geno1.push_back(false);
	 						  pCBSNP->geno2.push_back(true);
	 						  pCBSNP->aOrder.push_back(true);
						  }else if(elem1=="0"){
							  // homozygote allele2: false false
							  pCBSNP->geno1.push_back(true);
							  pCBSNP->geno2.push_back(true);
							  pCBSNP->aOrder.push_back(true);
						  }else
						  {
							  // homozygote major: false false
							  pCBSNP->geno1.push_back(true);
							  pCBSNP->geno2.push_back(false);
							  pCBSNP->aOrder.push_back(true);
						  }
							++i;
						 elem1="";
						// cout << ncol<<endl;
						if(zeilenEnde|| _ifs_eof(file))
						{
							// check if i has the equal to snp size
							if(vPed.size()!=i)
							{
								string msg ="Internal error, problem in Saving: \n"
											"There are only "+change_int_into_string(i)+" columns in the "+change_int_into_string(nrow)+"th row of the file.\n"
											"Please check if this line has information for each individual.\n"
										    "Otherwise: Hint: Size of vPed and index i  must  be equal. Please contact the software provider\n";
								error(msg);
							}
							if((vPed.size()!=  pCBSNP->geno1.size() )|| (  vPed.size()!=pCBSNP->geno2.size())){
							 string msg ="Internal error, problem in Saving: \n"
									  	  "\tSize of vPed and geno1 and geno2 must be same. Please contact the software provider\n";
							 error(msg);
							}


								ncol		=0;
								ccol		=0;
								break;

						}


					}
					genVec.push_back(pCBSNP);
					if(_ifs_eof(file))
						break;

				} //end of the line while(!zeilenEnde && !_ifs_eof(file))

		//-----------------------------------------------------------------------------------------------------------//
				//end of reading
		//-----------------------------------------------------------------------------------------------------------//
			if(_ifs_eof(file))
			{
				zeilenEnde=true;
				if(nrow==0)
				{
						string msg=" There are no elements in "+inFileName +" file.\n";
						error(msg);
					}
					else
					{
						printLIN("\t**->Total no of Individuals read: "+ change_int_into_string(vPed.size())+".\n");
						nrow=0;
					}
					nrow=0;
					break;
				}
	}

		}//end of file while
		_ifs_close(file);
		
		//-----------------------------------------------------------------------------------------------------------//

}

//read transposed genotype file 
void CRSNP::lese_transposed_genotype_file(vector<CBPED*>& vPed, vector<CBSNP*>& genVec, const string & inFileName)
{
	//printLIN("*->Reading file \""+ inFileName + "\": \n");

	IFS file(inFileName.c_str(),"r");
	if(!file)
	{
		string msg= "File \""+ inFileName+ "\" either does not exits or could not be opened.\n"	;
			error(msg);
	}
		// if file is okay then define more varaibles
		string elem					=""; // gelesen von ped File;
		string elem1				="";
		int ncol					=0;
		int ccol					=0;
		static int nrow				=0;
		bool zeilenEnde				=false;
		bool	assign_individ		=true;
		//unsigned int nPerson	=0;
		while(!_ifs_eof(file))
		{
				//address the end of the file
			if(_ifs_eof(file))
				{
					zeilenEnde=true;
					if(nrow==0)
					{
						string msg=" There are no elements in "+inFileName +" file.\n";
						error(msg);
					}
					else
					{
						cout << "Total no of Individual read: "<< (nrow-1)<< endl; // only for test
						nrow=0;
					}
					break;
				}
			zeilenEnde=false;
			// define  pointer of CBSNP
	//-----------------------------------------------------------------------------------------------------------//
			//reading first row
	//-----------------------------------------------------------------------------------------------------------//
		if(assign_individ)
		{
			//read first column
			ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
			elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
			if(_ifs_eof(file)) break; // break is okay because we do not want to read only one column
			if(zeilenEnde) break; // break is also okay because we do not want to read only one  column
			if(elem1.empty() )
				do
				{
					ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
					elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol,elem);
					if(zeilenEnde|| _ifs_eof(file))
						break;
				}while(elem1.empty());
			// If we read the first column and then there is the end of file then that is not correct
			//the follwoing line handels this mistake
			// there should be at least two columns with first the header that is why
			if( zeilenEnde || _ifs_eof(file) )
			{
				if(!elem1.empty())
				{
					string msg= " There is only one column in the first row of "+inFileName+" file. \n" ;
					error(msg);
				}
				else
					break; // only in the first column
			}
					//cout << elem1 <<endl;
			bool _b_tmp = (elem1=="rsid")||((elem1=="snpid")|| (elem1=="rsID")||(elem1=="SNPID")||(elem1=="snp_id")||(elem1=="snpID")||(elem1=="rs_id")||(elem1=="rs") ||(elem1=="SNP"));
			//cout <<boolalpha<<_b_tmp<<endl;
			if(!_b_tmp)
			{
				error("The first element of first column must be specified as the first column containing rsids.\n"
						"The first element of first column should have name like : \"rsid\""
						"\"snpid\",\"rsID\",\"snpID\" \n");
			}
			//  This case is true except last column
			if(zeilenEnde || _ifs_eof(file))
			{
				string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+inFileName+"\".\n";
				error(msg);
			}
					//starting reading second column of first row
					//each column gives the rsid or snpid
			while(!zeilenEnde && !_ifs_eof(file)) //reading start from second column to last column
					{

						//define pointer
						elem1="";
						elem="";
						ccol=0;
						ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
						elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);

						if(!elem1.empty())
						{
														//create SNP
							CBPED* pCBPED=new CBPED;
							if(pCBPED->indId=="")
							pCBPED->indId	=elem1;
							vPed.push_back(pCBPED);

						}

					}//end of 2nd column to last column
					// inform the total number of individuals
			printLIN("\t**->Total number of individuals found: "+change_int_into_string( (int)vPed.size() )+".\n ");
			assign_individ=false;
		}
	//-----------------------------------------------------------------------------------------------------------//
			//reading from second row
	//-----------------------------------------------------------------------------------------------------------//

		while(!_ifs_eof(file))
			{  // do all works till the end of line and then break

				//if(zeilenEnde || _ifs_eof(file))
				//read first column
				elem1="";
				elem="";
				ccol=0;
				ncol=0;
				ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
				elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
				if(_ifs_eof(file)) break; // break is okay because we do not want to read only one column
				if(zeilenEnde) break; // break is also okay because we do not want to read only one  column
				if(elem1.empty() )
					do
					{
						ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
						elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol,elem);
						if(zeilenEnde|| _ifs_eof(file))
							break;
					}while(elem1.empty());

				//this is true except for the last element of the line .
				if(zeilenEnde || _ifs_eof(file))
				{
								string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+inFileName+"\".\n";
								error(msg);
				}

				// read family id in the first column
				//cout << ccol <<", " <<ncol<<",elem: "<<elem1<<endl;
				CBSNP* pCBSNP 			=new CBSNP;
					
				if(ccol!=0 && ncol==1 ) // this condition must always be true ;
				{
					// The first column contains sample ids so save them
					pCBSNP->rsId			=elem1;
					pCBSNP->snpName			=elem1;
					
				}
				else
				{
					string msg= "problem in reading indiviudal ids "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+inFileName+"\". \n";
					error(msg);
				}

				// start of reading from the second column till end
				//CBSNP * pCBSNP =genVec[0];
				unsigned int i =0;
				while(true)
				{
					
					
					elem1="";
					elem="";
					ccol=0;
					if(elem1.empty() )
						do
						{
							ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
							elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol,elem);
							if(zeilenEnde|| _ifs_eof(file))
								break;
						}while(elem1.empty());
					// start saving the genotypes
					 //pCBSNP =genVec[i];
					 if(elem1=="2")
					 {
						 // homozygote allele1: false false
						 pCBSNP->geno1.push_back(false);
						 pCBSNP->geno2.push_back(false);
						 pCBSNP->aOrder.push_back(true);
					  }else if(elem1=="1"){
						  // heterozygote major: false false
						  pCBSNP->geno1.push_back(false);
 						  pCBSNP->geno2.push_back(true);
 						  pCBSNP->aOrder.push_back(true);
					  }else if(elem1=="0"){
						  // homozygote allele2: false false
						  pCBSNP->geno1.push_back(true);
						  pCBSNP->geno2.push_back(true);
						  pCBSNP->aOrder.push_back(true);
					  }else
					  {
						  // homozygote major: false false
						  pCBSNP->geno1.push_back(true);
						  pCBSNP->geno2.push_back(false);
						  pCBSNP->aOrder.push_back(true);
					  }
						++i;
					 elem1="";
					// cout << ncol<<endl;
					if(zeilenEnde|| _ifs_eof(file))
					{
						// check if i has the equal to snp size
						if(vPed.size()!=i)
						{
							string msg ="Internal error, problem in Saving: \n"
										"There are only "+change_int_into_string(i)+" columns in the "+change_int_into_string(nrow)+"th row of the file.\n"
										"Please check if this line has information for each individual.\n"
									    "Otherwise: Hint: Size of vPed and index i  must  be equal. Please contact the software provider\n";
							error(msg);
						}
						if((vPed.size()!=  pCBSNP->geno1.size() )|| (  vPed.size()!=pCBSNP->geno2.size())){
						 string msg ="Internal error, problem in Saving: \n"
								  	  "\tSize of vPed and geno1 and geno2 must be same. Please contact the software provider\n";
						 error(msg);
						}
				
						
							ncol		=0;
							ccol		=0;
							break;
						
					}


				}
				genVec.push_back(pCBSNP);
				if(_ifs_eof(file))
					break;

			} //end of the line while(!zeilenEnde && !_ifs_eof(file))

	//-----------------------------------------------------------------------------------------------------------//
			//end of reading
	//-----------------------------------------------------------------------------------------------------------//
		if(_ifs_eof(file))
		{
			zeilenEnde=true;
			if(nrow==0)
			{
				string msg=" There are no elements in "+inFileName +" file.\n";
				error(msg);
			}
			else
			{
				printLIN("\t**->Total no of Individuals read: "+ change_int_into_string(vPed.size())+".\n");
				nrow=0;
			}
			nrow=0;
			break;
		}

		}//end of file while
		_ifs_close(file);

}
