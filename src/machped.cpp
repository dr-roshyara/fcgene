/*
 * machped.cpp
 *
 *  Created on: 21.12.2011
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */
#include "machped.h"
	extern ofstream LIN;
	const unsigned  int CBPED:: coutWidth	=50;
	unsigned int CBPED::nMale				=0;
	unsigned int CBPED::nFemale				=0;
	unsigned int CBPED::nNosex				=0;
	unsigned int CBPED::nAffected			=0; // no of unaffected
	unsigned int CBPED::nUnaffected			=0; // no of affected
	unsigned int CBPED::nMiss_status		=0; //  undefined status
	bool CBPED::_given_cov_names			=	false;
	bool CBPED::_given_cov_types			=false;
	vector<string> CBPED::_cov_names_frm_argv; //this saves covariate names from argv
	vector<string> CBPED::_cov_types_frm_argv; //this saves covariate names from argv
	vector<string> CBPED::covariate_names; // this saves covariate names form file
	vector<CBPED*>CBPED::vPed; 
	vector<string> CBPED::save_line;
	//vector<CBPED*> CBPED::vPed(0);
	

	CBPED::~CBPED()
	{ 	
		/*if(vPed.size()>0)
		{
			for(vector<CBPED*>::iterator it=vPed.begin();it!=vPed.end();++it)
				delete *it;
			
		}
		vPed.clear();
		vPed.resize(0);
		*/
	}
	
CMPED::~CMPED()
	{	
		/*for(vector<CBPED*>::iterator it=vPed.begin();it!=vPed.end();++it)
			delete *it;
		
		vPed.clear();
		*/
		//cout << "\n I am a CMPED destructor\n";
		
	}
void CBPED::copy(const CBPED & pCBPED)
{

	 sex					=pCBPED.sex;
	 miss_sex				=pCBPED.miss_sex;
	 age 					=pCBPED.age;
	 pheno					=pCBPED.pheno;
	 miss_pheno				=pCBPED.miss_pheno;
	 ind_CR					=pCBPED.ind_CR;
	 quality				=pCBPED.quality;
	 indId					=pCBPED.indId;
	 famId					=pCBPED.famId;
	 patId					=pCBPED.patId;
	 matId					=pCBPED.matId;
	 groupLabel				=pCBPED.groupLabel;
	 save_line				=pCBPED.save_line;
	 covariate				=pCBPED.covariate;
	 //coutWidth				=pCBPED.coutWidth;
	 nMale					=pCBPED.nMale;
	 nFemale				=pCBPED.nFemale;
	 nNosex					=pCBPED.nNosex;
	 nAffected				=pCBPED.nAffected;
	 nUnaffected			=pCBPED.nUnaffected;
	 nMiss_status			=pCBPED.nMiss_status;
	 ngood_snps				=pCBPED.ngood_snps;
	 covariate_names		=pCBPED.covariate_names;
	 _cov_names_frm_argv	=pCBPED._cov_names_frm_argv;
	 _cov_types_frm_argv	=pCBPED._cov_types_frm_argv;
	 _given_cov_names		=pCBPED._given_cov_names;
	 _given_cov_types		=pCBPED._given_cov_types;
	 vPed					=pCBPED.vPed;
	 given_indiv_crate		=pCBPED.given_indiv_crate;

	 return ;}
CBPED::CBPED(const CBPED& pCBPED)
 {
	 copy(pCBPED);
 }
 CBPED&  CBPED::operator=(const CBPED& pCBPED)
{
	copy(pCBPED);
	return *this;
}

 //--------------------------------------------------------------------//
 //function for CBSNP
 //---------------------------------------------------------------------//

 void CBPED::genInfoLeser(const vector<CBSNP*>& genInfo, const string& first , const string& second, const int & ncol , const int& nrow, const int& cSNPN )
 {
 	bool fatal=false;
 	string EMSG ="\t**->probblem with reading genotypes in "+change_int_into_string(ncol-1)+"th column and "+change_int_into_string(ncol)+" th column of  line "+change_int_into_string(nrow+1)+".\n";
 	string errmsg="";
 	//
	//cout << "no of SNPs: "<<cSNPN <<endl;
	CBSNP* pCBSNP=genInfo[cSNPN-1];
 	if(first.length() != second.length()&& (first.length() !=1 && second.length()!=1) )
		{
			string rId	= pCBSNP->snpName;
			string msg=EMSG+ "SNP "+rId+" has genotypes [ " + first + " ] and ["+second+"],  each should be only one character long. \n";
			//delete pCBSNP;
			error(msg);
		}
		//cout <<  pCBSNP->snpName<< endl;
		//pCBSNP->geno1.push_back(false);

	// this is the case where both genotypes are non missing
 	if(first!=CBSNP::_miss_gvalue&& second!=CBSNP::_miss_gvalue)
 	{
 				if(first==second)
 				{
 					// add allele
 					if(pCBSNP->allele1==""&& pCBSNP->allele2!=first)
						pCBSNP->allele1=first;
 					else if (pCBSNP->allele2==""&& pCBSNP->allele1!=first)
						pCBSNP->allele2=first;
 					// add genotype
 					if(first==pCBSNP->allele1&& second==pCBSNP->allele1 )
 					{// this is the case where homozygote1 is denoted by 00
 						//cout << "reading homo1 \n";
 						//cout << first << " "<<second<<"\n";
 						pCBSNP->geno1.push_back(false);
 						pCBSNP->geno2.push_back(false);
						pCBSNP->aOrder.push_back(true);
 					}
 					else if(first==pCBSNP->allele2 && second==pCBSNP->allele2)
 					{// this is the case where homozygote 2  is denoted by 11
 						//cout << "reading homo2\n";// for test only
 						//cout << first << " "<<second<<"\n";// for test only
 						pCBSNP->geno1.push_back(true);
 						pCBSNP->geno2.push_back(true);
						pCBSNP->aOrder.push_back(true);
 					}
 					else
 					{
 						string wrongAllele="";
 						if(first!=pCBSNP->allele1 || first!=pCBSNP->allele2 )
 							wrongAllele=first;
 						else if(second!=pCBSNP->allele1 || second!=pCBSNP->allele2 )
 							wrongAllele=second;
 					   errmsg="SNP "+pCBSNP->snpName+" containing two alleles "+\
 						pCBSNP->allele1+ " and "+pCBSNP->allele2+ " has found a new allele type"+wrongAllele+ "  at individual "+\
 						change_int_into_string(nrow)+ ".\n ";
 					  fatal=true;
 					}
 				}
 				else if(first!=second) // save it as 01 .i.e. false true.
 				{// this is the case where heterozygote   is denoted by 01
 					//cout << "reading heterozygote\n";// for test only
					//cout << first << " "<<second<<"\n";// for test only
 					if(pCBSNP->allele1==""&&pCBSNP->allele2!=first )pCBSNP->allele1=first;
 					else if (pCBSNP->allele2==""&&pCBSNP->allele1!=first ) pCBSNP->allele2=first;
 					if (pCBSNP->allele1==""&& pCBSNP->allele2!=second ) pCBSNP->allele1=second;
 					if (pCBSNP->allele2==""&& pCBSNP->allele1!=second ) pCBSNP->allele2=second;
 					//now check if allele are correct in first and second
 					if(first==pCBSNP->allele1&& second==pCBSNP->allele2 )
 					{
 					 	pCBSNP->geno1.push_back(false);
 						pCBSNP->geno2.push_back(true);
						pCBSNP->aOrder.push_back(true);
 					}
 					else if(first==pCBSNP->allele2&& second==pCBSNP->allele1 )
 					 	 {
 							pCBSNP->geno1.push_back(false);
 							pCBSNP->geno2.push_back(true);
							pCBSNP->aOrder.push_back(false);
							// in this case order  given in original file
							// has been chnaged . So remember the order
							//pCBSNP->order[cSNPN]=false;
 					 	 }
 					else
 					{
 						string wrongAllele="";
 						if(first!=pCBSNP->allele1 || first!=pCBSNP->allele2 )
 						wrongAllele=first;
 						else if(second!=pCBSNP->allele1 || second!=pCBSNP->allele2 )
 						wrongAllele=second;
 						errmsg="SNP "+pCBSNP->snpName+" containing two alleles "+\
 								pCBSNP->allele1+ " and "+pCBSNP->allele2+ " has found a new allele type"+wrongAllele+ "  at individual "+\
 								change_int_into_string(nrow+1)+ ".\n ";
 						 fatal=true;

 					}


 				}
 	}
 	// This is the case where both genotypes are missing

	else if (first==CBSNP::_miss_gvalue || second==CBSNP::_miss_gvalue )
 		{ //In this case we use the code 10 to denote the missing genotypes
				//cout << "reading missing \n"; // for test only
				//cout << first << " "<<second<<"\n"; // for test only
				pCBSNP->geno1.push_back(true);
				pCBSNP->geno2.push_back(false);
				pCBSNP->aOrder.push_back(true);
 		}
 		// this is the case where non of above condition is true then error.
 	if (fatal)
 	 {
		//delete pCBSNP;
		printLIN( EMSG); //LIN<< EMSG; // okay in original .
 		error(errmsg);
 	}

	// check if size of geno1 and geno2 are equal they must be equal.

	if(pCBSNP->geno1.size()!= pCBSNP->geno2.size() )
	{
		//delete pCBSNP;
		error(EMSG);
	}

return ;}

 // writing output files
 void CBPED:: writeOutPlinkFiles(const vector<CBPED*>& pedInfo, const vector<CBSNP*> & genInfo, const string&outFileName)
 {
	 //----------------------------------#--------------------------------------#
	  	printLIN("*->Writing files in "+CBSNP::_out_format+" format: \n");
	  //----------------------------------#--------------------------------------#
	  	CBSNP::write_plinkMapFile(genInfo, outFileName);
	  	CBPED::write_plinkPedFile(pedInfo,genInfo, outFileName);
 }
 void CBPED:: write_plinkPedFile(const vector<CBPED*>& pedInfo, const vector<CBSNP*> & genInfo, const string&outFileName)
 {
	 //----------------------------------#--------------------------------------#
	  	// writing ped files.
	  //----------------------------------#--------------------------------------#
		 ofstream ofs;
		 ofs.clear();ofs.close();
	  	string pedFileName=outFileName+".ped";
	  	ofs.open(pedFileName.c_str());
	  	if(!ofs)
	  			error("out file "+outFileName+" does not exit.\n");
	  	ofs.clear();
		//cout <<pedInfo.size() <<endl;
		//cout <<genInfo.size() <<endl;
	  	for(unsigned int i=0; i<pedInfo.size();++i)
	  	{
	  		if(pedInfo[i]->quality)
	  		{
				if(pedInfo[i]->famId!="") ofs<<pedInfo[i]->famId<< " ";
				else ofs << "famid"<<i+1<< " ";
				if(pedInfo[i]->indId!="")ofs<<pedInfo[i]->indId << " ";
				else ofs<< "indv"<<i+1<<" ";
				if(pedInfo[i]->patId!="")ofs<<pedInfo[i]->patId << " ";
				else ofs<< 0<<" ";
				if(pedInfo[i]->matId!="") ofs <<pedInfo[i]->matId <<" ";
				else ofs<< 0 << " ";
				if(pedInfo[i]->sex) ofs<<"1"<< " "; // 1 is male
				else if(!pedInfo[i]->sex && !pedInfo[i]->miss_sex) ofs<<"2"<<" "; // 2 is female
				else if(pedInfo[i]->miss_sex) ofs<<"0"<< " "; // 0 is missing
				if(pedInfo[i]->pheno)
					 ofs<<"2"<< " ";
				 else if(!pedInfo[i]->pheno && !pedInfo[i]->miss_pheno)
						ofs<<"1"<< " ";
				else if(!pedInfo[i]->pheno && pedInfo[i]->miss_pheno)
				{
					if (CBSNP::_out_format =="haploview")
						ofs<<"0" << " ";
					else
					ofs<<"-9"<< " ";
				}
				//	cout << "test1"<<endl;
				for(unsigned int j=0;j<genInfo.size(); ++j)
				{
				//-----------------------------------------------------------------//
					CBSNP* pCBSNP=genInfo[j];
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
												ofs <<" "<<CBSNP::_miss_gvalue << " "<<CBSNP::_miss_gvalue;

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
								ofs <<" "<<CBSNP::_miss_gvalue << " "<<CBSNP::_miss_gvalue;
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
													ofs <<" "<<CBSNP::_miss_gvalue << " "<<CBSNP::_miss_gvalue;

								else
								{
									error("Problem in writing out "+change_int_into_string(j) +"th SNP.\n");


								}



						}
					}

				}
				//-----------------------------------------------------------------//
				//cout << "test2"<<endl;
				ofs<< "\n";
			} // end of quality loop
	  	}//end of pCBPED for loop
	  	printLIN("\t**->PLINK ped file has been written out and saved as \""+pedFileName+\
	  			"\".\n");
	  	ofs.clear();	 ofs.close();
	  	//----------------------------------#--------------------------------------#


 return ;}
 // for extrafiles
 void CBPED::pedinfoFileLeser(const vector<CBPED*> & pedInfo, const string& fileName)
 {
 	string r_msg="\n*->Reading  file \""+fileName+ "\":\n\n";
 	printLIN(r_msg);
 	string line			="";
 	rFile file;
 	file.close();
 	file.clear(); 	// check if it has been already opend.
 	static int extraMapFileLineZaehler	=0;			 // extra map file line zaehler
 	//ios_base::iostate i= file.rdstate();
 		file.open(fileName.c_str(),ios::in);
 		// this will check if each row has same number of columns.
 	if(!file)
 		{
 		file.setstate(ios_base::failbit);
 		file.clear();
 		error("your extra  file "+fileName+ " either does not exit or could not be opened.");
 		}
 		// check first header and  the number of columns
 		do
 		{
 		getline(file,line, '\n');
 		if(line[0]=='#') // leave the comments;
 			continue;
 		if(file.eofbit) break;
 		}while( line.empty());
 		//variable declaration need later .
 		string buffer="";
 		vector<string> header;
 		const unsigned int dvalue	=0;
 		unsigned int col_famid		=dvalue; // initialization of columns
 		unsigned int col_indid		=dvalue;
 		unsigned int col_patid		=dvalue;
 		unsigned int col_matid		=dvalue;
 		unsigned int col_sex		=dvalue;
 		unsigned int col_status		=dvalue;
 		unsigned int nCols			=dvalue;
 		bool found					=false;
 		//bool found1 				=false;
 		unsigned int nMan			=0;
 		unsigned int nWeib			=0;
 		unsigned int nTrans			=0;
 		//
 		unsigned int nAff			=0;
 		unsigned int nUnaff			=0;
 		unsigned int nMiss_state	=0;

 		stringstream line_parser(line);
 		if(line_parser.good())
 		{
 			while(line_parser>>buffer)
 			{
 				header.push_back(buffer);
 				//cout << " buffer " //only for test
 			}
 			nCols=header.size();
 			vector<string>sort_header=header;
 			vector<string>::iterator it;
 			sort(sort_header.begin(),sort_header.end());
 			//famid
 			found	= binary_search(sort_header.begin(), sort_header.end(),"famid");
 			//found1	= binary_search(sort_header.begin(), sort_header.end(),"pid");
 			if(found)
 				{
 					//famid
 					it=find(header.begin(),header.end(),"famid");
 					col_famid= (int)(distance(header.begin(),it)+1);
 				}
 			found	= binary_search(sort_header.begin(), sort_header.end(),"famid");
 			//indid
 			if(found)
 			{
 				//indid
 				it=find(header.begin(),header.end(),"indid");
 				col_indid= (int)(distance(header.begin(),it)+1);
 			}
 			//patid
 			found= binary_search(sort_header.begin(), sort_header.end(),"patid");
 			if(found)
 			{	// patid
 				it=find(header.begin(),header.end(),"patid");
 				col_patid= (int)(distance(header.begin(),it)+1);
 			}
 			//matid
 			found= binary_search(sort_header.begin(), sort_header.end(),"matid");
 			if(found)
 			{
 			//matid
 				it=find(header.begin(),header.end(),"matid");
 				col_matid= (int)(distance(header.begin(),it)+1);
 			}
 			//sex
 			found= binary_search(sort_header.begin(), sort_header.end(),"sex");
 			if(found)
 			{	// sex
 				it=find(header.begin(),header.end(),"sex");
 				col_sex= (int)(distance(header.begin(),it)+1);
 			}
 			//status
 			 found= binary_search(sort_header.begin(), sort_header.end(),"phenotype");
 			 if(found)
 			 {	// staus
 				 it=find(header.begin(),header.end(),"phenotype");
 				 col_status= (int)(distance(header.begin(),it)+1);
 			 }
 			//displaying for test



 		}
 		else
 		{
 			string msg= "Problem in parsing the first line of"+fileName+"!\n";
 			error(msg);

 		}
 		++extraMapFileLineZaehler;
 		vector<string> tokens;
 		vector<CBPED*>::const_iterator it_pCBPED=pedInfo.begin();

 		while(getline(file,line, '\n')) // until eof hit
 		{

 				if(line.length()==0)
 				continue;
 				if(line[0]=='#')
 				continue;
 				tokens.resize(0);
 				tokens.clear();
 				stringstream line_parser1(line);
 				//cout <<"\n line_parser1.good(): "<<line_parser1.good() <<endl; // only for test
 				if(line_parser1.good())
 				{
 					while(line_parser1>>buffer)
 						tokens.push_back(buffer);
 					//cout << buffer << " "; // only for test

 					if(tokens.size()!=nCols)
 					{
 						string msg= "Your extra  file \""+fileName+"\"  has "+change_int_into_string(nCols)+\
 								" columns in header line but in the line "+change_int_into_string(extraMapFileLineZaehler)+\
 								" there are "+change_int_into_string(tokens.size())+" columns.\n";
 						error(msg);
 					}


 				//famid
 					if(col_famid!=dvalue)
 					{
 						if((*it_pCBPED)->famId!="")
 						{
 							if((*it_pCBPED)->famId !=tokens[col_famid-1])
 							{
 								string msg= " family id  \""+(*it_pCBPED)->famId+"\" in ped file  does not match with family id \""+\
 										tokens[col_famid-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of  file\""+fileName+\
 										"\" \n.";
 								error(msg);
 							}
 						}
 							else
 								(*it_pCBPED)->famId =tokens[col_famid-1];
 					}

 					//indid

 					if(col_indid!=dvalue)
 					{
 						if((*it_pCBPED)->indId!="")
 						{
 							if((*it_pCBPED)->indId !=tokens[col_indid-1])
 							{
 								string msg= " individual id  \""+(*it_pCBPED)->indId+"\" in ped file  does not match with individual id \""+\
 										tokens[col_indid-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of  file\""+fileName+\
 										"\" \n.";
 								error(msg);
 							}

 						}
 						else
 							(*it_pCBPED)->indId =tokens[col_indid-1];
 					}

 					//patId
 					if(col_patid!=dvalue)
 					{
 						if((*it_pCBPED)->patId!="" && (*it_pCBPED)->matId!="0")
 						{
 							if((*it_pCBPED)->patId !=tokens[col_patid-1])
 							{
 								string msg= " Father Id \""+(*it_pCBPED)->patId+"\" in ped  file  does not match with father id \""+\
 										tokens[col_patid-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra  file\""+fileName+\
 										"\" \n.";
 								error(msg);
 							}
 						}
 						else
 							(*it_pCBPED)->patId =tokens[col_patid-1];
 					}
 					// mat id
 					if(col_matid!=dvalue)
 					{
 						if((*it_pCBPED)->matId!="" && (*it_pCBPED)->matId!="0")
 						{
 							if((*it_pCBPED)->matId !=tokens[col_matid-1])
 							{
 								string msg= " Mother id \""+(*it_pCBPED)->matId+
 										"\" in ped file  does not match with mother id \""+\
 										tokens[col_matid-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra  file\""+fileName+\
 										"\" \n.";
 								error(msg);
 							}
 						}
 						else
 						(*it_pCBPED)->matId =tokens[col_matid-1];
 					}

 					// sex
 					 if(col_sex!=dvalue)
 					 {
 						 string code =tokens[col_sex-1];
 						 if(!(*it_pCBPED)->sex && (*it_pCBPED)->miss_sex) // missing case
 						 {
 							 //do something
 							if(code=="M" || code=="1")
 							{
 								(*it_pCBPED)->sex				=true;
 								(*it_pCBPED)->miss_sex			=false;
 								++nMan; //local code
 								++nMale;
 								if(nNosex!=0)
 									--nNosex;
 							}
 							else if(code=="F"|| code=="2")
 							{
 								++nWeib; // local code
 								++nFemale; // local code
 								if(nNosex!=0)
 									--nNosex;
 								(*it_pCBPED)->sex				=false;
 								(*it_pCBPED)->miss_sex			=false;

 							}
 							else
 							{
 								(*it_pCBPED)->sex				=false;
 								(*it_pCBPED)->miss_sex			=true;
 									++nTrans;
 							}
 													// It is not necessary to assign if elem=f because sex is then already assigned false
 						 }
 						 else
 						 {
 							 //prüfen

 							 // case1
 							 	 //cout << "sex : "<< boolalpha << (*it_pCBPED)->sex <<endl;
 								 //cout << "aorb :"<<(code!="1" &&  code!="M") << endl;

 							  if((*it_pCBPED)->sex && (code!="1" && code!="M"))
 							 {
 								 string msg= " Sex status in ped file is male but it  does not match with sex status \""+\
 										 tokens[col_sex-1]+"\" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra  file\""+fileName+\
 										 "\" \n.";
 								 error(msg);
 							 }
 							 // case2
 							 if(!(*it_pCBPED)->sex && (code!="2" && code!="F"))
 							 {
 								 string msg= " Sex status  in ped file is female but it  does not match with sex status \""+\
 										 tokens[col_matid-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra  file\""+fileName+\
 							 										 "\" \n.";
 								 error(msg);
 							 }

 						 }
 					 }

 					// status
 					if(col_status!=dvalue)
 					{
 						string code =tokens[col_status-1];
 					 	 //cout << "phenotype code: " <<code <<endl; // only for test
 						if((*it_pCBPED)->miss_pheno) // missing case
 						{
 							//do something

 							if(code=="aff" || code=="2")
 							{
 								(*it_pCBPED)->pheno				=true;
 								(*it_pCBPED)->miss_pheno		=false;
 								++nAff;//local code
 								++nAffected;
 								--nMiss_status; //one missing will be reduced.
 							}
 							else if(code=="F"|| code=="1")
 							{
 								++nUnaff; //local code
 								++nUnaffected;
 								--nMiss_status;
 								(*it_pCBPED)->pheno				=false;
 								(*it_pCBPED)->miss_pheno		=false;
 							}
 							else
 							{
 								(*it_pCBPED)->pheno				=false;
 								(*it_pCBPED)->miss_pheno		=true;
 	 							++nMiss_state; //
 							}
 							// It is not necessary to assign if elem=f because pheno is then already assigned false
 						}
 						else
 						{
 							//prüfen
 							// case1
 							if((*it_pCBPED)->pheno && (code!="2" && code!="aff"))
 							{
 								string msg= " Phenotype  in ped file is affected but it  does not match with pheno pheno \""+\
 										tokens[col_matid-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra  file\""+fileName+\
 										"\" \n.";
 								error(msg);
 							}
 					 							 // case2
 							if(!(*it_pCBPED)->pheno && (code!="1" && code!="unaff"))
 							{
 								string msg= " Phenotype  in ped file is unaffected but it  does not match with phenotype \""+\
 										tokens[col_matid-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra  file\""+fileName+\
 										"\" \n.";
 								error(msg);
 							}

 						}
 					}


 				// end of loop
 				if(file.eof()) break;
 				}
 				else
 				{
 					printLIN( "Line parser not good while reading extra  info file\""+fileName+"\" .\n");
 					break;
 				}
 				if(it_pCBPED==(pedInfo.end()-1)) break;
 				else
 				++it_pCBPED;

 				++extraMapFileLineZaehler;

 		}
 		// check if parsing fails
 		// if parsing is good then go further
 printWithSpace ("\t**->Total number of successfully updated Males:",change_int_into_string(nMan)+"\n ",coutWidth );
 printWithSpace ("\t**->Total number of successfully updated Females:", change_int_into_string(nWeib)+"\n ",coutWidth );
 printWithSpace ("\t**->Total number of successfully updated missing sex:",change_int_into_string(nTrans)+"\n\n ",coutWidth );
 printWithSpace ("\t**->Total number of successfully updated cases: ",change_int_into_string(nAff)+"\n ",coutWidth );
 printWithSpace ("\t**->Total number of successfully updated controls:",change_int_into_string(nUnaff)+"\n ",coutWidth);
 printWithSpace ("\t**->Total number of successfully updated missing status:",change_int_into_string(nMiss_state)+"\n\n",coutWidth );

 file.close();
 return ;}

 void CBPED::write_pedinfoFiles(const vector<CBPED*>& pedInfo, const string& outFileName)
 {
	 //----------------------------------#--------------------------------------#
	 	 	// writing extra ped info files.
	 	 //----------------------------------#--------------------------------------#
	 	 	 ofstream ofs;
	 	 	string extraPedFileName=outFileName+"_pedinfo.txt";
	 		ofs.clear();	 ofs.close();
	 	 	ofs.open(extraPedFileName.c_str());
	 	 	if(!ofs)
	 	 			error("out file "+outFileName+" does not exit or can not be written out.\n");
	 	 	//
	 	 	ofs<<"famid" <<" ";
	 	 	ofs<<"indid" <<" ";
	 	  	ofs<<"patid" <<" ";
	 		ofs<<"matid" <<" ";
	 	 	ofs<<"sex"	 <<" ";
	 	 	ofs<<"phenotype"<<"\n";
	 	  	for(unsigned int i=0; i<pedInfo.size();++i)
	 	 	{
	 	 		CBPED* pCBPED =pedInfo[i];
	 	 		if(pCBPED->quality)
	 	 		{
					//#####
					if(pCBPED->famId!="") ofs<<pCBPED->famId<< " ";
					else ofs << "famid"<<i<< " ";
					if(pCBPED->indId!="")ofs<<pCBPED->indId << " ";
					else ofs<< "indv"<<i<<" ";
					if(pCBPED->patId!="")ofs<<pCBPED->patId << " ";
					else ofs<< 0<<" ";
					if(pCBPED->matId!="") ofs <<pCBPED->matId <<" ";
					else ofs<< 0 << " ";
					if(pCBPED->sex) ofs<<"1"<< " "; // 1 is male
					else if(!pCBPED->sex && !pCBPED->miss_sex) ofs<<"2"<<" "; // 2 is female
					else if(pCBPED->miss_sex) ofs<<"0"<< " "; // 0 is missing
					if(pCBPED->pheno)
						 ofs<<"2"<< " ";
					 else if(!pCBPED->pheno && !pCBPED->miss_pheno)
							ofs<<"1"<< " ";
					else if(!pCBPED->pheno && pCBPED->miss_pheno)
						ofs<<"-9"<< " ";
					 ofs <<"\n";
						//########
	 	 		}	//end of if quality loop
	 	 	}//end of for loop

	 		ofs.clear();	 ofs.close();
	 	printLIN("\t**->One extra \"pedinfo.txt\" file has  been written out and saved as \""+extraPedFileName+\
	 		 			"\".\n");
	  return ;

 return ;}

void CBPED::write_impute_gensFile(const vector<CBSNP*>& genInfo,const vector<CBPED*>&pedInfo,
		const string& outFileName,
		const bool _given_extra_snpinfo
		)
{
	 //----------------------------------#--------------------------------------#
		  	printLIN("*->Writing files in "+CBSNP::_out_format+" format: \n");

		//----------------------------------#--------------------------------------#
	 	// writing gens files.
	  	//----------------------------------#--------------------------------------#
	 	ofstream ofs;
	 	 ofs.clear();ofs.close();
	 	 ofstream ofs1;
	 	 ofs1.clear(); ofs1.close();
			CBPED* pCBPED =pedInfo[0];
			CBSNP* pCBSNP =genInfo[0];
	 string gensFileName=outFileName+".gens";
	 if(CBSNP::_out_format=="snptest")
		gensFileName=outFileName+".gen";
	 string strandFileName=outFileName+".strand.txt";
	 string commandFileName=outFileName+".commands.txt";
	 	ofs.open(gensFileName.c_str());
	 	if(!ofs)
	 			error("out file "+outFileName+" does not exit or can not be written out.\n");
	 	ofs.clear();
	 	// ofs1
	 		ofs1.open(strandFileName.c_str());
	 	 	if(!ofs1)
	 	 		error("out file "+outFileName+" does not exit.\n");
	 		 	ofs1.clear();
	 		 //collect bp
	 		 	vector<int> temp_vec_bp;
	 		 	vector<string> temp_vec_rsid;
	 		for(unsigned int j=0;j<genInfo.size(); ++j)
	 		{
	 			 pCBSNP=genInfo[j];
	 			//cout <<"(pCBSNP->quality): "<< boolalpha << pCBSNP->quality <<endl;
	 			if(pCBSNP->quality)
	 			{
	 				// write strand file
					if(pCBSNP->bp!=-1){
						ofs1 << pCBSNP->bp<< " ";
					}else{
						string str_err ="To work with IMPUTE and SNPTEST  programs, base pair positions are necessary. \nPlease update them first then only you can convert them.\n";
						error(str_err);
					}
					if(pCBSNP->coding_strand==0)
					{
						ofs1 << "+\n";
					}else if(pCBSNP->coding_strand==1){
						ofs1 << "+\n";
					}
					else if(pCBSNP->coding_strand==2){
						ofs1 << "-\n";
					}
					if(pCBSNP->rsId!="")
					{
						temp_vec_rsid.push_back(pCBSNP->rsId);
					}
					else if(pCBSNP->snpName!="") 
					{
						temp_vec_rsid.push_back(pCBSNP->snpName);
					}
						
					if(pCBSNP->geno1.size()!=pedInfo.size() || pCBSNP->geno2.size()!=pedInfo.size() )
					{
						cout <<"geno1.size() : "<<pCBSNP->geno1.size() <<"geno2.size():  "<<pCBSNP->geno1.size();
						for(unsigned int i=0;i<pCBSNP->geno1.size();++i)
							ofs<< pCBSNP->geno1[i] << " " <<pCBSNP->geno2[i] << " ";
							ofs.close();
							error("problem in writing out "+change_int_into_string(j) +"th SNP. geno1 and geno2 have no equal sizes\n");
					}
					else
					{
						if(pCBSNP->snpName!="") ofs << pCBSNP->snpName<< " ";
						else if(pCBSNP->rsId!="") ofs << pCBSNP->rsId<< " ";
						if(pCBSNP->rsId!="") ofs << pCBSNP->rsId<< " ";
							 else if(pCBSNP->snpName!="") ofs << pCBSNP->snpName<< " ";
						//cout <<pCBSNP->bp << " ";
						if(pCBSNP->bp!=(-1))
						{
							ofs<< pCBSNP->bp<< " ";
							temp_vec_bp.push_back(pCBSNP->bp);
						}
						else
						{
							if(!(j))
								printLIN("---------------------------\n"
										"WARNING!!!\n"
										"Your input files do not contain Basepair positions and  also you haven't updated them with --snpinfo option.\n"
										"It may be a problem if you want to impute this data with IMPUTE or use it in SNPTEST.\n"
										"We recommend you to generate files together with the \"--snpinfo\" option.\n"
										"---------------------------\n"
								);
														ofs << j<< " ";
							temp_vec_bp.push_back(j);

						}
						if(pCBSNP->allele1!="") ofs << pCBSNP->allele1 << " ";
						else if(pCBSNP->allele2!="") ofs << pCBSNP->allele2<< " ";
						else
						{
							ofs << "A" << " ";
						}
						if(pCBSNP->allele2!="") ofs << pCBSNP->allele2;
						  else if(pCBSNP->allele1!="") ofs << pCBSNP->allele1;
						 else
						 {
							 ofs << "B" ;
						 }
						// write geno prob if they are given.

						for(unsigned int i=0; i<pedInfo.size();++i)
						{
							pCBPED =pedInfo[i];
							if(pCBPED->quality)
							{
								//case1 if genotype probabilities are given
								if((pCBSNP->pgeno1.size()==pedInfo.size()) &&(pCBSNP->pgeno2.size()==pedInfo.size())&&(pCBSNP->pgeno3.size()==pedInfo.size()))
								{
									ofs <<" "<< pCBSNP->pgeno1[i] << " "<< pCBSNP->pgeno2[i] << " "<< pCBSNP->pgeno3[i];
								}
								//case 2 if genotype probabilities are not given
								else
								{
									if( !pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if homozygote1  00
										ofs <<" "<< "1" <<" "<< "0" <<" "<< "0";
									else if( pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if homozygote2 11
										ofs <<" "<< "0" <<" "<< "0" <<" "<< "1";
									else if( !pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if heterozygote
										ofs <<" "<< "0" <<" "<< "1" <<" "<< "0";
									else if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) //if missing
										ofs <<" "<<atoi(CBSNP::_miss_gvalue.c_str()) << " "<<CBSNP::_miss_gvalue << " "<<CBSNP::_miss_gvalue;
									else
									{
										error("problem in writing out "+change_int_into_string(j) +"th SNP.\n");
									}
									//	ofs <<"TEst\n";
								}
						} //end of if pedind quality

						}//end of for pedinfo loop

					}
				ofs<< "\n";
				}
	 		}

	 	ofs.clear();	 ofs.close();
	 	ofs1.clear();	 ofs1.close();
	 	printLIN( "\t**->Impute gens file has been written out and saved as \""+gensFileName+"\". \n");
	 	printLIN( "\t**->Impute strand file has been written out and saved as \""+strandFileName+"\". \n");

	 	//---------------------------------------------------------------------//
		ofs.open(commandFileName.c_str());
		 	if(!ofs)
		 		error("out file "+outFileName+" does not exits.\n");
		 	ofs.clear();
		 	bool sifo_tfvec =temp_vec_bp.size()<CBSNP::sifo_bp.size();
		 	ofs<<"Templates of IMPUTE commands:\n";
			if(sifo_tfvec&&_given_extra_snpinfo)
		 	{
				ofs<<"It seems like that you are interested to impute your genotyped data using a large reference panel. ";
				ofs<<"Therefore templates for the commands of  both \"imputation without using reference panel\" and  \"imputation using references\"  are given below.\n";
		 		
			}
			
				vector<int>temp_vec_bp1 =temp_vec_bp;
				vector<int>lower_range_vec;
				vector<int>upper_range_vec;
				vector<float>interval_size_MB;
				sort(temp_vec_bp1.begin(),temp_vec_bp1.end());
			
					
			{ // impute commmands for without using references.
				
				//Finding SNP Position which geq 5 MB
				// define necessary variables
					ofs<<"\n-------------------------------------------------------------------\n";
					ofs<<"Templates of IMPUTE commands for imputation without using references:\n";
					ofs<<"-------------------------------------------------------------------\n";
				fix_ints_4_impute_command(temp_vec_bp1, lower_range_vec,upper_range_vec, interval_size_MB);	
			
				if(lower_range_vec.size()==1)
				{
					ofs<<"The impute command written below is just an example to show how you can impute your genotype data with IMPUTE. You may have to  modify it ";
					ofs<<"at necessary places (especially the names of the input data) so that the command can be executed with IMPUTE properly.\n";
					ofs<<"\n------------\n";
					ofs<<" Command: \n";
					ofs<<"------------\n";
					ofs<<"./impute2 \\ \n";
					ofs<<"-m genetic_map.txt\\ \n";
					//ofs<<"-h hapmap3_r2_b36/hapmap3.r2.b36.all.chr22.haps\\ \n";
					//ofs<<"-l hapmap3_r2_b36/hapmap3.r2.b36.all.chr22.legend\\ \n";
					ofs<<"-g  " << gensFileName <<" \\  \n";
					ofs <<"-strand_g  " << strandFileName<<"\\ \n";
					ofs<< "-int"<<" "<< lower_range_vec[0] <<"  " << upper_range_vec[0]<<"\\ \n";
					ofs<<"-Ne 11418  -call_thresh 0.9 -pgs\\ \n";
					ofs<< "-o "+ outFileName+".impute2 \n";
				}
				else
				{
					ofs<<"Your genotype Data has a region of larger than 7Mb.IMPUTE software providers recommend to split the genotyping region into samll chunks for imputation purpose. ";
					ofs<<"These chunks can be  be imputed in parallel on multiple computer processors.\n";
					ofs<<"The following table  describes the lower are upper boudaries (in base pair position) of the region in which imputation should be carried out.\n";
					ofs<<"------------------------------------------------------------------------\n";
					ofs<<"Table containing lower and upper boundaries of chunks: \n";
					ofs<<"------------------------------------------------------------------------\n";
					ofs <<"chunk No:\t"<<"lower boundary\t"<< "Upper boundary\t"<<"Size of Chunk in MB\n";
					for(unsigned int i=0;i<lower_range_vec.size();++i)
							{

								ofs << "chunk "<<i+1<<":\t"<< lower_range_vec[i]<< "\t"<<upper_range_vec[i]<<"\t   "<<interval_size_MB[i] << endl;
							}
					ofs<<"------------------------------------------------------------------------\n\n";
					ofs<<"The impute commands written below are just the templates to show  how you can impute your genotype data with IMPUTE. You can modify them ";
					ofs<< "at necessary places (especially the input data names) so that the commands can be executed with IMPUTE properly. Moreover the upper and lower base pair positions given with \"-int\" command, \n";
					ofs<< "are based on the genotyped data (NOT on the reference panel). If you prefer to impute with the reference panel, then  the perl script or commands given  in the section of imputation using reference panel, should be used.\n";
					//perl command 
					{	
						ofs<<"\n------------------------------------------------------------\n";
						ofs<<"perl scripts to give impute commands in Linux cluster: \n";
						ofs<<"-----------------------------------------------------------\n";
						ofs << "#!/usr/bin/perl -w\n";
						ofs<< "my @lowerBound=(";
						 for(unsigned int i=0;i<(lower_range_vec.size()-1);++i)
						 ofs<<lower_range_vec[i] << ", ";
						 ofs<<lower_range_vec[lower_range_vec.size()-1];
						 ofs<< ");\n";
						 ofs<< "my @upperBound= (";
						 //upper bound
						 for(unsigned int i=0;i<(upper_range_vec.size()-1);++i)
						 ofs<<upper_range_vec[i] << ", ";
						 ofs<<upper_range_vec[upper_range_vec.size()-1];
						 ofs<< ");\n";
						 ofs << "my $nChunks = "<<upper_range_vec.size()<<";\n";
						 ofs<<"my @segs= (0 .. $nChunks-1);\n";
						 ofs<< "foreach $i (@segs)\n"; 
						 ofs <<"{\n";
							ofs<<"\tmy $actChunk=$i+1;\n";
							ofs<<"\t system (\n";
								ofs<<"\t\t\"./impute2 \". \n";
								ofs<<"\t\t\"-m genetic_map.txt \". \n";
								ofs<<"\t\t\"-g  " << gensFileName <<" \".\n";
								ofs <<"\t\t\"-strand_g  " << strandFileName<<" \". \n";
								ofs<< "\t\t\"-int  $lowerBound[$i]  $upperBound[$i] \". \n";
								ofs<<"\t\t\"-Ne 11418  -call_thresh 0.9 -pgs \". \n";
								ofs<< "\t\t\"-o "+ outFileName+".impute2_chunk_$actChunk\"\n";
							ofs<<"\t);\n";
							ofs<< "\t print(\"\\nImputation of chunk $actChunk  is finished.\\n\\n\");\n";
						ofs<<"}\n";
					}
					ofs<<"------------------------------------------------------------\n";
					ofs<<"\nWe have also written all the individual commands one by one as below.\n";
					for(unsigned int i=0;i<lower_range_vec.size();++i)
					{
						ofs<<"------------\n";
						ofs<<" Command "<<i+1<<": \n";
						ofs<<"------------\n";
						ofs<<"./impute2 \\ \n";
						ofs<<"-m genetic_map.txt\\ \n";
						//ofs<<"-h hapmap3_r2_b36/hapmap3.r2.b36.all.chr22.haps\\ \n";
						//ofs<<"-l hapmap3_r2_b36/hapmap3.r2.b36.all.chr22.legend\\ \n";
						ofs<<"-g  " << gensFileName <<" \\  \n";
						ofs <<"-strand_g  " << strandFileName<<"\\ \n";
						ofs<< "-int"<<" "<< lower_range_vec[i] <<"  " << upper_range_vec[i]<<"\\ \n";
						ofs<<"-Ne 11418  -call_thresh 0.9 -pgs\\ \n";
						ofs<< "-o "+ outFileName+".impute2_chunk_"+change_int_into_string(i+1)+"\n";
					}
					 //---------------------------------------------------------------------//
						//find the rs id
					ofs<<"\n If you prefer to impute the whole region at one time then you can use the following command.\n";
					ofs<<"---------------------------------------------------------\n";
					ofs<<"Command for imputing whole region at a time: \n";
					ofs<<"---------------------------------------------------------\n";
					ofs<<"./impute2 \\ \n";
					ofs<<"-m hapmap3_r2_b36/genetic_map_chr22_combined_b36.txt\\ \n";
					//ofs<<"-h hapmap3_r2_b36/hapmap3.r2.b36.all.chr22.haps\\ \n";
					//ofs<<"-l hapmap3_r2_b36/hapmap3.r2.b36.all.chr22.legend\\ \n";
					ofs<<"-g  " << gensFileName <<" \\  \n";
					ofs <<"-strand_g  " << strandFileName<<"\\ \n";
					ofs<< "-int"<<" "<< lower_range_vec[0] <<"  " << upper_range_vec[upper_range_vec.size()-1]<<"\\ \n";
					ofs<<"-Ne 11418  -call_thresh 0.9 -pgs -allow_large_regions\\ \n";
					ofs<< "-o "+ outFileName+".impute2\n";
				
				 }
				
			}
			//--------------------------------------------------------------------//
			//-------------------------------------------------------------------------//
			temp_vec_bp1.clear();
			lower_range_vec.clear();
			upper_range_vec.clear();
			interval_size_MB.clear();
				
			//second  case for hapmap imputation
			if(sifo_tfvec&&_given_extra_snpinfo)
		 	{
				ofs<<"\n-----------------------------------------------------------------\n";
				ofs<<"Templates of IMPUTE commands for imputation  using references panel:\n";
				ofs<<"-----------------------------------------------------------------\n";
				temp_vec_bp1 =CBSNP::sifo_bp;
		 		sort(temp_vec_bp1.begin(),temp_vec_bp1.end());
		 		fix_ints_4_impute_command(temp_vec_bp1, lower_range_vec,upper_range_vec, interval_size_MB);	
				if(lower_range_vec.size()==1)
				{
					//
					ofs<<"The impute commands written below is  just a template to show  how you can impute your genotype data with IMPUTE. You may have to modify it ";
					ofs<<" at necessary places ((especially the input data names)) so that the command can be executed with IMPUTE2 properly.\n";
					ofs<<"------------\n";
					ofs<<" Command: \n";
					ofs<<"------------\n";
					ofs<<"./impute2 \\ \n";
					ofs<<"-m hapmap3_r2_b36/genetic_map_chr_combined_b36.txt\\ \n";
					ofs<<"-h hapmap3_r2_b36/hapmap3.r2.b36.all.chr.haps\\ \n";
					ofs<<"-l hapmap3_r2_b36/hapmap3.r2.b36.all.chr.legend\\ \n";
					ofs<<"-g  " << gensFileName <<" \\  \n";
					ofs <<"-strand_g  " << strandFileName<<"\\ \n";
					ofs<< "-int"<<" "<< lower_range_vec[0] <<"  " << upper_range_vec[0]<<"\\ \n";
					ofs<<"-Ne 11418  -call_thresh 0.9 -pgs \\ \n";
					ofs<< "-o "+ outFileName+".impute2 \n";
				}
				else
				{
					// whole genome imputation
					ofs<<"Your genotype Data has a region of larger than 7Mb.IMPUTE software providers recommend to split the genotyping region into samll chunks for imputation purpose. ";
					ofs<<"These chunks can be  be imputed in parallel on multiple computer processors.\n";
					ofs<<"The following table  describes the lower are upper boudaries (in base pair position) of the region in which imputation should be carried out.\n";
					ofs<<"Moreover, the upper and lower limits which should be given with the command \"-int\" are based on the legend file you have given  with \"--snpinfo\" command. Fcgene assumes that the legend file is based on your reference panel.\n";
					//-----------------------------------------//
					ofs<<"------------------------------------------------------------------------\n";
					ofs<<"Table containing lower and upper boundaries of chunks: \n";
					ofs<<"------------------------------------------------------------------------\n";
					ofs <<"chunk No:\t"<<"lower boundary\t"<< "Upper boundary\t"<<"Size of Chunk in MB\n";
					for(unsigned int i=0;i<lower_range_vec.size();++i)
							{

								ofs << "chunk "<<i+1<<":\t"<< lower_range_vec[i]<< "\t"<<upper_range_vec[i]<<"\t   "<< interval_size_MB[i] << endl;
							}
					ofs<<"------------------------------------------------------------------------\n\n";
					ofs<<"The impute commands written below are just the templates to show  how you can impute your genotype data with IMPUTE. You can modify them ";
					ofs<< "at necessary places (especially the names of input data) so that the commands can be executed with IMPUTE properly.\n";
					//perl command 
					{	
						ofs<<"------------------------------------------------------------\n";
						ofs<<"perl scripts to give impute commands in Linux cluster: \n";
						ofs<<"------------------------------------------------------------\n";
						ofs << "#!/usr/bin/perl -w\n";
						ofs<< "my @lowerBound=(";
						 for(unsigned int i=0;i<(lower_range_vec.size()-1);++i)
						 ofs<<lower_range_vec[i] << ", ";
						 ofs<<lower_range_vec[lower_range_vec.size()-1];
						 ofs<< ");\n";
						 ofs<< "my @upperBound= (";
						 //upper bound
						 for(unsigned int i=0;i<(upper_range_vec.size()-1);++i)
						 ofs<<upper_range_vec[i] << ", ";
						 ofs<<upper_range_vec[upper_range_vec.size()-1];
						 ofs<< ");\n";
						 ofs << "my $nChunks = "<<upper_range_vec.size()<<";\n";
						 ofs<<"my @segs = (0 .. $nChunks-1);\n";
						 ofs<< "foreach $i (@segs)\n"; 
						 ofs <<"{\n";
							ofs<<"\tmy $actChunk=$i+1;\n";
							ofs<<"\tsystem (\n";
								ofs<<"\t\t\"./impute2 \". \n";
								ofs<<"\t\t\"-m hapmap3_r2_b36/genetic_map_chr_combined_b36.txt \". \n";
								ofs<<"\t\t\"-h hapmap3_r2_b36/hapmap3.r2.b36.all.chr.haps \".\n";
								ofs<<"\t\t\"-l hapmap3_r2_b36/hapmap3.r2.b36.all.chr.legend \". \n";
								ofs<<"\t\t\"-g  " << gensFileName <<" \".\n";
								ofs <<"\t\t\"-strand_g  " << strandFileName<<" \". \n";
								ofs<< "\t\t\"-int  $lowerBound[$i]  $upperBound[$i] \". \n";
								ofs<<"\t\t\"-Ne 11418  -call_thresh 0.9 -pgs \". \n";
								ofs<< "\t\t\"-o "+ outFileName+".impute2_chunk_$actChunk\"\n";
							ofs<<"\t);\n";
							ofs<< "\tprint(\"\\nImputation of chunk $actChunk  is finished.\\n\\n\");\n";
						ofs<<"}\n";
					}
					ofs<<"-------------------------------------------------------------------\n";
					ofs<<"Individual commands  for each chunks are given below.\n";	
				
					for(unsigned int i=0;i<lower_range_vec.size();++i)
					{
						ofs<<"------------\n";
						ofs<<" Command "<<i+1<<": \n";
						ofs<<"------------\n";
						ofs<<"./impute2 \\ \n";
						ofs<<"-m hapmap3_r2_b36/genetic_map_chr_combined_b36.txt\\ \n";
						ofs<<"-h hapmap3_r2_b36/hapmap3.r2.b36.all.chr.haps\\ \n";
						ofs<<"-l hapmap3_r2_b36/hapmap3.r2.b36.all.chr.legend\\ \n";
						ofs<<"-g  " << gensFileName <<" \\  \n";
						ofs <<"-strand_g  " << strandFileName<<"\\ \n";
						ofs<< "-int"<<" "<< lower_range_vec[i] <<"  " << upper_range_vec[i]<<"\\ \n";
						ofs<<"-Ne 11418  -call_thresh 0.9 -pgs \\ \n";
						ofs<< "-o "+ outFileName+".impute2_chunk_"+change_int_into_string(i+1)+"\n";
					 }
					 //---------------------------------------------------------------------//
						//find the rs id
					ofs<<"\n If you prefer to impute the whole region at one time then you can use the following command.\n";
					ofs<<"----------------------------------------------------------\n";
					ofs<<"Command for imputing whole region at a time: \n";
					ofs<<"----------------------------------------------------------\n";
					ofs<<"./impute2 \\ \n";
					ofs<<"-m hapmap3_r2_b36/genetic_map_chr22_combined_b36.txt\\ \n";
					ofs<<"-h hapmap3_r2_b36/hapmap3.r2.b36.all.chr22.haps\\ \n";
					ofs<<"-l hapmap3_r2_b36/hapmap3.r2.b36.all.chr22.legend\\ \n";
					ofs<<"-g  " << gensFileName <<" \\  \n";
					ofs <<"-strand_g  " << strandFileName<<"\\ \n";
					ofs<< "-int"<<" "<< lower_range_vec[0] <<"  " << upper_range_vec[upper_range_vec.size()-1]<<"\\ \n";
					ofs<<"-Ne 11418  -call_thresh 0.9 -pgs -allow_large_regions\\ \n";
					ofs<< "-o "+ outFileName+".impute2\n";
				 }
				


		 	//-------------------------------------------------------------------------//
		 	}
		// writing command fiel
	 	ofs.clear(); ofs.close();
	 	printLIN( "\t**->A file containing impute command(s)  has been written out and saved as \""+commandFileName+"\". The command(s) can be used to impute your genotype data. \n");
return ;}

void CBPED::write_beagleGenotypeFile(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName)
{
	//----------------------------------#--------------------------------------#
	printLIN("*->Writing files in "+CBSNP::_out_format+" format: \n");

	//----------------------------------#--------------------------------------#
	// writing bgl file.
	//----------------------------------#--------------------------------------#
	ofstream ofs;
	ofs.clear();ofs.close();
	//
	ofstream ofs1;
	ofs1.clear();ofs1.close();
	
	CBPED* pCBPED =pedInfo[0];
	CBSNP* pCBSNP	=genInfo[0];
	string bglFileName=outFileName+".bgl";
	string markerFileName=outFileName+".snps";
	
	ofs.open(bglFileName.c_str());
	if(!ofs)
		error("out file "+outFileName+" does not exit or can not be written out.\n");
	ofs.clear();
	
	//for marker file 
	ofs1.open(markerFileName.c_str());
	if(!ofs1)
		error("out file "+outFileName+" does not exit or can not be written out.\n");
	ofs1.clear();
	
	
	ofs<< "I" <<" "<< "id"<< " ";
	for(unsigned int i=0; i<pedInfo.size();++i)
	{
		pCBPED =pedInfo[i];
		if(pCBPED->quality)
			ofs <<pCBPED->indId<< " "<<pCBPED->indId<< " ";
	}
	ofs<<"\n";
	//end of first line
	ofs<<"A"<<"  "<< "phenotype" << " ";
	for(unsigned int i=0; i<pedInfo.size();++i)
	{
		pCBPED =pedInfo[i];
		if(pCBPED->quality)
		{
			if(pCBPED->pheno)
				ofs <<2<< " "<<2<< " ";
			if(!pCBPED->pheno && !pCBPED->miss_pheno)
				ofs <<1<< " "<<1<< " ";
			if(pCBPED->miss_pheno)
				ofs <<0<< " "<<0<< " ";
		}
	}
	ofs<<"\n";
	//end of second line.
		 	//writing genotypes
	for(unsigned int i=0;i<genInfo.size();i++)
	{

		pCBSNP	=genInfo[i];
		if(pCBSNP->quality)
		{
			if(pCBSNP->geno1.size()!= pCBSNP->geno2.size() )
			{
				string msg= " Problem in writing line"+change_int_into_string(i)+"! SNP \""+pCBSNP->snpName+"\" was not saved correctly while reading. Different number of individual alleles\n";
				error(msg);
			}
			ofs<< "M"<<" ";
			if(pCBSNP->rsId!="")
			{
				ofs<< pCBSNP->rsId<<" " ;
				ofs1<< pCBSNP->rsId<<" " ;
			}	
			else{
			
				ofs<< pCBSNP->snpName<<" " ;
				ofs1<< pCBSNP->snpName<<" " ;
			}	
			//for the marker list 
			ofs1<<pCBSNP->bp <<" "<<pCBSNP->allele1 << " "<<pCBSNP->allele2<<"\n";
			ofs1.close(); 
			ofs1.clear();
			for(unsigned int j=0; j<pCBSNP->geno1.size();++j)
			{
				pCBPED=pedInfo[j];
				if(pCBPED->quality)
				{
					//---------------------------------------------------//
					if(pCBSNP->allele1!=""&& pCBSNP->allele2!="")
					{
						if(!pCBSNP->geno1[j] && !pCBSNP->geno2[j])
							ofs <<pCBSNP->allele1  <<" " <<pCBSNP->allele1<<" ";
						else if(pCBSNP->geno1[j] && pCBSNP->geno2[j])
							ofs <<pCBSNP->allele2  <<" " <<pCBSNP->allele2<<" ";
						else if(!pCBSNP->geno1[j] && pCBSNP->geno2[j])
						{
							if(pCBSNP->aOrder[j])
								ofs <<pCBSNP->allele1  <<" " <<pCBSNP->allele2<<" ";
							else
								ofs <<pCBSNP->allele2  <<" " <<pCBSNP->allele1<<" ";
						}
						else
							ofs <<"?" <<" " <<"?"<<" ";
					}
					else if (pCBSNP->allele1!=""&& pCBSNP->allele2=="")
					{
						if(!pCBSNP->geno1[j] && !pCBSNP->geno2[j])
							ofs<< pCBSNP->allele1 <<" " <<pCBSNP->allele1<< " "; //<<",  ";
						else if( pCBSNP->geno1[j] &&  !pCBSNP->geno2[j]  ) // if missing
							ofs <<" "<<"?"<< " "<<"?"<< " ";
						else
						{
							string prmsg= pCBSNP->rsId +" can have only homozygous major allele SNP because no other allele is detected for it.\n";
							error(prmsg)	;
						}
					}
					//---------------------------------------------------//
					else
					{
						if(!pCBSNP->geno1[j] && !pCBSNP->geno2[j])
							ofs <<"A"  <<" " <<"A"<<" ";
						else if(pCBSNP->geno1[j] && pCBSNP->geno2[j])
							ofs <<"B" <<" " <<"B"<<" ";
						else if(!pCBSNP->geno1[j] && pCBSNP->geno2[j])
						{
							if(pCBSNP->aOrder[j])
								ofs <<"A" <<" " <<"B"<<" ";
							else
								ofs <<"B"  <<" " <<"A"<<" ";
						}
						else
							ofs <<"?" <<" " <<"?"<<" ";
					}
							//---------------------------------------------------//
				} //end of if quality pedinfo
			} //end of for loop j
			ofs << "\n";
		}

		 	}
		 ofs.clear(); ofs.close();

		 string msg="*->Genotype data for Beagle has been saved as \""+bglFileName+".\".\n";
		 printLIN(msg);
return;	}
//writing bimbam output geno file
//--------------------------------------------------//
void CBPED::write_bimbam_geno_file(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName)
{
	//----------------------------------#--------------------------------------#
	printLIN("*->Writing files in "+CBSNP::_out_format+" format: \n");

	//----------------------------------#--------------------------------------#
	// writing bimbam geno.txt file.
	//----------------------------------#--------------------------------------#
		 	ofstream ofs;
		 	 ofs.clear();ofs.close();
		 	CBPED* pCBPED =pedInfo[0];
		 	CBSNP* pCBSNP	=genInfo[0];

		  string bglFileName=outFileName+".geno.txt";
		 	ofs.open(bglFileName.c_str());
		 	if(!ofs)
		 			error("out file "+outFileName+" does not exit or can not be written out.\n");
		 	ofs.clear();

		 	ofs<< pedInfo.size()<< "\n";
		 	ofs<< genInfo[0]->ngood_snps <<"\n";
		 	//ofs<< genInfo.size()<< "\n";
		 	ofs<< "IND"; // <<",  ";
		 	for(unsigned int i=0; i<pedInfo.size();++i)
		 	{
		 		ofs<<", ";
		 		 pCBPED =pedInfo[i];
		 		if(pCBPED->quality)
		 		{
					if(pCBPED->indId=="")
					{
						string indNewName="indv"+change_int_into_string(i)+"";
						ofs<< indNewName; //<< ",  ";
					}
					else
					ofs <<pCBPED->indId;
				}
		 	}
		 	ofs<<"\n";
		 	//end of first line
		 	//writing genotypes
		 	for(unsigned int i=0;i<genInfo.size();i++)
		 	{
		 		 pCBSNP	=genInfo[i];
		 		if(pCBSNP->quality)
		 		{
					if(pCBSNP->geno1.size()!= pCBSNP->geno2.size() )
					{
						string msg= " Problem in writing line"+change_int_into_string(i)+"! SNP \""+pCBSNP->snpName+"\" was not saved correctly while reading. Different number of individual alleles\n";
						error(msg);
					}
					//writing snp id or rsid
					if(pCBSNP->rsId!="") ofs <<pCBSNP->rsId; //<< ",  ";
					else if(pCBSNP->snpName!="") ofs <<pCBSNP->snpName; // << ",  ";
					else if(pCBSNP->rsId==""&& pCBSNP->snpName=="")
					{
						string snpNewName="snp"+change_int_into_string(i)+"";
						ofs<< snpNewName; //<< ",  ";
					}

					for(unsigned int j=0; j<pCBSNP->geno1.size();++j)
					{
						pCBPED =pedInfo[j];
						if(pCBPED->quality)
						{
							ofs <<", ";
							//-----------------------------------------------------------//
							if(pCBSNP->allele1!=""&& pCBSNP->allele2!="")
							{
								if(!pCBSNP->geno1[j] && !pCBSNP->geno2[j]) // both false homo1 case
									ofs<< pCBSNP->allele1  <<pCBSNP->allele1; //<<",  ";
								else if(pCBSNP->geno1[j] && pCBSNP->geno2[j]) // both true homo2 case
									ofs <<pCBSNP->allele2 <<pCBSNP->allele2; // <<",  ";
								else if(!pCBSNP->geno1[j] && pCBSNP->geno2[j]) //  first false second true hetero case
								{
									if(pCBSNP->aOrder[j])
										ofs <<pCBSNP->allele1 <<pCBSNP->allele2;// <<",  ";
									else
										ofs <<pCBSNP->allele2  <<pCBSNP->allele1; //<<",  ";

								}
								else // missing
									ofs <<"??"; //<<",  "; // missing
							}else if (pCBSNP->allele1!=""&& pCBSNP->allele2=="")
							{
								if(!pCBSNP->geno1[j] && !pCBSNP->geno2[j])
									ofs<< pCBSNP->allele1  <<pCBSNP->allele1; //<<",  ";
								else if(pCBSNP->geno1[j] &&  !pCBSNP->geno2[j]  ) // if missing
									ofs <<"??"; //missing
								else
								{
										string prmsg= pCBSNP->rsId +" can have only homozygoze major SNP because no other allele is detected for it.\n";
										error(prmsg);
								}
							}
							//----------------------------------------------------//
							//---------------------------------------------------//
							else{
								if(!pCBSNP->geno1[j] && !pCBSNP->geno2[j])
									ofs <<"AA";
								else if(pCBSNP->geno1[j] && pCBSNP->geno2[j])
									ofs <<"BB";
								else if(!pCBSNP->geno1[j] && pCBSNP->geno2[j])
								{
									if(pCBSNP->aOrder[j])
										ofs <<"AB";
									else
										ofs <<"BA";
								}
								else
									ofs <<"??";
								}
							//---------------------------------------------------//
						} //endof if pedinfo quality loop
					} //end of j for loop
					ofs << "\n";
		 		}
		 	}
		 ofs.clear(); ofs.close();

		 string msg="\t**->Genotype data for bimbam has been saved as \""+bglFileName+".\".\n";
		 printLIN(msg);
		 //----------------------------------#--------------------------------------#
		 // writing bimbam pheno.txt file.
		 //----------------------------------#--------------------------------------#
		 string bimbam_pheno_fileName=outFileName+".pheno.txt";
		 ofs.open(bimbam_pheno_fileName.c_str());
		 	if(!ofs)
		 		error("out file "+outFileName+" does not exits.\n");
		 	ofs.clear();
		 //writing pedinfo
			for(unsigned int i=0; i<pedInfo.size();++i)
			{
				pCBPED =pedInfo[i];
				if(pCBPED->quality)
				{
					if(pCBPED->pheno)
						ofs <<0<<"\n"; // case
					if(!pCBPED->pheno && !pCBPED->miss_pheno)
						ofs <<1<< "\n";
					if(pCBPED->miss_pheno)
						ofs <<"NA"<<"\n";
				}
			}
		 //closing the file
		 ofs.clear(); ofs.close();
		  msg="\t**->Phenotype data for bimbam has been saved as \""+bimbam_pheno_fileName+".\".\n";
		 printLIN(msg);


return;	}
//--------------------------------------------------//
void CBPED::write_mach_files(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName)
{
		unsigned const int nIndv 	=pedInfo.size();
		unsigned const int nSnp		=genInfo.size();

	//----------------------------------#--------------------------------------#
		 	printLIN("*->Writing files in MaCH format: \n");
		   //----------------------------------#--------------------------------------#
			  	ofstream ofs;
			   ofs.clear();ofs.close();
		 //----------------------------------#--------------------------------------#
			   // writing mach dat file.
		 //----------------------------------#--------------------------------------#
			   string datFileName = outFileName+".dat";
			   ofs.open(datFileName.c_str(),ios::out);
				if(!ofs)
					error("out file "+outFileName+" does not exit.\n");
			   for(unsigned int i=0; i<nSnp;++i)
			   {
			    	if(genInfo[i]->quality)
			    		ofs<< "M" << " "<<genInfo[i]->rsId <<"\n";
			   }
			    	ofs.clear();ofs.close();
			    	printLIN("\t **->MaCH dat file has been written out and saved as  \""+\
			    			datFileName+"\".\n");
			    ofs.clear();ofs.close();

			    //cout<< pedInfo.size();

		//----------------------------------#--------------------------------------#
		 	// writing ped file.
		 //----------------------------------#--------------------------------------#

			  	string pedFileName=outFileName+".ped";
			  	ofs.open(pedFileName.c_str());
			  	if(!ofs)
			  			error("out file "+outFileName+" does not exit.\n");
			  	ofs.clear();
			// cout << pedInfo.size();

			  	if(nIndv==0)
			  	{
			  		//cout << "pedInfo has no elements"<<endl;
			  		//first change genotype probability into geno1 geno2
			  		// cout << "debug: "<< genInfo.size()<< genInfo[0]->geno1.size()<<endl;
			  	 		for(unsigned int i=0; i<(genInfo[0]->geno1.size());++i)
			  		{
			  			//writing ped Info.
			  			ofs << "fam"+change_int_into_string(i)+"" << " "; 	// family id
			  			ofs<< "indv"+change_int_into_string(i)+"" << " "; 	// inv id
			  			ofs<< 0  << " "; 									// pat id
			  			ofs<< 0  << " "; 									//mat  id
			  			ofs<< 0  << "	 "; 								// sex
			  			for(unsigned int j=0; j<nSnp; ++j)
			  			{	CBSNP* pCSNP=genInfo[j];
			  				if(pCSNP->quality)
			  				{
									//-------------------------------------------------------//
								if(pCSNP->allele1!=""&& pCSNP->allele2!="")
								{
									if( !pCSNP->geno1[i] &&  !pCSNP->geno2[i]  ) // if homozygote1  00
										ofs <<" "<< pCSNP->allele1 << " "<<pCSNP->allele1;
									else if(pCSNP->geno1[i] &&  pCSNP->geno2[i]  ) // if homozygote2 11
										ofs <<" " << pCSNP->allele2 << " "<<pCSNP->allele2;
									else if( !pCSNP->geno1[i] &&  pCSNP->geno2[i]  ) // if heterozygote
									{
										//ofs << " "<<pCSNP->allele1 << " "<<pCSNP->allele2;
										if(pCSNP->aOrder[i])
											ofs <<" "<< pCSNP->allele1 << " " <<pCSNP->allele2;
										if(!pCSNP->aOrder[i])
											ofs<<" " <<pCSNP->allele2 << " " <<pCSNP->allele1;
									}
									else if( pCSNP->geno1[i] &&  !pCSNP->geno2[i]  ) // if missing
										ofs <<" "<<CBSNP::_miss_gvalue << " "<<CBSNP::_miss_gvalue;
									else
									{
										error("Problem in writing out "+change_int_into_string(j) +"th SNP.\n");
									}
								}
								else if (pCSNP->allele1!=""&& pCSNP->allele2=="")
								{
										if(!pCSNP->geno1[i] && !pCSNP->geno2[i])
											ofs<<" "<< pCSNP->allele1<< " " <<pCSNP->allele1; //<<",  ";
										else if( pCSNP->geno1[i] &&  !pCSNP->geno2[i]  ) // if missing
											ofs <<" "<<CBSNP::_miss_gvalue << " "<<CBSNP::_miss_gvalue;
										else 
										{
												cout << pCSNP->allele1 << " " <<pCSNP->allele2 << endl;
												string prmsg= pCSNP->rsId +" can have only homozygoze major everywhere because no other allele is detected for it.\n";
												error(prmsg)	;
										}
								}
								//---------------------------------------------------//
								else{
									if(!pCSNP->geno1[i] && !pCSNP->geno2[i])
										ofs <<" "<< "A"  <<" " <<"A";
									else if(pCSNP->geno1[i] && pCSNP->geno2[i])
										ofs <<" "<<"B" <<" " <<"B";
									else if(!pCSNP->geno1[i] && pCSNP->geno2[i])
									{
										if(pCSNP->aOrder[i])
											ofs<<" " <<"A" <<" " <<"B";
										else
											ofs<<" " <<"B"  <<" " <<"A";
									}
									else
										ofs <<" "<<CBSNP::_miss_gvalue  <<" " <<CBSNP::_miss_gvalue;
									}

								//---------------------------------------------------//

							}
			  			}
			  			ofs << "\n";
			  		}

			  }

			  // here comes if pedinfo is already given.
			  	else
			  	{
			  		for(unsigned int i=0; i<nIndv;++i)
			  		{
						CBPED* ppedInfo = pedInfo[i];
						if(ppedInfo->quality)
						{
							if(ppedInfo->famId!="") ofs<<ppedInfo->famId<< " ";
							else
							{
								string newid="fam"+change_int_into_string(i+1);
								ofs << newid << " ";
							}
							if(ppedInfo->indId!="") ofs<<ppedInfo->indId<< " ";
							else
							{
								string newid="indv"+change_int_into_string(i+1);
								ofs << newid << " ";
							}
							if(ppedInfo->patId!="") ofs<<ppedInfo->patId<< " ";
							else
							{
								string newid="0";
								ofs << newid << " ";
							}
							if(ppedInfo->matId!="") ofs<<ppedInfo->matId<< " ";
							else
							{
								string newid="0";
								ofs << newid << " ";
							}
						//if(!ppedInfo->sex && !ppedInfo->miss_sex)
						//{
						//	string newid=0;
						//	ofs << newid << " ";
						//}
						//{
							if(ppedInfo->sex) ofs<<"m"<< " ";
							else if(!ppedInfo->sex && !pedInfo[i]->miss_sex) ofs<<"f"<<" ";
							else if(ppedInfo->miss_sex) ofs<<"0"<< " "; // 0 is missing
						//}


						for(unsigned int j=0;j<nSnp; ++j)
						{
							CBSNP* pCSNP=genInfo[j];
							if(pCSNP->quality)
							{
								if(pCSNP->geno1.size()!=nIndv || pCSNP->geno2.size()!=nIndv)
								{
									cout <<" pedInfo.size(): "<<  pedInfo.size()<<",  "<< "geno1.size() : "<<pCSNP->geno1.size()<<",  " <<"geno2.size():  "<<pCSNP->geno1.size()<< "\n";
													//for(unsigned int i=0;i<pCSNP->geno1.size();++i)
													//	ofs<< pCSNP->geno1[i] << " " <<pCSNP->geno2[i] << " ";
									ofs.close();
									error("Problem in writing out "+change_int_into_string(j) +"th SNP. While reading time, probably the genotypes were saved not correctly.\n");
								}
								//---------------------------------------------------//
								//	cout<< "i: "<<i <<", j: "<<j << ","<< pCSNP->geno1.size()<<","<< pCSNP->geno2.size() <<endl;

								if(pCSNP->allele1!=""&& pCSNP->allele2!="")
								{
									if( !pCSNP->geno1[i] &&  !pCSNP->geno2[i]  ) // if homozygote1  00
										ofs <<" "<< pCSNP->allele1 << " "<<pCSNP->allele1;
									else if( pCSNP->geno1[i] &&  pCSNP->geno2[i]  ) // if homozygote2 11
										ofs <<" " << pCSNP->allele2 << " "<<pCSNP->allele2;
									else if( !pCSNP->geno1[i] &&  pCSNP->geno2[i]  ) // if heterozygote
									{
										if(pCSNP->aOrder[i])
											ofs <<" "<< pCSNP->allele1 << " " <<pCSNP->allele2;
										if(!pCSNP->aOrder[i])
											ofs<<" " <<pCSNP->allele2 << " " <<pCSNP->allele1;
									}
									else if( pCSNP->geno1[i] &&  !pCSNP->geno2[i]  ) // if missing
										ofs <<" "<<CBSNP::_miss_gvalue << " "<<CBSNP::_miss_gvalue;
									else
										{
											error("Problem in writing out "+change_int_into_string(j) +"th SNP.\n");
										}

								}

								else if (pCSNP->allele1!=""&& pCSNP->allele2=="")
								{
									if(!pCSNP->geno1[i] && !pCSNP->geno2[i])
										ofs<<" "<< pCSNP->allele1<<" " <<pCSNP->allele1; //<<",  ";
									else if( pCSNP->geno1[i] &&  !pCSNP->geno2[i]  ) // if missing
										ofs <<" "<<CBSNP::_miss_gvalue << " "<<CBSNP::_miss_gvalue;
									else
									{
											cout << pCSNP->geno1[i] <<  " " << pCSNP->geno2[i] <<endl;
											cout << pCSNP->allele1 << " " <<pCSNP->allele2 << endl;
											string prmsg= pCSNP->rsId +" can have only major homozygoze  SNP because no other allele is detected for it.\n";
											error(prmsg)	;
									}
								}
								//---------------------------------------------------//
								else{
										if(!pCSNP->geno1[i] && !pCSNP->geno2[i])
											ofs <<" "<< "A"  <<" " <<"A";
										else if(pCSNP->geno1[i] && pCSNP->geno2[i])
											ofs <<" "<<"B" <<" " <<"B";
										else if(!pCSNP->geno1[i] && pCSNP->geno2[i])
										{
											if(pCSNP->aOrder[i])
												ofs<<" " <<"A" <<" " <<"B";
											else
												ofs<<" " <<"B"  <<" " <<"A";
										}
										else
											ofs <<" "<<CBSNP::_miss_gvalue  <<" " <<CBSNP::_miss_gvalue;
									}
								//---------------------------------------------------//

							}

						}

						ofs<< "\n";
						}//end of if (ppedInfo->quality)
			  		} //end of ppedInfo for loop
			  	}

			  	printLIN("\t**->MaCH ped file has been written out and saved as \""+pedFileName+\
			  			"\".\n");
			  	ofs.clear();	 ofs.close();
			  	//----------------------------------#--------------------------------------#





return;}
//write mach references
//------------------------------
//--------------------------------------------------//
void CMPED::write_mach_ref_files(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName)
{
		unsigned const int nIndv 	=pedInfo.size();
		unsigned const int nSnp		=genInfo.size();
		//cout <<nSnp <<endl;//test

	//----------------------------------#--------------------------------------#
		// 	printLIN("*->Writing files in "+CBSNP::_out_format+" format: \n");
		   //----------------------------------#--------------------------------------#
			  	ofstream ofs;
			   ofs.clear();ofs.close();
		 //----------------------------------#--------------------------------------#
			   // writing mach dat file.
		 //----------------------------------#--------------------------------------#
			   string datFileName = outFileName+"_machRef.snps";
			   ofs.open(datFileName.c_str(),ios::out);
				if(!ofs)
					error("out file "+outFileName+" does not exit.\n");
			   for(unsigned int i=0; i<nSnp;++i)
			   {
			    	if(genInfo[i]->quality)
			    		ofs<< genInfo[i]->rsId <<"\n";
			   }
			    	ofs.clear();ofs.close();
			    	printLIN("\t **->MACH reference SNP  file has been written out and saved as  \""+\
			    			datFileName+"\".\n");
			    ofs.clear();ofs.close();
		//----------------------------------#--------------------------------------#
		 // writing ped file.
		 //----------------------------------#--------------------------------------#

			  	string pedFileName=outFileName+"_machRef.haps";
			  	ofs.open(pedFileName.c_str());
			  	if(!ofs)
			  			error("out file "+outFileName+" does not exit.\n");
			  	ofs.clear();

			  		for(unsigned int i=0; i<nIndv;++i)
			  		{
						CBPED* ppedInfo = pedInfo[i];
						if(ppedInfo->quality)
						{

							if(ppedInfo->indId!="") ofs<<ppedInfo->indId<<"_A"<< " "<<ppedInfo->indId<<"_A"<<" ";
							else
							{
								string newid="indv"+change_int_into_string(i+1);
								ofs << newid << "_A"<<" "<<newid << "_A"<< " ";
							}


						for(unsigned int j=0;j<nSnp; ++j)
						{
							CBSNP* pCSNP=genInfo[j];
							//cout<< boolalpha<< pCSNP->quality<< endl;
							if(pCSNP->quality)
							{
								if(pCSNP->geno1.size()!=nIndv || pCSNP->geno2.size()!=nIndv)
								{
									cout <<" pedInfo.size(): "<<  pedInfo.size()<<",  "<< "geno1.size() : "<<pCSNP->geno1.size()<<",  " <<"geno2.size():  "<<pCSNP->geno1.size()<< "\n";
													//for(unsigned int i=0;i<pCSNP->geno1.size();++i)
													//	ofs<< pCSNP->geno1[i] << " " <<pCSNP->geno2[i] << " ";
									ofs.close();
									error("Problem in writing out "+change_int_into_string(j) +"th SNP. While reading time, probably the genotypes were saved not correctly.\n");
								}
								//---------------------------------------------------//
								//cout << pCSNP->allele1 << ", and  "<<pCSNP->allele2 <<endl; // test
								if(pCSNP->allele1!=""&& pCSNP->allele2!="")
								{
									if( !pCSNP->geno1[i] &&  !pCSNP->geno2[i]  ) // if homozygote1  00
										ofs<< pCSNP->allele1; //<< " "<<pCSNP->allele1;
									else if( pCSNP->geno1[i] &&  pCSNP->geno2[i]  ) // if homozygote2 11
										ofs << pCSNP->allele2 ; // << " "<<pCSNP->allele2;
									else if( !pCSNP->geno1[i] &&  pCSNP->geno2[i]  ) // if heterozygote
									{
										if(pCSNP->aOrder[i])
											ofs << pCSNP->allele1; // << " " <<pCSNP->allele2;
										if(!pCSNP->aOrder[i])
											ofs<<pCSNP->allele2; // << " " <<pCSNP->allele1;
									}
									else if( pCSNP->geno1[i] &&  !pCSNP->geno2[i]  ) // if missing
										ofs<<"0"; //<< " "<<CBSNP::_miss_gvalue;
									else
										{
											error("Problem in writing out "+change_int_into_string(j) +"th SNP.\n");
										}


								}


							}
						}
					ofs<< "\n";

					// til here we have written the first haplotype obtained form impute references.
					//writing the second haplotype of impute reference
					// for this: we need to write id
					if(ppedInfo->indId!="") ofs<<ppedInfo->indId<<"_B"<< " "<<ppedInfo->indId<<"_B"<<" ";
					else
					{
						string newid="indv"+change_int_into_string(i+1);
						ofs << newid << "_B"<<" "<<newid << "_B"<<" ";
					}

					for(unsigned int j=0;j<nSnp; ++j)
					{
						CBSNP* pCSNP=genInfo[j];
						if(pCSNP->quality)
						{
							if(pCSNP->geno1.size()!=nIndv || pCSNP->geno2.size()!=nIndv)
							{
								cout <<" pedInfo.size(): "<<  pedInfo.size()<<",  "<< "geno1.size() : "<<pCSNP->geno1.size()<<",  " <<"geno2.size():  "<<pCSNP->geno1.size()<< "\n";
								//for(unsigned int i=0;i<pCSNP->geno1.size();++i)
								//	ofs<< pCSNP->geno1[i] << " " <<pCSNP->geno2[i] << " ";
								ofs.close();
								error("Problem in writing out "+change_int_into_string(j) +"th SNP. While reading time, probably the genotypes were saved not correctly.\n");
							}
							//---------------------------------------------------//
							if(pCSNP->allele1!=""&& pCSNP->allele2!="")
							{
								if( !pCSNP->geno1[i] &&  !pCSNP->geno2[i]  ) // if homozygote1  00
									ofs<< pCSNP->allele1; //<< " "<<pCSNP->allele1;
								else if( pCSNP->geno1[i] &&  pCSNP->geno2[i]  ) // if homozygote2 11
									ofs<< pCSNP->allele2 ; // << " "<<pCSNP->allele2;
								else if( !pCSNP->geno1[i] &&  pCSNP->geno2[i]  ) // if heterozygote
								{
									if(pCSNP->aOrder[i])
										ofs<< pCSNP->allele2; // << " " <<pCSNP->allele2;
									if(!pCSNP->aOrder[i])
										ofs<<pCSNP->allele1; // << " " <<pCSNP->allele1;
								}
								else if( pCSNP->geno1[i] &&  !pCSNP->geno2[i]  ) // if missing
									ofs<<"0"; //<< " "<<CBSNP::_miss_gvalue;
								else
								{
									error("Problem in writing out "+change_int_into_string(j) +"th SNP.\n");
								}

							}


						}
					}
					ofs<< "\n";

					//----------------------------------------------
					//end of writing haplotye for each person.
					}//end of if (ppedInfo->quality)
			  		} //end of ppedInfo for loop

			  	printLIN("\t**->MACH reference hap  file has been written out and saved as \""+pedFileName+\
			  			"\".\n");
			  	ofs.clear();	 ofs.close();
			  	//----------------------------------#--------------------------------------#





return;}
//write mach references

//-----------------------------------------------------//
//new funcitons
bool  CBPED::lese_pedfile_zeile(IFS &pf, vector<string>& pZeile){
	string strng(""); //  zum spalte speichern
	bool ende =false;
	bool multi_whChar=false; // this will detect multi white space characters and avoid them.
	long int spaltenZahl=0;
	while(true)
	{
		if(_ifs_eof(pf))
		break;
		char myChar=_ifs_getc(pf);

		if((myChar=='\t')||(myChar==' ')|| (myChar=='\n')|| (myChar=='\r')|| _ifs_eof(pf))
		{
			if((myChar=='\n')||(myChar=='\r') || _ifs_eof(pf))
			{
				if(!ende && !multi_whChar && !strng.empty()&& strng!=""&& strng!=" \n")
				{
					pZeile.push_back(strng);
					spaltenZahl++;
					//cout << "spaltenZahl am ende der Zeile: "<< spaltenZahl <<": "<<strng <<endl;
					spaltenZahl=0;
					ende=true;
					//for(unsigned int i=0;i<pZeile.size();++i)
					//cout << pZeile[i]<< " ";
					return ende;
				}
				else ende=true;
					return ende;
			}
			if((myChar=='\t')||(myChar==' '))
			{
				if(!ende&& !multi_whChar && !strng.empty()&& strng!="" && strng!=" \n" )
				{
					pZeile.push_back(strng);
					spaltenZahl++;

					//cout << "spaltenZahl am ende des spaltes : "<< spaltenZahl <<": "<< strng <<endl; // for test only
				}
				multi_whChar=true;
				strng="";
				//continue;

			}


		}
		else if (myChar=='#')
		{

			ende=true;
	     // Ignore rest of line and advance to next line
			//fscanf (pf, "%*[^\n]");
			//(void) fgetc (pf);
			while(myChar!='\n'&& myChar!='\r'&& !_ifs_eof(pf))
						myChar =_ifs_getc(pf);
			//if(_ifs_eof(pf) || (myChar=='\r') || (myChar=='\n'))
			if(!strng.empty())
			{
				pZeile.push_back(strng);
				++spaltenZahl;
			}

		}
		else
		{
		strng+=myChar;
		multi_whChar=false;
		if(strng.empty())
			break;
		if(strng=="")
			break;
		}

	}
return ende;}
// This function is to read genotype  AA or A/A and to save them.
void CBPED::lese_pedfile_genInfo(const vector<CBSNP*>&genInfo, const vector<string>& zeilenVec,const int akt_zeile)
{
	string first	="";
	string second	="";
	string third	="";
	bool fatal=false;
	string EMSG ="\t**->probblem with reading line "+change_int_into_string(akt_zeile)+":\n";
	string errmsg="";
	//cout << "test" <<endl;
	//cout << genInfo.size() << " zeilen vec: " << zeilenVec.size()<<endl;
	if(zeilenVec.size()==genInfo.size())
	{
		for(unsigned int i=0;i<zeilenVec.size();++i)
		{

			CBSNP * pCBSNP=genInfo[i];
			third=zeilenVec[i];
			//cout << third << endl; // only for test
			//cout << "/"<<endl;
			if(third.length() !=3) // || third[2]!="/")
			{
				errmsg="Problem in reading  SNP "+ pCBSNP->snpName+" in line"+change_int_into_string(akt_zeile)+". The SNP has genotype [ " + first + " ] and ["+second+"],  each should be only one characters long\n";
				error(errmsg);

			}
			// if okay now parse the string
			first=third[0];
			second =third[2];
			//cout << " " << first << "  " <<second <<" "<<endl; // for test only
			if(first!=CBSNP::_miss_gvalue&& second!=CBSNP::_miss_gvalue)
					{
						if(first==second)
						{
							// add allele
							if(pCBSNP->allele1==""&& pCBSNP->allele2!=first)pCBSNP->allele1=first;
							else if (pCBSNP->allele2==""&& pCBSNP->allele1!=first) pCBSNP->allele2=first;

							// add genotype
							if(first==pCBSNP->allele1&& second==pCBSNP->allele1 )
							{// this is the case where homozygote1 is denoted by 00
								//cout << "reading homo1 \n";
								//cout << first << " "<<second<<"\n";
								pCBSNP->geno1.push_back(false);
								pCBSNP->geno2.push_back(false);
								pCBSNP->aOrder.push_back(true);
							}
							else if(first==pCBSNP->allele2 && second==pCBSNP->allele2)
							{// this is the case where homozygote 2  is denoted by 11
								//cout << "reading homo2\n";// for test only
								//cout << first << " "<<second<<"\n";// for test only
								pCBSNP->geno1.push_back(true);
								pCBSNP->geno2.push_back(true);
								pCBSNP->aOrder.push_back(true);
							}
							else
							{
								string wrongAllele="";
								if(first!=pCBSNP->allele1 && first!=pCBSNP->allele2 )
									wrongAllele=first;
								else if(second!=pCBSNP->allele1 && second!=pCBSNP->allele2 )
									wrongAllele=second;
							   errmsg="SNP "+pCBSNP->snpName+" containing two alleles "+\
								pCBSNP->allele1 + " and "+pCBSNP->allele2+ " has found a new allele type  "+wrongAllele+ "  at individual "+\
								change_int_into_string(akt_zeile)+ ".\n ";
							  fatal=true;
							}
						}
						else if(first!=second) // save it as 01 .i.e. false true.
						{// this is the case where heterozygote   is denoted by 01
							//cout << "reading heterozygote\n";// for test only
							//cout << first << " "<<second<<"\n";// for test only
							if(pCBSNP->allele1==""&&pCBSNP->allele2!=first )pCBSNP->allele1=first;
							else if (pCBSNP->allele2==""&&pCBSNP->allele1!=first ) pCBSNP->allele2=first;
							if (pCBSNP->allele1==""&& pCBSNP->allele2!=second ) pCBSNP->allele1=second;
							if (pCBSNP->allele2==""&& pCBSNP->allele1!=second ) pCBSNP->allele2=second;
							//now check if allele are correct in first and second
							if(first==pCBSNP->allele1&& second==pCBSNP->allele2 )
							{
							 	pCBSNP->geno1.push_back(false);
								pCBSNP->geno2.push_back(true);
								pCBSNP->aOrder.push_back(true);
							}
							else if(first==pCBSNP->allele2&& second==pCBSNP->allele1 )
							 	 {
									pCBSNP->geno1.push_back(false);
									pCBSNP->geno2.push_back(true);
									pCBSNP->aOrder.push_back(false);
							 	 }
							else
							{
								cout<< first<< " " << second << " \n";
								cout << pCBSNP->allele1 << " "<< pCBSNP->allele2 << "\n";
								string wrongAllele="";
								if((first!=pCBSNP->allele1 ) && first!=pCBSNP->allele2 )
								wrongAllele=first;
								else if(second!=pCBSNP->allele1 && second!=pCBSNP->allele2 )
								wrongAllele=second;
								errmsg="SNP "+pCBSNP->snpName+" containing two alleles "+\
										pCBSNP->allele1+ " and "+pCBSNP->allele2+ " has found a new allele type "+wrongAllele+ "  at individual "+\
										change_int_into_string(akt_zeile)+ ".\n ";
								 fatal=true;

							}


						}
					}
				// This is the case where both genotypes are missing
				else if (first==CBSNP::_miss_gvalue || second==CBSNP::_miss_gvalue )
				{ //In this case we use the code 10 to denote the missing genotypes
						//cout << "reading missing \n"; // for test only
						//cout << first << " "<<second<<"\n"; // for test only
						pCBSNP->geno1.push_back(true);
						pCBSNP->geno2.push_back(false);
						pCBSNP->aOrder.push_back(true);
				}
				// this is the case where non of above condition is true then error.
				 if (fatal)
				 {
					 cerr<< EMSG;
					 LIN << EMSG;
					 error(errmsg);

				 }
				// check if size of geno1 and geno2 are equal they must be equal.
					if(pCBSNP->geno1.size()!= pCBSNP->geno2.size() )
						error(EMSG);


		 //-------------------------------------------------------------#
		}// end of first for loop
		//cout << "It seems like that you have  genotypes not in the form \"A  B\", but in the form \"A/B\". Please let us know if you want to convert \"A/B\" types of genotypes. We will update the program for you to convert such type of genotypes.\n ";
		//cout << "email:nab.raj.roshyara@imise.uni-leipzig.de \n";
		//exit(1);
		// this is the case where snps info is given as  AB
	}
	else if((zeilenVec.size()/2)==genInfo.size())
	{ // this is the case where snp info is given as: A B
		for(unsigned int i=1;i<=zeilenVec.size()/2;++i)
		{
			first		=zeilenVec[2*i-2];
			second		=zeilenVec[2*i-1];
		//	cout << " " << first << "  " <<second ; // for test only
				CBSNP * pCBSNP=genInfo[i-1];
			if(first.length() != second.length()&& (first.length() !=1 && second.length()!=1) )
				error("Problem in reading  SNP "+ pCBSNP->snpName+" in line"+change_int_into_string(akt_zeile)+". The SNP has genotype [ " + first + " ] and ["+second+"],  each should be only one characters long\n");
			// this is the case where both genotypes are non missing
			//cout << CBSNP::_miss_gvalue<<endl;
			if(first!=CBSNP::_miss_gvalue&& second!=CBSNP::_miss_gvalue)
			{
				if(first==second)
				{
					// add allele
					if(pCBSNP->allele1==""&& pCBSNP->allele2!=first)pCBSNP->allele1=first;
					else if (pCBSNP->allele2==""&& pCBSNP->allele1!=first) pCBSNP->allele2=first;

					// add genotype
					if(first==pCBSNP->allele1&& second==pCBSNP->allele1 )
					{// this is the case where homozygote1 is denoted by 00
						//cout << "reading homo1 \n";
						//cout << first << " "<<second<<"\n";
						pCBSNP->geno1.push_back(false);
						pCBSNP->geno2.push_back(false);
						pCBSNP->aOrder.push_back(true);
					}
					else if(first==pCBSNP->allele2 && second==pCBSNP->allele2)
					{// this is the case where homozygote 2  is denoted by 11
						//cout << "reading homo2\n";// for test only
						//cout << first << " "<<second<<"\n";// for test only
						pCBSNP->geno1.push_back(true);
						pCBSNP->geno2.push_back(true);
						pCBSNP->aOrder.push_back(true);
					}
					else
					{
						string wrongAllele="";
						if(first!=pCBSNP->allele1 || first!=pCBSNP->allele2 )
							wrongAllele=first;
						else if(second!=pCBSNP->allele1 || second!=pCBSNP->allele2 )
							wrongAllele=second;
					   errmsg="SNP "+pCBSNP->snpName+" containing two alleles "+\
						pCBSNP->allele1+ " and "+pCBSNP->allele2+ " has found a new allele type"+wrongAllele+ "  at individual "+\
						change_int_into_string(akt_zeile)+ ".\n ";
					  fatal=true;
					}
				}
				else if(first!=second) // save it as 01 .i.e. false true.
				{// this is the case where heterozygote   is denoted by 01
					//cout << "reading heterozygote\n";// for test only
					//cout << first << " "<<second<<"\n";// for test only
					if(pCBSNP->allele1==""&&pCBSNP->allele2!=first )pCBSNP->allele1=first;
					else if (pCBSNP->allele2==""&&pCBSNP->allele1!=first ) pCBSNP->allele2=first;
					if (pCBSNP->allele1==""&& pCBSNP->allele2!=second ) pCBSNP->allele1=second;
					if (pCBSNP->allele2==""&& pCBSNP->allele1!=second ) pCBSNP->allele2=second;
					//now check if allele are correct in first and second
					if(first==pCBSNP->allele1&& second==pCBSNP->allele2 )
					{
					 	pCBSNP->geno1.push_back(false);
						pCBSNP->geno2.push_back(true);
						pCBSNP->aOrder.push_back(true);
					}
					else if(first==pCBSNP->allele2&& second==pCBSNP->allele1 )
					 	 {
							pCBSNP->geno1.push_back(false);
							pCBSNP->geno2.push_back(true);
							pCBSNP->aOrder.push_back(false);
					 	 }
					else
					{
						string wrongAllele="";
						if(first!=pCBSNP->allele1 || first!=pCBSNP->allele2 )
						wrongAllele=first;
						else if(second!=pCBSNP->allele1 || second!=pCBSNP->allele2 )
						wrongAllele=second;
						errmsg="SNP "+pCBSNP->snpName+" containing two alleles "+\
								pCBSNP->allele1+ " and "+pCBSNP->allele2+ " has found a new allele type"+wrongAllele+ "  at individual "+\
								change_int_into_string(akt_zeile)+ ".\n ";
						 fatal=true;

					}


				}
			}
		// This is the case where both genotypes are missing
		else if (first==CBSNP::_miss_gvalue || second==CBSNP::_miss_gvalue )
		{ //In this case we use the code 10 to denote the missing genotypes
				//cout << "reading missing \n"; // for test only
				//cout << first << " "<<second<<"\n"; // for test only
				pCBSNP->geno1.push_back(true);
				pCBSNP->geno2.push_back(false);
				pCBSNP->aOrder.push_back(true);
		}
		// this is the case where non of above condition is true then error.
		 if (fatal)
		 {
			 cerr<< EMSG;
			 LIN<< EMSG;
			 error(errmsg);
		 }
		// check if size of geno1 and geno2 are equal they must be equal.
			if(pCBSNP->geno1.size()!= pCBSNP->geno2.size() )
				error(EMSG);
			/*if((int)pCBSNP->geno1.size()!= akt_zeile)
			{
				cout <<"genotypes: "<< first << " " <<second<< "\n";// only for test
				cout << "geno1.size()"<<pCBSNP->geno1.size()<<endl;
				LIN <<"geno1.size()"<<pCBSNP->geno1.size()<< " genotypes: "<< first << " " <<second<< "\n";// only for test
				//
				cout <<"sizeof zeilenVec: "<<zeilenVec.size()<< "\n";// only for test
				LIN <<"sizeof zeilenVec: "<<zeilenVec.size()<< "\n";// only for test
				//
				cout <<"sizeof genInfo: "<<genInfo.size()<< "\n";// only for test
				LIN <<"sizeof genInfo: "<<genInfo.size()<< "\n";// only for test
				//for (unsigned int j=1;j<zeilenVec.size();++j)
				//	LIN<<" "<< zeilenVec[j] ;
				//error(EMSG);
				break;
			}
			*/


		}
	}

	return ;}
CBPED*  CBPED::lese_pedfile_pedInfo(const vector<string>& zeilenVec, const long int& snpAnzahl, const int& akt_zeile)
{
		string pheno_value	="";
		string sex_value	="";	
		long int spaltenZahl=zeilenVec.size();
		// check if  vector<CBSNP*> genInfo has  size equal to snpAnzahl.
		// this declaration will produce pointer to CBPED classes
		CBPED* pCBPED=new CBPED;
		// check if there are correct pheno infos and genotype infos. plink has either 6 phenotype infos.
		int value = spaltenZahl-6; // 6 is number columns of phenotype infos
		if(CMSNP::_is_code_mach)
			value=spaltenZahl-5;
		//cout << value <<endl;
			//spaltenZahl = 6 + nSNPs
			//cout <<"value: "<<value<<endl;
			//cout << boolalpha<< (value==snpAnzahl)<<endl;
			//cout << boolalpha<< (value/2==snpAnzahl)<<endl;
		bool tfValue =((value==snpAnzahl) || (value/2==snpAnzahl));

		if(!tfValue)	// Fall mit 6 Zeilen.
		{
			string msg="\t**->Problem on the " + change_int_into_string(akt_zeile)+ "th row.\n";
			printLIN(msg);
			msg="There must be 5 columns  for phenotype infos and  either exactly "+change_int_into_string(snpAnzahl)+\
			" columns  with \"A/B\" type of genotype Infos or exactly "+change_int_into_string(2*snpAnzahl)+" columns \"A   B\" type of  genotype infos.\n";
			delete pCBPED;
			 error(msg);

		}else
		{
			pCBPED->famId		=	zeilenVec[0];
			pCBPED->indId		=	zeilenVec[1];
			pCBPED->patId		=	zeilenVec[2];
			pCBPED->matId		=	zeilenVec[3];
			//the following code is to save sex. if sexcode=="1" then it   male is saved as:  CBPED::sex=true;
			// if sexcode=="2" then it is female and we save it as:  sex=false and miss_sex=false
			// if sexcode=="missing" then CBPED::sex=false and CBPED::sex_miss=true;
			 //cout << "sex: "<< zeilenVec[4]<<" "<<endl;
			 // for sex_info
			 sex_value= zeilenVec[4];
			 lese_sex_info(pCBPED, sex_value);


			//for phenotype
			 //for phenotype
			 pheno_value=zeilenVec[5];
			 lese_phenotype_info(pCBPED, pheno_value);
			
		}

 return pCBPED;}
//this is plink type, I have to delete it.
void CBPED::lese_pedfile_all(vector<CBPED*> &vPedInfo,	const vector<CBSNP*>&genInfo, vector<string>&zeilenVector, const string& pedFileName)
 {
	string pr_msg="\n*-> Reading  ped file \""+pedFileName+ "\":\n";
	printLIN(pr_msg);
	 const long int snpAnzahl=genInfo.size();
 	//cout << "snpAnzahl: "<<snpAnzahl <<endl;
 	// we need CBPED*
 	CBPED* pCBPED;
 	IFS  file(pedFileName.c_str(),"r");
 	if(!file)
 		error("file does not exits´or could not open");
 	//vector<string> zeilenVector;
 	string elem;
 	string elem1;
 	//
 	//int nMan=0;
 	//int nWeib=0;
 	//int nTrans=0;
 	while(true)
 	{
 		if(_ifs_eof(file))		break;
 		static int akt_zeile=0;	 // aktuelle Zeile
 		// lese Zeile
 		bool tfValue=CBPED::lese_pedfile_zeile(file, zeilenVector);

 		int spaltenZahl=zeilenVector.size();
 		if(spaltenZahl==0)
 			continue;
 		//cout << "\n zeilenVector.size(): "<<zeilenVector.size() <<endl; // for test only
 		if(tfValue) akt_zeile++;
 		//string msg="\t**->Total number of columns in the "+change_int_into_string(akt_zeile)+ "th row is : "+change_int_into_string(spaltenZahl)+".\n";
 		//cout <<msg;
 		//LIN <<msg;
 		// lese spalten
 		pCBPED= lese_pedfile_pedInfo(zeilenVector,snpAnzahl, akt_zeile);
 		// add the pCBPED into vector
 		vPedInfo.push_back(pCBPED);
 		// now erase the first 6 columns of pedInfo from ZeilenVector
 	 	zeilenVector.erase(zeilenVector.begin(),zeilenVector.begin()+5);
 		// read genotpye now
 			lese_pedfile_genInfo(genInfo,zeilenVector, akt_zeile);
 		//end des loops.
 		spaltenZahl=0;
 		zeilenVector.clear();
 		zeilenVector.resize(0);
 	}

   _ifs_close(file);
// 	cout << vPedInfo.size()<<endl; // for test only
// check if there correct number of genotpyes saved in each CBSNP.
 	for (unsigned int i=0; i<genInfo.size();++i)
 		{
 			CBSNP* pCBSNP=genInfo[i];

 			bool tfValue=false;
 			if(pCBSNP->geno1.size()!=vPedInfo.size()&& pCBSNP->geno2.size()!=vPedInfo.size())
 			{
 				tfValue=true;
 				cout << "geno1.size(): " << pCBSNP->geno1.size()<<endl;
 				cout << "geno2.size(): " << pCBSNP->geno2.size()<<endl;
 				cout << "vPedInfo.size(): " << vPedInfo.size()<<endl;
 				cout << "allele1: " << pCBSNP->allele1<<endl;
 				cout << "allele2: " << pCBSNP->allele2<<endl;
 			}
 			if(tfValue)
 			{
 				//delete	pCBSNP;
 				error("problem in saving genotypes in  "+ change_int_into_string(i) +"th SNP with Name "+pCBSNP->snpName+".\n");

 			//for(unsigned int j=0; j<pCBSNP->geno1.size();++j)
 			//	cout << pCBSNP->geno1[j]<< " "<< pCBSNP->geno2[j]<< " \n";
 				//cout << endl;
 			}
 		}

 	string pmsg="**->Total number of successfully read individuals: "+	change_int_into_string(vPedInfo.size())+ "\n";
 	//printLIN(pmsg);
 	cout <<"\t"<< setw(coutWidth);
 	printLIN(pmsg);
 }
//alternative
 void CBPED::lese_pedfile( vector<CBPED*>& pedInfo, const vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName)
 {
	 printLIN("*->Reading file \""+fileName+"\": \n");
	 //CBPED* pCBPED;
	 string myString					="";
	 bool ende							= false;
	 unsigned int ccol					=0;
	 unsigned int ncol					=0;
	 unsigned int nrow					=0;
	 unsigned int nLines				=0;
	 unsigned int max_ncol				=0;
	 bool def_max_ncol_once				=true;
	 unsigned long int	snpAnzahl		=0;
	string	pheno_value					="";
	string	sex_value					="";
	 try
	 {
		IFS pf(fileName.c_str(),"r");
		while(true)
			 {
				if(_ifs_eof(pf)) break;
				while(!_ifs_eof(pf))
				{
					ccol=CBSNP::stringLeser(pf, myString, ende);
					ncol+=ccol;
					 if(ccol!=0)
					{
						if(!myString.empty())
							cZeile.push_back(myString);
						else
						{
							cout << "error";
							exit(1);	
						}
					}
					if(ende|| _ifs_eof(pf))
					{
					 //ende=false;
						++nLines;
						if(def_max_ncol_once &&(ncol!=0))
						{
							max_ncol=ncol;
							def_max_ncol_once=false;
							//cout << "max_ncol: "<<max_ncol <<" " <<"ncol: "<<ncol <<endl; // only test; // only test
						}
						if((ncol==max_ncol)&& (max_ncol!=0))
							++nrow;
						myString="";
						//cout <<"-----------------\n";//test
						break;
					}
					if(ccol==0)
					 continue;	
					// this loop is broken only if end  is found.
				}//end of inner while loop 
				 if(cZeile.size()!=0)
			 {
				 snpAnzahl=genInfo.size();
				 CBPED* pCBPED= lese_pedfile_pedInfo(cZeile,snpAnzahl, nrow);
				 pedInfo.push_back(pCBPED);
				 //cout << "nrow: "<< nrow<< " pedInfo.size()"<< pedInfo.size()<<endl;
				 if(nrow!=pedInfo.size())
				 {
					 int sz=pedInfo.size()-1;
					 string msg= "Problem in saving genotypes of person "+pedInfo[sz]->indId+" in the "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
					 error(msg);
				 }
				 if(max_ncol!=ncol)
				 {
					 string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
					 error(msg);
				 }
			 // now erase the first 6 columns of pedInfo from ZeilenVector
				 cZeile.erase(cZeile.begin(),cZeile.begin()+5);
				 lese_pedfile_genInfo(genInfo, cZeile,nLines);
			 }
			 else continue;	
			 // assignment of genotypes

			//at the end des loops.
			 ncol=0;
			 cZeile.clear();
			 cZeile.resize(0);
			 if(_ifs_eof(pf))break;
			 } // end of outer while loop 
		_ifs_close(pf);
		nrow=0; // total rows
		nLines=0; // total lines
		//display_ped_summary();
		string pmsg="\t**->Total number of successfully read individuals: ";
		//printLIN(pmsg);
		printWithSpace(pmsg,  change_int_into_string(pedInfo.size()),coutWidth);
	 }catch(bad_alloc& memoryAllocationException)
	 {
		 cout << "error";

		 bool def_max_ncol_once				=true;
		
		IFS pf(fileName.c_str(),"r");

	 while(!_ifs_eof(pf))
	 {
		//read first column

		 if(_ifs_eof(pf)) break;
		 ccol=CBSNP::stringLeser(pf, myString, ende);
		 ncol+=ccol;
		// cout << "ccol: "<<ccol <<" "<< ncol<<" "<<myString<<  " "<< boolalpha<< ende<<", "<<"nLines: "<< nLines <<" nrow "<< nrow <<endl;
		 if(ende)
		 {
			 ende=false;
			 ++nLines;
			 if(def_max_ncol_once &&(ncol!=0))
			 {
				max_ncol=ncol;
				def_max_ncol_once=false;
			//	cout << "max_ncol: "<<max_ncol <<" ";
			//	cout << "ncol: "<<ncol <<endl;
			}
			 if((ncol==max_ncol)&& (max_ncol!=0))
				 ++nrow;

			 if(!myString.empty())
			 {
				 string msg= " There is only one column in the "+change_int_into_string(nLines)+" th row of "+fileName+" file. \n" ;
				 error(msg);
			  }

			 ncol=0;
			 continue;
	 		//cout <<"-----------------\n";//test

		 }
		 if(ccol==0)
			 continue;
		 //defining pointer of class CBPED
		 CBPED* pCBPED= new CBPED;
		 //
		 if(ncol==1&& !myString.empty() && pCBPED->famId=="")
			 pCBPED->famId=myString;
		 else
		  {
			  delete pCBPED;
			  string msg= "problem in saving family Id from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nLines+1)+ " line of   pedfile \""+fileName+"\". \n";
			  error(msg);
		  }
		 //read second column
		 ccol=CBSNP::stringLeser(pf, myString, ende);
		 ncol=ncol+ccol;
		// cout << "ccol: "<<ccol <<" "<< ncol<<" "<<myString<<  " "<< boolalpha<< ende<<", "<<"nLines: "<< nLines <<" nrow "<< nrow <<endl;
		 if((!ende && !_ifs_eof(pf))&& ncol==2&& !myString.empty() && pCBPED->indId=="")
			 pCBPED->indId=myString;
		 else
		  {
			  delete pCBPED;
			  string msg= "problem in saving individual Id from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nLines+1)+ " line of   pedfile \""+fileName+"\". \n";
			  error(msg);
		  }
		// third column
		ccol=CBSNP::stringLeser(pf, myString, ende);
		ncol=ncol+ccol;
		// cout << "ccol: "<<ccol <<" "<< ncol<<" "<<myString<<  " "<< boolalpha<< ende<<", "<<"nLines: "<< nLines <<" nrow "<< nrow <<endl;
		if((!_ifs_eof(pf)&&!ende)&& (ncol==3 && !myString.empty()) && pCBPED->patId=="")
			pCBPED->patId=myString;
		else
		{
			delete pCBPED;
			string msg= "problem in saving paternal Id from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nLines+1)+ " line of   pedfile \""+fileName+"\". \n";
			error(msg);
		}
		// fourth column
		ccol=CBSNP::stringLeser(pf, myString, ende);
		ncol=ncol+ccol;
		// cout << "ccol: "<<ccol <<" "<< ncol<<" "<<myString<<  " "<< boolalpha<< ende<<", "<<"nLines: "<< nLines <<" nrow "<< nrow <<endl;
		 if((!_ifs_eof(pf)&&!ende) && (ncol==4 && !myString.empty()) && pCBPED->matId=="")
			pCBPED->matId=myString;
			else
			{
				delete pCBPED;
				string msg= "problem in maternal Id from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nLines+1)+ " line of   pedfile \""+fileName+"\". \n";
				error(msg);
			}

		// fifth column
		ccol=CBSNP::stringLeser(pf, myString, ende);
		ncol=ncol+ccol;
		// cout << "ccol: "<<ccol <<" "<< ncol<<" "<<myString<<  " "<< boolalpha<< ende<<", "<<"nLines: "<< nLines <<" nrow "<< nrow <<endl;
		 if((!_ifs_eof(pf)&&!ende) && (ncol==5 && !myString.empty())&&(!pCBPED->sex&& pCBPED->miss_sex) )
		 {
				 if(myString=="M")
				 {
					 pCBPED->sex			=true;
					 pCBPED->miss_sex		=false;
					 ++pCBPED->nMale;
				 }
				 else if(myString=="F")
				 {
					 pCBPED->sex			=false;
					 pCBPED->miss_sex		=false;
					 ++pCBPED->nFemale;
				 }
				 else
				  ++pCBPED->nNosex;
		 }
		 else
		 {
			 delete pCBPED;
			 string msg= "problem in saving sex status Id from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nLines+1)+ " line of   pedfile \""+fileName+"\". \n";
			 error(msg);
		 }
		lese_phenotype_info(pCBPED, pheno_value);
		pedInfo.push_back(pCBPED);
		//sart of genotype reading
		unsigned long int ncSNP	=0; //no of current snp
		unsigned long int nSNPs	=genInfo.size();
		string first 			=""; // geno1
		string second			=""	;// geno2
	// genotype reading
		while(1)
		{
			//first term
			//cout << "\nGENOTYPES:\n";
			ccol=0;
			if(_ifs_eof(pf)) break;
			if(ende)
			{
				ende=false;
				if(!myString.empty())
				{
					string msg= " There is only one genotypes in the "+change_int_into_string(nLines+1)+" th row of "+fileName+" file. \n" ;
					error(msg);
				}
				++nLines;
				if(def_max_ncol_once &&(ncol!=0))
				{
					max_ncol=ncol;
					def_max_ncol_once=false;
					//cout << "max_ncol: "<<max_ncol <<" ";
					//cout << "ncol: "<<ncol <<endl;
				}
				if((ncol==max_ncol)&& (max_ncol!=0))
					++nrow;

				if(ncSNP!=nSNPs|| ncol!=(5+2*ncSNP)|| max_ncol!=ncol)
				{
					string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
					error(msg);
				}
				if(nrow!=pedInfo.size())
				{
					int sz=pedInfo.size()-1;
					string msg= "Problem in saving genotypes of person "+pedInfo[sz]->indId+" in the "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
					error(msg);
				}

				ncol=0;
				//cout <<"-----------------\n";//test
				break;
			}
		if(true)// single genotypes ie.e  one genotype in one column
		{
			if(ncSNP==nSNPs) break;
			++ncSNP;
			ccol=CBSNP::stringLeser(pf, myString, ende);
			ncol+=ccol;
			//cout << "ccol: "<<ccol <<" "<< ncol<<" "<<myString<<  " "<< boolalpha<< ende<<", "<< "nLines: "<< nLines <<" nrow "<< nrow <<endl;
			if((!ende && !_ifs_eof(pf))&& ncol==(2*ncSNP+4) && !myString.empty())
				first=myString;
			else
			{
				string msg= "Problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nLines+1)+"th row of   file \""+fileName+"\".\n";
				error(msg);

			}
			// second genotype
			ccol=CBSNP::stringLeser(pf, myString, ende);
			ncol+=ccol;
			//cout << "ccol: "<<ccol <<" "<< ncol<<" "<<myString<<  " "<< boolalpha<< ende<<", "<<"nLines: "<< nLines+1 <<" nrow "<< nrow <<endl;
			if(ncol==(2*ncSNP+5) && !myString.empty())
			{
				 second= myString;
				 myString="";
			}
			else
			{
				string msg= "Problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nLines+1)+"th row of   file \""+fileName+"\".\n";
				error(msg);
			}

			genInfoLeser(genInfo, first , second, ncol ,nrow,ncSNP);

			if(ende|| _ifs_eof(pf))
			{
				ende=false;
				++nLines;

				if(def_max_ncol_once &&(ncol!=0))
				{
					max_ncol=ncol;
						def_max_ncol_once=false;
						cout << "max_ncol: "<<max_ncol <<" "; // only test
						cout << "ncol: "<<ncol <<endl; // only test
				}
				if((ncol==max_ncol)&& (max_ncol!=0))
					++nrow;
				myString="";
				//
			//	cout << "nSNPs: "<<nSNPs<<endl;
			//	cout << "nrow: "<<nrow<<endl;
			//	cout << "size of pedInfo: "<<pedInfo.size()<<endl;

				if(ncSNP!=nSNPs|| ncol!=(5+2*ncSNP)|| max_ncol!=ncol)
				{
					string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
					error(msg);
				}
				if(nrow!=pedInfo.size())
				{
					int sz=pedInfo.size()-1;
					string msg= "Problem in saving genotypes of person "+pedInfo[sz]->indId+" in the "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
					error(msg);
				}

				ncol=0;
				//cout <<"-----------------\n";//test
				break;
			}


		}//end of first type of genotypes reading
		else if (false)
		{
			// double gentoypes

		}
		else
		{
			string msg= "could not read genotypes!";
			error(msg);

		}



		if(_ifs_eof(pf)) break;
		}// end of while loop first

		 if(_ifs_eof(pf)) break;

	 }//end of while big loop
	 
	 _ifs_close(pf);
	 }
	 
return; }
 //-------------------------------------------------------------//
 //helper functions 
 //-------------------------------------------------------------//
 void CBPED::lese_sex_info( CBPED* const pCBPED, const string sexValue,const bool update_sex)
 {
	 if(update_sex)
	{
		bool is_male 			= (pCBPED->sex && !pCBPED->miss_sex); // male
		bool is_female 			= (!pCBPED->sex && !pCBPED->miss_sex); // female
		//bool _miss_sex 			= (!pCBPED->sex && pCBPED->miss_pheno); // miss sex 
		if(is_male)
			--CBPED::nMale;
		else if(is_female)
			--CBPED::nFemale;
		else 
			--CBPED::nNosex;
			
	} 
	 if(sexValue!="")
	 {
		 if(sexValue=="M" ||sexValue=="1"||sexValue=="m" )
				{//male
					pCBPED->sex			=true;
					pCBPED->miss_sex		=false;
					++CBPED::nMale;

				}
				else if (sexValue=="F" ||sexValue== "2"||sexValue=="f")
				{ //female
					pCBPED->sex=false;
					pCBPED->miss_sex=false;
					++CBPED::nFemale;

				}
				else
				{	// missing
					pCBPED->sex=false;
					pCBPED->miss_sex=true;
					++CBPED::nNosex;
				}
	 }
	 else
	 {
		 if(pCBPED->sex)
			 ++CBPED::nMale;
		 else if(!pCBPED->sex && !pCBPED->miss_sex)
			 ++CBPED::nFemale;
		 else if(pCBPED->miss_sex)
			 ++CBPED::nNosex;
	 }

 return ;}
bool CBPED::lese_genotype_info(CBSNP*pCBSNP, const string first,  const string second,const  int akt_zeile, string&ERR_MSG)
{
	bool fatal =false;
	if((first!=CBSNP::_miss_gvalue&& first!="?")  && (second!=CBSNP::_miss_gvalue&& second!="?"))
	{

		if(first==second)
		{
			// add allele
			if(pCBSNP->allele1==""&& pCBSNP->allele2!=first)pCBSNP->allele1=first;
			else if (pCBSNP->allele2==""&& pCBSNP->allele1!=first) pCBSNP->allele2=first;
			// add genotype
			if(first==pCBSNP->allele1&& second==pCBSNP->allele1 )
			{// this is the case where homozygote1 is denoted by 00
				//cout << "reading homo1 \n";
				//cout << first << " "<<second<<"\n";
				pCBSNP->geno1.push_back(false);
				pCBSNP->geno2.push_back(false);
				pCBSNP->aOrder.push_back(true);
			}
			else if(first==pCBSNP->allele2 && second==pCBSNP->allele2)
			{// this is the case where homozygote 2  is denoted by 11
				//cout << "reading homo2\n";// for test only
				//cout << first << " "<<second<<"\n";// for test only
				pCBSNP->geno1.push_back(true);
				pCBSNP->geno2.push_back(true);
				pCBSNP->aOrder.push_back(true);
			}
			else
			{
				string wrongAllele="";
				if(first!=pCBSNP->allele1 && first!=pCBSNP->allele2 )
					wrongAllele=first;
				else if(second!=pCBSNP->allele1 && second!=pCBSNP->allele2 )
					wrongAllele=second;
				ERR_MSG="SNP "+pCBSNP->snpName+" containing two alleles "+\
						pCBSNP->allele1+ " and "+pCBSNP->allele2+ " has found a new allele type"+wrongAllele+ "  at individual "+\
						change_int_into_string(akt_zeile)+ ".\n ";
				fatal=true;
			}
		}
		else if(first!=second) // save it as 01 .i.e. false true.
					{// this is the case where heterozygote   is denoted by 01
						//cout << "reading heterozygote\n";// for test only
						//cout << first << " "<<second<<"\n";// for test only
						if(pCBSNP->allele1==""&&pCBSNP->allele2!=first )pCBSNP->allele1=first;
						else if (pCBSNP->allele2==""&&pCBSNP->allele1!=first ) pCBSNP->allele2=first;
						if (pCBSNP->allele1==""&& pCBSNP->allele2!=second ) pCBSNP->allele1=second;
						if (pCBSNP->allele2==""&& pCBSNP->allele1!=second ) pCBSNP->allele2=second;
						//now check if allele are correct in first and second
						if(first==pCBSNP->allele1&& second==pCBSNP->allele2 )
						{
						 	pCBSNP->geno1.push_back(false);
							pCBSNP->geno2.push_back(true);
							pCBSNP->aOrder.push_back(true);
						}
						else if(first==pCBSNP->allele2&& second==pCBSNP->allele1 )
						 	 {
								pCBSNP->geno1.push_back(false);
								pCBSNP->geno2.push_back(true);
								pCBSNP->aOrder.push_back(false);
						 	 }
						else
						{
							string wrongAllele="";
							if(first!=pCBSNP->allele1&& first!=pCBSNP->allele2 )
							wrongAllele=first;
							else if(second!=pCBSNP->allele1 && second!=pCBSNP->allele2 )
							wrongAllele=second;
							ERR_MSG="SNP "+pCBSNP->snpName+" containing two alleles "+\
									pCBSNP->allele1+ " and "+pCBSNP->allele2+ " has found a new allele type"+wrongAllele+ "  at individual "+\
									change_int_into_string(akt_zeile)+ ".\n ";
							 fatal=true;

						}
					}
		}// end of if (!CBSNP::_miss_gvalue)
			// This is the case where both genotypes are missing
		else if ((first==CBSNP::_miss_gvalue || first=="?") || (second==CBSNP::_miss_gvalue || second=="?"))
		{ //In this case we use the code 10 to denote the missing genotypes
			//cout << "reading missing \n"; // for test only
			//cout << first << " "<<second<<"\n"; // for test only
			pCBSNP->geno1.push_back(true);
			pCBSNP->geno2.push_back(false);
			pCBSNP->aOrder.push_back(true);
		}
			// this is the case where non of above condition is true then error.



return fatal;}
//xxxxxxxxxxxxxxxxxxx
void CBPED::add_sex_n_phenotype_info(CBPED* const pCBPED, const bool update_sex,const bool update_pheno)
{
	//cout << nMale << ", nFemale:  "<< nFemale << " \n";
	if(update_sex)
		{
			bool is_male 			= (pCBPED->sex && !pCBPED->miss_sex); // male
			bool is_female 			= (!pCBPED->sex && !pCBPED->miss_sex); // female
			//bool _miss_sex 			= (!pCBPED->sex && pCBPED->miss_pheno); // miss sex
		//	cout << is_male << is_female <<endl;
			if(is_male)
				++CBPED::nMale;
			else if(is_female)
				++CBPED::nFemale;
			else
				++CBPED::nNosex;

		}
	 if(update_pheno)
		{
			bool aff 			= (pCBPED->pheno && !pCBPED->miss_pheno); // affected
			bool unaff 			= (!pCBPED->pheno && !pCBPED->miss_pheno); // unaffected
			//bool _miss_status 	= (pCBPED->miss_pheno); // miss status
			if(aff)
				++CBPED::nAffected;
			else if(unaff)
				++CBPED::nUnaffected;
			else
				++CBPED::nMiss_status;

		}

}
//xxxxxxxxxxxxxxxxxxx

void CBPED::display_ped_summary(void)
{
	 // count affected and unnaffected
		cout <<left;
		cout <<"*->Individual summary:\n";
		// male and female Geschiste
		cout <<"\t"<< setw(coutWidth)<< 	"**->Total female: " << CBPED::nFemale<<"\n";
		cout <<"\t"<<setw(coutWidth)<< 	"**->Total male: " << CBPED::nMale<< "\n";
		cout <<"\t"<<setw(coutWidth)<< 	"**->Total undefined individual: " <<CBPED::nNosex<< "\n";
		cout.flush()<<endl;
		//write to fcgene.out.txt file
		LIN << left;
		LIN <<						 	"\n*->Individual summary:\n";
		LIN <<"\t"<< setw(coutWidth)<< 		"**->Total female: " <<  CBPED::nFemale<<"\n";
		LIN <<"\t"<<setw(coutWidth)<< 			"**->Total male: " <<  CBPED::nMale<< "\n";
		LIN <<"\t"<<setw(coutWidth)<< 			"**->Total undefined individual: " << CBPED::nNosex<< "\n";
		LIN.flush()<<endl;


		// for affected and unaffected.
		cout <<"\t"<< setw(coutWidth)<< 	"**->Total cases: " <<  CBPED::nAffected<<"\n";
		cout <<"\t"<<setw(coutWidth)<< 	"**->Total controls: " <<  CBPED::nUnaffected<< "\n";
		cout <<"\t"<<setw(coutWidth)<< 	"**->Total undefined individuals: " << CBPED::nMiss_status<< "\n";
		cout.flush()<<endl;
		//write to fcgene.out.txt file
		LIN <<"\t"<< setw(coutWidth)<< 	"**->Total cases: " <<  CBPED::nAffected<<"\n";
		LIN<<"\t"<<setw(coutWidth)<< 	"**->Total controls: " <<CBPED::nUnaffected<< "\n";
		LIN<<"\t"<<setw(coutWidth)<< 	"**->Total undefined individuals: " << CBPED::nMiss_status<< "\n";
		LIN.flush()<<endl;

return ;}
 void CBPED::lese_phenotype_info(CBPED*const pCBPED, const string pheno, const bool update_pheno )
 {
 	 if(update_pheno)
	{
		bool aff 			= (pCBPED->pheno && !pCBPED->miss_pheno); // affected 
		bool unaff 			= (!pCBPED->pheno && !pCBPED->miss_pheno); // unaffected 
		//bool _miss_status 	= (pCBPED->miss_pheno); // miss status
		if(aff)
			--CBPED::nAffected;
		else if(unaff)
			--CBPED::nUnaffected;
		else
			--CBPED::nMiss_status;
			
	}	

	
	if(pheno!="")
 	{
 		if(pheno=="2" || pheno=="aff")
 		{ // affected
 			pCBPED->pheno		=true;
 			pCBPED->miss_pheno	=false;
 			++CBPED::nAffected;
 		}
 		else if(pheno=="1"|| pheno=="unaff")
 		{ // unaffected
 			pCBPED->miss_pheno	=false;
 			pCBPED->pheno		=false;
 			++CBPED::nUnaffected;
 		}
 		else
 		{ // missing
 			pCBPED->miss_pheno	=true;
 			pCBPED->pheno		=false;
 			++CBPED::nMiss_status;
 		}
 	}
 	else
 	{
 		if(pCBPED->miss_pheno)
 			++CBPED::nMiss_status;
 		else if(pCBPED->pheno && !pCBPED->miss_pheno )
 			++CBPED::nAffected;
 		else if(!pCBPED->pheno && !pCBPED->miss_pheno )
 			++CBPED::nAffected;
 	}
return ;}
//-------------------------------------------------------------//
// computation related functions 
//-------------------------------------------------------------//
 void CBPED::calculate_callrate(const vector<CBPED*>& pedInfo,const vector<CBSNP*>& genInfo)
 {
	 // missings are denoted by TRUE FALSE
	 int 	p_nMiss				=0;
	 const unsigned int pTotal	=pedInfo.size();
	 const unsigned int sTotal	=genInfo.size();
	 const unsigned int pgTotal= CBPED::no_of_good_indivs(pedInfo);
	 const unsigned int sgTotal= CBSNP::no_of_good_snps(genInfo);
	 CBPED* pCBPED 				=pedInfo[0];
	 CBSNP* pCBSNP				=genInfo[0];
	 vector<int>s_nMiss(sTotal,0);
	// cout << "\n indivs:"<<pTotal << " "<< "nsps: "<< sTotal <<endl;
	 for(unsigned int i =0;i<pTotal;++i)
	 {
		 pCBPED =pedInfo[i];
		 if(pCBPED->quality)
		 {
				 for(unsigned int j=0;j<sTotal;++j)
			 {
				 pCBSNP				=genInfo[j];
				 if(pCBSNP->quality)
				 {
					if((pCBSNP->geno1[i])&& !(pCBSNP->geno2[i]))
					{
						 ++p_nMiss;
						s_nMiss[j]= (s_nMiss[j])+1;
					}
				 }
				 //++pCBSNP;

			 }
		 }
		 //calculate individual call rate
		 pCBSNP=genInfo[0];
		// cout << "p_nMiss: "<<p_nMiss ; //<< endl;
		//cout <<" "<< (1.0*(sTotal-p_nMiss)/sTotal)<< "  "<<endl; //debug
		pCBPED->ind_CR=(1.0-(float)(p_nMiss)/sgTotal);
		//++pCBPED;
		p_nMiss=0;
	 }

	 pCBPED	=pedInfo[0];
	 //calculate snpwise call rate
	 for(unsigned int j=0;j<sTotal;++j)
	 {
		 pCBSNP=genInfo[j];
		 pCBSNP->snp_CR=(1.0-s_nMiss[j]/pgTotal);
		// ++pCBSNP;

	 }
		 pCBSNP=genInfo[0];

 return; }
//------------------------------------//

 void CBPED::genoProb_to_geno_converter(const vector<CBSNP*>& genInfo,const float& thresh)
 {
 	for(unsigned int j=0;j<genInfo.size(); ++j)
 		  		{
 		  			CBSNP* pCSNP=genInfo[j];

 		  				unsigned int size=pCSNP->pgeno1.size();
 		  				for(unsigned int i=0; i<size;++i)
 		  				{
 		  					float pg1	=pCSNP->pgeno1[i];
 		  					float pg2	=pCSNP->pgeno2[i];
 		  					float pg3	=pCSNP->pgeno3[i];
 		  					float maxv	=max(pg1,pg2);
 		  					maxv 		= max(maxv,pg3);
 		  					if(maxv==0)
 		  					{
 		  						// missing case
 		  						pCSNP->geno1.push_back(true);
 		  						pCSNP->geno2.push_back(false);
 		  						pCSNP->aOrder.push_back(true);

 		  					}
 		  					else if(maxv<=thresh)
 		  					{

 		  						pCSNP->geno1.push_back(true);
 		  						pCSNP->geno2.push_back(false);
 		  						pCSNP->aOrder.push_back(true);
 		  					}
 		  					else if(maxv==pg1)
 		  					{
 		  						// case homo1
 		  						pCSNP->geno1.push_back(false);
 		  						pCSNP->geno2.push_back(false);
 		  						pCSNP->aOrder.push_back(true);
 		  					}
 		  					else if(maxv==pg2)
 		  					{
 		  						// case hetero
 		  						pCSNP->geno1.push_back(false);
 		  						pCSNP->geno2.push_back(true);
 		  						pCSNP->aOrder.push_back(true);

 		  					}
 		  					else if (maxv==pg3)
 		  					{
 		  						// case homo2
 		  						pCSNP->geno1.push_back(true);
 		  						pCSNP->geno2.push_back(true);
 		  						pCSNP->aOrder.push_back(true);
 		  					}
 		  					else
 		  					{
 		  						string msg= "problem in converting impute genotype format to ped genotype format\n";
 		  						error(msg);

 		  					}


 		  			   }
 		  		}
 return ; }

 //--------------------------------------//
void CBPED::update_pedInfo_given_gProb(vector<CBPED*>&pedInfo,const vector<CBSNP*>genInfo, const bool update_sex, const bool update_pheno)
	{
		const unsigned int sz=CBSNP::no_of_persons;
		
		for(unsigned int i=0;i<sz;++i)
		{
			if(pedInfo.size()==i)
			{
				CBPED* pCBPED =new CBPED;
				pedInfo.push_back(pCBPED);
			}
			CBPED* pCBPED=pedInfo[i];
			 // count affected and unnaffected
			if(pCBPED->sex) ++nMale;
			if(!pCBPED->sex&& !pCBPED->miss_sex) ++nFemale;
			if(!pCBPED->sex&& pCBPED->miss_sex) ++nNosex;
			if(pCBPED->miss_pheno)
				++nMiss_status;
			if(pCBPED->pheno && !pCBPED->miss_pheno )
				++nAffected;
			if(!pCBPED->pheno && !pCBPED->miss_pheno)
			++nAffected;

		}

return;}

//--------------------------------------------------------------------//
 //function for CMSNP
 //---------------------------------------------------------------------//

void CMPED::pedFileleser(const string& pedFileName,	 vector<CBPED*> &vPedInfo, const vector<CBSNP*>& genInfo)
{

	printLIN("*->Reading file named "+ pedFileName+ ":\n");
	IFS file(pedFileName.c_str(),"r");
	if(!file)
	{
		string msg= pedFileName+ " either does not exits or could not be opend.\n"	;
		error(msg);
	}
	// if file is okay then define more varaibles
	string elem(""); // gelesen von ped File;
	string elem1("");
	 int ncol				 =0;
	 int ccol				=0;
	static int nrow			=0;
	bool zeilenEnde			=false;

	while(!_ifs_eof(file))
	{
		if(_ifs_eof(file))
		{
			zeilenEnde=true;
			if(nrow==0)
			{
				string msg=" There are no elements in "+pedFileName +" file.\n";
				error(msg);
			}
			else
			{
				cout << "total no of Individuals read: "<< nrow << endl; // only for test
				nrow=0;
				cout << "size of vector vPedInfo: "<<vPedInfo.size()<<endl; // only for test
			}
			break;
		}
		zeilenEnde=false;
		// define  pointer of CBPED

		while(!zeilenEnde && !_ifs_eof(file))
		{ // do all works till the end of line and then break
			if(zeilenEnde) break;
			if(_ifs_eof(file)) break;
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

			if( zeilenEnde || _ifs_eof(file) )
			{
				if(!elem1.empty())
				{
					string msg= " There is only one column in the first row of "+pedFileName+" file. \n" ;
					error(msg);
				}
				else
				break; // only in the first column
			}
			//define pointer
			CBPED* pCBPED=new CBPED;
			// read family id
			if(!elem1.empty())
			{
				if(ccol!=0 && ncol==1 ) // this condition must always be true ;
				{
					if(pCBPED->famId=="")
					pCBPED->famId		=elem1;
					else
					{
						string msg= "problem in saving family Id from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
						error(msg);
					}
				}
				else
				{
						string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
						error(msg);
				}
			}
			else
			{
				string msg= "String is empty! problem in reading  family Id  from pedfile.\n";
				error(msg);
			}
			// for the second column
			if(zeilenEnde || _ifs_eof(file))
			{
				string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
				 error(msg);
			}
			// start of rading second column
			elem1="";
			elem="";
			ccol=0;
			ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
			elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
			if(elem1.empty())
			{
				string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
				error(msg);
			}
			else
			{
				if(ccol!=0 && ncol==2 ) // this condition must always be true ;
				{
					if(pCBPED->indId=="")
						pCBPED->indId		=elem1;
					else
					{
						string msg= "problem in saving Individual Id from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
						error(msg);
					}
				}
				else
				{
						string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
						error(msg);
				}

			//   This case is true except last column
				if(zeilenEnde || _ifs_eof(file))
				{		string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow +1)+"th row"+ "of   file \""+pedFileName+"\".\n";
						error(msg);
				}
			}
			// start of rading third column
			elem1="";
			elem="";
			ccol=0;
			ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
			elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
			if(elem1.empty())
			{
				string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
				error(msg);
			}
			else
			{
				if(ccol!=0 && ncol==3 ) // this condition must always be true ;
				{
					if(pCBPED->patId=="")
					pCBPED->patId		=elem1;
					else
					{
					string msg= "problem in saving family Id from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
					error(msg);
					}
				}
				else
				{
					string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
					error(msg);
				}

					//   This case is true except last column
				if(zeilenEnde || _ifs_eof(file))
				{		string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
						error(msg);
				}


			}
			// start of rading fourth column
			elem1="";
			elem="";
			ccol=0;
			ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
			elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
			if(elem1.empty())
				{
					string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
					error(msg);
				}
			else
			{
				if(ccol!=0 && ncol==4 ) // this condition must always be true ;
				{
					if(pCBPED->matId=="")
						pCBPED->matId		=elem1;
						//if(true)
						//	cout << "element value: "<< elem1 <<endl; // test
						else
						{
							string msg= "problem in saving Mother Id from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
							error(msg);
						}
					}
				else
				{
					string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
					error(msg);
				}

					//   This case is true except last column
				if(zeilenEnde || _ifs_eof(file))
				{
					string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
					error(msg);
				}
			}
				// start of rading fifth column
			elem1="";
			elem="";
			ccol=0;
			ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
			elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
			if(elem1.empty())
			{
				string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
				error(msg);
			}
			else
			{
				if(ccol!=0 && ncol==5 ) // this condition must always be true ;
				{
					if(!pCBPED->sex&& pCBPED->miss_sex)
					{
						lese_sex_info(pCBPED, elem);
						// It is not necessary to assign if elem=f because sex is then already assigned false
					}
					else
					{
						string msg= "problem in saving sex status from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
						error(msg);
					}
				}
				else
				{
					string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
					error(msg);
				}

			//   This case is true except last column
				if(zeilenEnde || _ifs_eof(file))
				{		string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
					error(msg);
				}
			}

			//since phenotype status is not given in mach, put all of them in
			//undefined categorie
			if(pCBPED->miss_pheno)
				++nMiss_status;
			if(pCBPED->pheno && !pCBPED->miss_pheno )
				++nAffected;
			if(!pCBPED->pheno && !pCBPED->miss_pheno )
				++nAffected;


				vPedInfo.push_back(pCBPED);
				//cout <<"manlich: "<< nManlich << "  "<<nWeiblich << " "<< nTranse<<endl;
				// end of phenotype reading.
				//-------------------------------------------------//
				//sart of genotype reading
				unsigned long int snpN	=0;
				unsigned long int nSNPs			=genInfo.size();
				string first 			=""; // geno1
				string second			=""	;// geno2

			while(1) // start of genotype reading loop
			{
					elem1="";
					elem="";
					ccol=0;
					if(snpN==nSNPs) break;
					++snpN;
					if(true) // true  single genotypes
					{
					//-----------------------------------------------//
					// start of reading genotypes
					//reading first genotype
						ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
						elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
						if(elem1.empty())
						{
							string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row of   file \""+pedFileName+"\".\n";
							error(msg);
						}
						else
						{
							if(ccol!=0 && ncol==(int)(2*snpN+4) ) // this condition must always be true ;
							{
								first=elem1;
							}
							else
							{
								string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
								error(msg);
							}

							//   This case is true except last column
							if(zeilenEnde || _ifs_eof(file))
							{
								string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
								error(msg);
							}
						}
						// reading second genotype
						ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
						elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
						if(elem1.empty())
						{
							string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
							error(msg);
						}
						else
						{
							if(ccol!=0 && ncol==(int)(2*snpN+5) ) // this condition must always be true ;
							{
								second=elem1;
							}
							else
							{
								string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   pedfile \""+pedFileName+"\". \n";
								error(msg);
							}

						//   This case is true except last column
							//cout << snpN << " " <<nSNPs << endl; // test
							if(snpN!=nSNPs)
							{
								//cout << snpN << " " <<nSNPs << endl; // test
								if(zeilenEnde || _ifs_eof(file))
								{
									string msg= "sNot enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row of   file \""+pedFileName+"\".\n";
									error(msg);
								}
							}

						}

						genInfoLeser(genInfo, first , second, ncol ,nrow,snpN);
						//cout << "snpN : "<<snpN <<endl; // only for test
						//CBSNP * pCBSNP=genInfo[snpN-1]; // only for test
						//cout << pCBSNP->snpName <<endl; // only for test

					//-----------------------------------------------//
					}
					else if(false) // compund genotypes
					{
						cout << "false\n ";
						exit(1);
					}
					else
					{// neither single genotypes nor compound genotypes
						string msg= "Problemm in  reading genotypes \n";
						error(msg);
					}
					if(snpN==nSNPs) break;
					if(zeilenEnde) break;
			} // end of genotyping while loop
			snpN		=0;
			//--------------------------------------------------//
			// end of genotype reading



		//---------------------------end of first loop -----------#
		if(zeilenEnde|| _ifs_eof(file))
			break;
		}
		//---------------------------end of first loop -----------#
		ncol		=0;
		ccol		=0;

		if(_ifs_eof(file))
		{
			zeilenEnde=true;
			if(nrow==0)
			{
				string msg=" There are no elements in "+pedFileName +" file.\n";
				error(msg);
			}
			else
			{
				int indSize=vPedInfo.size();
				string msg= "\t**->total no of individuals read successfully: "+ change_int_into_string(indSize)+" \n";
				printLIN(msg);
				nrow=0;

			}
			nrow=0;
			break;
		}


	}

	_ifs_close(file);

	// update sex and  phenotype

			/*for(unsigned int i=0;i<vPedInfo.size();i++)
			{
				CBPED* pCBPED=vPedInfo[i];
				cout << pCBPED->famId << " ";
				cout <<pCBPED->indId <<" "<< pCBPED->patId <<" ";
				cout << pCBPED->matId<< " ";
				if(pCBPED->sex) cout << "M" <<" ";
				else
				if(!pCBPED->sex&& !pCBPED->miss_sex) cout << "F" <<" ";
				if(!pCBPED->sex&& pCBPED->miss_sex) cout << "0"<<" ";
				for(unsigned int j=0;j<genInfo.size();j++)
				{
					CBSNP* pCBSNP=genInfo[j];
					cout << boolalpha<< pCBSNP->aOrder[i]<< " ";
					if(!pCBSNP->geno1[i] && !pCBSNP->geno2[i])
					cout << pCBSNP->allele1 << " " <<pCBSNP->allele1 << " ";
					// homo2
					if(pCBSNP->geno1[i] && pCBSNP->geno2[i])
					cout << pCBSNP->allele2 << " " <<pCBSNP->allele2 << " ";
					// hetero
					if(!pCBSNP->geno1[i] && pCBSNP->geno2[i])
					{
						//cout << boolalpha<< pCBSNP->aOrder[i]<< " ";
						if(pCBSNP->aOrder[i])
							cout << pCBSNP->allele1 << " " <<pCBSNP->allele2 << " ";
						if(!pCBSNP->aOrder[i])
						cout << pCBSNP->allele2 << " " <<pCBSNP->allele1 << " ";
					}

					if(pCBSNP->geno1[i] && !pCBSNP->geno2[i])
						cout << CBSNP::_miss_gvalue << " " <<CBSNP::_miss_gvalue << " ";
				}
				cout <<endl;
			}

			*/

return ;}
//mlgeno file (mach output file leser)
 // display function
 // mach mlgeno file leser
 //alternative
  void CMPED::lese_mach_mlProbFile(vector<CBPED*>& pedInfo,  const vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName, const float & thresh)
  {
 	 printLIN("*->Reading file \""+fileName+"\": \n");
 	 string myString					="";
 	 string temp_value					="";
 	 string emsg						="";
 	 bool ende							= false;
 	 unsigned int ccol					=0;
 	 unsigned int ncol					=0;
 	 unsigned int nrow					=0;
 	 unsigned int nLines				=0;
 	 unsigned int max_ncol				=0;
 	 bool def_max_ncol_once				=true;
 	 unsigned long int	snpAnzahl		=0;
 	 float first_prob					=0.0;
 	 float second_prob					=0.0;
 	 float third_prob					=0.0;
 	string	pheno_value					="";
 	string	sex_value					="";
 	 try
 	 {
 		IFS pf (fileName.c_str(),"r");
		if(!pf)
		{
			string msg= fileName+ " either does not exits or could not be opend.\n"	;
			error(msg);
		}

 		while(true)
 			 {
 				if(_ifs_eof(pf))break;
 				while(!_ifs_eof(pf))
 				{
 					ccol=CBSNP::stringLeser(pf, myString, ende);
 					ncol+=ccol;
 					 if(ccol!=0)
 					{
 						if(!myString.empty())
 							cZeile.push_back(myString);
 						else
 						{
 							cout << "error";
 							exit(1);
 						}
 					}
 					if(ende|| _ifs_eof(pf))
 					{
 					 //ende=false;
 						++nLines;
 						if(def_max_ncol_once &&(ncol!=0))
 						{
 							max_ncol=ncol;
 							def_max_ncol_once=false;
 							//cout << "max_ncol: "<<max_ncol <<" " <<"ncol: "<<ncol <<endl; // only test; // only test
 						}
 						if((ncol==max_ncol)&& (max_ncol!=0))
 							++nrow;
 						myString="";
 						//cout <<"-----------------\n";//test
 						break;
 					}
 					if(ccol==0)
 					 continue;
 					// this loop is broken only if end  is found.
 				}//end of inner while loop
 				 if(cZeile.size()!=0)
 			 {
 				 snpAnzahl=genInfo.size();
 				 //cout << snpAnzahl << endl; //debug
 				// CBPED* pCBPED= lese_pedfile_pedInfo(cZeile,snpAnzahl, nrow);
 				//                lese_pedfile_pedInfo(zeilenVec,  snpAnzahl, akt_zeile)

 				 //-----------------------------------------------------------------//
 				string pheno_value	="";
 				string sex_value	="";
 				long int spaltenZahl=cZeile.size();
 				// check if  vector<CBSNP*> genInfo has  size equal to snpAnzahl.
 				// this declaration will produce pointer to CBPED classes

 				unsigned long int value = spaltenZahl-2; //
 					//cout <<"value: "<<value<<endl;
 					//cout << boolalpha<< (value==snpAnzahl)<<endl;
 					//cout << boolalpha<< (value/2==snpAnzahl)<<endl;
 				bool tfValue =(value==2*snpAnzahl) ;
 				if(tfValue)
 				{
 					CBPED* pCBPED=new CBPED;
 					myString 			=cZeile[0];
 					int temp_indx		=myString.rfind("->");
 					pCBPED->famId		=myString.substr(0,temp_indx);
 					temp_value			=myString.substr(temp_indx+2,myString.size());
 					pCBPED->indId		=temp_value;
 					pCBPED->patId		="";
 					pCBPED->matId		="";
 					//for sex
 					lese_sex_info(pCBPED, sex_value);
 					//for phenotype
 					lese_phenotype_info(pCBPED, pheno_value);
 					pedInfo.push_back(pCBPED);
 					//cout << "test" <<endl; //debug
 					myString			="";

 				}else
 				{

 					emsg="There must be 2 columns  for phenotype infos and   exactly "+change_int_into_string(2*snpAnzahl)+\
 							" columns  with  genotype probabilities for AA and AB.\n";
 					//delete pCBPED;
 					error(emsg);
 				}
 				 //-------------------------------------------------------------//
 				 //cout << "nrow: "<< nrow<< " pedInfo.size()"<< pedInfo.size()<<endl;
 				 if(nrow!=pedInfo.size())
 				 {
 					 int sz=pedInfo.size()-1;
 					 emsg= "Problem in saving genotypes of person "+pedInfo[sz]->indId+" in the "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
 					 error(emsg);
 				 }


 				 if(max_ncol!=ncol)
 				 {
 					 emsg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
 					 error(emsg);
 				 }
 			 // now erase the first 2 columns of pedInfo from ZeilenVector
 				 cZeile.erase(cZeile.begin(),cZeile.begin()+2);
 				// cout << cZeile.size()<<endl;
 				//saving genotype probabilities.
 					 //CBSNP* pCBSNP = new CBSNP;
 			 				//genInfo.push_back(pCBSNP);
 				 for(unsigned int i=0;i<snpAnzahl; ++i)
 				 {
 					//cout << "test" <<endl; //debug
 					 CBSNP* pCBSNP=genInfo[i];
 					 if(i==1) pCBSNP->no_of_persons=pedInfo.size(); // save total no of persons also in snp class.
 					 first_prob  =atof( cZeile[2*i].c_str()); // first is homo1
 					 second_prob = atof(cZeile[2*i+1].c_str()); // hetero
 					 third_prob = 1.0 -(first_prob+ second_prob); // homo2
 					if(third_prob<0.0)
 					{
 						third_prob			=0.0;
 						float temp_summe 	=first_prob+second_prob+third_prob;
 						first_prob 			=first_prob/temp_summe;
 						second_prob 		=second_prob/temp_summe;
 					}
 					float temp_float= (first_prob+second_prob);
 					//cout << first_prob << " "<< second_prob <<" and sum:  "<< temp_float << endl; //debug
 					//cout << pCBSNP->allele1 << " "<< pCBSNP->allele2 <<endl; //debug
 					if((temp_float>=0.0) &&(temp_float<=1.1) && pCBSNP->allele1!=""&& pCBSNP->allele2!="")
 					{
 						 pCBSNP->pgeno1.push_back(first_prob);
 						 pCBSNP->pgeno2.push_back(second_prob);
 						 pCBSNP->pgeno3.push_back(third_prob);
 					 }
 					 else if( (temp_float>1.000001)  )
 					 {
 						 emsg="Sum of genotype distribution should not be greater than 1. But this is the case in SNP "+pCBSNP->snpName+" and individual " +temp_value+" \n";
 						 error(emsg);
 					 }
 					 else
 					 {
 						 emsg=" first and second allele of  "+pCBSNP->snpName+"  is not given. please update alleles information first. \n";
 						 error(emsg);

 					 }

 				 }


 			 }
 			 else continue;
 			 // assignment of genotypes
 			//	 for(unsigned int i=0;i<cZeile.size();++i) // debug
 			//	 cout << cZeile[i] << "  " ; // debug
 			//	 cout <<endl; // debug
 			//-----------------------------------------//
 			//	 cout << "test" <<endl;
 			//-------------------------------------------//
 			//lese_pedfile_genInfo(genInfo, cZeile,nLines);
 			//at the end des loops.
 			 ncol=0;
 			 cZeile.clear();
 			 cZeile.resize(0);
 			 if(_ifs_eof(pf))break;
 	} // end of outer while loop
 		_ifs_close(pf);
 		nrow=0; // total rows
 		nLines=0; // total lines
 		//update_pedInfo_given_gProb(pedInfo,genInfo);
 		genoProb_to_geno_converter(genInfo, thresh);
 		//string pmsg="\t**->Total number of successfully read SNPs: ";
 		//printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
 		//CBPED::display_ped_summary();
 		string pmsg="\t**->Total number of successfully read individuals: ";
 		printWithSpace(pmsg,  change_int_into_string(pedInfo.size()),coutWidth);
 		//-------------------------------------------//


 	 }catch(bad_alloc& memoryAllocationException)
 	 {
 		 cout<<"\n Please inform the programmer that you could not read mlgeno file. Probbaly he will give you a solution. \n ";
 		 error( "\n bad _alloc()  in reading mach mlgeno file \n");
  	 }

 return; }
// lese mach mlprob file
    //alternative
      void CMPED::lese_mach_geno_mlGenoFile( vector<CBPED*>& pedInfo, const vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName)
      {
     	 printLIN("*->Reading file \""+fileName+"\": \n");
     	 //CBPED* pCBPED;
     	 string myString					="";
     	 string temp_value					="";
     	 bool ende							= false;
     	 unsigned int ccol					=0;
     	 unsigned int ncol					=0;
     	 unsigned int nrow					=0;
     	 unsigned int nLines				=0;
     	 unsigned int max_ncol				=0;
     	 bool def_max_ncol_once				=true;
     	 unsigned long int	snpAnzahl		=0;
     	string	pheno_value					="";
     	string	sex_value					="";
     	 try
     	 {
     		IFS pf (fileName.c_str(),"r");
			if(!pf)
			{
				string msg= fileName+ " either does not exits or could not be opend.\n"	;
				error(msg);
			}

     		while(true)
     			 {
     				if(_ifs_eof(pf))break;
     				while(!_ifs_eof(pf))
     				{
     					ccol=CBSNP::stringLeser(pf, myString, ende);
     					ncol+=ccol;
     					 if(ccol!=0)
     					{
     						if(!myString.empty())
     							cZeile.push_back(myString);
     						else
     						{
     							cout << "error";
     							exit(1);
     						}
     					}
     					if(ende|| _ifs_eof(pf))
     					{
     					 //ende=false;
     						++nLines;
     						if(def_max_ncol_once &&(ncol!=0))
     						{
     							max_ncol=ncol;
     							def_max_ncol_once=false;
     							//cout << "max_ncol: "<<max_ncol <<" " <<"ncol: "<<ncol <<endl; // only test; // only test
     						}
     						if((ncol==max_ncol)&& (max_ncol!=0))
     							++nrow;
     						myString="";
     						//cout <<"-----------------\n";//test
     						break;
     					}
     					if(ccol==0)
     					 continue;
     					// this loop is broken only if end  is found.
     				}//end of inner while loop
     				 if(cZeile.size()!=0)
     			 {
     				 snpAnzahl=genInfo.size();
     				// cout << snpAnzahl <<endl;

     				 string pheno_value	="";
     				string sex_value	="";
     				long int spaltenZahl=cZeile.size();
     				// check if  vector<CBSNP*> genInfo has  size equal to snpAnzahl.
     				// this declaration will produce pointer to CBPED classes
     				//	CBPED* pCBPED=new CBPED;
     				unsigned long int value = spaltenZahl-2; //
     					//cout <<"value: "<<value<<endl;
     					//cout << boolalpha<< (value==snpAnzahl)<<endl;
     					//cout << boolalpha<< (value/2==snpAnzahl)<<endl;
     				bool tfValue =(value==snpAnzahl) ;
     				if(!tfValue)	// Fall mit 6 Zeilen.
     				{
     					string msg="\t**->Problem on the " + change_int_into_string(nrow)+ "th row.\n";
     					printLIN(msg);
     					msg="There must be 2 columns  for phenotype infos and   exactly "+change_int_into_string(snpAnzahl)+\
     					" columns  with \"A/B\" type of genotype Infos.\n";
     					//delete pCBPED;
     					 error(msg);

     				}else
     				{
     					CBPED* pCBPED=new CBPED;
     					myString 			=cZeile[0];
     					int temp_indx		=myString.rfind("->");
     					pCBPED->famId		=myString.substr(0,temp_indx);
     					temp_value			=myString.substr(temp_indx+2,myString.size());
     					pCBPED->indId		=temp_value;
     					pCBPED->patId		="";
     					pCBPED->matId		="";
     					//for sex
     					lese_sex_info(pCBPED, sex_value);
     					//for phenotype
     					lese_phenotype_info(pCBPED, pheno_value);
     					pedInfo.push_back(pCBPED);
     					//cout << "test" <<endl; //debug
     					myString			="";
     				}
     				 //-------------------------------------------------------------//

     				 //cout << "nrow: "<< nrow<< " pedInfo.size()"<< pedInfo.size()<<endl;
     				 if(nrow!=pedInfo.size())
     				 {
     					 int sz=pedInfo.size()-1;
     					 string msg= "Problem in saving genotypes of person "+pedInfo[sz]->indId+" in the "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
     					 error(msg);
     				 }
     				 if(max_ncol!=ncol)
     				 {
     					 string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
     					 error(msg);
     				 }
     			 // now erase the first 2 columns of pedInfo from ZeilenVector
     				 cZeile.erase(cZeile.begin(),cZeile.begin()+2);

     			 }
     			 else continue;
     			 // assignment of genotypes
     			//	 for(unsigned int i=0;i<cZeile.size();++i) // debug
     			//	 cout << cZeile[i] << "  " ; // debug
     			//	 cout <<endl; // debug
     			lese_pedfile_genInfo(genInfo, cZeile,nLines);
     			//at the end des loops.
     			 ncol=0;
     			 cZeile.clear();
     			 cZeile.resize(0);
     			 if(_ifs_eof(pf))break;
     			 } // end of outer while loop
     		_ifs_close(pf);
     		nrow=0; // total rows
     		nLines=0; // total lines
     		//string pmsg="\t**->Total number of successfully read SNPs: ";
     		//printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
     		string pmsg="\t**->Total number of successfully read individuals: ";
     		//display_ped_summary();
     		printWithSpace(pmsg,  change_int_into_string(pedInfo.size()),coutWidth);
     	 }catch(bad_alloc& memoryAllocationException)
     	 {
     		 cout<<"\n Please inform the programmer that you could not read mlgeno file. Probbaly he will give you a solution. \n ";
     		 error( "\n bad _alloc()  in reading mach mlgeno file \n");
      	 }

     return; }
    // lese mach mlprob file

    //end

      //

    //---------------------------------------------//


// to determine impute  ints in impute comamnds 
void fix_ints_4_impute_command(const vector<int>& bp, vector<int>&lower_range_vec, vector<int>&upper_range_vec,vector<float> & interval_size_MB){

	const float I_S 		=5.0; // interval_size;
	const float MAX_I_S		=7.0;// maiximum interval size
	const float   MB_UNIT 	=(1.0e6); // mega base unit 
	int nGroups				=0;
	const float TOTAL_SIZE	=(bp[bp.size()-1]-bp[0])/(MB_UNIT);
	//cout << "Total size:"<< (bp[bp.size()-1]-bp[0])<<" "<<TOTAL_SIZE <<endl;
	lower_range_vec.push_back(bp[0]);
	//if(TOTAL_SIZE>7.0)
	
	{
	
		//float curr_value = bp[0];
		//cout << "curr_value: "<< curr_value<<endl;
		//double x =19.598998; 
		//	double y= 5.0;
		const int dividend=static_cast<int>(TOTAL_SIZE/I_S);
		// cout <<"dividend: "<< dividend << endl;
		 float rem= TOTAL_SIZE-static_cast<float>(dividend)*I_S;
		// cout <<"remainder: " << rem<< endl;
		if(I_S+rem/dividend<=MAX_I_S)
			 nGroups =dividend;
		else
			 nGroups =dividend+1;
		//cout << "nGroups: "<<nGroups<<endl;
		const float eachSize = TOTAL_SIZE/(static_cast<float>( nGroups));
		// cout << "each size: "<< 	eachSize <<endl;
		 int upper_bp =0;
		for(int i=0;i<nGroups;i++)
		{
			if(i>0)
				lower_range_vec.push_back(((upper_range_vec[i-1]+1)));
			//cout << "lower.range_vec.size() "<<lower_range_vec.size() <<endl;
			upper_bp= lower_range_vec[i]+(int)round(eachSize*(MB_UNIT));
			//cout << "upper_bp = "<< upper_bp<<endl;
			if(i<(nGroups-1))
			{
				upper_range_vec.push_back(upper_bp);
				interval_size_MB.push_back(eachSize);
			}else if(i==(nGroups-1))
			{ 	
				upper_range_vec.push_back(bp[(bp.size()-1)]);
				//cout << upper_range_vec[i]<< "\t" <<lower_range_vec[i]<<endl;
				float ci_s =(upper_range_vec[i]-lower_range_vec[i])/(MB_UNIT);
				//cout <<"ci_s:\t"<< ci_s<<endl;
				interval_size_MB.push_back(ci_s);
			}
		
	
		}		
		
	
	//
		/* test
		for(unsigned int i =0;i<lower_range_vec.size();i++)
		{
			cout << "lower_range_vec: "<< lower_range_vec[i]<<"\t" ;
		cout << upper_range_vec[i]<<"\t" ;
		 cout << interval_size_MB[i]<<endl;

		}
		*/
		
	}

}	
