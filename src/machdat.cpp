/*
 * machdat.cpp
 *
 *  Created on: Dec 29, 2011
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */
#include "machdat.h"
#include <stdio.h>
#include <stdlib.h>
//#include "../helper.h"
//#include "../bpar.h"

 int CBSNP::col_width=55;
 extern ofstream LIN;
// extern CGENERAL GENERAL(Basicparam);
 //assign previously defiend variables which come very often

// vector<CBSNP*> CBSNP::genInfo;
ostream& operator<<(ostream &out, const CBSNP &pCBSNP)
	  {
	  	    return out << "chr No: "<< pCBSNP.nchr<< " "<< "SNP ID: "<< pCBSNP.snpName\
	  	    		<< " "<< "rs ID:"<< pCBSNP.rsId\
	               << "  cm_pos: " << pCBSNP.cm_pos<< "  bp: "<< pCBSNP.bp;
	  }
// plink way
CBSNP::~CBSNP()
	{
		/*
		if(genInfo.size()>0)
		{
			for (unsigned int i=0;i<genInfo.size();++i)
			{
				delete genInfo[i];
				//genInfo[i]=NULL;
			}
			genInfo.clear();
			genInfo.resize(0);
		}
		*/
		
	}

CMSNP::~CMSNP()
{
	/*
	if(genInfo.size()>0)
	{
		for (unsigned int i=0;i<genInfo.size();++i)
		{
			delete genInfo[i];
		}
		genInfo.clear();
	}
	*/
 	//cout << "\n Good bye from CMSNP destructor.\n";
 	//cout << "\n Good bye from CBSNP destructor.\n";
}

unsigned int CBSNP::no_of_persons=0;
unsigned int CBSNP::ngood_snps=0;	// no of good quality snps.
 vector<string> CBSNP::sifo_snpid; // snp info rsid
  vector<string> CBSNP::sifo_rsid; // snp info rsid
  vector<string> CBSNP:: sifo_nchr; // chr no.
  vector<int>    CBSNP::sifo_bp; //snp info position
  vector<double> CBSNP::sifo_cm_pos; // snp info cm_pos
  vector<string> CBSNP:: sifo_allele1;
 vector<string> CBSNP:: sifo_allele2;
  //bool CBSNP::given_pgeno		=false;
  bool CMSNP::_is_code_mach		=false;
  bool CBSNP::_read_extramap	=false;
  bool CBSNP::_read_extraped	=false;
  string 	CBSNP::_miss_gvalue	="0";
  string 	CBSNP::_out_format	="";

string CBSNP::stringAssigntoRightPlace(const bool rowEnd, const int ccol, \
int & nrow, int& ncol,const string& stringValue)
{
	//rowEnd= end of row.
	// ccol= current columns either one or zero.
	// nrow = no of row coming
	// ncol= current column will be determined from ncol+=col
	// n'th col
	string cvalue("");
	//while(cvalue.empty())
	//{
		if(rowEnd)
		{
			if(ncol!=0)
			{
				if(ccol!=0)
				{
					ncol+=ccol;
					cvalue=stringValue;
					//for test
					/*
					cout << "after assignment: \n";
					cout << ncol<<"th col : "<< "ccol: " << ccol ;
					cout << " zeilenEnde : " << boolalpha<< rowEnd;
					cout << "\t string: "<< cvalue<<endl;
					*/
				///++nrow;
				}
				//else
				//cout << "Total columns: "<< ncol;
				++nrow;
				//cout << "\t end of " << nrow <<"th row."<<endl;

			}

			//ncol=0; // this will be normally done outside of the function
			if(ccol==0 && stringValue.empty())
			{

			}
		}
		else if (ccol==0 )
		{
			if (!stringValue.empty() && ncol!=0)
			{
			cout << " problem with reading the end of  file (EOF)! "<<endl;
			//nthCol=ncol;
			//ncol=0;
			exit(1);
			}
			else if (ccol==0 && ncol==0 )
			{
			 cout << " problem with reading the end of  file (EOF)! "<<endl;

			}
		}
		else
		{

			ncol+=ccol;
			cvalue=stringValue;
			// for test
			/*
			cout << "after assignment: \n";
			cout << ncol<<"th col : "<< "ccol: " << ccol ;
			cout << " zeilenEnde : " << boolalpha<< rowEnd;
			cout << "\t string: "<< cvalue<<endl;
			*/
		}
	//if(cvalue!="") break;
	//}

return cvalue; }
int CBSNP::stringLeser(IFS& pf , string& sValue, bool& rowEnd)
{
	sValue=""; // value o f string in each column
	int nstring=0; // no of columns in a row
	char myChar='\0';
	bool multi_whChar=false; // this will detect multi white space characters and avoid them.
	bool 	ende=false; // end of line
while(true)
{
	if(_ifs_eof(pf))
		break;
	if(nstring==1) break;
		myChar=_ifs_getc(pf);
	  char _endchar='\r\n';
	if(myChar=='\n'|| myChar=='\r'|| myChar=='\t'||  myChar==' '|| _ifs_eof(pf)||(myChar==_endchar) )
	{
		//case1
		//cout << "myChar: "<< myChar<<" "<<(myChar=='\t')<< " "  << (myChar=='\n')<<" "<<(myChar=='\r') <<" "<< _ifs_eof(pf)<<endl; //debug
		if((myChar=='\n')||(myChar=='\r') || _ifs_eof(pf)||(myChar==_endchar))
		{
			myChar='\0';
			ende=true;
			if(!sValue.empty())
			{

				++nstring;
				break;

			}
			else break;
		}
		//case 2
		if((myChar=='\t')||(myChar==' '))
		{
			myChar='\0';
			if(!ende&& !multi_whChar && !sValue.empty()&& sValue!=" \n" )
			{
				++nstring;
				break;
			}
			else
			{
				multi_whChar=true;
				continue;
			}

		}
		
	}
	//case 2
	else if(myChar=='#')
	{
		ende=true;
		 // Ignore rest of line and advance to next line
		while(myChar!='\n'&& myChar!='\r'&& !_ifs_eof(pf))
			myChar =_ifs_getc(pf);
		//fscanf(pf,"%*[^\n]");
		//(void) _ifs_getc(pf);
		if(!sValue.empty())
			++nstring;
		myChar='\0';
		break;
	}
	else
	{ // normal case
		sValue+=myChar;
		multi_whChar=false;
		if(sValue.empty()) break;

	}

}
rowEnd=ende;

return nstring;
}

// void copy function
void CBSNP::copy(const CBSNP& pCBSNP)
{
	nchr			=pCBSNP.nchr;
	snpName			=pCBSNP.snpName;
	rsId			=pCBSNP.rsId;

	allele1			=pCBSNP.allele1;
	allele2			=pCBSNP.allele2;

	coding_strand	=pCBSNP.coding_strand;

	min_allele		=pCBSNP.min_allele;
	maj_allele		=pCBSNP.maj_allele;

	dose_allele1	=pCBSNP.dose_allele1;
	dose_allele2	=pCBSNP.dose_allele2;
	nchrobs			=pCBSNP.nchrobs;
	change_0123		=pCBSNP.change_0123;
	change_dose		=pCBSNP.change_dose;
	quality			=pCBSNP.quality;
	no_of_persons	=pCBSNP.no_of_persons;
	ngood_snps		=pCBSNP.ngood_snps;
	freq			=pCBSNP.freq;
	pvalue_hwe_exact=pCBSNP.pvalue_hwe_exact;
	cm_pos			=pCBSNP.cm_pos;
	bp				=pCBSNP.bp;
	nonMiss			=pCBSNP.nonMiss;
	snp_CR			=pCBSNP.snp_CR;

	geno1			=pCBSNP.geno1;
	geno2			=pCBSNP.geno2;

	aOrder			=pCBSNP.aOrder;
	strand			=pCBSNP.strand;

	pgeno1			=pCBSNP.pgeno1;
	pgeno2			=pCBSNP.pgeno2;
	pgeno3			=pCBSNP.pgeno3;

	given_pgeno		=pCBSNP.given_pgeno;
	geno_0123		=pCBSNP.geno_0123;
	geno_dose		=pCBSNP.geno_dose;
	//given_geno_0123	=pCBSNP.given_geno_0123;
 	//given_geno_dose	=pCBSNP.given_geno_dose;

	//GENOTYPE_MISSING =pCBSNP.GENOTYPE_MISSING;
	sifo_snpid		=pCBSNP.sifo_snpid;
	sifo_rsid		=pCBSNP.sifo_rsid;
	sifo_nchr		=pCBSNP.sifo_nchr;
	sifo_bp			=pCBSNP.sifo_bp;
	sifo_cm_pos		=pCBSNP.sifo_cm_pos;
	sifo_allele1	=pCBSNP.sifo_allele1;
	sifo_allele2	=pCBSNP.sifo_allele2;
	_read_extramap	=pCBSNP._read_extramap;
	_read_extraped	=pCBSNP._read_extraped;
	_miss_gvalue	=pCBSNP._miss_gvalue;
	_out_format		=pCBSNP._out_format;
	// next
	maf				=pCBSNP.maf;
	rsq				=pCBSNP.rsq;
	given_both_crate=pCBSNP.given_both_crate;
	given_snp_crate	=pCBSNP.given_snp_crate;
	given_indiv_crate=pCBSNP.given_indiv_crate;
	given_snp_hwe	=pCBSNP.given_snp_hwe;
	col_width		=pCBSNP.col_width;


return ;}
CBSNP::CBSNP(const CBSNP& pCBSNP){copy(pCBSNP); }
CBSNP& CBSNP::operator =(const CBSNP& pCBSNP){
	copy(pCBSNP);

return *this;}
void CMSNP::mapFileLeser(const string&  fileName,  vector<CBSNP*>& genInfo)
{
	//cout << "-----------------------start of dat file: --------------------------------\n";
	printLIN("*->Reading file named "+ fileName+ ": \n");
	IFS file(fileName.c_str(),"r");
		if(!file)
		{
			string msg= fileName+ " either does not exits or could not be opened.\n"	;
			error(msg);
		}
		// if file is okay then define more varaibles
		string elem(""); // gelesen von ped File;
		string elem1("");
		int ncol(0);
		int ccol(0);
		static int nrow=0;
		bool zeilenEnde=false;

		while(!_ifs_eof(file))
	{


		if(_ifs_eof(file))
		{
			zeilenEnde=true;
			if(nrow==0)
			{
				string msg=" There are no elements in the "+fileName +" file.\n";
				error(msg);
			}
			else
			{
				long int snpSize=genInfo.size();
				string msg= "\t**->total no of SNPs read successfully: "+ change_int_into_string(snpSize)+" \n";
				printLIN(msg);	nrow=0;
				//cout << "size of vector genInfo: "<<genInfo.size()<<endl; // only for test
			}
			nrow=0;
			break;
		}

		zeilenEnde=false;
		// read string
			while(!zeilenEnde && !_ifs_eof(file))
		{ // do all works till the end of line and then break
		  //
			if(zeilenEnde) break;
			if(_ifs_eof(file)) break;
			ccol=stringLeser(file , elem, zeilenEnde);
			elem1=stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
			/*//if(ccol==0) continue;
			cout << "#------#\n";
			cout <<" elm : " <<elem;
			cout <<"    elem1:  "  <<elem1;
			cout << "   ccol:  "<<	ccol;
			cout << "   ncol:  "<<	ncol;
			cout << "   nrow:  "<<	nrow;
			cout << "   zeilenende:  "<<	zeilenEnde;
			cout <<endl;
			*/
			if(_ifs_eof(file)) break; // break is okay because we do not want to read only one column
			if(zeilenEnde) break; // break is also okay because we do not want to read only one  column
			if(elem1.empty() )
			do
			{
				ccol=stringLeser(file , elem, zeilenEnde);
				elem1=stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol,elem);
				/*	//if(ccol==0) continue;
				cout << "#------#\n";
				cout <<" elm : " <<elem;
				cout <<"   elem1:  "  <<elem1;
				cout << "   ccol:  "<<	ccol;
				cout << "   ncol:  "<<	ncol;
				cout << "   nrow:  "<<	nrow;
				cout << "   zeilenende:  "<<	zeilenEnde;
				cout <<endl;
				*/
				if(zeilenEnde|| _ifs_eof(file))
				break;
				//ccol=0;
			}while(elem1.empty());
			if(( zeilenEnde || elem1.empty()) || _ifs_eof(file) )
			{
				if(!elem1.empty())
				{
					string msg= " There is only one column in the first row of "+fileName+" file. \n" ;
					error(msg);
				}
				else
				break; // only in the first column
			}
			// define  pointer of CBSNP
			CBSNP* pCBSNP= new CBSNP;

			if(!elem1.empty())
			{
				if(ccol!=0 && ncol==1 ) // this condition must always be true ;
				{
					//cout << "ncol ->1st:->>>>>>>>>>>>>>>>>>> " << ncol << "elem: "<< elem1<<endl;
					elem1="";
					ccol=0;
				}
				else
					{
						delete pCBSNP;
						string msg= "problem in reading mapfile and discarding Ms of mach dat file\n";
						error(msg);
					}


			}
			else
			{
				delete pCBSNP;
				string msg= "problem in reading  first column in the "+change_int_into_string(nrow+1)+ "of   file \""+fileName+"\".\n";
				error(msg);
			}
			// for the second column
			if(zeilenEnde || _ifs_eof(file))
			{	delete pCBSNP;
				string msg= "Not enough columns in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+fileName+"\".\n";
				error(msg);
			}
			elem="";
			elem1="";
			ccol=0;
			ccol=stringLeser(file , elem, zeilenEnde);
			elem1=stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
			/*cout << "#------#\n";
			cout <<" elm : " <<elem;
			cout <<"   elem1:  "  <<elem1;
			cout << "   ccol:  "<<	ccol;
			cout << "   ncol:  "<<	ncol;
			cout << "   nrow:  "<<	nrow;
			cout << "   zeilenende:  "<<	zeilenEnde;
			cout <<endl;
			*/
			//if(zeilenEnde || _ifs_eof(file))
			//{  this case is true  only at  the middle of the line.
			//	string msg= "Not enough columns in the line of "+change_int_into_string(nrow+1)+"th row. \n";
			//	error(msg);
			//}

			if(elem1.empty())
			{
				delete pCBSNP;
				string msg= "Problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+fileName+"\".\n";
				error(msg);
			}
				//if(zeilenEnde|| _ifs_eof(file))
				//{ //the other cases are true only in the middle of the line
					//string msg= "Not enough columns in the line of "+change_int_into_string(nrow+1)+"th row. \n";
					//error(msg);
				//}
			else
			{
				if(ccol!=0 && ncol==2 ) // this condition must always be true ;
				{
					pCBSNP->snpName	=elem;
					if(pCBSNP->rsId=="")
						pCBSNP->rsId=elem;
					//pCBSNP->rsId	=elem;
					//cout << "ncol ->2ndst:->>>>>>>>>>>>>>>>>>> " << ncol << " \t elem: "<< elem<<endl;
					elem1="";
					ccol=0;
				}
				else if(zeilenEnde) // at the end of line ncol is zero again .
				{
					pCBSNP->snpName=elem;
					//pCBSNP->rsId	=elem;
					//cout << "ncol ->2ndst:->>>>>>>>>>>>>>>>>>> " << ncol << " \t elem: "<< elem<<endl;
					elem1="";
					ccol=0;

				}
				else{
					delete pCBSNP;
					string msg= "problem in reading  second column of mapfile and discarding Ms of mach dat file\n";
					error(msg);
					}
			//   This case is true except last column
				//if(zeilenEnde || _ifs_eof(file))
				//{	string msg= "Not enough columns in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+pedFileName+"\".\n";
				//error(msg);
				//}


			}
			if(pCBSNP->snpName=="")
				delete pCBSNP;
			else
				genInfo.push_back(pCBSNP);

		//--------------------------------------------------------------#
		// here must be the end of line
			if(ncol>2|| ncol<2)
			{
				string msg="";
				if(ncol<2 )
					 msg= "Not enough columns in the line of "+change_int_into_string(nrow+1)+"th row. \n";
				else
					 msg= "More than 2 columns in the line of "+change_int_into_string(nrow+1)+"th row. \n";
				error(msg);
			}
		// if there is end of line then break;
			if(zeilenEnde|| _ifs_eof(file))
				break;

		}


			ncol		=0;
			ccol		=0;

		if(_ifs_eof(file))
		{
			zeilenEnde=true;
			if(nrow==0)
			{
				string msg=" There are no elements in  "+fileName +" file.\n";
				error(msg);
			}
			else
			{
				//long int snpSize=genInfo.size();
				//string msg= "\t**->total no of SNPs read successfully: "+ change_int_into_string(snpSize)+" \n";
				//printLIN(msg);
				nrow=0;
				//cout << "size of vector genInfo: "<<genInfo.size()<<endl; // only for test
			}
			nrow=0;
			break;
		}

	}

	//end now
		//genInfo.push_back(pCBSNP);
		_ifs_close(file);
	 //display_snp_summary(genInfo);
		//for(vector<CBSNP*>::iterator it=genInfo.begin(); it!=genInfo.end();++it)
		//cout << "SNPs read: \n";
		//for(unsigned int i=0; i<genInfo.size();++i)
		//cout << genInfo[i]->snpName << "  " <<endl;
		//cout << endl;

return ;}

//function to convert genotypes into 0, 1,2 and 3 according as  minor allele
//needs geno_123 vec and calculate freq.
void CBSNP::convert_snp_genotype_into_0123(CBSNP* const pCBSNP)
{
	/* sometimes major allele homozygote as 0 hetero as 1, and minor allele homo as
	 * 2. but here I haven't done it.
	 * I have saved as allele1 homo as 0 and allele2  homo as 2 and hetero as 1
	 * one should change the notation according as necessary
	 */
	if(pCBSNP->change_0123)
	{
			pCBSNP->change_0123=false;
			const unsigned int _nPerson=pCBSNP->geno1.size();
			if(!pCBSNP->geno_0123.size())
			{
				//writing genotypes
				if(pCBSNP->geno1.size()!=pCBSNP->geno2.size())
				{
					cout <<"geno1.size() = "<<pCBSNP->geno1.size() <<",and  geno2.size()=  "<<pCBSNP->geno1.size()<<endl;
					//for(unsigned int i=0;i<pCBSNP->geno1.size();++i)
					//	ofs<< pCBSNP->geno1[i] << " " <<pCBSNP->geno2[i] << " ";
					error("Problem in writing out SNP with rsid \""+pCBSNP->rsId +"\". geno1 and geno2 do not have equal size.\n");
				}
				else if(pCBSNP->allele1!="" || pCBSNP->allele2!="") // if alleles are not empty
				{
					for(unsigned int i=0;i<_nPerson;++i)
					{
						if( !pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if homozygote1  00
							pCBSNP->geno_0123.push_back(0);
						else if( pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if homozygote2 11
							pCBSNP->geno_0123.push_back(2);
						else if( !pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if heterozygote
							pCBSNP->geno_0123.push_back(1);
						else if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if missing
							pCBSNP->geno_0123.push_back(3);
						else
						{
							printLIN("Problem in saving  SNP with rsid \""+pCBSNP->rsId +"\"  and "+change_int_into_string(i)+"th individual :\n");
							error("Please contact the software developer. \n");
						}
									//ofs << "\n"; //debug
					}
				}
				else{
						printLIN("It seems like that your SNP is completely missing with no information on alleles."
								"It can also be an internal error! Hints: alleles are empty."+
								pCBSNP->rsId +" "+change_int_into_string(pCBSNP->bp)+" \n");
						for(unsigned int i=0;i<_nPerson;++i)
						{
							if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if missing
								pCBSNP->geno_0123.push_back(3);
							else
							{
								printLIN("Problem in saving  SNP with rsid \""+pCBSNP->rsId +"\"  and "+change_int_into_string(i)+"th individual :\n");
								//error("Please contact the software developer. \n");

								//printLIN("first allele: "+ pCBSNP->allele1+"  second allele : "+ pCBSNP->allele2 +" . \n");
								error("Please  check if this SNP is completely missing, if yes, give allele info with --snpinfo option.\n"
										"If not, then contact the software provider.\n");
							}
						}

				}
				// find out if maf is correct
				//calculate_maf( pCBSNP);
				//cout<<"allelle freq: "<< pCBSNP->freq<< " allele1: "<< pCBSNP->allele1<< " \n"; //debug
				//for(unsigned int jj=0;jj< (pCBSNP->geno_0123.size()); ++jj ) //debug
				//	cout <<pCBSNP->geno_0123[jj]<< ",  "; //debug
				//cout <<endl; //debug
			}

	}

}//end of function
//xxxxxxxxxxxxxxxxx
void CBSNP::convert_snp_genotype_into_dose(CBSNP* const pCBSNP,const string &ref_allele)
{

	float _temp_dose=0.0;
	if(pCBSNP->change_dose)
	{
			pCBSNP->change_dose=false;
			const unsigned int _nPerson=pCBSNP->geno1.size();
			const bool _temp_probs_given =(_nPerson==pCBSNP->pgeno1.size())&&
										  (_nPerson==pCBSNP->pgeno2.size())&&
										  (_nPerson==pCBSNP->pgeno3.size());
			//cout << " geno.size():"<<pCBSNP->geno1.size()<<endl;
			if(!pCBSNP->geno_dose.size()&& _temp_probs_given)
			{

				// find out if maf is correct
					if(pCBSNP->quality)
					{
							if(pCBSNP->change_0123)
							// this function also calculates maf
						 convert_snp_genotype_into_0123(pCBSNP);
						//ref_allele is the count of the allele that means the dose_allele2
						//cout <<"ref_allele:"<<ref_allele <<endl;
						else if((ref_allele=="") ||ref_allele=="minor-allele") //this is for general case
						{
							pCBSNP->dose_allele1= pCBSNP->maj_allele; 	//
							pCBSNP->dose_allele2= pCBSNP->min_allele; 	//
						}
						else if((ref_allele=="allele2")) //this is for general case
						{
							pCBSNP->dose_allele1= pCBSNP->allele1; //
							pCBSNP->dose_allele2= pCBSNP->allele2; //
						}
						//for ref allele is major  allele;
						else if(ref_allele=="allele1") //this is for general case
						{
							pCBSNP->dose_allele1= pCBSNP->allele2; //
							pCBSNP->dose_allele2= pCBSNP->allele1; //
						}
						//for ref allele is major  allele;
						else if(ref_allele=="major-allele") //this is for general case
						{
								pCBSNP->dose_allele1= pCBSNP->min_allele; //
								pCBSNP->dose_allele2= pCBSNP->maj_allele; //
						}
						//for ref allele is major  allele;
						//cout <<"ref_allele:"<<ref_allele <<endl;
						//	cout <<"pCBSNP->allele1: "<< pCBSNP->allele1 << ", pCBSNP->allele2: "<< pCBSNP->allele2 <<endl; //debug
						//cout <<"pCBSNP->dose_allele1: "<< pCBSNP->dose_allele1 << ", pCBSNP->dose_allele2: "<< pCBSNP->dose_allele2 <<endl; //debug
						//cout<<"\n--------------\n";
						if(pCBSNP->dose_allele1==pCBSNP->allele1)
						{
							for(unsigned int i=0; i<_nPerson;++i)
							{
								if((pCBSNP->pgeno1[i]+pCBSNP->pgeno2[i]+pCBSNP->pgeno3[i])==0.0)
										pCBSNP->geno_dose.push_back(-1.0);
								else
								{
									_temp_dose =0.0*pCBSNP->pgeno1[i]+1.0*pCBSNP->pgeno2[i]+2.0*pCBSNP->pgeno3[i];
									pCBSNP->geno_dose.push_back(_temp_dose);
									//cout <<pCBSNP->pgeno1[i]<<" "<<pCBSNP->pgeno2[i]<<" " <<pCBSNP->pgeno3[i]<< " "<< _temp_dose<<endl ;
									_temp_dose=0.0;
								}
							}
							//_temp_dose =

						}
						else if(pCBSNP->dose_allele1==pCBSNP->allele2)
						{
							for(unsigned int i=0; i<_nPerson;++i)
							{
								if((pCBSNP->pgeno1[i]+pCBSNP->pgeno2[i]+pCBSNP->pgeno3[i])==0.0)
									pCBSNP->geno_dose.push_back(-1.0);
								else
								{
									_temp_dose =0.0*pCBSNP->pgeno3[i]+1.0*pCBSNP->pgeno2[i]+2.0*pCBSNP->pgeno1[i];
									pCBSNP->geno_dose.push_back(_temp_dose);
									//cout <<pCBSNP->pgeno1[i]<<" "<<pCBSNP->pgeno2[i]<<" " <<pCBSNP->pgeno3[i]<< " "<< _temp_dose ;
									_temp_dose=0.0;
								}
							}
						}
						else
						{
							printLIN("Something is wrong  in saving genotypes.\n"
									  "Dose reference allele must be equal to either allele1 or allele2  of the SNP:"+pCBSNP->rsId+".\n"
									  "But this is not the case here. Please contact the software provider.\n  ");
							error("Problem in converting genotypes into dose value");
						}

						//cout<<"allelle freq: "<< pCBSNP->freq<< " allele1: "<< pCBSNP->allele1<< " \n"; //debug
						//for(unsigned int jj=0;jj< (pCBSNP->geno_0123.size()); ++jj ) //debug
						//	cout <<pCBSNP->geno_0123[jj]<< ",  "; //debug
						//cout <<endl; //debug
				}
			}

			//writing genotypes
			else //if(!_temp_probs_given) //if(pCBSNP->pgeno1.size()!=pCBSNP->pgegeno2.size())
			{
				cout <<"pgeno1.size() = "<<pCBSNP->pgeno1.size() <<", pgeno2.size() = "<<pCBSNP->pgeno2.size()<< " and  pgeno3.size()=  "<<pCBSNP->pgeno3.size()<<endl;
				//for(unsigned int i=0;i<pCBSNP->geno1.size();++i)
				//	ofs<< pCBSNP->geno1[i] << " " <<pCBSNP->geno2[i] << " ";
				error("Problem in writing out SNP with rsid \""+pCBSNP->rsId +"\". genotype probabilities are not saved.\n");
			}

	}


}//end of function
//xxxxxxxxxxxxxxxxx
void CBSNP::convert_hardcalls_into_genoprob(CBSNP* const pCBSNP)
{
	const int _sz_pgeno =pCBSNP->pgeno1.size();

	if(!_sz_pgeno&& !pCBSNP->given_pgeno)
	{
		 pCBSNP->given_pgeno=true;
		for(unsigned int ii =0;ii<pCBSNP->geno1.size();++ii)
		{
			if( !pCBSNP->geno1[ii] &&  !pCBSNP->geno2[ii]  ) // if homozygote1  00
			{
				 pCBSNP->pgeno1.push_back(1.0);
				 pCBSNP->pgeno2.push_back(0.0);
				 pCBSNP->pgeno3.push_back(0.0);
			}
			else if( pCBSNP->geno1[ii] &&  pCBSNP->geno2[ii]  ) // if homozygote2 11
			{
				 pCBSNP->pgeno1.push_back(0.0);
				 pCBSNP->pgeno2.push_back(0.0);
				 pCBSNP->pgeno3.push_back(1.0);
			}
			else if( !pCBSNP->geno1[ii] &&  pCBSNP->geno2[ii]  ) // if heterozygote
			{
				 pCBSNP->pgeno1.push_back(0.0);
				 pCBSNP->pgeno2.push_back(1.0);
				 pCBSNP->pgeno3.push_back(0.0);
			}
			else if( pCBSNP->geno1[ii] &&  !pCBSNP->geno2[ii]  ) //if missing
			{
				 pCBSNP->pgeno1.push_back(0.0);
				 pCBSNP->pgeno2.push_back(0.0);
				 pCBSNP->pgeno3.push_back(0.0);
			}

		}//end of for
	}//end of if
}
//xxxxxxxxxxxxxxxxx
void CBSNP::convert_genoprob_into_hardcalls(CBSNP* const pCSNP, const float &thresh)
{
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
			 string msg= "Problem in converting probabilities of genotype format to hard call.\n";
			 error(msg);
		 }

	 }

return ; }

// calcuate freq //needs geno_123 vector .

//snp callrate
void CBSNP::calculate_snp_callrate(CBSNP* const _pCBSNP )
{

	if(!_pCBSNP->given_snp_crate)
	{
		convert_snp_genotype_into_0123(_pCBSNP);
		//cout<< count( _pCBSNP->geno_0123.begin(),_pCBSNP->geno_0123.end(),0)<<endl;
		const int THREES	=count( _pCBSNP->geno_0123.begin(),_pCBSNP->geno_0123.end(),3);
		//	cout <<THREES <<endl;
		_pCBSNP->snp_CR =1.0-(float)(1.0*THREES/float(_pCBSNP->geno_0123.size()));
	}
	_pCBSNP->given_snp_crate=true;
 }
//ind call rate

// hardy
/*
// This function implements an exact SNP test of Hardy-Weinberg
// Equilibrium as described in Wigginton, JE, Cutler, DJ, and
// Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
// Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton
*/

double CBSNP:: SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
{

	  if (obs_hom1 + obs_hom2 + obs_hets == 0 ) return 1;

	  if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
			{
			printf("FATAL ERROR - SNP-HWE: Current genotype configuration (%d  %d %d ) includes a"
				   " negative count", obs_hets, obs_hom1, obs_hom2);
			exit(EXIT_FAILURE);
			}
	  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

	  int rare_copies = 2 * obs_homr + obs_hets;
	  int genotypes   = obs_hets + obs_homc + obs_homr;

	  double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
	  if (het_probs == NULL)
		error("Internal error: SNP-HWE: Unable to allocate array" );

	  int i;
	  for (i = 0; i <= rare_copies; i++)
		het_probs[i] = 0.0;

	  /* start at midpoint */
	  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

	  /* check to ensure that midpoint and rare alleles have same parity */
	  if ((rare_copies & 1) ^ (mid & 1))
		mid++;

	  int curr_hets = mid;
	  int curr_homr = (rare_copies - mid) / 2;
	  int curr_homc = genotypes - curr_hets - curr_homr;

	  het_probs[mid] = 1.0;
	  double sum = het_probs[mid];
	  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
		{
		  het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
		/ (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
		  sum += het_probs[curr_hets - 2];

		  /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
		  curr_homr++;
		  curr_homc++;
		}

	  curr_hets = mid;
	  curr_homr = (rare_copies - mid) / 2;
	  curr_homc = genotypes - curr_hets - curr_homr;
	  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
		{
		  het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
		/((curr_hets + 2.0) * (curr_hets + 1.0));
		  sum += het_probs[curr_hets + 2];

		  /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
		  curr_homr--;
		  curr_homc--;
		}

	  for (i = 0; i <= rare_copies; i++)
		het_probs[i] /= sum;

	  /* alternate p-value calculation for p_hi/p_lo
	   double p_hi = het_probs[obs_hets];
	   for (i = obs_hets + 1; i <= rare_copies; i++)
		 p_hi += het_probs[i];

	   double p_lo = het_probs[obs_hets];
	   for (i = obs_hets - 1; i >= 0; i--)
		  p_lo += het_probs[i];

	   double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
	  */

	  double p_hwe = 0.0;
	  /*  p-value calculation for p_hwe  */
	  for (i = 0; i <= rare_copies; i++)
		{
		  if (het_probs[i] > het_probs[obs_hets])
		continue;
		  p_hwe += het_probs[i];
		}

	  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

	  free(het_probs);

	  return p_hwe;
}


void CBSNP::schreibe_snp_callrate(const vector<CBSNP*>&genInfo,const string & outFileName)
{
	CBSNP* pCBSNP =genInfo[0];
	string crFileName =outFileName+"_snp.crate";
	ofstream ofs;
	ofs.clear(); ofs.close();

	ofs.open(crFileName.c_str(),ios::out);
	if(!ofs)
		error("out file "+outFileName+" does not exit or can not be written out. Please check the filename and file path.\n");

	ofs <<"CHR"<<" "<< "SNP"<<" "<<"CRATE" ;
	ofs<<"\n";
	for(unsigned int j=0; j<genInfo.size();++j)
	{
		pCBSNP =genInfo[j];
		calculate_snp_callrate(pCBSNP);
		ofs<<pCBSNP->nchr <<" "<<pCBSNP->rsId <<" "<< pCBSNP->snp_CR;
		ofs<<"\n"; //end of line

	}
	pCBSNP= genInfo[0];
	ofs.clear();ofs.close();
	printLIN("\t**->File containing SNP callrate has been written out and saved as \""+\
			crFileName+"\". \n");

}//end of schreibe_snp_callrate

void CBSNP::schreibe_snp_hwe (const vector<CBSNP*>&genInfo,const string & outFileName){
	CBSNP* pCBSNP =genInfo[0];
	string hweFileName =outFileName+"_snp.hwe";
	ofstream ofs;
	ofs.clear(); ofs.close();
	ofs.open(hweFileName.c_str(),ios::out);
	if(!ofs)
		error("out file "+outFileName+" does not exit or can not be written out. Please check the filename and file path.\n");
	ofs <<"CHR"<<" "<< "SNP"<<" "<<"PVALUE_EXACT" ;
	ofs<<"\n";
	for(unsigned int j=0; j<genInfo.size();++j)
		{

			pCBSNP =genInfo[j];
			//calculate_snp_hwe(pCBSNP); // this will be done later
			if(pCBSNP->quality)
			{
				ofs<<pCBSNP->nchr <<" "<<pCBSNP->rsId <<" "<< pCBSNP->pvalue_hwe_exact;
				ofs<<"\n"; //end of line
			}
		}
	// make sure that hardy wont be calculated again.
	pCBSNP= genInfo[0];
	ofs.clear();ofs.close();
	printLIN("\t**->File containing SNP callrate has been written out and saved as \""+\
			hweFileName+"\". \n");



}

void CBSNP::write_SNP_frequency_file(const vector<CBSNP*>& genInfo, const string& outFileName){
	ofstream ofs;
	ofs.clear();ofs.close();
	string freqFileName=outFileName+"_snp.frq";
	ofs.open(freqFileName.c_str(),ios::out);
	if(!ofs)
		error("out file "+outFileName+" does not exit or can not be written out. Please check the filename and file path.\n");
	ofs <<"CHR"<< " "<< "SNP"<<" "<< "A1"<<" "<< "A2"<<" "<< "MAF" << " "<<"NCHROBS\n" ;
	CBSNP* pCBSNP= genInfo[0];
	for(unsigned int i=0; i<genInfo.size();++i)
	{
		pCBSNP= genInfo[i];
		//convert_snp_genotype_into_0123(pCBSNP);
		//calculate_maf(pCBSNP);
		//cout<<boolalpha<< pCBSNP->quality<<" ";
		if(pCBSNP->quality ) // default value of quality is true;
		{
			ofs<< pCBSNP->nchr<<" ";
			if(pCBSNP->rsId!="")
				ofs<<pCBSNP->rsId<<" ";
			else
				ofs<<pCBSNP->snpName<<" ";
			if(pCBSNP->min_allele!="" )
				ofs<< pCBSNP->min_allele <<" ";
			else ofs <<"0"<< " ";
			if(pCBSNP->maj_allele!="")
				ofs<<pCBSNP->maj_allele<<" ";
			else
				ofs<<"0"<<" ";
			ofs<<pCBSNP->maf <<" ";
			ofs<<pCBSNP->nchrobs;
			ofs<<"\n"; //end of line
		}
	}
	pCBSNP= genInfo[0];
	//end of writing
	ofs.clear();ofs.close();
	printLIN("\t**->File containing allele frequency has been written out and saved as \""+\
			freqFileName+"\". \n");
}//end of function write freq.

//#####################################################
void CBSNP::write_plinkMapFile(const vector<CBSNP*>& genInfo, const string& outFileName)
{
	ofstream ofs;
	ofs.clear();ofs.close();
  	string datFileName=outFileName+".map";
  	ofs.open(datFileName.c_str(),ios::out);
  	if(!ofs)
  		error("out file "+outFileName+" does not exit or can not be written out.\n");
  	for(unsigned int i=0; i<genInfo.size();++i)
  	{
  		CBSNP* pgenInfo= genInfo[i];
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
				ofs<< "0.0"<< " ";
			}
			if(pgenInfo->bp!=-1  ) ofs<< pgenInfo->bp<< "\n";
			else
			{
				ofs << i << "\n";

			}
		}
  	}
  	ofs.clear();ofs.close();
  	printLIN("\t **->PLINK map file has been written out and saved as  \""+\
	  			datFileName+"\".\n");

return;}
//-------------------------------------------------------------------//
// the following function will be used for haploview info file
void CBSNP::write_rsid_basePairposition(const vector<CBSNP*>genInfo, const string& outFileName)
{
	//this function is for different formats
	ofstream ofs;
	ofs.clear();ofs.close();
	string datFileName="";
	if(CBSNP::_out_format=="haploview") // haploview has an end name *.info
		datFileName=outFileName+".info";
	else
		datFileName=outFileName+".map";
	ofs.open(datFileName.c_str(),ios::out);
	  	if(!ofs)
	  		error("out file "+outFileName+" does not exit or can not be written out.\n");
	  	for(unsigned int i=0; i<genInfo.size();++i)
	  	{
	  		CBSNP* pgenInfo= genInfo[i];
	  		if(pgenInfo->quality)
	  		{
				//ofs<< pgenInfo->nchr << " ";
				if(pgenInfo->rsId!="") ofs <<pgenInfo->rsId<< " ";
				else if(pgenInfo->snpName!="") ofs <<pgenInfo->snpName<< " ";
				else if(pgenInfo->rsId==""&& pgenInfo->snpName=="")
				{
					string snpNewName="snp"+change_int_into_string(i)+"";
					ofs<< snpNewName << " ";
				}

				if(pgenInfo->bp!=-1  ) ofs<< pgenInfo->bp<< "\n";
				else
				{
					ofs << i << "\n";

				}
			}
	  	}
	  	ofs.clear();ofs.close();


}
// allele info leser 
void CBSNP::allele_info_leser(const string& a1, const string& a2,   CBSNP* const it_pCBSNP, const int& nrow, const string& fileName)
{
	//This function  
	if((*it_pCBSNP).allele1!="" && (*it_pCBSNP).allele1!="0")
	{
		//---------------------------------------------------//
		/*
		if((*it_pCBSNP).allele1 !=a1)
		{
			bool temp_bool1=(*it_pCBSNP).allele1 ==a2;
			bool temp_bool2= (*it_pCBSNP).allele2 ==a1;
			bool temp_bool3=(*it_pCBSNP).allele2 ==a2;
			bool temp_bool4= ((*it_pCBSNP).allele2=="");
			if(temp_bool1 & (temp_bool2||temp_bool3 || temp_bool4))
			{
				//cout << "change of allele! allele 1 is allele2 and allele 2 will be allele1 \n"; //ONly for test
				for(unsigned int i=0; i<(*it_pCBSNP).geno1.size();++i)
				{
					//homo1.homo2
					if(!((*it_pCBSNP).geno1[i])&& !((*it_pCBSNP).geno2[i]))
					{
						((*it_pCBSNP).geno1[i])	=true;
						((*it_pCBSNP).geno2[i])	=true;
					}
										//homo2 .homo1
										else if(((*it_pCBSNP).geno1[i])&& ((*it_pCBSNP).geno2[i]))
										{
											(*it_pCBSNP).geno1[i]=false;
											(*it_pCBSNP).geno2[i]=false;
										}
										else if(!((*it_pCBSNP).geno1[i])&& ((*it_pCBSNP).geno2[i]))
										{
											if((*it_pCBSNP).aOrder[i])
												(*it_pCBSNP).aOrder[i]=false;
											if(!((*it_pCBSNP).aOrder[i]))
												(*it_pCBSNP).aOrder[i]=true;
										}
										// change probability genotypes also
										if((*it_pCBSNP).pgeno1.size()>0 && (*it_pCBSNP).pgeno3.size()>0)
										{
											//cout << no_of_persons <<endl;
											for(unsigned int j=0;j<(*it_pCBSNP).pgeno1.size(); ++j)
												swap((*it_pCBSNP).pgeno1[j],(*it_pCBSNP).pgeno3[j]);
												  	 //swap will exchange prob(AA) and prob(BB).prob(AB) remains the same
										}

									}
									// change the allele also
									(*it_pCBSNP).allele1 =a1;
									(*it_pCBSNP).allele2 =a2;
									//only for test
									//cout <<"allele1: "<<(*it_pCBSNP).allele1 <<endl ;//debug
									//cout <<"allele2: "<<(*it_pCBSNP).allele2 <<endl ;//debug


								}
								else
								{
									cout <<"Allele in original files: " << 	(*it_pCBSNP).allele1 << " " <<(*it_pCBSNP).allele2 <<endl;
									cout <<"Allele in snpinfo files: " << 	a1 << " " <<a2 <<endl;
									string msg= " Allele1 \""+(*it_pCBSNP).allele1+\
											"\" in datfile  does not match with the allele \""+\
											a1+" given in "+change_int_into_string(nrow)+"th line of  SNP info file\""+fileName+\
											"\" \n.";
									error(msg);
								}
		}
		*/						//-----------------------------------------//
	}
	//else
	(*it_pCBSNP).allele1 =a1; // this condition never comes normally
	 cout << (*it_pCBSNP).allele1 <<"  allele 1: \n" ;

return ;} 

//extra file lesers.
void CBSNP::extraMapFileLeser(const string& fileName, const vector<CBSNP*> & genInfo)
{
	// The following for loop is just for test 
	/*
	 for(unsigned int i=0;i<genInfo.size();++i)
	 {
		cout << genInfo[i]->rsId <<" \t" << genInfo[i]->allele1 << " "<< genInfo[i]->allele2 << " \n";
	 
	 }
	*/
	string r_msg="\n*->Reading SNP info  file \""+fileName+ "\":\n";
	printLIN(r_msg);
	string line;
	rFile file;	 file.close(); file.clear(); 	// check if it has been already opend.
	static int extraMapFileLineZaehler	=0;			 // extra map file line zaehler
	//ios_base::iostate i= file.rdstate();
		file.open(fileName.c_str(),ios::in);
		// this will check if each row has same number of columns.
	if(!file)
		{
			file.setstate(ios_base::failbit);
			file.clear();
			error("your extra  file "+fileName+ " either does not exit or could not be opend.");
		}
	// check first header and  the number of columns
		do
		{
			getline(file,line, '\n');
			if(line[0]=='#') // leave the comments;
			continue;
			if(file.eofbit) break;
		}while( line.empty());
		string buffer="";
		vector<string> header;
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
		// this part is only for only mach hapmap data 
		if(header.size()==1 && (header[0]!="rsid" || header[0]!="rs" ||header[0]!="snpid"||header[0]!="rsId"||header[0]!="snpId"||header[0]!="rsID"||header[0]!="snpID"))
		{//cout << "This part is only for mach hapmap \n";
			vector<CBSNP*>::const_iterator it_pCBSNP =genInfo.begin();
			if((*it_pCBSNP)->rsId=="")
				(*it_pCBSNP)->rsId=header[0];
			else if((*it_pCBSNP)->rsId!=header[0])
			{
				cerr << "Problem in reading "<<++extraMapFileLineZaehler<<"th line : \n";
				LIN<<"Problem in reading "<<++extraMapFileLineZaehler<<"th line : \n";
				string pr_msg="rsid "+(*it_pCBSNP)->rsId+", which is already given , does not match with  rsid "+header[0]+".\n";
				error(pr_msg);
			} 
			//cout << (*it_pCBSNP)->rsId	<<" ";
			++it_pCBSNP;
			vector<string> tokens;
			
			while(getline(file,line, '\n')) // until eof hit
			{
			//---------------------//
				++extraMapFileLineZaehler;
				//cout <<"line no: "<< extraMapFileLineZaehler<<endl; //debug
				if(line.length()==0)
				continue;
				if(line[0]=='#')
				continue;
				//if(!((*it_pCBSNP)->quality))
				//	continue;
				tokens.resize(0);
				tokens.clear();
				stringstream line_parser1(line);
				//cout <<"\n line_parser1.good(): "<<line_parser1.good() <<endl; // only for test
				if(line_parser1.good())
				{
					while(line_parser1>>buffer)
					{
						tokens.push_back(buffer);
						//cout << buffer << " "; // only for test
					}
					if(tokens.size()!=1)
					{
						string msg= "Your extra map file \""+fileName+"\"  has 1 column in first line but in the line "+change_int_into_string(extraMapFileLineZaehler)+\
								" there are "+change_int_into_string(tokens.size())+" columns.\n";
						error(msg);
					}
					if((*it_pCBSNP)->rsId=="")
						(*it_pCBSNP)->rsId=tokens[0];
					else if((*it_pCBSNP)->rsId!=tokens[0])
					{
						cerr << "Problem in reading "<<++extraMapFileLineZaehler<<"th line : \n";
						LIN<<"Problem in reading "<<++extraMapFileLineZaehler<<"th line : \n";
						string pr_msg="rsid "+(*it_pCBSNP)->rsId+", which is already given , does not match with  rsid "+tokens[0]+".\n";
						error(pr_msg);
					}	
					
				}
				//cout << (*it_pCBSNP)->rsId	<<" ";
				++it_pCBSNP;
				if(it_pCBSNP==genInfo.end()) break;
				
			}
			if(it_pCBSNP!=genInfo.end())
			{
				cerr <<"Total rsid in "+fileName+": "<< (int)(it_pCBSNP-genInfo.begin())<< endl;
				cerr << "Total SNPs read from pedigree file: "<< genInfo.size() << endl;
				string msg= "Your  file \""+fileName+"\"  has not enough rsid Names given one in one column.\n";
						error(msg);
			
			}	
			header.clear();
			tokens.clear();
			// here is the end of hapmap snpinfo file !!
			
		}
		else
		{
			//cout<<"header size: "<< header.size() <<endl;
			const unsigned int dvalue	=0;
			unsigned int col_nchr		=dvalue; // initialization of columns
			unsigned int col_rsid		=dvalue;
			unsigned int col_snpid		=dvalue;
			unsigned int col_bp			=dvalue;
			unsigned int col_cm_pos		=dvalue;
			unsigned int nCols			=dvalue;
			unsigned int col_allele1	=dvalue;
			unsigned int col_allele2	=dvalue;
			string 	 temp_col_value			="";
			unsigned int temp_indx			=0;
			static 	int  no_of_update_snps	=0;
			//bool	temp_bool				=false;
			nCols=header.size();
			vector<string>sort_header=header;
			vector<string>::iterator it;
			sort(sort_header.begin(),sort_header.end());
				//rsid
			bool found0	= binary_search(sort_header.begin(), sort_header.end(),"rs");
			bool found0_1	= binary_search(sort_header.begin(), sort_header.end(),"id");
			bool found0_2	= binary_search(sort_header.begin(), sort_header.end(),"ID");
			bool found	= binary_search(sort_header.begin(), sort_header.end(),"rsid");
			bool found1	= binary_search(sort_header.begin(), sort_header.end(),"snpid");
			bool found2	= binary_search(sort_header.begin(), sort_header.end(),"rsID");
			bool found3	= binary_search(sort_header.begin(), sort_header.end(),"snpID");
			//cout << found2 <<endl;
				if(found0||found0_1||found0_2||found|| found1||found2||found3)
				{
					if(found0)
					{
						//rsid
						it=find(header.begin(),header.end(),"rs");
						col_rsid= (int)(distance(header.begin(),it)+1);
					}else 
					if(found0_1)
					{
						//rsid
						it=find(header.begin(),header.end(),"id");
						col_rsid= (int)(distance(header.begin(),it)+1);
					}else
					if(found0_2)
					{
						//rsid
						it=find(header.begin(),header.end(),"ID");
						col_rsid= (int)(distance(header.begin(),it)+1);
					}else
					if(found)

					{
						//rsid
						it=find(header.begin(),header.end(),"rsid");
						col_rsid= (int)(distance(header.begin(),it)+1);
					}
					else if(found1)
					{
						it=find(header.begin(),header.end(),"snpid");
						col_snpid= (int)(distance(header.begin(),it)+1);
					}
					else if(found2)
					{
						it=find(header.begin(),header.end(),"rsID");
						col_rsid= (int)(distance(header.begin(),it)+1);
					}
					else if(found3)
					{
						it=find(header.begin(),header.end(),"snpID");
						col_snpid= (int)(distance(header.begin(),it)+1);
					}

				}
				
				//else
				//{
				//	string msg="your extra map file \""+fileName+"\"  neither  possesses acolumn containing rsid nor a column containing snpid.\n";
				//	error(msg);
				//}
				found= binary_search(sort_header.begin(), sort_header.end(),"nchr");
				//cout <<" nchr found?? " <<found <<endl; // debug nchr found ??
				if(found)
				{	// nchr
					it=find(header.begin(),header.end(),"nchr");
					col_nchr= (int)(distance(header.begin(),it)+1);
				}
				found= binary_search(sort_header.begin(), sort_header.end(),"bp");
				found1= binary_search(sort_header.begin(), sort_header.end(),"position");
				//cout << "position "<< found1 <<endl; //debug 
				if(found||found1)
				{
				// bp
					if(found)
					{
						// bp
						it=find(header.begin(),header.end(),"bp");
						col_bp= (int)(distance(header.begin(),it)+1);
					}
					else if(found1)
					{
						it=find(header.begin(),header.end(),"position");
						col_bp= (int)(distance(header.begin(),it)+1);
					}
				}
				found= binary_search(sort_header.begin(), sort_header.end(),"cm_pos");
				if(found)
				{	// cm_pos
					it=find(header.begin(),header.end(),"cm_pos");
					col_cm_pos= (int)(distance(header.begin(),it)+1);
				}
				//allele1
				 found= binary_search(sort_header.begin(), sort_header.end(),"allele1");
				 found1= binary_search(sort_header.begin(), sort_header.end(),"a0");
				 found2= binary_search(sort_header.begin(), sort_header.end(),"Allele1");
				 if(found||found1||found2)
				 {
					 if(found)
					{	// allele1
								it=find(header.begin(),header.end(),"allele1");
								col_allele1= (int)(distance(header.begin(),it)+1);
					}
					else if(found1)
					 {	// allele1
						 it=find(header.begin(),header.end(),"a0");
						 col_allele1= (int)(distance(header.begin(),it)+1);
					 }
					else  if(found2)
					 {	// allele1
						 it=find(header.begin(),header.end(),"Allele1");
						 col_allele1= (int)(distance(header.begin(),it)+1);
					 }
				 }
							//allele2
				 found= binary_search(sort_header.begin(), sort_header.end(),"allele2");
				 found1= binary_search(sort_header.begin(), sort_header.end(),"a1");
				 found2= binary_search(sort_header.begin(), sort_header.end(),"Allele2");
				 if(found||found1||found2)
				 {
					 if(found)
					 {	// allele2
						 it=find(header.begin(),header.end(),"allele2");
						 col_allele2= (int)(distance(header.begin(),it)+1);
					 }
					 else if(found1)
					 {	// allele2
						 it=find(header.begin(),header.end(),"a1");
						 col_allele2= (int)(distance(header.begin(),it)+1);
					 }
					else if(found2)
					 {	// allele2
						 it=find(header.begin(),header.end(),"Allele2");
						 col_allele2= (int)(distance(header.begin(),it)+1);
					 }

				 }


				 if((col_allele1!=dvalue&& col_allele2==dvalue)|| (col_allele1==dvalue &&col_allele2!=dvalue))
				 {
					 string msg= "It is not allowed to give only one allele information in snpinfo file \""+\
								 fileName+"\". \n";
					 error(msg);

				 }
					
				//displaying for test
					//cout << "col_rsId: "<< col_rsid << " "<< "col_snpId: "<< col_snpid << " "<< "col_nchr: "<< col_nchr << " "<< "col_bp: "<< col_bp << " " << "col_cm_pos: "<< col_cm_pos << " ";
				
			
			
		
			vector<string> tokens;
			vector<CBSNP*>::const_iterator it_pCBSNP=genInfo.begin();
			//-----------------------------------------------------------//
			vector<string> temp_vec_rsid;
			vector<string> temp_vec_snpid;
			vector<long int>	temp_vec_bp;
			vector<string>::iterator temp_it_string =temp_vec_rsid.end();

			for(unsigned int i=0;i<genInfo.size();++i)
			{
				CBSNP* it_pCBSNP =genInfo[i];
				temp_vec_rsid.push_back(it_pCBSNP->rsId);
				temp_vec_snpid.push_back(it_pCBSNP->snpName);
				temp_vec_bp.push_back(it_pCBSNP->bp);

			}
			//-------------------------------------------------------------------//
			r_msg="\t**->Updating SNP information from \""+fileName+ "\": \n";
			printLIN(r_msg);
			while(getline(file,line, '\n')) // until eof hit
			{
			//---------------------//
				++extraMapFileLineZaehler;
				//cout <<"line no: "<< extraMapFileLineZaehler<<endl; //debug
				if(line.length()==0)
				continue;
				if(line[0]=='#')
				continue;
				//if(!((*it_pCBSNP)->quality))
				//	continue;
				tokens.resize(0);
				tokens.clear();
				stringstream line_parser1(line);
				//cout <<"\n line_parser1.good(): "<<line_parser1.good() <<endl; // only for test
				if(line_parser1.good())
				{
					while(line_parser1>>buffer)
					{
						tokens.push_back(buffer);
						//cout << buffer << " "; // only for test
					}
					if(tokens.size()!=nCols)
					{
						string msg= "Your extra map file \""+fileName+"\"  has "+change_int_into_string(nCols)+\
								" columns in header line but in the line "+change_int_into_string(extraMapFileLineZaehler)+\
								" there are "+change_int_into_string(tokens.size())+" columns.\n";
						error(msg);
					}
					// add all values of tokens to the vectors constructed for extra snp info file.
					//cout <<"col_rsid: "<< col_rsid<<" "<<(col_rsid!=dvalue) <<endl ;
					 if(col_rsid!=dvalue)
						 sifo_rsid.push_back(tokens[col_rsid-1]);
					 if(col_snpid!=dvalue)
						 sifo_snpid.push_back(tokens[col_snpid-1]);
					 if(col_nchr!=dvalue)
						 sifo_nchr.push_back(tokens[col_nchr-1]);
					 if(col_bp!=dvalue)
						 sifo_bp.push_back(atoi((tokens[col_bp-1]).c_str()));
					 if(col_cm_pos!=dvalue)
						 sifo_cm_pos.push_back(atof(tokens[col_cm_pos-1].c_str()));
					 if( col_allele1!=dvalue)
						 sifo_allele1.push_back(tokens[col_allele1-1]);
					 if( col_allele2!=dvalue)
						 sifo_allele2.push_back(tokens[col_allele2-1]);
					 //cout <<sifo_rsid.size()<< " " <<sifo_snpid.size() <<" " << " "<<sifo_nchr.size()<<" "<< sifo_bp.size()<<" "<<sifo_cm_pos.size()<<" "<<sifo_allele1.size()<<" "<<sifo_allele2.size();// debug
					 //cout << "\n"; // debug
					 //rsid
					// first  check if rsid or snpid is given. This is a necessary command
					if(col_rsid!=dvalue)
					{
						temp_col_value	=tokens[col_rsid-1];
						temp_it_string	= find(temp_vec_rsid.begin(), temp_vec_rsid.end(),temp_col_value);
						//cout << *temp_it_string<<endl; // debug
						if(temp_it_string!=temp_vec_rsid.end())
						{
							temp_indx	  	=(int) (temp_it_string-temp_vec_rsid.begin());
							it_pCBSNP		=genInfo.begin()+temp_indx;
						}
						else
						{
							temp_col_value="";
							continue; // start the while loop from the begning.
						// check if rs id match with other or not
						}
					}
					else if(col_snpid!=dvalue)
					{
						temp_col_value	=tokens[col_snpid-1];
						temp_it_string	= find(temp_vec_snpid.begin(), temp_vec_snpid.end(),temp_col_value);
						//cout << *temp_it_string<<endl; // debug
						if(temp_it_string!=temp_vec_snpid.end())
						{
							temp_indx	  	=(int)(temp_it_string-temp_vec_snpid.begin());
							it_pCBSNP		=genInfo.begin()+temp_indx;
						}
						else
						{
							continue; // start the while loop again
							temp_col_value="";
						}

					}
					else
					{
						string msg= "To update SNP information  your file snpinfo file\""+fileName+\
								"\" should contain at least a column either containing rs IDs or SNP IDs \n.";
						error(msg);
					}
					//
					//----------------------------------------------------//
					//Now check if the column of base pair position is given and if the base
					//pair positions match with each others
					// check bp first
					//cout << "test1"<<endl; //debug
					//cout << boolalpha << (it_pCBSNP!=genInfo.end()) <<endl; //debug
					if( (it_pCBSNP!=genInfo.end()) )
					{
						no_of_update_snps++;
						//bp 
						if(col_bp!=dvalue)
						{
							temp_col_value= tokens[col_bp-1];
							if((*it_pCBSNP)->bp!=(-1) && (*it_pCBSNP)->bp!=0)
							{
								if((*it_pCBSNP)->bp !=atoi(temp_col_value.c_str()))
								{
									string msg= " Base pair position \""+change_int_into_string((*it_pCBSNP)->bp)+
											"\" in the original file  does not match with base pair position \""+\
											tokens[col_bp-1]+" given in snpinfo file "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra map file\""+fileName+\
											"\" \n.";
									error(msg);
								}
							}
							else
							(*it_pCBSNP)->bp =atoi(temp_col_value.c_str());
							temp_col_value="";

						}
						//cout << "test2"<<endl; //debug
						// now rsId
						if(col_rsid!=dvalue)
						{
							if((*it_pCBSNP)->rsId!="")
							{
								if((*it_pCBSNP)->rsId !=tokens[col_rsid-1])
								{
									string msg= "rs id  \""+(*it_pCBSNP)->rsId+"\" in original file  does not match with rs id \""+\
											tokens[col_rsid-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra map file\""+fileName+\
												"\" \n.";
										error(msg);

								}

							}
							else
								(*it_pCBSNP)->rsId =tokens[col_rsid-1];
						}
						//
						//cout << "test3"<<endl; //debug
						//snpid
						if(col_snpid!=dvalue)
						{
							if(((*it_pCBSNP)->snpName!="")&&((*it_pCBSNP)->snpName!="---") &&((*it_pCBSNP)->snpName!="."))
							{
								if((*it_pCBSNP)->snpName !=tokens[col_snpid-1])
								{
									if((*it_pCBSNP)->snpName !=tokens[col_rsid-1])
									{
										string msg= " snp id  \""+(*it_pCBSNP)->snpName+"\" in datfile  does not match with snp id \""+\
												tokens[col_snpid-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra map file\""+fileName+\
												"\" \n.";
										error(msg);
									}
									else
										(*it_pCBSNP)->snpName =tokens[col_snpid-1];
								}
							}
							else
								(*it_pCBSNP)->snpName =tokens[col_snpid-1];
						}
						//nchr
						//cout << "test4"<<endl; //debug
						if(col_nchr!=dvalue)
						{
							if(((*it_pCBSNP)->nchr!="0")&&((*it_pCBSNP)->nchr!="."))
							{
								if((*it_pCBSNP)->nchr !=tokens[col_nchr-1])
								{
									string msg= " chromosome \""+(*it_pCBSNP)->nchr+"\" in datfile  does not match with chromsome \""+\
											tokens[col_nchr-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra map file\""+fileName+\
											"\" \n.";
									error(msg);
								}
							}
							else
								(*it_pCBSNP)->nchr =tokens[col_nchr-1];
						}
								//cm_pos
						if( col_cm_pos!=dvalue)
						{
							if((*it_pCBSNP)->cm_pos!=-1.0 && (*it_pCBSNP)->cm_pos!=0.0)
							{
								if((*it_pCBSNP)->cm_pos !=atof(tokens[col_cm_pos-1].c_str()))
								{
									string msg= " Recombination rate \""+change_double_into_string((*it_pCBSNP)->cm_pos)+\
											"\" in datfile  does not match with the recombination rate \""+\
											tokens[col_cm_pos-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra map file\""+fileName+\
											"\" \n.";
									error(msg);
								}
							}
							else
								(*it_pCBSNP)->cm_pos =atof(tokens[col_cm_pos-1].c_str());
						}
						//cout << "test5"<<endl; //debug
								//allele1
							//cout << col_allele1 << " "<<col_allele2 <<endl;	//test 
						if(tokens[col_allele1-1] ==tokens[col_allele2-1]){
							tokens[col_allele2-1]="";
						}
						if(col_allele1!=dvalue && col_allele2!=dvalue)
						{
							//--------------------------------------------------------------------------------------------------------------------------------------------------------------// 
							bool cond_1=(*it_pCBSNP)->allele1!=""; // && (*it_pCBSNP)->allele1!="0";
							//bool cond_2=(*it_pCBSNP)->allele2!="" && (*it_pCBSNP)->allele2!="0";
							//cout << "tmpb2: "<< tmpb2<<endl; // test 
							// cout << cond_1 << endl; // test 	
							 if(cond_1)
							 {
								//1
								bool tmpb1=((*it_pCBSNP)->allele1 ==tokens[col_allele1-1]) && ((*it_pCBSNP)->allele2 ==tokens[col_allele2-1]);
								bool tmpb1_1=((*it_pCBSNP)->allele1 ==tokens[col_allele1-1]) && ((tokens[col_allele2-1]=="" ) ||(tokens[col_allele2-1]=="." ) );
								bool tmpb1_2=(tokens[col_allele1-1]=="") && ((*it_pCBSNP)->allele2 ==tokens[col_allele2-1]);
								bool tmpb2=((*it_pCBSNP)->allele1 ==tokens[col_allele2-1]) && ((*it_pCBSNP)->allele2 ==tokens[col_allele1-1]);

								//3
								bool tmpb3=((*it_pCBSNP)->allele1 ==tokens[col_allele1-1])&& ((*it_pCBSNP)->allele2 =="");//||(*it_pCBSNP)->allele2 =="0");
								//4
								bool tmpb4=((*it_pCBSNP)->allele1 ==tokens[col_allele2-1]) && (((*it_pCBSNP)->allele2 ==""));// || ((*it_pCBSNP)->allele2 =="0"));
								// this condition does not come!! 
								bool tmpb5=((*it_pCBSNP)->allele2 ==tokens[col_allele1-1]) && (((*it_pCBSNP)->allele1 =="") || ((*it_pCBSNP)->allele1 =="0") || ((*it_pCBSNP)->allele1 =="."));
								bool tmpb6=((*it_pCBSNP)->allele2 ==tokens[col_allele2-1]) && (((*it_pCBSNP)->allele1 =="") || ((*it_pCBSNP)->allele1 =="0") || ((*it_pCBSNP)->allele1 =="."));
								// this condition does not also come!! 
								//bool tmpb7=((*it_pCBSNP)->allele1 ==tokens[col_allele2-1]) && (( tokens[col_allele2-1] =="" ) || (tokens[col_allele2-1] =="0" ) );	
								// this condition is ok nothing to change 
								bool tmpb8=((*it_pCBSNP)->allele1 ==tokens[col_allele2-1]) && ((tokens[col_allele1-1] =="" )	|| (tokens[col_allele1-1] ==".") ) ;// || (tokens[col_allele1-1] =="0" ) );
								bool tmpb9=((*it_pCBSNP)->allele2 ==tokens[col_allele1-1]) && ((tokens[col_allele2-1] =="" )	|| (tokens[col_allele1-1] ==".") );// || (tokens[col_allele2-1] =="0" ) );
								
								//cout << "tmpb2: "<< tmpb2<<endl; // test 
								if(tmpb1|| tmpb1_1|| tmpb1_2)
								{
									//do nothing 
								}
								else if(tmpb2 ||tmpb8 || tmpb9 || tmpb5)
								{
									// chnage whole allele 1 and allele 2 with each other . 
									//cout << "change of allele! allele 1 is allele2 and allele 2 will be allele1 \n"; //ONly for test
										for(unsigned int i=0; i<(*it_pCBSNP)->geno1.size();++i)
										{
											//homo1->homo2
											if(!((*it_pCBSNP)->geno1[i])&& !((*it_pCBSNP)->geno2[i]))
											{
												(*it_pCBSNP)->geno1[i]=true;
												(*it_pCBSNP)->geno2[i]=true;
											}
											//homo2 ->homo1
											else if(((*it_pCBSNP)->geno1[i])&& ((*it_pCBSNP)->geno2[i]))
											{
												(*it_pCBSNP)->geno1[i]=false;
												(*it_pCBSNP)->geno2[i]=false;
											}
											//geno1 =F geno2 =T => heterozygote . here only the order should be changed. 
											else if(!((*it_pCBSNP)->geno1[i])&& ((*it_pCBSNP)->geno2[i]))
											{
												if((*it_pCBSNP)->aOrder[i])
													(*it_pCBSNP)->aOrder[i]=false;
												if(!((*it_pCBSNP)->aOrder[i]))
													(*it_pCBSNP)->aOrder[i]=true;
											}
											// change probability genotypes also
											if((*it_pCBSNP)->pgeno1.size()>0 && (*it_pCBSNP)->pgeno3.size()>0)
											{
												//cout << no_of_persons <<endl;
												for(unsigned int j=0;j<(*it_pCBSNP)->pgeno1.size(); ++j)
													swap((*it_pCBSNP)->pgeno1[j],(*it_pCBSNP)->pgeno3[j]);
														 //swap will exchange prob(AA) and prob(BB).prob(AB) remains the same
											}

										}
										// change the allele also
										string _tmp_a_1=(*it_pCBSNP)->allele1;
										string _tmp_a_2=(*it_pCBSNP)->allele2;
										if(tmpb2 || tmpb5)
										{
											(*it_pCBSNP)->allele1 =tokens[col_allele1-1];
											(*it_pCBSNP)->allele2 =tokens[col_allele2-1];
										}else
										if(tmpb8){

												string _tmp_a_2=(*it_pCBSNP)->allele2;
												(*it_pCBSNP)->allele2 =tokens[col_allele2-1];
												(*it_pCBSNP)->allele1 =_tmp_a_2;

										}
										else if(tmpb9)
										{
											string _tmp_a_1=(*it_pCBSNP)->allele1;
											(*it_pCBSNP)->allele1 =tokens[col_allele1-1];
											(*it_pCBSNP)->allele2 =_tmp_a_1;
										}
										//only for test
													
								}else
									if(tmpb6)
								{
									//this condition in which saved allele1 is empty but allele2 major allele match with the alelle 2 of
									//snpinfo file so just assign allele1 to allele1 of snpinfo nothing more than that
									(*it_pCBSNP)->allele1 =tokens[col_allele1-1];
								}
								else if(tmpb3)
								{
									// it assumes that there is no allele 2 in the data so check if there is not the case and 
									// if this is the case throw error . if not just assign allele2 =second allele2 
									for(unsigned int i=0; i<(*it_pCBSNP)->geno1.size();++i)
									{
											//if geno1=T and geno2=T or geno1 =F and geno2=T  (i.e. case homo2 and hetero 
											if ((((*it_pCBSNP)->geno1[i])&& ((*it_pCBSNP)->geno2[i] ))|| (!((*it_pCBSNP)->geno1[i])&& ((*it_pCBSNP)->geno2[i])))
											{
												cout <<"Allele in original files: " << 	(*it_pCBSNP)->allele1 << " " <<(*it_pCBSNP)->allele2 <<endl;
												cout <<"Allele in snpinfo files: " << 	tokens[col_allele1-1] << " " <<tokens[col_allele2-1] <<endl;
												string msg= " Allele2 \""+(*it_pCBSNP)->allele2+\
												"\" in datfile  is empty but there are Heterozygous 2  types and heterozygoous genotypes exisit in SNP  \""+\
												(*it_pCBSNP)->rsId +"\" \n";
												error(msg);
												
											} 
											
									}
									(*it_pCBSNP)->allele2 =tokens[col_allele2-1];
								}
								
								else if(tmpb4)
								{
								
									//if geno1=T and geno2=T or geno1 =F and geno2=T  (i.e. case homo2 and hetero 
										for(unsigned int i=0; i<(*it_pCBSNP)->geno1.size();++i)
										{	
											//if geno1=T and geno2=T or geno1 =F and geno2=T  (i.e. case homo2 and hetero 
											if ((((*it_pCBSNP)->geno1[i])&& ((*it_pCBSNP)->geno2[i] ) )|| (!((*it_pCBSNP)->geno1[i])&& ((*it_pCBSNP)->geno2[i])))
											{
												cout <<"Allele in original files: " << 	(*it_pCBSNP)->allele1 << " " <<(*it_pCBSNP)->allele2 <<endl;
												cout <<"Allele in snpinfo files: " << 	tokens[col_allele1-1] << " " <<tokens[col_allele2-1] <<endl;
												string msg= " Allele2 \""+(*it_pCBSNP)->allele2+\
												"\" in datfile  is empty but there are Homozygous 2   and heterozygoous genotypes  in SNP  \""+(*it_pCBSNP)->rsId +"\" \n";
												error(msg);
												
											}
											else
											{
												//homo1->homo2
												if(!((*it_pCBSNP)->geno1[i])&& !((*it_pCBSNP)->geno2[i]))
												{
													(*it_pCBSNP)->geno1[i]=true;
													(*it_pCBSNP)->geno2[i]=true;
												}
													if((*it_pCBSNP)->pgeno1.size()>0 && (*it_pCBSNP)->pgeno3.size()>0)
												{
													//cout << no_of_persons <<endl;
													for(unsigned int j=0;j<(*it_pCBSNP)->pgeno1.size(); ++j)
														swap((*it_pCBSNP)->pgeno1[j],(*it_pCBSNP)->pgeno3[j]);
															 //swap will exchange prob(AA) and prob(BB).prob(AB) remains the same
												}

											}
											
											
										}
										(*it_pCBSNP)->allele1 =tokens[col_allele1-1];
										(*it_pCBSNP)->allele2 =tokens[col_allele2-1];	
								}
								else {
								
									cout <<"Allele in original files: " << 	(*it_pCBSNP)->allele1 << " " <<(*it_pCBSNP)->allele1 <<endl;
									cout <<"Allele in snpinfo files: " << 	tokens[col_allele1-1] << " " <<tokens[col_allele2-1] <<endl;
									string msg= " Allele1 \""+(*it_pCBSNP)->allele1+ "\" or  Allele2 \""+\
									(*it_pCBSNP)->allele2+ "\" in datfile  does not match with the allele1 \""+\
									tokens[col_allele1-1]+"\" allele2 \""+tokens[col_allele2-1]+\
									"\"	given in "+change_int_into_string(extraMapFileLineZaehler)+\
									"th line of  SNP info file\""+fileName+"\" \n.";
									error(msg);
								
								
								}
								//
									 
							 
							 }
							 //
							 else 
							 {
								bool tmp_b1=((*it_pCBSNP)->allele2 ==tokens[col_allele1-1]) && (((*it_pCBSNP)->allele1 ==""));// || ((*it_pCBSNP)->allele1 =="0"));
								bool tmp_b2=((*it_pCBSNP)->allele1 =="") && ((*it_pCBSNP)->allele2 =="");// ||(*it_pCBSNP)->allele2 =="0")  ;

							   if(tmp_b1)
							   {
									//if geno1=F and geno2=F or geno1 =F and geno2=T  (i.e. case homo1 and hetero 
									for(unsigned int i=0; i<(*it_pCBSNP)->geno1.size();++i)
									{
											//if geno1=T and geno2=T or geno1 =F and geno2=T  (i.e. case homo2 and hetero 
											if ( (!((*it_pCBSNP)->geno1[i])&& !((*it_pCBSNP)->geno2[i] ) ) || (!((*it_pCBSNP)->geno1[i])&& ((*it_pCBSNP)->geno2[i])))
											{
												cout <<"Allele in original files: " << 	(*it_pCBSNP)->allele1 << " " <<(*it_pCBSNP)->allele2 <<endl;
												cout <<"Allele in snpinfo files: " << 	tokens[col_allele1-1] << " " <<tokens[col_allele2-1] <<endl;
												string msg= " Allele1 \""+(*it_pCBSNP)->allele1+\
												"\" in datfile  is empty but there are Homozygous 1   and heterozygoous  genotypes  in SNP  \""+(*it_pCBSNP)->rsId +"\" \n";
												error(msg);
												
											}
											else
											{
												//homo1->homo2
												if(((*it_pCBSNP)->geno1[i])&& ((*it_pCBSNP)->geno2[i]))
												{
													(*it_pCBSNP)->geno1[i]=false;
													(*it_pCBSNP)->geno2[i]=false;
												}
													if((*it_pCBSNP)->pgeno1.size()>0 && (*it_pCBSNP)->pgeno3.size()>0)
												{
													//cout << no_of_persons <<endl;
													for(unsigned int j=0;j<(*it_pCBSNP)->pgeno1.size(); ++j)
														swap((*it_pCBSNP)->pgeno1[j],(*it_pCBSNP)->pgeno3[j]);
														//swap will exchange prob(AA) and prob(BB).prob(AB) remains the same
												}

											}
											
											
									}
									(*it_pCBSNP)->allele1 =tokens[col_allele1-1];
									(*it_pCBSNP)->allele2 =tokens[col_allele2-1];	   
							   
							   }
							   else if(tmp_b2)
							   {
									
								(*it_pCBSNP)->allele1 =tokens[col_allele1-1];
								(*it_pCBSNP)->allele2 =tokens[col_allele2-1];	   
								
							   }
							   else
							   {
									cout <<"rsid: "<< (*it_pCBSNP)->rsId<< endl;
									cout <<"Allele in original files: " << 	(*it_pCBSNP)->allele1 << " " <<(*it_pCBSNP)->allele1 <<endl;
									cout <<"Allele in snpinfo files: " << 	tokens[col_allele1-1] << " " <<tokens[col_allele2-1] <<endl;
									string msg= " Allele1 \""+(*it_pCBSNP)->allele1+ "\" or  Allele2 \""+(*it_pCBSNP)->allele2+\
									"\" in datfile  does not match with the allele1 \""+tokens[col_allele1-1]+\
									"\"  and allele2 \""+tokens[col_allele2-1]+"\"	given in "+\
									change_int_into_string(extraMapFileLineZaehler)+"th line of  SNP info file\""+fileName+"\" \n.";
									error(msg);
								
							   
							   }
									
							}
							//--------------------------------------------------------------------------------------------------------------------------------------------------------------// 
							
						}

					//-----------------------------------------------------------------//
				  }// end of it_pCBSNP!=genInfo.end()
					//else continue;

				} // end of if line_praser1.good();
				else
				{
					printLIN( "Line parser not good while reading extra snp info file.\n");
					break;
				}

				it_pCBSNP=genInfo.end();

			} // end of while loop

			
		// for last test 
		/*
		cout << "After update:\n";
		for(unsigned int i=0;i<genInfo.size();++i)
		 {
			cout << genInfo[i]->rsId <<" \t" << genInfo[i]->allele1 << " "<< genInfo[i]->allele2 << " \n";
		 
		 }
		*/	
		//-----------------------------------------------------------------------//
		// here is for impute hapmap file 
		//-----------------------------------------------------------------------//
		//cout << temp_vec_rsid.size()<< " "<< sifo_rsid.size() <<endl; 
		found0	= (temp_vec_rsid.size()==sifo_rsid.size());
		found1 	=(temp_vec_snpid.size()==sifo_snpid.size()); 	
		found2  =(temp_vec_bp.size()==sifo_bp.size());
		found3 	=(found0||found1 );
		//cout << found2 <<endl; //debug 
		if(no_of_update_snps==0)
		{
			if( found3)
			{
				//cout<< "impute hapmap \n"; // debug 
				for(unsigned int i=0;i<genInfo.size();++i)
				{
					CBSNP* it_pCBSNP =genInfo[i];
					if(it_pCBSNP->rsId =="" && found0)
						it_pCBSNP->rsId=sifo_rsid[i];
					if((it_pCBSNP->snpName =="") && found1)
						it_pCBSNP->snpName	=sifo_snpid[i];
					if((it_pCBSNP->nchr =="0" ) && (sifo_nchr.size()==genInfo.size()))
						it_pCBSNP->nchr=sifo_nchr[i];
					if((it_pCBSNP->bp ==(-1)) && found2)
						it_pCBSNP->bp=sifo_bp[i];
					if(it_pCBSNP->cm_pos==(-1.0)&& (sifo_cm_pos.size()==genInfo.size()))
						 it_pCBSNP->cm_pos=sifo_cm_pos[i];
					if(it_pCBSNP->allele1==""&& (sifo_allele1.size()==genInfo.size()))
						 it_pCBSNP->allele1=sifo_allele1[i];
					if(it_pCBSNP->allele2==""&& (sifo_allele2.size()==genInfo.size()))
						 it_pCBSNP->allele2=sifo_allele2[i];
					++no_of_update_snps;
				}
			}
		}	
		//-----------------------------------------------------------------------//
		//end of impute hapmap file 
		//-----------------------------------------------------------------------//
		printLIN("\t Updated  elements:\n");
				if(col_snpid!=dvalue ||col_rsid!=dvalue ||col_nchr!=dvalue || col_bp!=dvalue||col_cm_pos!=dvalue ||col_allele1!=dvalue||col_allele2!=dvalue)
				{
					if(col_rsid!=dvalue)
					{
						printLIN("\t rsid \n");
					}
					if(col_snpid!=dvalue)
					{
						printLIN("\t snpid \n");
					}
					if(col_nchr!=dvalue)
					{
						printLIN("\t chromosome number\n");
					}
					if(col_bp!=dvalue)
					{
						printLIN("\t basepair position\n");
					}
					if(col_cm_pos!=dvalue)
					{
						printLIN("\t basepair's centimorgen value \n");
					}
					if(col_allele1!=dvalue && col_allele2!=dvalue)
					{
						printLIN("\t allele information \n");
					}
					printLIN("\t Information for all total "+ change_int_into_string(no_of_update_snps)+" SNPs have been updated.\n");
					
				}else
				{
					printLIN("\t Nothing:\n");
					printLIN("\t Information for all total "+ change_int_into_string(no_of_update_snps)+" SNPs have been read but not updated.\n \t So please check if you have given correct names in the header of your snpinfo file.\n");
				
				}
		 } // closing of if header.size()!=1 
		
		file.close();

return ;}
void CBSNP::lese_exclude_snps_file(const string& fileName, const vector<CBSNP*> & genInfo){
	printLIN("\t**->Excluding SNPs mentioned in file \""+fileName+"\": \n");
	string line ;
	string buffer;
	rFile file;
	file.close(); file.clear();
	file.open(fileName.c_str(),ios::in);
	static int _tmp_nLine=0;
	static int _tmp_nFound=0;
	if(!file)
	{
		file.setstate(ios_base::failbit);
		file.clear();
		error("file \""+fileName+"\" does not exist or could not be opend. Please make sure that you have given the correct filepath and file name.\n");
	}
	vector<string> tmp_vec_rsid;
	vector<string> tmp_vec_snpid;
	CBSNP* pCBSNP =genInfo[0];
	vector<string>::iterator it;
	int snp_index =0;

	bool _found =false;
	for(unsigned int j=0;j<genInfo.size();++j)
	{
		pCBSNP =genInfo[j];
		tmp_vec_rsid.push_back(pCBSNP->rsId);
		tmp_vec_snpid.push_back(pCBSNP->snpName);
	}
	//getline(file,line);
	//cout <<line <<endl;
	while(getline(file,line))
	{
		if(line.size()>100)
		{
			file.clear();file.close();
			error("Are you sure you are feeding a list of SNPs with \"--exclude\" command?\n");
		}
		if(line.empty())
			continue;
		if(line[0]=='#')
			continue;
		//if(file.eofbit)
		//	break;
		//cout<< line <<endl;
		stringstream line_parser(line);
		if(line_parser.good())
			line_parser>>buffer;
		//cout << buffer <<endl;
		it =find(tmp_vec_rsid.begin(),tmp_vec_rsid.end(),buffer);
		if(it!=tmp_vec_rsid.end())
		{
				_found =true;
				snp_index =(int)(it-tmp_vec_rsid.begin());
				genInfo[snp_index]->quality =false;
				++_tmp_nLine;
		}
		else
		{
			it =find(tmp_vec_snpid.begin(),tmp_vec_snpid.end(),buffer);
			if(it!=tmp_vec_snpid.end())
			{
				_found =true;
				snp_index =(int)(it-tmp_vec_snpid.begin());
				genInfo[snp_index]->quality =false;
				++_tmp_nLine;
			}else if(buffer!="")
			{
				++_tmp_nFound;
				if(_tmp_nFound<20)
				 printLIN("\t...SNP named \""+ buffer+"\" not found in the data \n");

			}


		}
		_found=false;
		line="";
	}//end of while getline loop;
	 if(_tmp_nFound>20)
		 printLIN("\tThere are more than 20 SNPs not found in the data\n");
	printLIN("\t**->Total number of excluded SNPs: "+change_int_into_string(_tmp_nLine)+".\n");
	//closing file
	file.clear(); file.close();
return ;}
//xxxxxxxxxxxxxxxxxxxxx
void CBSNP::write_snpinfoFile( const vector<CBSNP*>& genInfo,const string& outFileName)
 {
 	//----------------------------------#--------------------------------------#
  	// writing extra mapinfo file
  	 //----------------------------------#--------------------------------------#
 		ofstream ofs;
  	 	 string extraMapInfoName=outFileName+"_snpinfo.txt";
  	 	 	 ofs.clear();ofs.close();
  	 	  	ofs.open(extraMapInfoName.c_str(),ios::out);
  	 	 	if(!ofs)
  	 	 		error("out file "+extraMapInfoName+" does not exits.\n");
  	 	 	//if(genInfo[1]->nchr!="")
  	 	 	  	ofs<< "nchr" << " ";
  	 	 	 //if(genInfo[1]->snpName!="")
  	 	 	 		ofs<< "snpid" << " ";
  	 	 	 //if(genInfo[1]->rsId!="")
  	 	 	 	 	ofs<< "rsid" << " ";
  	 	 	 //if(genInfo[1]->bp!=0)
  	 	 	 		ofs<<"bp"<< " ";
  	 	 	 //if(genInfo[1]->cm_pos!=0)
  	 	 	  	ofs<<"cm_pos"<< " ";
  	 	 	 	ofs<<"allele1"<< " ";
  	 	 	 	ofs<<"allele2"<< "\n";
  	 	 for(unsigned int i=0; i<genInfo.size();++i)
  	 	 	{
  	 		 	 CBSNP* pCBSNP	= genInfo[i];
  	 		 	 if(pCBSNP->quality)
  	 		 	 {
					//if(genInfo[i]->nchr!="")
						ofs<<pCBSNP->nchr<<" ";
					//if(pCBSNP->snpName!="")
					 if(pCBSNP->snpName!="") ofs <<pCBSNP->snpName<< " ";
					  else if(pCBSNP->rsId!="") ofs <<pCBSNP->rsId<< " ";
					//if(pCBSNP->rsId!="")
					 if(pCBSNP->rsId!="") ofs << pCBSNP->rsId<< " ";
						 else if(pCBSNP->snpName!="") ofs <<pCBSNP->snpName<<" ";
					if(pCBSNP->bp!=-1)
						ofs<<pCBSNP->bp<< " ";
					else
						ofs<<0<<" ";
					if(pCBSNP->cm_pos!=(-1.0))
						ofs<<pCBSNP->cm_pos<< " ";
					else ofs <<0.0 << " ";
					//
					if(pCBSNP->allele1!="")
						ofs<< pCBSNP->allele1<< " ";
					else if(pCBSNP->allele2!="")
						ofs<< pCBSNP->allele2<< " ";
					else ofs <<"0"<< " ";
					//
					if(pCBSNP->allele2!="")
						ofs<< pCBSNP->allele2<< "\n";
					else if(pCBSNP->allele1!="")
						ofs<< pCBSNP->allele1<< "\n";
					else ofs <<"0"<< "\n";

				} //end of if quality loop
			}

  	 	 	ofs.clear();ofs.close();
 	 	printLIN("\t**->One extra \"snpinfo.txt\" file has been written out and saved as \""+\
  	 	 			extraMapInfoName+"\". \n");


 return ;}
// writee bimbam snp position file
void CBSNP::write_snplist_file( const vector<CBSNP*>& genInfo,const string& outFileName)
 {
 	//----------------------------------#--------------------------------------#
  	// writing snp list
  	 //----------------------------------#--------------------------------------#
	ofstream ofs;
	string extraMapInfoName=outFileName+"_snplist.txt";
	if(CBSNP::_out_format=="mach"||CBSNP::_out_format=="minimac")
		extraMapInfoName=outFileName+".snps";
	ofs.clear();ofs.close();
	ofs.open(extraMapInfoName.c_str(),ios::out);
	if(!ofs)
  	 	 		error("out file "+extraMapInfoName+" does not exits.\n");
	CBSNP* pCBSNP	= genInfo[0];
	for(unsigned int i=0; i<genInfo.size();++i)
	{
		pCBSNP	= genInfo[i];
		if(pCBSNP->quality)
		{
			if(pCBSNP->rsId!="") ofs << pCBSNP->rsId;
			else if(pCBSNP->snpName!="") ofs <<pCBSNP->snpName;
			ofs<<"\n"; //end of line

		}
	}

  	 	 	ofs.clear();ofs.close();
 	 	printLIN("\t **->snplist file has been written out and saved as \""+\
  	 	 			extraMapInfoName+"\". \n");


 return ;}
// writee bimbam snp position file

void CBSNP::write_bimbam_snp_pos_file(const vector<CBSNP*>& genInfo, const string& outFileName)
{
	ofstream ofs;
	ofs.clear();ofs.close();
  	string datFileName=outFileName+".pos.txt";
  	ofs.open(datFileName.c_str(),ios::out);
  	if(!ofs)
  		error("out file "+outFileName+" does not exits.\n");
  	for(unsigned int i=0; i<genInfo.size();++i)
  	{
  		CBSNP* pgenInfo= genInfo[i];
  		if(pgenInfo->quality)
  		{
			if(pgenInfo->rsId!="") ofs <<pgenInfo->rsId<< ",  ";
			else if(pgenInfo->snpName!="") ofs <<pgenInfo->snpName<< ",  ";
			else if(pgenInfo->rsId==""&& pgenInfo->snpName=="")
			{
				string snpNewName="snp"+change_int_into_string(i)+"";
				ofs<< snpNewName << ",  ";
			}
			if(pgenInfo->bp!=-1  ) ofs<< pgenInfo->bp<< "\n";
			else
			{
				ofs << i << "\n";

			}
  	  }
  	}
  	ofs.clear();ofs.close();
  	printLIN("\t**->SNP location file has been written out and saved as  \""+\
	  			datFileName+"\".\n");

return;}
//
void CMSNP::lese_mach_info_mlinfo_file(const string& fileName, vector<CBSNP*>&genInfo, const float &rsq_value, const float &maf_value)
 {
 	//-------------------------------------------------------//
 	printLIN("*->Reading file \""+fileName+"\": \n");
 	string myString					="";
 	string temp_value					="";
 	string emsg						="";
 	bool ende							= false;
 	unsigned int ccol					=0;
 	unsigned int ncol					=0;
 	unsigned int nrow					=0;
 	unsigned int nLines					=0;
 	unsigned int max_ncol				=0;
 	bool def_max_ncol_once				=true;
 	//float first_prob					=0.0;
 	//float second_prob					=0.0;
 	//float third_prob					=0.0;
 	//string	pheno_value					="";
 	//string	sex_value					="";
 	vector<string>cZeile;
 	 	 try
 	 	 {
 	 		IFS pf(fileName.c_str(),"r");
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

 	 				 //-----------------------------------------------------------------//
 	 				 int spaltenZahl=cZeile.size();
 	 				 bool tfValue =(spaltenZahl==7) ;
 	 				if(!tfValue)	// Fall mit 6 Zeilen.
 	 				{
 	 					string msg="\t**->Problem on the " + change_int_into_string(nrow)+ "th row.\n";
 	 					printLIN(msg);
 	 					msg="There must be 7 columns in each line of file: "+fileName+" but there  are "+change_int_into_string(spaltenZahl)+" columns.\n";
 	 					//delete pCBPED;
 	 					 error(msg);

 	 				}else
 	 				{
 	 					if(nrow==1)
 	 					{
 	 						tfValue =(cZeile[0]!="SNP" || cZeile[1]!="Al1" ||cZeile[2]!="Al2"|| cZeile[3]!="Freq1") ||(cZeile[4]!="MAF" ||cZeile[5]!="Quality" ||cZeile[6]!="Rsq");
 	 						if(tfValue)
 	 						{
 	 							emsg =" First line of  file \""+fileName+\
 	 									"\" should contain SNP     \"Al1     Al2     Freq1   MAF     Quality Rsq\"  as column specifiers.\n";
 	 							error(emsg);
 	 						}

 	 					}
 	 					else
 	 					{
 							CBSNP* pCBSNP 	=new CBSNP;
 							//myString 		=cZeile[0];
 							pCBSNP->snpName	=cZeile[0];
 							pCBSNP->allele1 =cZeile[1];
 							pCBSNP->allele2 =cZeile[2];
 							pCBSNP->freq	=atof(cZeile[3].c_str());
 							pCBSNP->maf		=atof(cZeile[4].c_str());
 							pCBSNP->rsq		=(atof(cZeile[6].c_str()));
 							pCBSNP->quality =(pCBSNP->rsq>=rsq_value) && (pCBSNP->maf>=maf_value);
 							//cout << pCBSNP->allele1 << " "<< pCBSNP->allele2 <<endl;
 							genInfo.push_back(pCBSNP);

 						}
 	 				}

 	 			 }
 	 			 ncol=0;
 	 			 cZeile.clear();
 	 			 cZeile.resize(0);
 	 			 if(_ifs_eof(pf))break;
 	 	} // end of outer while loop
 	 		_ifs_close(pf);
 	 		nrow=0; // total rows
 	 		nLines=0; // total lines
 	 		//display_snp_summary(genInfo);

 	 	 }catch(bad_alloc& memoryAllocationException)
 	 	 {
 	 		 cout<<"\n Please inform the programmer that you could not read mlinfo file. Probbaly he will give you a solution. \n ";
 	 		 error( "\n bad _alloc()  in reading mach mlgeno file \n");
 	  	 }


 return;}

 //
void CBSNP::display_snp_summary(const vector<CBSNP*>&genInfo)
{
	int good_snps=0;
	int bad_snps=0;

	for(vector<CBSNP*>::const_iterator it=genInfo.begin(); it<genInfo.end(); ++it )
	{

		if((*it)->quality)
			++good_snps;
		//cout << boolalpha << (*it)->quality << " "; //debug
	}
	bad_snps=(genInfo.size()-good_snps);
	ngood_snps=good_snps;
	if(bad_snps>0)
	{
		cout <<left;
		cout <<"\n*->SNP summary:\n";
		// male and female Geschiste
		cout <<"\t"<< setw(coutWidth)<<	"**->Total good Quality SNPs: " << good_snps<<"\n";
		cout <<"\t"<< setw(coutWidth)<<	"**->Total bad Quality SNPs: " << bad_snps<<"\n";
		cout <<"\t"<< setw(coutWidth)<< 	"**->SNPs marked as bad Quality SNPs will not be considered for Format converting process.\n";
		cout.flush()<<endl;
		//
		LIN <<left;
		LIN <<"\n*->SNP summary:\n";
		// male and female Geschiste
		LIN <<"\t"<< setw(coutWidth)<< 	"**->Total good Quality SNPs: " << good_snps<<"\n";
		LIN <<"\t"<< setw(coutWidth)<< 	"**->Total bad Quality SNPs: " << bad_snps<<"\n";
		LIN <<"\t"<< setw(coutWidth)<< 	"**->SNPs marked as bad Quality SNPs will not be considered any more for format converting process.\n";
		LIN.flush()<<endl;
	}
	else{

		cout <<left;
		cout <<"\n*->SNP summary:\n";
		cout <<"\t"<< setw(coutWidth)<<	"**->Total number of successfully read SNPs: " << good_snps<<"\n";
		cout.flush()<<endl;
		LIN <<"\t"<< setw(coutWidth)<< "**->Total number of successfully read: " << good_snps<<"\n";
		LIN.flush()<<endl;
	}


return; }

// read mach_ref_snps file

void CMSNP::read_mach_ref_snpsFile(vector<CBSNP*>& genInfo,const string& snpsFileName)
{
	// here comes the read extra file like structre.
	string r_msg = "\n*->Reading mach  reference snps file:\""+snpsFileName+"\":\n";
	printLIN(r_msg);
	string line="";
	string buffer="";
	rFile file; file.close(); file.clear();
	static int extraMapFileLineZaehler =0;
	file.open(snpsFileName.c_str(),ios::in);
	if(!file)
	{
		file.setstate(ios_base::failbit);
		file.clear();
		error("Your snps file \""+snpsFileName+"\" eitehr does not exit or can't be opened. Have you given the correct path and name of the file?\n");

	}

	vector<string> tokens;
	while(getline(file,line, '\n')) // until eof hit
	{
		//---------------------//
		++extraMapFileLineZaehler;
		//cout <<"line no: "<< extraMapFileLineZaehler<<endl; //debug
		//cout<< line.length()<<endl;
		if(line.length()==0)
			continue;
		if(line[0]=='#')
			continue;
		//if(!((*it_pCBSNP)->quality))
		//	continue;
		tokens.resize(0);
		tokens.clear();
		stringstream line_parser1(line);
		//cout <<boolalpha<<"\n line_parser1.good(): "<<line_parser1.good() <<endl; // only for test
		if(line_parser1.good())
		{
			while(line_parser1>>buffer)
			{
				tokens.push_back(buffer);
				//cout << buffer << " "; // only for test
			}
			if(tokens.size()!=1)
			{
				string msg= "Your snps file \""+snpsFileName+"\"  must contain 1 column in the line "+change_int_into_string(extraMapFileLineZaehler)+\
						", however there is/are "+change_int_into_string(tokens.size())+" column/s.\n";
				error(msg);
			}
			else
			{
				CBSNP* pCBSNP		=new CBSNP;
				pCBSNP->rsId 		=tokens[0];
				//pCBSNP->bp 			=atoi(tokens[1].c_str());
				//pCBSNP->allele1		=tokens[2];
				//pCBSNP->allele2		=tokens[3];
				genInfo.push_back(pCBSNP);
				//cout <<" genInfo.size()"<< genInfo.size()<<endl; //test

			}

		}


	}
	//cout << "snps : "<< extraMapFileLineZaehler<<endl;

}//end of function


