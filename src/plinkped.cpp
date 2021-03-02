/*
 * pedgree.cpp
 *
 *  Created on: Dec 9, 2011
  *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.

 */
#include <bitset>
#include "plinkped.h"
#include "helper.h"
#include "gz_stream.h"
#include <iostream>
#include<bitset>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
extern ofstream LIN;
 const int cw=35; // width of cout
 CPPED::~CPPED()
 {
	 /*
	 for (unsigned int i=0;i<vPed.size();++i)
	 	{
	 		delete vPed[i];
	 	
	 	}
	 	vPed.clear();
	 	*/
		//cout << "\n Good bye from Cpedgree destructor.\n";
 }
  bool CPPED::_read_recodeA_type=false;
  bool CPPED::_read_recodeAD_type=false;

// functions
 bool  CPPED::zeilenLeser(IFS pf, vector<string>& pZeile){
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
			++spaltenZahl;
			ende=true;
	     // Ignore rest of line and advance to next line
			while(myChar!='\n'&& myChar!='\r'&& !_ifs_eof(pf))
				myChar =_ifs_getc(pf);
			//ifs_read(pf, "%*[^\n]");
			//(void) _ifs_getc (pf);
			//if(_ifs_eof(pf) || (myChar=='\r') || (myChar=='\n'))
			if(!strng.empty())
						pZeile.push_back(strng);

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
 //
  // spalten leser
 CBPED*  CPPED::pedInfoLeser(const vector<string>& zeilenVec, const long int& snpAnzahl, const int& akt_zeile)
 {
 		long int spaltenZahl=zeilenVec.size();
 		// check if  vector<CBSNP*> genInfo has  size equal to snpAnzahl.
 		// this declaration will produce pointer to CpCBPEDclasses
 		CBPED* pCBPED=new CBPED;
 		// check if there are correct pheno infos and genotype infos. plink has either 6 phenotype infos.
 		int value = spaltenZahl-6; // 6 is number columns of phenotype infos
 			//cout <<"value: "<<value<<endl;
 			//cout << boolalpha<< (value==snpAnzahl)<<endl;
 			//cout << boolalpha<< (value/2==snpAnzahl)<<endl;
 		bool tfValue =((value==snpAnzahl) || (value/2==snpAnzahl));

 		if(!tfValue)	// Fall mit 6 Zeilen.
 		{
 			string msg="**->Problem on the " + change_int_into_string(akt_zeile)+ "th row.\n";
 			cout <<"\t"<< setw(cw)<<msg;
 			msg="There must be 6 columns  for phenotype infos and  either exactly "+change_int_into_string(snpAnzahl)+\
 			" columns  with \"A/B\" type of genotype Infos or exactly "+change_int_into_string(2*snpAnzahl)+" columns \"A   B\" type of  genotype infos.\n";
 			cout <<"\t"<< setw(cw);
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
 			if(zeilenVec[4]=="1")
 			{	//male
 				pCBPED->sex			=true;
 				pCBPED->miss_sex		=false;

 			}
 			else if (zeilenVec[4]=="2" )
 			{ //female
 				pCBPED->sex=false;
 				pCBPED->miss_sex=false;

 			}
 			else
 			{	// missing
 				pCBPED->miss_sex=false;
 				pCBPED->miss_sex=true;
 			}
 			// for phenotype
 			if(zeilenVec[5]=="2")
 			{ // affected
 				pCBPED->pheno		=true;
 				pCBPED->miss_pheno	=false;
 			}
 			else if(zeilenVec[5]=="1")
 			{ // unaffected
 				pCBPED->miss_pheno	=false;
 				pCBPED->pheno		=false;
 			}
 			else
 			{ // missing
 				pCBPED->miss_pheno	=true;
 				pCBPED->pheno		=false;
 			}

 			if(pCBPED->sex)									++nMale;
 			else if(!pCBPED->sex && !pCBPED->miss_sex)	++nFemale;
 			else if(pCBPED->miss_sex)						++nNosex;
 			// about phenotye
 			if(pCBPED->pheno) 							++nAffected;
 			else if(!pCBPED->pheno && !pCBPED->miss_pheno)	++nUnaffected;
 			else if(pCBPED->miss_pheno) 						++nMiss_status;
 		}

  return pCBPED;}

 void CPPED::pedFileleser(const string& pedFileName,\
 	const vector<CBSNP*>&genInfo, vector<CBPED*> &vPedInfo)
 {
	string pr_msg="\n*-> Reading  ped file: \""+pedFileName+ "\":\n";
	printLIN(pr_msg);
	 const long int snpAnzahl=genInfo.size();
 	//cout << "snpAnzahl: "<<snpAnzahl <<endl;
 	// we need CBPED*
 	 CBPED* pCBPED;
 	//rFile::checkFileExists(pedFileName);
	//rFile ofs(pedFileName);

 	  IFS  file;
 	  file =_ifs_open(pedFileName.c_str(),"r");
 		if(!file)
 		error("file does not exits´or could not open");
 	vector<string> zeilenVector;
 	string sline="";
 	string elem;
 	string elem1;
 	//
 	//int nMan=0;
 	//int nWeib=0;
 	//int nTrans=0;
 	while(true)
 	{
 		if(_ifs_eof(file))		break;
 		//if(ofs.eof()) break;
 		static int akt_zeile=0;	 // aktuelle Zeile
 		// lese Zeile
 		//sline =ofs.readLine();
 		//zeilenVector =ofs.split_line_as_vector(sline);
 		//bool tfValue=(zeilenVector.size()>0);

 		bool tfValue=zeilenLeser(file, zeilenVector);
 		int spaltenZahl=zeilenVector.size();
 		if(spaltenZahl==0)
 			continue;
 		//cout << "\n zeilenVector.size(): "<<zeilenVector.size() <<endl; // for test only
 		if(tfValue) akt_zeile++;
 		//string msg="\t**->Total number of columns in the "+change_int_into_string(akt_zeile)+ "th row is : "+change_int_into_string(spaltenZahl)+".\n";
 		//cout <<msg;
 		//LIN <<msg;
 		// lese spalten
 		pCBPED= pedInfoLeser(zeilenVector,snpAnzahl, akt_zeile);
 		// add the pCBPED into vector
 		vPedInfo.push_back(pCBPED);
 		// now erase the first 6 columns of pedInfo from ZeilenVector
 	 	zeilenVector.erase(zeilenVector.begin(),zeilenVector.begin()+6);
 		// read genotpye now
 			//genInfoLeser(zeilenVector,genInfo, akt_zeile);
 			lese_pedfile_genInfo(genInfo,zeilenVector,akt_zeile);
 		//end des loops.
 		spaltenZahl=0;
 		zeilenVector.clear();
 		zeilenVector.resize(0);
 		sline="";
 	}

   _ifs_close(file);
  // ofs.close(); ofs.clear();
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
 	cout <<"\t"<< setw(cw);
 	printLIN(pmsg);
 }
//read plink-dosage format file
 void CPPED::lese_plink_dosage_file(const vector<CBPED*>&pedInfo, vector<CBSNP*>&genInfo,vector<string>&cZeile,const float& thresh,const string&dosageFileName)
 {
	// cout<<" HI I am inside the  lese plink_dosage file now\n";//debug
	 printLIN("*->Reading file named \""+dosageFileName+"\": \n");
	IFS file(dosageFileName.c_str(),"r");
	if(!file){
		string msg= dosageFileName+" either does not exits or could not be opend.\n Plase make sure that you have given the correct file path and name.\n";
		error(msg);
	}
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
	bool 	def_nSample					=true;
	bool	given_mapfile				=false;
	bool 	check_mapfile				=true;
	unsigned long int nSample			=0;
	try
	{


		if(check_mapfile)
		{
			if(genInfo.size()!=0)
			{
				snpAnzahl =genInfo.size();
				given_mapfile=true;
			}
		}
		check_mapfile=false;
		//if(given_mapfile)
		//{
			vector<CBSNP*>::const_iterator  ppCBSNP =genInfo.begin();
		//}
		while(true)
		{
			if(_ifs_eof(file))break;

			while(!_ifs_eof(file))
			{
				ccol=CBSNP::stringLeser(file,myString,ende);
				ncol+=ccol;
				if(ccol!=0){
					if(!myString.empty())
						cZeile.push_back(myString);

					else{
							cout << " error in reading "<<dosageFileName<<" \n";
							exit(1);
					}
				}
				if(ende|| _ifs_eof(file))
				{
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
			//-----------------------------//
			}//end of _ifs_eof(file)

			if(cZeile.size()!=0)
			{
				//first row contains identifier.
				if(nrow==1)
				{
					bool _temp_bool =(cZeile[0]!="marker"&&cZeile[0]!="SNP")||(cZeile[1]!="A1"||cZeile[2]!="A2");
					if(_temp_bool)
					{
						 emsg ="First line of plink dosage format file \""+dosageFileName+\
								 "\" should start with the word \"SNP\", the second column with the word \"A1\" and the thrid  column with  the word \"A2\" .\n";
						 error(emsg);
					}
					else{
						cZeile.erase(cZeile.begin(),cZeile.begin()+3);
						if(def_nSample){
							nSample =cZeile.size()/2;
							def_nSample=false;
						}
						if(nSample!=pedInfo.size()){
							error("Number of individuals in fam file is not the same as the number of individuals in dosage file.\n");
						}
						CBPED* pCBPED = pedInfo[0];
						bool _temp_bool =true;
						 for(unsigned int i=0; i<nSample; ++i)
						 {

							// CBPED* pCBPED = new CBPED;
							 pCBPED = pedInfo[i];
							 _temp_bool= (pCBPED->famId==cZeile[2*i])&& (pCBPED->indId==cZeile[2*i+1]);
							if(!_temp_bool)
							{
								printLIN("family id or individual id in fam file and dosage file in"+change_int_into_string(i+1)+"th individual does not match. \n");
							}

						 }

					//----
					}//end of else
					//

					//check the  columns third and onwards

				}// end of if nrow==1
				// from 2nd row.
				//cout <<pedInfo.size()<<" ";
				//####################################
				else
				{
					const string rsName=cZeile[0];
					const string allele_A=cZeile[1];
					const string allele_B=cZeile[2];
					  cZeile.erase(cZeile.begin(),cZeile.begin()+3);
					  // now the line should have 2* times nSNPs elements
					  if(cZeile.size()==(2*nSample))
					  {
						  if(given_mapfile)
						  {
							  (*ppCBSNP)->no_of_persons=nSample;
							  if((*ppCBSNP)->allele1=="")
								  (*ppCBSNP)->allele1=allele_A;
							  if((*ppCBSNP)->allele2==""&& (allele_A!=allele_B))
								  (*ppCBSNP)->allele2=allele_B;
							  float _temp_summe=0.0;
							  for(unsigned int i=0; i<nSample; ++i)
							  {
								  first_prob  =atof(cZeile[2*i].c_str());
								  second_prob =atof(cZeile[2*i+1].c_str());
								  third_prob =abs(1.0-(first_prob+second_prob));
								  _temp_summe=(first_prob+second_prob+third_prob);
								  first_prob 	= first_prob/_temp_summe;
								  second_prob	 = second_prob/_temp_summe;
								  third_prob	 = third_prob/_temp_summe;
								  (*ppCBSNP)->pgeno1.push_back(first_prob);
								  (*ppCBSNP)->pgeno2.push_back(second_prob);
								  (*ppCBSNP)->pgeno3.push_back(third_prob);
							  }
							 if(nrow-1<snpAnzahl)
								 ++ppCBSNP;
							 else break;
						}
						  else
						  {
							  CBSNP* pCBSNP = new CBSNP;
							  if(pCBSNP->rsId=="")
								  pCBSNP->rsId=rsName;
							  if(pCBSNP->snpName=="")
								  pCBSNP->snpName=rsName;
							  if(pCBSNP->allele1=="")
								  pCBSNP->allele1=allele_A;
							  if(pCBSNP->allele2==""&&(allele_A!=allele_B))
								  pCBSNP->allele2=allele_B;
							  genInfo.push_back(pCBSNP);
							  if(genInfo.size()==(nrow-1)) // -1 means first row was indiv info.
							  {
								  CBSNP* pCBSNP=genInfo[genInfo.size()-1];
								  pCBSNP->no_of_persons=nSample;
								  for(unsigned int i=0; i<nSample; ++i)
								  {
									  first_prob  =atof(cZeile[2*i].c_str());
									  second_prob =atof(cZeile[2*i+1].c_str());
									  third_prob =abs(1.0-(first_prob+second_prob));
									  pCBSNP->pgeno1.push_back(first_prob);
									  pCBSNP->pgeno2.push_back(second_prob);
									  pCBSNP->pgeno3.push_back(third_prob);
								  }
							  }
							  else
							  {

								  emsg =" problem in saving genotypes line in file \""+dosageFileName+\
										  "\".\n";
								  error(emsg);
							  }


						  }
					  }
					 else
					 {
						 emsg =" problem in saving genotypes line in file \""+dosageFileName+\
								 "\". There are not enough columns for genotype data.\n";
						 error(emsg);
					 }

				  }//end of genotype saving else loop

				 //CBPED* pCBPED= lese_pedfile_pedInfo(cZeile,snpAnzahl, nrow);
				 //pedInfo.push_back(pCBPED);
				 //cout << "nrow: "<< nrow<< " pedInfo.size()"<< pedInfo.size()<<endl;
				  if(max_ncol!=ncol)
				 {
					 string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+dosageFileName+"\".\n";
					 error(msg);
				 }


			//##############################
			}//end of if cZeile.size()!=0
			else continue;
			 ncol=0;
			 cZeile.clear();
			 cZeile.resize(0);
		 	if(_ifs_eof(file))break;

		}//while(true) loop
		//---------------------------//
		_ifs_close(file);
 		nrow=0; // total rows
 		nLines=0; // total lines
 		//update_pedInfo_given_gProb(pedInfo,genInfo);
 		genoProb_to_geno_converter(genInfo, thresh);
		string pmsg="\t**->Total number of successfully read SNPs: ";
		printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
	}catch(bad_alloc& memoryAllocationException)
	{
		cout<<"\n Please inform the programmer that you could not read mlgeno file. Probbaly he will give you a solution. \n ";
		error( "\n bad _alloc()  in reading plink dosage file"+dosageFileName +"\n");
 	 }

return;}
void CPPED::lese_plink_raw_file(vector<CBPED*>&pedInfo, vector<CBSNP*>&genInfo,vector<string>&cZeile,const string&rawFileName)
{
	//cout << "Hi! I am inside lese plink_raw file .\n";//debug
	printLIN("*->Reading file named \""+rawFileName+"\": \n");
	IFS file(rawFileName.c_str(),"r");
	if(!file){
		string msg= rawFileName+" either does not exits or could not be opend.\n Plase make sure that you have given the correct file path and name.\n";
		error(msg);
	}
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
		// unsigned long int	snpAnzahl		=0;
		string	pheno_value					="";
		string	sex_value					="";
		bool 	def_nSample					=true;
		//bool	given_mapfile				=false;
		//bool 	check_mapfile				=true;
		//unsigned long int nSample			=0;
		unsigned long int _nSNPs			=0;
		string 	_min_allele					="";
		string  _cur_rsid					="";
		string _temp_string 				="";
		int 	_temp_int					=-1;
		try
		{


			//vector<CBSNP*>::const_iterator  ppCBSNP =genInfo.begin();
			// we are going to save genotypes into geno_0123. So  the
			// identifier "changed_0123" (it recognizes whether the vector is full) must be true now.

			while(true)
			{
				if(_ifs_eof(file))break;

				while(!_ifs_eof(file))
				{
					ccol=CBSNP::stringLeser(file,myString,ende);
					ncol+=ccol;
					if(ccol!=0)
					{
						if(!myString.empty())
							cZeile.push_back(myString);

						else{
								cout << " error in reading "<<rawFileName<<" \n";
								exit(1);
						}
					}
					//cout << ende << endl;
					if(ende|| _ifs_eof(file))
					{
						++nLines;
						if(def_max_ncol_once &&(ncol!=0))
						{
							max_ncol=ncol;
							def_max_ncol_once=false;
							//cout << "max_ncol: "<<max_ncol <<" " <<"ncol: "<<ncol <<endl; // only test; // only test
						}
						//cout << ncol << " "<<max_ncol << endl;//debug
						if((ncol==max_ncol)&& (max_ncol!=0))
							++nrow;
						myString="";
						//cout <<"-----------------\n";//test
						break;
					}
					if(ccol==0)
						continue;
					// this loop is broken only if end  is found.
				//-----------------------------//
				}//end of _ifs_eof(file)

				if(cZeile.size()!=0)
				{
					//first row contains identifier.
					if(nrow==1)
					{

						bool _temp_bool =(cZeile[0]!="FID"||cZeile[1]!="IID")||(cZeile[2]!="PAT"||cZeile[3]!="MAT"||cZeile[4]!="SEX"||cZeile[5]!="PHENOTYPE");
						//cout <<"norw: "<<nrow<< " _temp_bool: "<<boolalpha << _temp_bool << " \n";
						if(_temp_bool)
						{
							_temp_bool=false;
							emsg ="First 6 columns of first line of plink raw format file \""+rawFileName+\
									 "\" should contain words:"
									 "\"FID\",\"IID\",\"PAT\",\"MAT\",\"SEX\" and \"PHENOTYPE\" .\n";
							 error(emsg);
						}
						else{
								cZeile.erase(cZeile.begin(),cZeile.begin()+6);
							if(def_nSample)
							{
								//cout << "_read_recodeA_type: "<<_read_recodeA_type <<" "<<cZeile.size() <<endl;
								if(_read_recodeAD_type)
									_nSNPs =cZeile.size()/2;
								else if(_read_recodeA_type)
									_nSNPs =cZeile.size();
								def_nSample=false;
							}
							// cout << _nSNPs<<endl;
							CBSNP* pCBSNP = genInfo[0];
							for(unsigned int j=0;j<_nSNPs;++j)
							{
								pCBSNP			=genInfo[j];
								if(_read_recodeAD_type)
									_temp_string =cZeile[2*j];
								else
									_temp_string =cZeile[j];
								unsigned found 	=_temp_string.find_last_of("_");
								_cur_rsid 		=_temp_string.substr(0,found);
								//cout <<_cur_rsid << endl;
								_min_allele		=_temp_string.substr(found+1);
							//	cout << _min_allele<<endl;
								if(_cur_rsid!=pCBSNP->rsId&& _cur_rsid!=pCBSNP->snpName)
								{
									printLIN("rsid in map file:"+ pCBSNP->rsId+"\n" );
									printLIN("rsid in raw  file:"+ _cur_rsid+"\n ");
									error("rsid both in map and raw file does not match for the "+change_int_into_string(j+1)+"th SNP.\n");

								}
								if(_min_allele==pCBSNP->allele1||_min_allele==pCBSNP->allele2)
								{
									//(_min_allele!="0") //
									pCBSNP->min_allele=_min_allele;
									//pCBSNP->allele2=_min_allele;
									// we do not know allele1 yet. but allele2 is minor allele.
									// if we update alleles with --snpinfo then:
									// 1.
								if(_min_allele==pCBSNP->allele1)
									pCBSNP->maj_allele=pCBSNP->allele2;
								else
									pCBSNP->maj_allele=pCBSNP->allele1;
								}
								else if(_min_allele=="0")
								{
									//if min_allele does not exits .i.e. if the SNP is monomorf
									pCBSNP->maj_allele =pCBSNP->allele1;
									pCBSNP->min_allele =pCBSNP->allele2;
									if(pCBSNP->allele1!=pCBSNP->allele2)
									{
										// This is a complicated case: It is a very difficult situation.
										// Throw warning : say first allele is  assumed as major and second allele as minor.
										// If you are not satisfied with this case, then you should change the allele position at this snp first.
										string _warn_message = "\t---------------------------------------------------------------------------\n"
												"\tWARNINGS:\n"
												"\tFile \""+rawFileName+"\" has no minor allele at SNP \""+pCBSNP->rsId+"\",\n"
												"\tBut two different alleles are already updated at this SNP with \"--snpinfo \" option. \n"
												"\tAt such a situation, it is difficult to make a decision about major and minor allele coding.\n"
												"\tHowever, we assume that the first allele is major and the second is minor.\n"
												"\tIf this is not the case, please exchnage the allele order in the file given with \"--snpinfo \" option. \n"
												"\t---------------------------------------------------------------------------\n";
										printLIN(_warn_message);
									}


								}
								else
								{
									printLIN("rsid in map file:"+ pCBSNP->rsId+"\n" );
									printLIN("rsid in raw  file:"+ _cur_rsid+"\n ");
									printLIN("alleles in map file: "+ pCBSNP->allele1+" and  "+pCBSNP->allele2+".  \n" );
									error(_min_allele+" given after \"_\" symbol in "+_temp_string+" must match with the alles given in map file. But this is not the case here in "+change_int_into_string(j+1)+"th SNP.\n");

								}
							}
							//here we will check if the header is correct
							// this is just to make sure that rsid in map file and headers
							//are correct.

						}

					}// end of if nrow==1
					// from 2nd row.
					//cout <<pedInfo.size()<<" ";
					//####################################
					else
					{
						//check first if there are enough snps in raw file
						//that means  no of snps in map file must be equal to no of snps in raw file.
						//cout <<"_nSNPs: "<< _nSNPs << " "<<genInfo.size()<<endl;
						if(_nSNPs!=genInfo.size())
							error("Total SNPs in map file and raw file must be equal to each other. But this is not the case in"+change_int_into_string(nrow)+"th row. \n");
							//cout <<"here: "<<cZeile.size()<<"-> "<<_nSNPs <<",row: "<<nrow<< endl;
						CBPED* pCBPED= lese_pedfile_pedInfo(cZeile,_nSNPs, nrow);
						pedInfo.push_back(pCBPED);
						//cout << "nrow: "<< nrow<< " pedInfo.size()"<< pedInfo.size()<<endl;
						if(nrow-1!=pedInfo.size())
						{
							int sz=pedInfo.size()-1;
							string msg= "Problem in saving genotypes of person "+pedInfo[sz]->indId+" in the "+change_int_into_string(nLines)+"th row of   file \""+rawFileName+"\".\n";
							error(msg);
						}
						if(max_ncol!=ncol)
						{
							string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+rawFileName+"\".\n";
							error(msg);
						}
						// now erase the first 6 columns of pedInfo from ZeilenVector
						cZeile.erase(cZeile.begin(),cZeile.begin()+6);
						//cout <<"czeile size: "<< cZeile.size() << " ";
						//cout <<"test2"<<endl;
						CBSNP* pCBSNP=genInfo[0];
						for(unsigned int j=0;j<_nSNPs;++j)
						{
							pCBSNP =genInfo[j];
							if(_read_recodeAD_type)
								_temp_string =cZeile[2*j];
							else
								_temp_string =cZeile[j];
							if (cZeile[2*j]=="NA")
							{
								pCBSNP->geno_0123.push_back(3);
								// for missing
								pCBSNP->geno1.push_back(true);
								pCBSNP->geno2.push_back(false);
								pCBSNP->aOrder.push_back(true);

							}
							else
							{
								_temp_int=atoi(cZeile[2*j].c_str());
								if((_temp_int==0))
								{
									// 0 means here major homozygote
									// if maj_allele == allele1 then this should be saved as homozygote major
									// else as homo minor
									//homozygote major
									pCBSNP->geno_0123.push_back(_temp_int);
									if(pCBSNP->allele1==pCBSNP->maj_allele)
									{
										pCBSNP->geno1.push_back(false);
										pCBSNP->geno2.push_back(false);
										pCBSNP->aOrder.push_back(true);
									}
									else { // this is the case where allele2 = maj_allele
										pCBSNP->geno1.push_back(true);
										pCBSNP->geno2.push_back(true);
										pCBSNP->aOrder.push_back(true);
									}
								}
								else if((_temp_int==1))
								{
									//heterozygote
									pCBSNP->geno_0123.push_back(_temp_int);
									pCBSNP->geno1.push_back(false);
									pCBSNP->geno2.push_back(true);
									pCBSNP->aOrder.push_back(true);
								}
								else if((_temp_int==2))
								{
									//homozygote minor
									pCBSNP->geno_0123.push_back(_temp_int);
									if(pCBSNP->allele1==pCBSNP->min_allele)
									{//this is the case where allele 1 and min allele
										pCBSNP->geno1.push_back(false);
										pCBSNP->geno2.push_back(false);
										pCBSNP->aOrder.push_back(true);
									}
									else { // this is the case where allele2 = min_allele
										pCBSNP->geno1.push_back(true);
										pCBSNP->geno2.push_back(true);
										pCBSNP->aOrder.push_back(true);
									}
								}
								_temp_int=-1;
							}

						}
					}
					 //CBPED* pCBPED= lese_pedfile_pedInfo(cZeile,snpAnzahl, nrow);
					 //pedInfo.push_back(pCBPED);
					// cout << "nrow: "<< nrow<< " pedInfo.size()"<< pedInfo.size()<<endl;
					  if(max_ncol!=ncol)
					 {
						 string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+rawFileName+"\".\n";
						 error(msg);
					 }

					  cZeile.clear();
					  cZeile.resize(0);

				//##############################
				}//end of if cZeile.size()!=0
				else continue;
				 ncol=0;
				 cZeile.clear();
				 cZeile.resize(0);
			 	if(_ifs_eof(file))break;

			}//while(true) loop
			//---------------------------//
			_ifs_close(file);
	 		nrow=0; // total rows
	 		nLines=0; // total lines
	 		// update maf and allele freq
	 		CBPED::calculate_maf(genInfo,pedInfo);
	 		//update_pedInfo_given_gProb(pedInfo,genInfo);
	 		//genoProb_to_geno_converter(genInfo, thresh);
	 		//--------------------------------------------//
	 		//to check if saved correctly
	 		//cout <<"test23"<<endl;
	 	 	for (unsigned int i=0; i<genInfo.size();++i)
	 	 		{
	 	 			 CBSNP* pCBSNP=genInfo[i];

	 	 			bool tfValue=false;
	 	 			if(pCBSNP->geno1.size()!=pedInfo.size()&& pCBSNP->geno2.size()!=pedInfo.size())
	 	 			{
	 	 				tfValue=true;
	 	 				cout << "geno1.size(): " << pCBSNP->geno1.size()<<endl;
	 	 				cout << "geno2.size(): " << pCBSNP->geno2.size()<<endl;
	 	 				cout << "vPedInfo.size(): " << pedInfo.size()<<endl;
	 	 				cout << "allele1: " << pCBSNP->allele1<<endl;
	 	 				cout << "allele2: " << pCBSNP->allele2<<endl;
	 	 			}
	 	 			if(tfValue)
	 	 			{
	 	 				//delete	pCBSNP;
	 	 				error("Problem in saving genotypes in  "+ change_int_into_string(i) +"th SNP with Name "+pCBSNP->snpName+".\n");

	 	 			//for(unsigned int j=0; j<pCBSNP->geno1.size();++j)
	 	 			//	cout << pCBSNP->geno1[j]<< " "<< pCBSNP->geno2[j]<< " \n";
	 	 				//cout << endl;
	 	 			}
	 	 		}

	 	 	string pmsg="**->Total number of successfully read individuals: "+	change_int_into_string(pedInfo.size())+ "\n";
	 	 	//printLIN(pmsg);
	 	 	cout <<"\t"<< setw(cw);
	 	 	printLIN(pmsg);

	 		//-------------------------------------------//
			 pmsg="\t**->Total number of successfully read SNPs: ";
			printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
		}catch(bad_alloc& memoryAllocationException)
		{
			cout<<"\n Please inform the programmer that you could not read  file. Probbaly he will give you a solution. \n ";
			error( "\n bad _alloc()  in reading plink dosage file"+rawFileName +"\n");
	 	 }


	//###########################################################
return ;}
 //fam file
void CPPED::lese_plink_fam_file( vector<CBPED*>&pedInfo,const string&famFileName)
{
	//cout <<"Hi! I am from plink's fam file leser!\n ";
	static long int countIndiv=0;
	string pr_msg="*->Reading fam file \""+famFileName+"\": \n";
	printLIN(pr_msg);
	string 	line 		="";
	string 	buffer 		="";
	string 	sex_value	="";
	string 	pheno_value	="";
	//unsigned int nCols  =0;
	vector<string> tokens;
	rFile file; file.close();file.clear();
	file.open(famFileName.c_str(),ios::in);
	if(!file){
		file.setstate(ios_base::failbit);
		file.clear();
		error("fam file \""+famFileName+"\" either does not exit or can not be opened.\n");
	}
	while(getline(file,line,'\n'))
		 {
			if(line=="") // this blank line
				continue;
			if(line[0]=='#') // this is comment line
				continue;
			stringstream line_parser(line);

			if(line_parser.good())
			{
				 ++countIndiv;
				 while(line_parser>>buffer)
					tokens.push_back(buffer);
					if(tokens.size()==6)
					{
						CBPED* pCBPED=new CBPED;
						pCBPED->famId		=	tokens[0];
						pCBPED->indId		=	tokens[1];
						pCBPED->patId		=	tokens[2];
						pCBPED->matId		=	tokens[3];
						sex_value			= tokens[4];
						lese_sex_info(pCBPED, sex_value);
						//for phenotype
						pheno_value			=tokens[5];
						lese_phenotype_info(pCBPED, pheno_value);
						pedInfo.push_back(pCBPED);
					}
					else{
						error("The " + change_int_into_string( countIndiv)+"th row (excluding comments and blank lines)"+
								"of your fam file does not contain 6 columns. Please delete or correct this row.\n");
					}
			}
			tokens.clear();
			tokens.resize(0);
		 }//end of while getline
	file.close();file.clear();
	printLIN("\t**->Total no of individuals: "+change_int_into_string((int)pedInfo.size())+ " \n");
}//end of function

 //alternative

 // converting functions
  void CPPED::writeOutMachFiles(const vector<CBPED*>& pedInfo ,const vector<CBSNP*> &genInfo, const string&outFileName)
 {
 //----------------------------------#--------------------------------------#
 	printLIN("*->Writing files in "+CBSNP::_out_format+" format: \n");
 //----------------------------------#--------------------------------------#
 	ofstream ofs;
 	 ofs.clear();ofs.close();
 	string datFileName=outFileName+".dat";
 	ofs.open(datFileName.c_str(),ios::out);
 	if(!ofs)
 		error("out file "+outFileName+" does not exits.\n");
 	for(unsigned int i=0; i<genInfo.size();++i)
 		ofs<< "M" << " "<<genInfo[i]->snpName <<"\n";
 	ofs.clear();ofs.close();

 	printLIN("\t **->MACH map file has been written out and saved as  \""+\
 			datFileName+"\".\n");


//----------------------------------#--------------------------------------#
 	// writing ped files.
 //----------------------------------#--------------------------------------#

 	string pedFileName=outFileName+".ped";
 	ofs.open(pedFileName.c_str());
 	if(!ofs)
 			error("out file "+outFileName+" does not exits.\n");
 	ofs.clear();

 	for(unsigned int i=0; i<pedInfo.size();++i)
 	{
 		CBPED* pCBPED =pedInfo[i];
 		ofs<<pCBPED->famId<< " " <<pCBPED->indId << " "<<pCBPED->matId <<" "<<pCBPED->patId <<" ";
 		if(pCBPED->sex) ofs<<"M";
 		else if(!pCBPED->sex && !pCBPED->miss_sex) ofs<<"F";
 		else if(pCBPED->miss_sex) ofs<<"0"; // 0 is missing
 		for(unsigned int j=0;j<genInfo.size(); ++j)
 		{
 			CBSNP* pCBSNP=genInfo[j];
 			if(pCBSNP->geno1.size()!=pedInfo.size() || pCBSNP->geno2.size()!=pedInfo.size() )
 			{
 				cout <<"geno1.size() : "<<pCBSNP->geno1.size() <<"geno2.size():  "<<pCBSNP->geno1.size();
 				for(unsigned int i=0;i<pCBSNP->geno1.size();++i)
 					ofs<< pCBSNP->geno1[i] << " " <<pCBSNP->geno2[i] << " ";
 				ofs.close();
 				error("Problem in writing out "+change_int_into_string(j) +"th SNP. geno1 and geno2 has no equal size\n");
 			}
 			else
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


 		}
 		ofs<< "\n";
 	}
 	printLIN("\t**->MACH ped file has been written out and saved as \""+pedFileName+\
 			"\".\n");
 	ofs.clear();	 ofs.close();


 return;}
 void CPPED::writeOutImputeFiles(const vector<CBPED*>& pedInfo ,const vector<CBSNP*> &genInfo, const string&outFileName)
 {
	 //----------------------------------#--------------------------------------#
	  	printLIN("*->Writing files in "+CBSNP::_out_format+" format: \n");

	//----------------------------------#--------------------------------------#
 	// writing gens files.
  	//----------------------------------#--------------------------------------#
 	ofstream ofs;
 	 ofs.clear();ofs.close();

 string pedFileName=outFileName+".gens";
 	ofs.open(pedFileName.c_str());
 	if(!ofs)
 			error("out file "+outFileName+" does not exits.\n");
 	ofs.clear();


 		for(unsigned int j=0;j<genInfo.size(); ++j)
 		{
 			CBSNP* pCBSNP=genInfo[j];


 			if(pCBSNP->geno1.size()!=pedInfo.size() || pCBSNP->geno2.size()!=pedInfo.size() )
 			{
 				cout <<"geno1.size() : "<<pCBSNP->geno1.size() <<"geno2.size():  "<<pCBSNP->geno1.size();
 				for(unsigned int i=0;i<pCBSNP->geno1.size();++i)
 					ofs<< pCBSNP->geno1[i] << " " <<pCBSNP->geno2[i] << " ";
 					ofs.close();
 					error("problem in writing out "+change_int_into_string(j) +"th SNP. geno1 and geno2 has no equal size\n");
 			}
 			else
 			{
 				if(pCBSNP->snpName!="") ofs << pCBSNP->snpName<< " ";
 				else if(pCBSNP->rsId!="") ofs << pCBSNP->rsId<< " ";
 				if(pCBSNP->rsId!="") ofs << pCBSNP->rsId<< " ";
 				 	 else if(pCBSNP->snpName!="") ofs << pCBSNP->snpName<< " ";
 				if(pCBSNP->bp!=(-1)) ofs<< pCBSNP->bp<< " ";
 				else
 				{
 					ofs << j<< " ";

 				}
 				if(pCBSNP->allele1!="") ofs << pCBSNP->allele1 << " ";
 				else if(pCBSNP->allele2!="") ofs << pCBSNP->allele2<< " ";
 				else
 				{
 					ofs << 0 << " ";
 				}
 				if(pCBSNP->allele2!="") ofs << pCBSNP->allele2;
 				  else if(pCBSNP->allele1!="") ofs << pCBSNP->allele1;
 				 else
 				 {
 					 ofs << 0 ;
 				 }
 				for(unsigned int i=0; i<pedInfo.size();++i)
 				{
 					if( !pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if homozygote1  00
 						ofs <<" "<< "1" <<" "<< "0" <<" "<< "0";
 					else if( pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if homozygote2 11
 						ofs <<" "<< "0" <<" "<< "0" <<" "<< "1";
 					else if( !pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if heterozygote
 						ofs <<" "<< "0" <<" "<< "1" <<" "<< "0";
 					else if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) //if missing
 						ofs <<" "<<CBSNP::_miss_gvalue << " "<<CBSNP::_miss_gvalue << " "<<CBSNP::_miss_gvalue;
 					else
 					{
 					error("problem in writing out "+change_int_into_string(j) +"th SNP.\n");
 					}

 				}
 			}
 		ofs<< "\n";
 		}
 	printLIN( "\t**->Impute gens file has been written out and saved as \""+pedFileName+"\". \n");
 	ofs.clear();	 ofs.close();


 return;}
 void CPPED::write_pedinfoFiles(vector<CBPED*> pedInfo, const string& outFileName)
 {
	//----------------------------------#--------------------------------------#
	 	// writing extra ped info files.
	 //----------------------------------#--------------------------------------#
	 	 ofstream ofs;
	 	string extraPedFileName=outFileName+"_pedinfo.txt";
		ofs.clear();	 ofs.close();
	 	ofs.open(extraPedFileName.c_str());
	 	if(!ofs)
	 			error("out file "+outFileName+" does not exits.\n");
	 	//
	 	ofs<<"famid" <<" ";
	 	ofs<<"indid" <<" ";
	 	ofs<<"matid" <<" ";
	 	ofs<<"patid" <<" ";
	 	ofs<<"sex"	 <<" ";
	 	ofs<<"phenotype"<<"\n";
	  	for(unsigned int i=0; i<pedInfo.size();++i)
	 	{
	 		CBPED* pCBPED =pedInfo[i];
	  		ofs<<pCBPED->famId<< " " <<pCBPED->indId << " "<<pCBPED->patId <<" "<<pCBPED->matId <<" ";
	 		if(pCBPED->sex) ofs<<"1" <<" ";
	 		else if(!pCBPED->sex && !pCBPED->miss_sex) ofs<<"2"<< " ";
	 		else if(pCBPED->miss_sex) ofs<<"0" <<" "; // 0 is missing
	 		if(pCBPED->pheno)
	 			ofs<<"2"<< "\n";
	 		else if(!pCBPED->pheno && !pCBPED->miss_pheno)
	 			ofs<<"1"<< "\n";
	 		else if(!pCBPED->pheno && pCBPED->miss_pheno)
	 		 	 			ofs<<"0"<< "\n";
	 	}

		ofs.clear();	 ofs.close();
	printLIN("\t**->One extra pedinfo.txt file has also been written out and saved as \""+extraPedFileName+\
		 			"\".\n");
 return ;}
 void CPPED::write_plink_raw_files(const vector<CBPED*>&pedInfo,const vector<CBSNP*>&genInfo, const string&outFileName)
 {
	 //----------------------------------#--------------------------------------#
	 	printLIN("*->Writing files in "+CBSNP::_out_format+" format: \n");
	 //----------------------------------#--------------------------------------#
	 //writing map file
	 CBSNP::write_plinkMapFile(genInfo, outFileName);
	 //
	 //writing raw file:
	 ofstream ofs; ofs.clear(); ofs.close();
	  string rawFileName =outFileName+".raw";
	  ofs.open(rawFileName.c_str());
	  CBPED* pCBPED =pedInfo[0];
	  CBSNP* pCBSNP =genInfo[0];
	  const bool _recodeAD=(CBSNP::_out_format=="plink-recodeAD");
	    ofs << "FID"<< " "<<"IID"<< " " <<"PAT"<< " "<<"MAT"<< " "<<"SEX"<< " ";
	  	  ofs <<"PHENOTYPE";
	  	for(unsigned int j =0; j<genInfo.size(); ++j)
	  	{
	  		pCBSNP =genInfo[j];
	  		//cout << "minallele: "<< pCBSNP->min_allele <<" ";
	  		if(pCBSNP->change_0123)
	  			CBSNP::convert_snp_genotype_into_0123(pCBSNP);
	  		//	CBSNP::convert_snp_genotype_into_0123(pCBSNP);
	  		//pCBSNP->change_0123=false;
	  		if(pCBSNP->quality)
	  		{
	  			if(pCBSNP->rsId!="")
	  				ofs<<" "<<pCBSNP->rsId;
	  			else
	  				ofs<<" "<<pCBSNP->snpName;
	  		//	cout <<  (pCBSNP->min_allele)<< endl;
	  			if(pCBSNP->min_allele!="")
	  				ofs<<"_"<<pCBSNP->min_allele<< " "; // minor allele
	  			else
	  				ofs<<"_"<<"0"<< " "; // minor allele
				if(_recodeAD && pCBSNP->rsId!="")
					ofs<<pCBSNP->rsId<<"_"<<"HET";
				else if(_recodeAD && pCBSNP->snpName!="")
					ofs<<pCBSNP->snpName<<"_"<<"HET";
	  		}//end of if quality loop
	  	} //end of for loop
	  	ofs<<"\n";//end of line
	 pCBSNP =genInfo[0];
	  for(unsigned int i=0; i<pedInfo.size();++i)
	  {
		  pCBPED =pedInfo[i];
		  if(pCBPED->quality)
		  {
			  if(pCBPED->famId!="") ofs<<pCBPED->famId<< " ";
			  else ofs << "famid"<<i+1<< " ";
			  if(pCBPED->indId!="")ofs<<pCBPED->indId << " ";
			  else ofs<< "indv"<<i+1<<" ";
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
			  {
				  ofs<<"-9"<< " ";
			  }
			  // here we have to write the genotypes according as minor and major allele coding 0, 1,2 and 3 3 will be NA
			  if(_recodeAD)
			  {
				  for(unsigned int j =0; j<genInfo.size(); ++j)
				  {
					  pCBSNP =genInfo[j];
					  if(pCBSNP->geno_0123[i]==3)
						  ofs<< " "<<"NA"<< " "<<"NA";

					  else
					  {
						 //--------------------
						  // check if firs allele is major or minor . If major then
						  // it is at it is . if the first allele is minor then we should
						  // exchange the order of major and minor when it comes
						  // for homozygotes in geno_0123[i]
						  if (pCBSNP->geno_0123[i]==1 || (pCBSNP->freq>0.5) )
							  ofs <<" "<< (pCBSNP->geno_0123[i]);
							else if(pCBSNP->freq<=0.5 )
							{
								if (pCBSNP->geno_0123[i]==0)
									ofs << " "<< "2";
								else if(pCBSNP->geno_0123[i]==2)
									ofs << " "<< "0";

							}
							 //--------------------
						 // ofs<<" "<< pCBSNP->geno_0123[i];
						  if(pCBSNP->geno_0123[i]==1)
							  ofs<<" "<<1;
						  else
							  ofs<<" "<<0;
					  }
				  }

			  }
			  else
			  {
				  for(unsigned int j =0; j<genInfo.size(); ++j)
				  {
					  pCBSNP =genInfo[j];
					  if(pCBSNP->geno_0123[i]==3)
						  ofs<< " "<<"NA";
					  else
					  {
						 // ofs<<" "<< pCBSNP->geno_0123[i];
						  //--------------------
						  // check if firs allele is major or minor . If major then
						  // it is at it is . if the first allele is minor then we should
						  // exchange the order of major and minor when it comes
						  // for homozygotes in geno_0123[i]
						  if (pCBSNP->geno_0123[i]==1 || (pCBSNP->freq>0.5) )
							  ofs <<" "<< (pCBSNP->geno_0123[i]);
						  else if(pCBSNP->freq<=0.5 )
						  {
							  if (pCBSNP->geno_0123[i]==0)
								  ofs << " "<< "2";
							  else if(pCBSNP->geno_0123[i]==2)
								  ofs << " "<< "0";

						  }
						  //--------------------
					  }
				  }

			  }
			  //phenotype
			  ofs<<"\n"; //end of the line;
			  // ++ppedInfo; //increament of pointer to individual
		  } // end of if qualtiy loop
	  }//end for pCBPED loop

	  //  ofs <<"test"; // debug

	  ofs.clear();
	  ofs.clear(); ofs.close();
	 printLIN("\t**->plink-raw formatted file has been written out and saved as \""+rawFileName+\
	  		 			"\".\n");
}
//xxxxxxxxxxxxxxxxxxxxx
 void CPPED::write_plink_raw_dose_files(const vector<CBPED*>&pedInfo,const vector<CBSNP*>&genInfo, const string&outFileName, const string & ref_allele)
  {
 	 //----------------------------------#--------------------------------------#
 	 	printLIN("*->Writing files in "+CBSNP::_out_format+" format: \n");
 	 //----------------------------------#--------------------------------------#
 	 //	cout << "test"<<endl;
 	 //writing map file
 	 CBSNP::write_plinkMapFile(genInfo, outFileName);
 	 //
 	 //writing raw file:
 	 ofstream ofs; ofs.clear(); ofs.close();
 	  string rawFileName =outFileName+"_dose.raw";
 	  ofs.open(rawFileName.c_str());
 	  CBPED* pCBPED =pedInfo[0];
 	  CBSNP* pCBSNP =genInfo[0];
 	  //const bool _recodeAD=(CBSNP::_out_format=="recodeAD-dose");
 	    ofs << "FID"<< " "<<"IID"<< " " <<"PAT"<< " "<<"MAT"<< " "<<"SEX"<< " ";
 	  	  ofs <<"PHENOTYPE";
 	  	  //cout << "genInfo.size() "<< genInfo.size()<<endl;
 	  	for(unsigned int j =0; j<genInfo.size(); ++j)
 	  	{
 	  		pCBSNP =genInfo[j];
 			CBSNP::convert_hardcalls_into_genoprob(pCBSNP);
 	  		if(pCBSNP->change_dose)
 	  			CBSNP::convert_snp_genotype_into_dose(pCBSNP,ref_allele);
 	  		//cout <<"pCBSNP->geno_dose.size() "<< pCBSNP->geno_dose.size() << endl;
 	  		if(pCBSNP->quality)
 	  		{
 				if(pCBSNP->rsId!="")
 	  			ofs<<" "<<pCBSNP->rsId;
 				else
 					ofs<<" "<<pCBSNP->snpName;
 				if(pCBSNP->dose_allele1!="")
 					ofs<<"_"<<pCBSNP->dose_allele1<< " "; // minor allele
 				else
 					ofs<<"_"<<"0"<< " "; // minor allele
 				//if(_recodeAD)
 				// ofs<<pCBSNP->rsId<<"_"<<"HET";
 	  		}//end of if quality loop
 	  	} //end of for loop
 	  	ofs<<"\n";//end of line
 	 pCBSNP =genInfo[0];

 	// cout << "\n pedInfo.size(): "<< pedInfo.size()<<endl;
 	  for(unsigned int i=0; i<pedInfo.size();++i)
 	  {
 		  pCBPED =pedInfo[i];
 		  //cout << pCBSNP->geno_dose.size()<<endl;
 		  if(pCBPED->quality)
 		  {
 			  if(pCBPED->famId!="") ofs<<pCBPED->famId<< " ";
 			  else ofs << "famid"<<i+1<< " ";
 			  if(pCBPED->indId!="")ofs<<pCBPED->indId << " ";
 			  else ofs<< "indv"<<i+1<<" ";
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
 			  {
 				  ofs<<"-9"<< " ";
 			  }
 			  // here we have to write the genotypes according as minor and major allele coding 0, 1,2 and 3 3 will be NA
 			  for(unsigned int j =0; j<genInfo.size(); ++j)
 			  {
 				  pCBSNP =genInfo[j];
 				  if(!(pCBSNP->geno_dose[i]+1.0))
 					  ofs<< " "<<"NA";
 				  else
 					  ofs<<" "<< pCBSNP->geno_dose[i];
 			  }
  			  //phenotype
 			  ofs<<"\n"; //end of the line;
 			  // ++ppedInfo; //increament of pointer to individual
 		  } // end of if qualtiy loop
 	  }//end for pCBPED loop

 	  //  ofs <<"test"; // debug

 	  ofs.clear();
 	  ofs.clear(); ofs.close();
 	 printLIN("\t**->plink-raw (recodeA-dose) formatted file has been written out and saved as \""+rawFileName+\
 	  		 			"\".\n");
 }
 //xxxxxxxxxxxxxxxxxxxxx


 void CPPED::displayIndividualSummary(const vector<CBPED*>& vPedInfo)
{
	 //
	cout <<left;
	cout <<"\n*->Individual summary:\n";
	// male and female Geschiste
	cout <<"\t"<< setw(coutWidth)<< 	"**->Total female: " << nFemale<<"\n";
	cout <<"\t"<<setw(coutWidth)<< 	"**->Total male: " << nMale<< "\n";
	cout <<"\t"<<setw(coutWidth)<< 	"**->Total undefined individual: " <<nNosex<< "\n";
	cout.flush()<<endl;
	//write to fcgene.out.txt file
	LIN << left;
	LIN <<						 	"\n*->Individual summary:\n";
	LIN <<"\t"<< setw(coutWidth)<< 		"**->Total female: " << nFemale<<"\n";
	LIN <<"\t"<<setw(coutWidth)<< 			"**->Total male: " << nMale<< "\n";
	LIN <<"\t"<<setw(coutWidth)<< 			"**->Total undefined individual: " <<nNosex<< "\n";
	LIN.flush()<<endl;


	// for affected and unaffected.
	cout <<"\t"<< setw(coutWidth)<< 	"**->Total cases: " << nAffected<<"\n";
	cout <<"\t"<<setw(coutWidth)<< 	"**->Total controls: " << nUnaffected<< "\n";
	cout <<"\t"<<setw(coutWidth)<< 	"**->Total undefined individuals: " <<nMiss_status<< "\n";
	cout.flush()<<endl;
	//write to fcgene.out.txt file
	LIN <<"\t"<< setw(coutWidth)<< 	"**->Total cases: " << nAffected<<"\n";
	LIN<<"\t"<<setw(coutWidth)<< 	"**->Total controls: " << nUnaffected<< "\n";
	LIN<<"\t"<<setw(coutWidth)<< 	"**->Total undefined individuals: " <<nMiss_status<< "\n";
	LIN.flush()<<endl;

return;}
// next


bool CPPED::openBinaryFile(string s, ifstream & BIT)
{

  BIT.open(s.c_str(), ios::in | ios::binary);

  // 1) Check for magic number
  // 2) else check for 0.99 SNP/Ind coding
  // 3) else print warning that file is too old
  bool bfile_SNP_major = false;
  bool v1_bfile = true;

  char ch[1];
  BIT.read(ch,1);
  bitset<8> b;
  b = ch[0];	  
  

  // If v1.00 file format
  // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
  if (   ( b[2] && b[3] && b[5] && b[6] ) && 
       ! ( b[0] || b[1] || b[4] || b[7] )    )
    {

     // Next number
     BIT.read(ch,1);
     b = ch[0];	  
     if (   ( b[0] && b[1] && b[3] && b[4] ) && 
          ! ( b[2] || b[5] || b[6] || b[7] )    )
      {
        // Read SNP/Ind major coding
        BIT.read(ch,1);
        b = ch[0];	  
        if ( b[0] ) bfile_SNP_major = true;
        else bfile_SNP_major = false;

        if (bfile_SNP_major) 
  	  printLIN("\t**->Detected that binary PED file is v1.00 SNP-major mode\n");
        else
	  printLIN("**->\tDetected that binary PED file is v1.00 individual-major mode\n");

      } else v1_bfile = false;
      
    } else v1_bfile = false;


  // Reset file if < v1
  if ( ! v1_bfile ) 
   {
    printLIN("Warning, old BED file <v1.00 : will try to recover...\n");
    printLIN("  but you should --make-bed from PED )\n");
    BIT.close();
    BIT.clear();
    BIT.open(s.c_str(), ios::in | ios::binary);
    BIT.read(ch,1);
    b = ch[0];	  
  }

  // If 0.99 file format
  if ( (!v1_bfile) && ( b[1] || b[2] || b[3] || b[4] || b[5] || b[6] || b[7] ) )
    {
      printLIN("\n *** Possible problem: guessing that BED is < v0.99      *** \n");
      printLIN(" *** High chance of data corruption, spurious results    *** \n");
      printLIN(" *** Unles you are _sure_ this really is an old BED file *** \n");
      printLIN(" *** you should recreate PED -> BED                      *** \n\n");

      bfile_SNP_major = false;
      BIT.close(); 
      BIT.clear();
      BIT.open(s.c_str(), ios::in | ios::binary);
    }
  else if ( ! v1_bfile ) 
    {
      if ( b[0] ) bfile_SNP_major = true;
      else bfile_SNP_major = false;

      printLIN("Binary PED file is v0.99\n");
      
      if (bfile_SNP_major) 
	printLIN("Detected that binary PED file is in SNP-major mode\n");
      else
	printLIN("Detected that binary PED file is in individual-major mode\n");
    }

 return bfile_SNP_major;

}

// read binary data
void CPPED::lese_bim_file(vector<CBSNP*>&genVec, const string & fileName){
	printLIN("*->Reading bim file: "+fileName+":\n");
	static long int countIndiv=0;
	string 	line 		="";
	string 	buffer 		="";
	//unsigned int nCols  =0;
	vector<string> tokens;
	rFile file; file.close();file.clear();
	file.open(fileName.c_str(),ios::in);
	if(!file){
		file.setstate(ios_base::failbit);
		file.clear();
		error("File \""+fileName+"\" either does not exit or can not be opened.\n");
	}
	while(getline(file,line,'\n'))
	{
		if(line=="") // this blank line
			continue;
		if(line[0]=='#') // this is comment line
			continue;
		stringstream line_parser(line);
		if(line_parser.good())
		{
			++countIndiv;
			while(line_parser>>buffer)
				tokens.push_back(buffer);
			//cout << tokens.size()<<endl;
			if(tokens.size()==6)
			{
				CBSNP* pCBSNP=new CBSNP;
				pCBSNP->nchr		=tokens[0];
				pCBSNP->rsId		=tokens[1];
				pCBSNP->snpName		=tokens[1];
				pCBSNP->cm_pos		=atof(tokens[2].c_str());
				pCBSNP->bp			=atoi(tokens[3].c_str());
				if(tokens[4]!="0")
					pCBSNP->allele1	=tokens[4];
				pCBSNP->allele2		=tokens[5];
				genVec.push_back(pCBSNP);
			}
			else{
				error("The " + change_int_into_string( countIndiv)+"th row (excluding comments and blank lines)"+
						"of your fam file does not contain 6 columns. Please delete or correct this row.\n");
			}
		}
		tokens.clear();
		tokens.resize(0);
	}//end of while getline
	file.close();file.clear();
	printLIN("\t**->Total no of SNPs read: "+change_int_into_string((int)genVec.size())+ " \n");
//----------------------------------------------------------------------------//

}
void CPPED::lese_bed_file(vector<CBPED*>&pedVec,vector<CBSNP*>&genVec, const string & fileName)
{
	printLIN("*->Reading bed file: "+fileName+":\n");
    ifstream BIT;
    bool bfile_SNP_major = openBinaryFile(fileName, BIT);
    unsigned int nSNP=genVec.size();
    unsigned int nInd=pedVec.size();
      
      //////////////////////////////
      // Allocate space for SNPs
     // vector<CBSNP*> SNP ;
      if (bfile_SNP_major)
	{
    	  CBSNP * newlocus = genVec[0];
	  for (unsigned int i=0; i<nSNP; i++)
	    {
		  newlocus = genVec[i];
	      newlocus->geno1.resize( nInd);
	      newlocus->geno2.resize(nInd);
	      newlocus->aOrder.resize(nInd);

	     // SNP.push_back(newlocus);
	    }     
	}
      
     
      
      ///////////////////////////
      // SNP-major mode
      
      if (bfile_SNP_major)
	{
	  
	  CBSNP * snp;
	  // Outer loop for SNPs
	  unsigned int s=0;
	  while (s<nSNP) // for all SNPs
	    {
	      // Do we want to include this SNP?
		  snp = genVec[s ];
	      // Inner loop for individuals
		  int indx = 0;
	      int ss = nInd; //sample.size();
	      while ( indx < ss )
	      {
		  
	    	  bitset<8> b;
	    	  char ch[1];
	    	  BIT.read(ch,1);
	    	//  cout <<boolalpha<<(bool) BIT <<endl;
	    	  if (!BIT)
	    		  error("Problem with the BED file.. Probably the associated FAM/BIM files are not correct.\n");
	    	  b = ch[0];
	    	  //		}
	    	  int c=0;
		  while (c<7 && indx < ss ) 
		  {
		      if (snp)
			{
			  //  		      snp->one.push_back( b[c++] );
			  //  		      snp->two.push_back( b[c++] );
			  
			snp->geno1[indx] 	= b[c++];
			snp->geno2[indx] 	= b[c++];
			snp->aOrder[indx]	=true;
			 
			  
			}
		      else
			{
			  c+=2;
			}
		      // 		  ++person;		  
		      ++indx;
		    }
		  
		}	  
	      
	      // next SNP
	  s++;
	    }
	  // remove indloop later
	  	s=0;
	 	for(unsigned int j=0; j<genVec.size();++j)
	  {
          snp =genVec[j];

			 if(snp->allele1=="")
			 {
				 //I have assumed all the time that first allele is not zero
				 //now if there is no minor allele , then first allele can be zero.
				 // so I change second allele to first so that the  first is never zero.
				 for(unsigned int ii=0; ii<nInd; ++ii)
				 {

	  					if(!snp->geno1[ii]&&!snp->geno2[ii])
	  					{
	  						snp->geno1[ii]=true;
	  						snp->geno2[ii]=true;
	  					}
	  					else if(snp->geno1[ii]&&snp->geno2[ii])
	  					{
	  						snp->geno1[ii]=false;
	  						snp->geno2[ii]=false;
	  					}

				 }

			  snp->allele1=snp->allele2;
			  snp->allele2="";
			 }
	  }


	}
      else{
    	  string _msg ="File couldn't be read. It has no SNP Major mode.\n";
    	  error(_msg);
      }
      //test
      /*
      for(unsigned int i=0; i<nInd; ++i)
    	  {
    		for(unsigned int j=0; j<genVec.size();++j)
    		  {

    				CBSNP* snp =genVec[j];
    				cout<< snp->geno1[i]<<" "<<snp->geno2[i]<<"\t";
    		  }
    		cout<<endl;
    	  }
      */
     /*
      ////////////////////////////////////
      // Individual-major mode
      
      else
	{
	  
	  // Outer loop for individuals
	  vector<Individual*>::iterator person = sample.begin();
	  while ( person != sample.end() )
	    {
	      
	      // Inner loop for SNPs
	      int s=0;
	      while (s<locus.size()) // for all SNPs
		{
		  
		  char ch[1];
		  BIT.read(ch,1);
		  if (!BIT) error("Problem with the BED file... has the FAM/BIM file been changed?\n");
		  
		  bitset<8> b;
		  b = ch[0];	  
		  
		  int c=0;
		  
		  while (c<7 && s<locus.size())
		    {
		      if ( include[s] > -1 )
			{
			  (*person)->one[ include[s] ] = b[c++];
			  (*person)->two[ include[s] ] = b[c++];	      
			}
		      else
			{
			  c+=2;
			}
		      s++;
		    }	 	  
		}
	      
	      person++;
	    }
	  
	  // Set file mode
	  par::SNP_major = false;
	}
      
      // Check that we got what we expected
      
      char ch[1];
      BIT.read(ch,1);
      if (BIT) 
	error("Problem with the BED file... has the FAM/BIM file been changed?\n");
            
      BIT.clear();
      BIT.close();
      
      
    
  // Free any buffer memory used
  //   if ( par::fast_binary )
  //     memblock.clear();
  
  
  ////////////////////////////////////////
  // If need be, now prune the MAP file 
  // i.e. if --chr or --from/--to were used
  
  if ( (!par::plink) && (!par::run_chr==0) )
    {
      
      vector<Locus*> l0(0);
      for(int l=0; l < locus.size(); l++)
	{    
	  if ( !( l < par::run_start || l > par::run_end ) ) 
	    l0.push_back(locus[l]);
	  else
	    delete locus[l];
	}
      locus.clear();
      locus = l0;
    }
  */
}


void CPPED::write_BITFILE(const vector<CBSNP*>&genVec,
		const vector<CBPED*>&pedVec, const string &outFile)
{
	string oFile =outFile+".fam";
	printLIN( "*->Writing fam file \"" + oFile+"\":\n");
	ofstream ofs(oFile.c_str(), ios::out);
  // For each individual
	const CBPED* pCBPED =pedVec[0];
	const CBSNP* pCBSNP =genVec[0];
	 for(unsigned int i=0; i<(pedVec.size());++i)
	 {
		 pCBPED =pedVec[i];
		 if(pCBPED->quality)
		 {
			 //#####
			 if(pCBPED->famId!="") ofs<<pCBPED->famId<< " ";
			 else ofs << "famid"<<i<< " ";
			//
			 if(pCBPED->indId!="")ofs<<pCBPED->indId << " ";
			 else ofs<< "indv"<<i<<" ";
			 //
			 if(pCBPED->patId!="")ofs<<pCBPED->patId << " ";
			 else ofs<< 0<<" ";
			 //
			 if(pCBPED->matId!="") ofs <<pCBPED->matId <<" ";
			 else ofs<< 0 << " ";
			 //
			 if(pCBPED->sex) ofs<<"1"<< " "; // 1 is male
			 else if(!pCBPED->sex && !pCBPED->miss_sex) ofs<<"2"<<" "; // 2 is female
			 else if(pCBPED->miss_sex) ofs<<"0"<< " "; // 0 is missing
			 //
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
	 printLIN("\t**->fam file is saved as \""+ oFile+"\".\n");

//------------------------------------------------------------------------------//
	 //Before writing Bim file and bed file
	 // We have to calculate minor and major allele and save them as first and
	 // second allele
	 CBPED::calculate_maf(genVec,pedVec);

	 string bimFile=outFile+".bim";
	 printLIN("*->Writing bim file \""+bimFile+"\": \n");
	 ofs.clear();ofs.close();
	 ofs.open(bimFile.c_str(),ios::out);
	 if(!ofs)
		 error("out file "+bimFile+" does not exits.\n");
	 for(unsigned int i=0; i<genVec.size();++i)
	 {
		 CBSNP* pCBSNP	=genVec[i];
		 if(pCBSNP->quality)
		 {
			  ofs<<pCBSNP->nchr;
			 //if(pCBSNP->snpName!="")
			 if(pCBSNP->rsId!="") ofs <<"\t"<< pCBSNP->rsId;
			 else if(pCBSNP->snpName!="")
				 ofs <<"\t"<<pCBSNP->snpName;
			 else
				 ofs<<"\t"<<"snp"<<i;

			 if(pCBSNP->bp!=-1)
				 ofs <<"\t"<<pCBSNP->bp;
			 else
				 ofs <<"\t"<<0;
			 if(pCBSNP->cm_pos!=(-1.0))
				 ofs <<"\t"<<pCBSNP->cm_pos;
			 else ofs  <<"\t"<<0.0 ;
			 //
			 if(pCBSNP->min_allele!="")
				 ofs <<"\t"<< pCBSNP->min_allele;
			  else ofs  <<"\t"<<"0";
			 //
			 ofs <<"\t"<< pCBSNP->maj_allele<<"\n";

		 } //end of if quality loop
	 }

	 ofs.clear();ofs.close();
	 printLIN("\t**->bim file is saved as \""+ bimFile+"\". \n");
	 //------------------------------------------------------------------------------//



  //////////////////////////////////////
  // Save genotype data in BITFILE format
  string bedFile =outFile+".bed";
  printLIN("*->Writing genotype bitfile \""+bedFile+"\":\n");
  ofstream BIT;
  BIT.clear();  BIT.close();
  BIT.open(bedFile.c_str(), ios::out | ios::binary);
//start here
  bitset<8> b;
    char ch[1];


    // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file

    b.reset();
    b.set(2);  b.set(3);  b.set(5);  b.set(6);
    ch[0] = (char)b.to_ulong();
    BIT.write(ch,1);

    b.reset();
    b.set(0);  b.set(1);  b.set(3);  b.set(4);
    ch[0] = (char)b.to_ulong();
    BIT.write(ch,1);


    // BIT represents status of SNP-major (true) or Ind-major (false)

    b.reset();
   // if (par::out_SNP_major)
      b.set(0);
    ch[0] = (char)b.to_ulong();
    BIT.write(ch,1);

   const unsigned int nSNPs =genVec.size();
 //  const unsigned int nINDs =genVec.size();
       	 unsigned int j		=0;
       	 unsigned int ii		=0;

      // Outer loop over SNPs
      while (j<nSNPs )
	{
    	 //cout <<"hi from SNP: "<<j<<endl;
    	  pCBSNP =genVec[j];
    	  if(pCBSNP->quality)
    	  {
    		 // vector<bool>::const_iterator i1 = pCBSNP->geno1.begin();
    		 // vector<bool>::const_iterator i2 = pCBSNP->geno2.begin();
    		  bool first, second;
    		  // check if allele1 and minor allele are same or not . If not then we have to change it
    		  // but major comes often and it can be that allele1 is mostly equal to major allele
    		  //so better to change if allele1 is not equal to major allele and down in the loop
    		  //change  the bool values  in geno1[j] and geno2[j] of plink representation.
    		  //
    		  while( ii <pCBSNP->geno1.size())
    		  {
    			//  cout<<"hi from indiv: "<<ii<<endl;
    			  bitset<8> b;
    			 // cout << b <<endl;
    			  b.reset();
    			  int c=0;
    			  while (c<8 && (ii <pCBSNP->geno1.size()) )
    			  {
    				  first 	=pCBSNP->geno1[ii];
    				  second	 =pCBSNP->geno2[ii];

    				  if(pCBSNP->allele1==pCBSNP->min_allele)
    				  {
    					  if ( first ) b.set(c);
    				  	  c++;
    				  	  if ( second ) b.set(c);
    				  	  c++;
    				  }//this above is from plink but in my case this is true only when if
    				  else
    				  {
    					  //check more
    					//cout <<first <<", "<< second <<" ii: ";
    					if ( !first &&!second )
    						  b.set(c);
    					//  cout <<ii<<": "<<"c: "<<c<<", b: "<< b<<endl;
    					  if ( first &&!second )
    						  b.set(c);
       					 // cout <<ii<<": "<<"c: "<<c<<", b: "<< b<<endl;

    					  c++;
    					 if ( !second &&!first  )
    						  b.set(c);
      					 // cout <<ii<<": "<<"c: "<<c<<", b: "<< b<<endl;
      					  if (second &&!first  )
    						  b.set(c);
       					 // cout <<ii<<": "<<"c: "<<c<<", b: "<< b<<endl;
       					 // cout<<"\n------------------------------------\n";
      					  c++;
    				  }

    				  ++ii;
    			  }
    			  //-------------------------------------//
    			  char ch[1];
    			  ch[0] = (char)b.to_ulong();
    			  BIT.write(ch,1);

    		  }

    		  ii=0; // assign ii to zero again

    	  }
    	  j++;
	}
     j=0;

  BIT.close();  BIT.clear();

  printLIN("\t**->Genotype bit file is saved as \""+ bedFile+"\".\n");
  //------------------------------------------------------------------------------//
}




