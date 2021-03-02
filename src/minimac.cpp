/*
 * minimac.cpp
 *  Created on: 22.01.2013
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */
#include "minimac.h"
#include "helper.h"
extern ofstream LIN;
CMINIMACSNP::~CMINIMACSNP()
{

}

CMINIMACPED::~CMINIMACPED()
	{

	}

void CMINIMACSNP::lese_minimac_info_file(const string& fileName, vector<CBSNP*>&genInfo, const float& rsq_value, const float& maf_value){
	//cout <<"hi I am from lese_minimac_info_file\n ";
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
	 				 bool tfValue =(spaltenZahl==13) ;
	 				if(!tfValue)	// Fall mit 6 Zeilen.
	 				{
	 					string msg="\t**->Problem on the " + change_int_into_string(nrow)+ "th row.\n";
	 					printLIN(msg);
	 					msg="There must be 13 columns in each line of file: "+fileName+" but there  are "+change_int_into_string(spaltenZahl)+" columns.\n";
	 					//delete pCBPED;
	 					 error(msg);

	 				}else
	 				{
	 					if(nrow==1)
	 					{
	 						tfValue =(cZeile[0]!="SNP" || cZeile[1]!="Al1" ||cZeile[2]!="Al2"|| cZeile[3]!="Freq1") ||(cZeile[4]!="MAF" ||cZeile[5]!="AvgCall" ||cZeile[6]!="Rsq");
	 						bool tfvalue1=cZeile[7]!="Genotyped" ||cZeile[8]!="LooRsq"||cZeile[9]!="EmpR"||cZeile[10]!="EmpRsq"||cZeile[11]!="Dose1" ||cZeile[7]!="Dose2";
	 						if(tfValue&&tfvalue1)
	 						{
	 							emsg =" First line of  file \""+fileName+\
	 									"\" should contain \"SNP Al1 Al2 Freq1 MAF AvgCall Rsq Genotyped LooRsq EmpR EmpRsq Dose1 Dose2\"  as column specifiers.\n";
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

}
// next
void CMINIMACPED::lese_minimac_probFile( vector<CBPED*>& pedInfo, const vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName, const float & thresh){
	//cout <<"hi I am from lese_minimac_prob_file\n ";
	printLIN("*->Reading file \""+fileName+"\": \n");
	string myString					="";
	string temp_value				="";
	string emsg						="";
	bool ende						= false;
	unsigned int ccol				=0;
	unsigned int ncol				=0;
	unsigned int nrow				=0;
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
					else
					{
						if( (temp_float>1.000001)  )
						{
							emsg="Sum of genotype distribution should not be greater than 1. But this is the case in SNP "+pCBSNP->snpName+" and individual " +temp_value+" \n";
							error(emsg);
						}
						else
						{
							emsg=" first and second allele of  "+pCBSNP->snpName+"  is not given. please update alleles information first. \n";
							error(emsg);
						}

					}// end of else

				}//end of for loop


			} //if(cZeile.size()!=0)
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


