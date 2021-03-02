/*
 * readimpute.cpp
 *
 *  Created on: 03.01.2012
  *  Created on: Dec 29, 2011
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */

#include"readimpute.h"
extern ofstream LIN;
CISNP::~CISNP()
{
	/*for (unsigned int i=0;i<genInfo.size();++i)
	{
		delete genInfo[i];
		
	}
	genInfo.clear();
	*/
 	//cout << "\n Good bye from CISNP destructor.\n";
 	//cout << "\n Good bye from CBSNP destructor.\n";
}

CIPED::~CIPED()
	{

		/*for(vector<CBPED*>::iterator it=vPed.begin();it!=vPed.end();++it)
		{
			delete *it;
			
		}
		vPed.clear();
		*/
		//cout << "\n I am a CIPED destructor\n";
	}


void CISNP::gensFileLeser(vector<CBSNP*>& genInfo,const string& gensFileName)
{

	printLIN("*->Reading file \""+ gensFileName+ "\": \n");

	IFS file(gensFileName.c_str(),"r");
	if(!file)
	{
		string msg= "File \""+ gensFileName+ "\" either does not exits or could not be opend.\n"	;
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
	bool	assign_nInd			=true;

	//unsigned int nPerson	=0;

	while(!_ifs_eof(file))
	{
		if(_ifs_eof(file))
		{
			zeilenEnde=true;
			if(nrow==0)
			{
				string msg=" There are no elements in "+gensFileName +" file.\n";
				error(msg);
			}
			else
			{
				cout << "total no of SNPs read: "<< nrow << endl; // only for test
				nrow=0;
			}
			break;
		}
		zeilenEnde=false;
		// define  pointer of CBSNP

		while(!zeilenEnde && !_ifs_eof(file))
		{  // do all works till the end of line and then break
			if(zeilenEnde || _ifs_eof(file))
				break;
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
					string msg= " There is only one column in the first row of "+gensFileName+" file. \n" ;
					error(msg);
				}
				else
				break; // only in the first column
			}
			//define pointer
			CBSNP* pCBSNP=new CBSNP;
			// read family id
			if(!elem1.empty())
			{
				if(ccol!=0 && ncol==1 ) // this condition must always be true ;
				{
					if(pCBSNP->snpName=="")
					pCBSNP->snpName		=elem1;
					else
					{
						string msg= "problem in saving snpName from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
						error(msg);
					}
				}
				else
				{
						string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
						error(msg);
				}
			}
			else
			{
				string msg= "String is empty! problem in reading  SNP ID  from "+gensFileName+".\n";
				error(msg);
			}
			// for the second column
			if(zeilenEnde || _ifs_eof(file))
			{
				string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
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
				string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
				error(msg);
			}
			else
			{
				if(ccol!=0 && ncol==2 ) // this condition must always be true ;
				{
					if(pCBSNP->rsId=="")
						pCBSNP->rsId		=elem1;
					else
					{
						string msg= "problem in saving rs id from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
						error(msg);
					}
				}
				else
				{
						string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
						error(msg);
				}

			//   This case is true except last column
				if(zeilenEnde || _ifs_eof(file))
				{		string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow +1)+"th row"+ "of   file \""+gensFileName+"\".\n";
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
				string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
				error(msg);
			}
			else
			{
				if(ccol!=0 && ncol==3 ) // this condition must always be true ;
				{
					if(pCBSNP->bp==-1.0)
					pCBSNP->bp		=atoi(elem1.c_str());
					else
					{
					string msg= "problem in saving  base pair position from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
					error(msg);
					}
				}
				else
				{
					string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
					error(msg);
				}

					//   This case is true except last column
				if(zeilenEnde || _ifs_eof(file))
				{		string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
						error(msg);
				}


			}
			// start of reading fourth column
			elem1="";
			elem="";
			ccol=0;
			ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
			elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
			if(elem1.empty())
				{
					string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
					error(msg);
				}
			else
			{
				if(ccol!=0 && ncol==4 ) // this condition must always be true ;
				{
					//	cout << "pCBSNP->allele1: "<< pCBSNP->allele1 <<endl; //only for test
					if(pCBSNP->allele1=="")
					{
						pCBSNP->allele1		=elem1;
						//cout << "element value: "<< elem1 <<endl; // test
					}
					else
					{
							string msg= "problem in saving first allele from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
							error(msg);
					}
				}
				else
				{
					string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
					error(msg);
				}

					//   This case is true except last column
				if(zeilenEnde || _ifs_eof(file))
				{
					string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
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
				string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
				error(msg);
			}
			else
			{
				if(ccol!=0 && ncol==5 ) // this condition must always be true ;
				{
					if(pCBSNP->allele2=="")
					{
						pCBSNP->allele2=elem1;
					}
					else
					{
						string msg= "problem in saving second allele from "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
						error(msg);
					}
				}
				else
				{
					string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
					error(msg);
				}

			//   This case is true except last column
				if(zeilenEnde || _ifs_eof(file))
				{		string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
					error(msg);
				}
			}

				// end of snp info
				//-------------------------------------------------//
				//sart of genotype reading
				unsigned  int nPerson	=0;
				//unsigned long int nSNPs			=genInfo.size();
				string first 			=""; // geno1
				string second			=""	;// geno2
				string  third 			="";
				//ofstream test;
			while(true)
			{
				// start of genotype reading loop
				if(zeilenEnde) break;

					elem1="";
					elem="";
					ccol=0;
					//if(nPerson==nSNPs) break;
					++nPerson;

					//-----------------------------------------------//
					// start of reading genotypes
					//reading first genotype
						ccol=CBSNP::stringLeser(file , elem, zeilenEnde);
						elem1=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
						if(elem1.empty())
						{
							string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row of   file \""+gensFileName+"\".\n";
							error(msg);
						}
						else
						{
							if(ccol!=0 && ncol== (int)(3*(nPerson+1)) ) // this condition must always be true ;
							{
								//first=atof(elem1.c_str());
								pCBSNP->pgeno1.push_back(atof(elem1.c_str()));


							}
							else
							{
								string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
								error(msg);
							}

							//   This case is true except last column
							if(zeilenEnde || _ifs_eof(file))
							{
								string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
								error(msg);
							}
						}

						//reading second genotype
						elem1	="";
						elem	="";
						ccol	=0;
						ccol	=CBSNP::stringLeser(file , elem, zeilenEnde);
						elem1	=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
						if(elem1.empty())
						{
							string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row of   file \""+gensFileName+"\".\n";
							error(msg);
						}
						else
						{
							if(ccol!=0 && ncol==(int)(3*nPerson+4) ) // this condition must always be true ;
							{
								//first=elem1;
								pCBSNP->pgeno2.push_back(atof(elem1.c_str()));
							}
							else
							{
								string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
								error(msg);
							}

							//   This case is true except last column
							if(zeilenEnde || _ifs_eof(file))
							{
								string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
								error(msg);
							}
						}

						// reading third genotype
						elem1	="";
						elem	="";
						ccol	=0;
						ccol	=CBSNP::stringLeser(file , elem, zeilenEnde);
						elem1	=CBSNP::stringAssigntoRightPlace(zeilenEnde, ccol,  nrow,  ncol, elem);
						if(elem1.empty())
						{
							string msg= " problem in reading "+change_int_into_string(ncol)+"th column in the line of "+change_int_into_string(nrow+1)+"th row"+ "of   file \""+gensFileName+"\".\n";
							error(msg);
						}
						else
						{
							if(ccol!=0 && ncol==(int)(3*nPerson+5) ) // this condition must always be true ;
							{
								//third=elem1;
								pCBSNP->pgeno3.push_back(atof(elem1.c_str()));
							}
							else
							{
								string msg= "problem in reading "+change_int_into_string(ncol)+ "th column in the "+change_int_into_string(nrow+1)+ " line of   file \""+gensFileName+"\". \n";
								error(msg);
							}


						//   This case is true except last column
							//cout << nPerson << " " <<nSNPs << endl; // test
							//if(nPerson!=nSNPs)
							//{
								//cout << nPerson << " " <<nSNPs << endl; // test
								//if(zeilenEnde || _ifs_eof(file))
								//{
									//string msg= "sNot enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nrow+1)+"th row of   file \""+gensFileName+"\".\n";
									//error(msg);
									//break;
								//}
							//}

						}



					if(zeilenEnde) break;
			} // end of genotyping while loop
			//--------------------------------------------------//
			// end of genotype reading
				//cout << "nPerson : "<<nPerson <<endl; // only for test
					if(assign_nInd)
						nInd=nPerson;
						assign_nInd=false;
						//cout << nInd<<endl;
				if(nInd!=nPerson)
				{
					string msg;
					msg=	"problem with reading individuals on line "+ change_int_into_string(nrow) +\
							"th line. Only "+change_int_into_string(nPerson) +" individuals were read.\n";
					error(msg);
				}
				nPerson		=0;
				//cout << pCBSNP->snpName <<endl; // only for test
				genInfo.push_back(pCBSNP);

		if(!zeilenEnde)
		{
			error("problem  with reading the end of file in "+ change_int_into_string(nrow) + "th line. \n");
		}
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
				string msg=" There are no elements in "+gensFileName +" file.\n";
				error(msg);
			}
			else
			{
				printLIN("\t**->Total no of SNPs read: "+ change_int_into_string(genInfo.size())+" \n");
				nrow=0;
			}
			nrow=0;
			break;
		}


	}

	_ifs_close(file);

//cout <<"pedInfo.size()"<<pedInfo.size();

			/*	for(unsigned int j=0;j<genInfo.size();j++)
				{
					CBSNP* pCBSNP	=genInfo[j];
					cout << pCBSNP->snpName <<" " <<pCBSNP->rsId<<" "<< pCBSNP->bp << " " <<pCBSNP->allele1 <<" " <<pCBSNP->allele2 << " " ;
					unsigned int size=pCBSNP->pgeno1.size();
					for(unsigned int i=0; i<size;++i)
					{
						cout <<  pCBSNP->pgeno1[i] << " "<<  pCBSNP->pgeno2[i] << " "<<  pCBSNP->pgeno3[i] << " ";

					}
				cout <<endl;

				}

			*/


return ;}
//
void CISNP::lese_imputed_info_score_file(vector<CBSNP*>&genInfo,vector<string>&cZeile, const double& info_thresh, const double& maf_thresh, const string& fileName)
{
//-----------------------------------------------------------------//
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
    	 //unsigned long int	snpAnzahl		=0;
    	 string emsg						="";
    	 long int temp_bp					=0;
    	 double  temp_info					=0.0;
    	 double  temp_freq					=0.0;

    	 //int  count_low_info				=0;
    	 try
    	 {
    		IFS pf(fileName.c_str(),"r");
			if(!pf)
			{
				string msg= fileName+ " either does not exits or could not be opend.\n"	;
				error(msg);
			}

    		vector<CBSNP*>::iterator it_CBSNP=genInfo.begin();
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
    							cout << "error while reading file "+fileName+"\n";
    							exit(1);
    						}
							//cout << myString << " "; //debug 
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
    				// snpAnzahl=genInfo.size();
    				 int spaltenZahl=cZeile.size();
	 				 //-----------------------------------------------------------------//
					 //cout << spaltenZahl<<endl;
	 				 bool tfValue =(spaltenZahl<=12 )||(spaltenZahl>=7 ) ;
	 				 if(!tfValue)	// Fall mit 6 Zeilen.
	 				 {
	 					 emsg="\t**->Problem on the " + change_int_into_string(nrow)+ "th row.\n";
	 					printLIN(emsg);
	 					emsg="There must be 7 or 12 columns in each line of file: "+fileName+" but there  are "+change_int_into_string(spaltenZahl)+" columns.\n";
	 					//delete pCBPED;
	 					 error(emsg);

	 				}else
	 				{
	 					if(nrow==1)
	 					{
	 						tfValue =(cZeile[0]!="snp_id" || cZeile[1]!="rs_id" ||cZeile[2]!="position"|| cZeile[3]!="a0"|| cZeile[4]!="a1"|| cZeile[5]!="exp_freq_a1") ||(cZeile[6]!="info" ||cZeile[7]!="certainty"); // || cZeile[5]!="certainty type");
	 						for (int i =0;i<(int)(cZeile.size());i++){
								cout <<cZeile[i]<<endl;
							}
							if(tfValue)
	 						{
	 							 emsg =" First line of  file \""+fileName+\
	 									"\" should contain \" snp_id      rs_id position exp_freq_a1  info certainty \" etc as column specifiers.\n";
	 							error(emsg);
	 						}

	 					}else
	 					{
	 						if((*it_CBSNP)->snpName!="" && cZeile[0]!="")
	 						{
	 							if((*it_CBSNP)->snpName!=cZeile[0]){
		 							emsg ="on the"+change_int_into_string(nrow)+"th row of  file  \""+fileName+\
		 									"\" SNP Name"+cZeile[0]+ " does not match with SNP Name "+(*it_CBSNP)->snpName+ " contained in original file.\n";
		 							error(emsg);

	 							}
	 							else
	 								(*it_CBSNP)->snpName=cZeile[0];

	 						}

	 						// check for rs id
	 						if((*it_CBSNP)->rsId!="" && cZeile[1]!="")
	 						{
	 							if((*it_CBSNP)->rsId!=cZeile[1])
	 							{
		 							emsg ="on the"+change_int_into_string(nrow)+"th row of  file  \""+fileName+\
		 									"\" rsId Name"+cZeile[1]+ " does not match with rs id "+(*it_CBSNP)->rsId+ " contained in original file.\n";
		 							error(emsg);

	 							}
	 							else
	 								(*it_CBSNP)->rsId=cZeile[1];

	 						}
							// for bp
	 						// check for rs bp
	 						if((*it_CBSNP)->bp!=-1 && cZeile[2]!="")
	 						{
	 								temp_bp=atoi(cZeile[2].c_str());
	 							if((*it_CBSNP)->bp!=temp_bp)
	 							{
		 							emsg ="on the"+change_int_into_string(nrow)+"th row of  file  \""+fileName+\
		 									"\" Base pair position "+cZeile[2]+ " does not match with base pair position "+change_int_into_string((*it_CBSNP)->bp)+ "contained in original file.\n";
		 							error(emsg);

	 							}
	 							else
	 								(*it_CBSNP)->bp= temp_bp;

	 						}
	 						//
	 						//for freq1
	 						if(cZeile[5]!="")
	 						{
	 								temp_freq			=atof(cZeile[5].c_str());
	 								(*it_CBSNP)->maf	= min((1.0-temp_freq),temp_freq);
	 								(*it_CBSNP)->freq	=max((1.0-temp_freq),temp_freq);

	 						}
	 						//
	 						//for freq1
	 						if(cZeile[5]!="")
	 						{
	 							temp_info	=atof(cZeile[5].c_str());
	 							//cout <<temp_freq<< " "<< info_thresh << " "<<( temp_freq<=info_thresh)<< endl; // debug
	 							/*
	 							 * if(temp_freq<=info_thresh)
	 							{
	 								(*it_CBSNP)->quality=false;
	 								//++count_low_info;
	 							}
	 							*/

	 						}
	 						//
	 						//cout <<temp_info<< " "<< info_thresh << " "<<( temp_info<=info_thresh)<< endl; // debug
	 						//cout <<temp_freq<<" "<<(*it_CBSNP)->maf<< " "<< maf_thresh << " "<<( (*it_CBSNP)->maf<=maf_thresh)<< endl; // debug

 							if(temp_info<info_thresh|| (*it_CBSNP)->maf<maf_thresh)
 							{
 								(*it_CBSNP)->quality=false;
 								//cout <<temp_info<< " "<< info_thresh << " "<<( temp_info<=info_thresh)<< endl; // debug
 								//cout <<temp_freq<<" "<<(*it_CBSNP)->maf<< " "<< maf_thresh << " "<<( (*it_CBSNP)->maf<=maf_thresh)<< endl; // debug

 								//++count_low_info;
 							}


	 						++it_CBSNP;

    					}

	 				}
    			 }
    				else continue;

	 				 //-------------------------------------------------------------//
	 				 /*
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
    			*/
    			ncol=0;
    			 cZeile.clear();
    			 cZeile.resize(0);
    			 if(_ifs_eof(pf))break;
    			} // end of outer while loop
    		_ifs_close(pf);
    		nrow=0;
    		nLines=0;
    		//cout << "count_low_info: "<< count_low_info <<endl;
    		// total lines
    		//string pmsg="\t**->Total number of successfully read SNPs: ";
    		//printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
    		//pmsg="\t**->Total number of successfully read individuals: ";
    		//display_ped_summary();
    		//printWithSpace(pmsg,  change_int_into_string(pedInfo.size()),coutWidth);
    	 }catch(bad_alloc& memoryAllocationException)
    	 {
    		 cout<<"\n Please inform the programmer that you could not read mlgeno file. Probbaly he will give you a solution. \n ";
    		 error( "\n bad _alloc()  in reading mach mlgeno file \n");
     	 }



//-------------------------------------------------------------------//
return;}

void CIPED::pedFileLeser(vector<CBPED*>&pedInfo,const vector<CBSNP*>genInfo)
{
	const unsigned int sz=CBSNP::no_of_persons;
	for(unsigned int i=0;i<sz;++i)
	{
		CBPED* pCBPED =new CBPED;
		pedInfo.push_back(pCBPED);
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

}

void CIPED::imputegen_to_pedgen_converter(const vector<CBSNP*>& genInfo,const float& thresh)
{
	CBSNP* pCSNP=genInfo[0];
	unsigned int _size =0;
	 unsigned int _mycase =-1;
	for(unsigned int j=0;j<genInfo.size(); ++j)
	{
		 pCSNP	=genInfo[j];
		 _size	=pCSNP->pgeno1.size();
		// cout <<pCSNP->pgeno1.size()<<" "<< pCSNP->pgeno2.size() <<" "<<pCSNP->pgeno3.size()<<"\n";
		 if(_size==pCSNP->pgeno2.size() && _size==pCSNP->pgeno3.size())
		 {
			 for(unsigned int i=0; i<_size;++i)
			 {
				 float pg1	=pCSNP->pgeno1[i];
				 float pg2	=pCSNP->pgeno2[i];
				 float pg3	=pCSNP->pgeno3[i];
				 float maxv	=max(pg1,pg2);
	 			 maxv 		= max(maxv,pg3);
	 			 //cout <<"pg1:"<<pg1<< ", "<<pg2 <<"," <<pg3<<", "<<", maxv: "<< maxv <<", "<<thresh<<",: " ; //debug
	 			 if(maxv-pg1<1e-20) 				_mycase=0;
				 else if(maxv-pg2<1e-20)	 		_mycase=1;
				 else if(maxv-pg3<1e-20)			 _mycase=2;
	 			 if(maxv==0.0|| maxv<thresh)		_mycase=3;
	 			// cout <<"\n _mycase: "<< _mycase <<endl;
				 switch(_mycase)
				 {

				 case 0:

						 // case homo1
						 pCSNP->geno1.push_back(false);
						 pCSNP->geno2.push_back(false);
						 pCSNP->aOrder.push_back(true);
						 break;
				 case 1:

					 // case hetero
					 pCSNP->geno1.push_back(false);
					 pCSNP->geno2.push_back(true);
					 pCSNP->aOrder.push_back(true);
					 break;
				 case 2 :
					 // case homo2
					 pCSNP->geno1.push_back(true);
					 pCSNP->geno2.push_back(true);
					 pCSNP->aOrder.push_back(true);
					 break;
				 case 3:
					 // missing case
					 pCSNP->geno1.push_back(true);
					 pCSNP->geno2.push_back(false);
					 pCSNP->aOrder.push_back(true);
					 break;
				 default:
					// cout <<"->error: "<<boolalpha<< (maxv==pg3)<<endl;
					 string msg= "Problem in converting impute-formatted genotype-probabilites in hard genotypes. \n"
							 "Please contact the software provider.\n"
							 "function: imputegen_to_pedgen_converter\n";
					 error(msg);
					 break;
				 } //end of switch

			 } //end of for loop for pedinfo size();
		 }// if all genotypes have same size;
		 else {
			 string msg= "Problem in converting impute-formatted genotype-probabilites in  hard genotypes. \n"
					 "Please contact the software provider.\n"
					 "function: imputegen_to_pedgen_converter\n"
					 "Genotype probabilities have no equal size\n"	 ;
			 error(msg);

		 }

	} // for loop for geninfo.size();
return ; }


void CISNP::lese_imputegens_file(vector<CBSNP*>&genInfo,vector<string>&cZeile, const string fileName)
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
	 //unsigned long int	snpAnzahl		=0;
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
							//cout << myString << " ";
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
				if(cZeile.size()>4)
				{
					if(max_ncol!=ncol)
					{
						string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
						error(msg);
					}
					CBSNP* pCBSNP		=new CBSNP;
					pCBSNP->snpName		=cZeile[0];
					pCBSNP->rsId		=cZeile[1];
					pCBSNP->bp			=atoi(cZeile[2].c_str());
					pCBSNP->allele1		=cZeile[3];
					pCBSNP->allele2		=cZeile[4];
					//now erase the first 5 columns of pedInfo from ZeilenVector
					cZeile.erase(cZeile.begin(),cZeile.begin()+5);
					//cout << "cZeile size: "<< cZeile.size()<<endl;
					CBSNP::no_of_persons=(max_ncol-5)/3;
					lese_imputegens_file_genInfo( pCBSNP, cZeile,nrow);
					genInfo.push_back(pCBSNP);
				}
				else
				{
					string emsg	="At least 5 columns are necessary in each row of impute.gens type of file. However, there are less than 5 columns on the "+change_int_into_string(nrow)+"th row of file \""+fileName +"\". ";
					error(emsg);
				}

			}
			 else continue;
			 // assignment of genotypes
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
		string pmsg="\t**->Total number of successfully read individuals: ";
		//printLIN(pmsg);
		printWithSpace(pmsg,  change_int_into_string(CBSNP::no_of_persons),CBPED::coutWidth);
		//pmsg="\t**->Total number of successfully read SNPs: ";
		//printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
	}
	 catch(bad_alloc& memoryAllocationException)

	 {
		 cout << "reading with second type  reading art \n";
	
	 }

return; }
//helper functions
void CISNP::lese_imputegens_file_genInfo( CBSNP* const pCBSNP, const vector<string>& zeilenVec,const int akt_zeile)
{
	bool fatal=false;
	string EMSG ="\t**->probblem with reading line "+change_int_into_string(akt_zeile)+":\n";
	string errmsg="";
	double first		=0.0;
	double second	=0.0;
	double third		=0.0;
	if(zeilenVec.size()%3!=0)
	{
		fatal =true;
		errmsg = "There are not correct number of columns.";
	}
	else
	{
		float temp_summ 	=0.0;
		bool  temp_boole	=false;
		bool  _tmp_miss 	=false;
		for(unsigned int i=1;i<=zeilenVec.size()/3;++i)
		{
			first		=atof(zeilenVec[3*i-3].c_str());
			second		=atof(zeilenVec[3*i-2].c_str());
			third 		=atof(zeilenVec[3*i-1].c_str());
			//cout << " first: "<< first << ", "<< second << ", "<< third <<"\n";
			temp_summ	=first+second+third;
			temp_boole	=(first<0.0 ||first>1.0001) || (second<0.0 ||second>1.0001) ||(third<0.0 ||third>1.0001);
			_tmp_miss  	=(temp_summ==0.0);
			if((temp_summ!=1 || temp_boole)&& !_tmp_miss)
			{
				if(temp_summ!=1.0 && temp_summ<1.1 )
				{
					first 	= first/temp_summ;
					second	= second/temp_summ;
					third 	= third/temp_summ;
				}
				else
				{
					printLIN("The probabality distribution of SNP: "+pCBSNP->rsId+" is not correct. \n" );
					cout <<"Three probabilities: "<< first <<" "<< second <<" "<<third <<" \n";
					exit(1);
				}
			}
			//cout << " " << first << "  " <<second <<" "<<third <<endl; // for test only
			pCBSNP->pgeno1.push_back(first);
			pCBSNP->pgeno2.push_back(second);
			pCBSNP->pgeno3.push_back(third);

		}
		if(CBSNP::no_of_persons!=pCBSNP->pgeno1.size()||CBSNP::no_of_persons!=pCBSNP->pgeno2.size()||CBSNP::no_of_persons!=pCBSNP->pgeno2.size())
		{
			fatal=true;
			errmsg = "Problem with saving the genotypes.\n";
		}


	}
	//for(unsigned int i=0; i<pCBSNP->pgeno1.size();++i)
		//cout << pCBSNP->pgeno1[i]<< " "<< pCBSNP->pgeno2[i]<<" "<< pCBSNP->pgeno3[i]<< " \n";

	if(fatal)
	{
		printLIN(EMSG);
		error(errmsg);
	}
return ;

}

void CIPED::lese_hapmap_genotype_vec(CBSNP* const pCBSNP,const vector<string>& zeilenVec, const int& nLines)
{
	 //cout << zeilenVec.size()<<endl;
	 if(zeilenVec.size()%2!=0)
	 {
		cerr<<"Total columns in "+change_int_into_string(nLines)+"th column: " << zeilenVec.size()<<endl;
		LIN<<"Total columns in "+change_int_into_string(nLines)+"th column: " << zeilenVec.size()<<endl;
		string pr_msg="In each line of genotype file, therere should even number of columns so that each individual can have a pair of haplotypes.\n However this is not the case in "+\
		change_int_into_string(nLines)+"th line. ";
		delete pCBSNP;
		error(pr_msg);
	 
	 }
	 int value1 =0;
	 int value2 =0;
	for(unsigned int i=0;i<(zeilenVec.size()/2);++i)
	{
		value1 = atoi(zeilenVec[(2*i)].c_str());
		value2	= atoi(zeilenVec[(2*i+1)].c_str());
		//cout << value1 << value2 << endl;
		if(!value1&& !value2) // if 00
		{
			// homozygote 1 ;
			pCBSNP->geno1.push_back(false);
			pCBSNP->geno2.push_back(false);
			pCBSNP->aOrder.push_back(true);
		}
		//heterozygote 
		else if(!value1&& value2)
		{
		
			pCBSNP->geno1.push_back(false);
			pCBSNP->geno2.push_back(true);
			
			pCBSNP->aOrder.push_back(true);
		}		
		else if(value1&& !value2)
		{
			pCBSNP->geno1.push_back(false);
			pCBSNP->geno2.push_back(true);
			pCBSNP->aOrder.push_back(false);
		}
		//homozygote 2
		else if(value1&& value2) // if 11
		{
		
			pCBSNP->geno1.push_back(true);
			pCBSNP->geno2.push_back(true);
			pCBSNP->aOrder.push_back(true);
		}
		// missing
		else 
		{
		
			pCBSNP->geno1.push_back(true);
			pCBSNP->geno2.push_back(false);
			pCBSNP->aOrder.push_back(true);
		}
		
	}
	//cout<<"pCBSNP->geno1.size(): "<<pCBSNP->geno1.size()<< " "<< pCBSNP->geno1.size() <<endl; //debug
return ;}	
//lese impute hapamp file
// read mach hapmap file
void CIPED::lese_impute_hapmap_file(vector<CBPED*>&vPed,vector<CBSNP*>& genInfo,vector<string>& cZeile,const string& fileName)
{
	//cout << "test"<<endl;//debug 
	printLIN("*->Reading file \""+fileName+"\": \n");
	string myString					="";
	string temp_value				="";
	string pheno					="";
	bool 	ende 					=false;
	unsigned long int ccol			=0;
	unsigned long int ncol			=0;
	unsigned int 	 nrow 			=0;
	unsigned int nLines				=0;
	unsigned long int max_ncol		=0;
	static unsigned long int snpAnzhal		=0;
	bool def_max_ncol_once			=true;
	bool define_pCBPED				=true;
	
	try
	{
		IFS pf(fileName.c_str(),"r");
		if(!pf){
			string msg =fileName+ " either does not exits or could not be opend.\n";
			error(msg);		
		}
		while(true)
		{
			if(_ifs_eof(pf)) break;
			while(!_ifs_eof(pf))
			{
				ccol	=CBSNP::stringLeser(pf,myString,ende);
				ncol	=ncol+ccol;
				//cout << "ccol : "<< ccol << ",\t ncol: "<<ncol << ", ende: "<<ende <<" _ifs_eof: "<<_ifs_eof(pf)<< endl;
				if(ccol!=0)
				{
					
					if(!myString.empty())
						cZeile.push_back(myString);
						
					else
					{
						string pr_msg= "Error in reading "+change_int_into_string(nLines+1)+"th line of file \""+fileName+"\". \n";
						error(pr_msg);
					}
					//cout << "test1" << endl; // debug 
					if(ende || _ifs_eof(pf))
					{
						++nLines;
						if(def_max_ncol_once && (ccol!=0))
						{
							max_ncol=ncol;
							def_max_ncol_once =false;
						//	cout << "max_ncol: "<< max_ncol << endl;//debug 
						}
						if((ncol==max_ncol) && (max_ncol!=0))
							++nrow;
						myString="";
						break;	
							
					}	
				
				}
				if(ccol==0)continue;
			}//end of inner while loop while loop ;
			//cout <<"size: "<< cZeile.size() << " "<<endl;
			if(cZeile.size()!=0)
			{
				if(cZeile.size()!=max_ncol)
				{
					string pr_msg="Problem in reading  the  "+change_int_into_string(nLines) +"th  line: \n";
					printLIN(pr_msg);
					cerr <<"Number of columns in first line : "<< max_ncol << endl; 
					LIN <<"Number of columns in first line : "<< max_ncol << endl; 
					//
					cerr <<"Number of columns in current line : "<< cZeile.size() << endl; 
					LIN <<"Number of columns in current line : "<< cZeile.size() << endl; 
					pr_msg="These two lines should have the same number of columns as in the first line.\n" ;
					cerr_with_delete(vPed,genInfo,pr_msg);
				}
				//cout <<"max_ncol:"<< max_ncol <<endl; 
				if(define_pCBPED)
				{
					define_pCBPED=false;
					for(unsigned int i=0;i<(max_ncol/2);++i)
					{
						CBPED* pCBPED =new CBPED;
						lese_phenotype_info(pCBPED, pheno);
						vPed.push_back(pCBPED);
					}
				}
				//now read the genotypes one by one.
				//CBSNP* pCBSNP		=new CBSNP;
				//cout << "indiv: " << vPed.size()<< ", snp size: "<< genInfo.size()<<endl;
				//for(unsigned int i=0; i<genInfo.size();i++)
				if(snpAnzhal<genInfo.size())
				{

					CBSNP* pCBSNP = genInfo[snpAnzhal];
					lese_hapmap_genotype_vec(pCBSNP,cZeile, nLines);
				//	cout << pCBSNP->geno1.size()<< pCBSNP->geno2.size()<<endl;
							//genInfo.push_back(pCBSNP);
					++snpAnzhal;
				}

			}
			else 
			continue;
			ncol=0;
			cZeile.clear();
			cZeile.resize(0);
			if(_ifs_eof(pf))break;
			
			//break;
		}//end of first while loop
		_ifs_close(pf);
		nrow	=0;
		ccol	=0;
		nLines	=0;
		
	
	}catch(bad_alloc& memoryAllocationException)
	{
 		 cout<<"\n Please inform the programmer that you could not read mlgeno file. Probbaly he will give you a solution. \n ";
 		 error( "\n bad _alloc()  in reading mach hapamp file \n");
  	 }

	//cout << "test" <<endl; //debug  

return;}

// hapmap impute legend file 

void CISNP::read_impute_ref_legendfile(vector<CBSNP*>& genInfo,const string& legendFileName)
{
	// here comes the read extra file like structre.
	string r_msg = "*->Reading Impute reference legend file:\""+legendFileName+"\":\n";
	printLIN(r_msg);
	string line;
	rFile file; file.close(); file.clear();
	static int extraMapFileLineZaehler =0;
	file.open(legendFileName.c_str(),ios::in);
	if(!file){
		file.setstate(ios_base::failbit);
		file.clear();
		error("Your legend file \""+legendFileName+"\" eitehr does not exit or can't be opened. Have you given the correct path and name of the file?\n");

	}
	// check the first header and the number of column
	do
	{
		getline(file,line,'\n');
		if(line[0]=='#') continue;
		if(file.eofbit) break;

	}while(line.empty());
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

	if(header.size()!=4 || header[0]!="rsID" || header[1]!="position" || !(header[2]=="a0"&& header[3]=="a1") )
	{//cout << "This part is only for mach hapmap \n";

			string msg= "Problem in parsing the first line! The first line should contain \"rsID\", \"position\", \"a0\" and \"a1\" as header.\n";
			error(msg);

	}
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
						if(tokens.size()!=4)
						{
							string msg= "Your legend  file \""+legendFileName+"\"  must contain 4 columns in the line "+change_int_into_string(extraMapFileLineZaehler)+\
									", however there are "+change_int_into_string(tokens.size())+" columns.\n";
							error(msg);
						}
						else
						{
							CBSNP* pCBSNP		=new CBSNP;
							pCBSNP->rsId 		=tokens[0];
							pCBSNP->bp 			=atoi(tokens[1].c_str());
							pCBSNP->allele1		=tokens[2];
							pCBSNP->allele2		=tokens[3];
							genInfo.push_back(pCBSNP);
							//cout <<" genInfo.size()"<< genInfo.size()<<endl; //test

						}

					}


				}

}//end of function

//shapeit haps file
 void CISNP::lese_shapeit_haps_file(vector<CBPED*>&vPed,vector<CBSNP*>& genInfo,const string& fileName)
 {
	 // cout <<"hi from haps file\n";
	  //cout << "test"<<endl;//debug
	  	printLIN("*->Reading file \""+fileName+"\": \n");
	  	string myString					="";
	  	string temp_value				="";
	  	string pheno					="";
	  	bool 	ende 					=false;
	  	unsigned long int ccol			=0;
	  	unsigned long int ncol			=0;
	  	unsigned int 	 nrow 			=0;
	  	unsigned int nLines				=0;
	  	unsigned long int max_ncol		=0;
	  	static unsigned long int snpAnzhal		=0;
	  	bool def_max_ncol_once			=true;
	  	vector<string> cZeile;
	  	cZeile.clear(); cZeile.resize(0);
	  	try
	  	{
	  		IFS pf(fileName.c_str(),"r");
	  		if(!pf)
	  		{
	  			string msg =fileName+ " either does not exits or could not be opend.\n";
	  			error(msg);
	  		}
	  		while(true)
	  		{
	  			if(_ifs_eof(pf)) break;
	  			while(!_ifs_eof(pf))
	  			{
	  				ccol	=CBSNP::stringLeser(pf,myString,ende);
	  				ncol	=ncol+ccol;
	  				//cout << "ccol : "<< ccol << ",\t ncol: "<<ncol << ", ende: "<<ende <<" _ifs_eof: "<<_ifs_eof(pf)<< endl;
	  				if(ccol!=0)
	  				{

	  					if(!myString.empty())
	  						cZeile.push_back(myString);

	  					else
	  					{
	  						string pr_msg= "Error in reading "+change_int_into_string(nLines+1)+"th line of file \""+fileName+"\". \n";
	  						error(pr_msg);
	  					}
	  					//cout << "test1" << endl; // debug
	  					if(ende || _ifs_eof(pf))
	  					{
	  						++nLines;
	  						if(def_max_ncol_once && (ccol!=0))
	  						{
	  							max_ncol=ncol;
	  							def_max_ncol_once =false;
	  						//	cout << "max_ncol: "<< max_ncol << endl;//debug
	  						}
	  						if((ncol==max_ncol) && (max_ncol!=0))
	  							++nrow;
	  						myString="";
	  						break;

	  					}

	  				}
	  				if(ccol==0)continue;
	  			}//end of inner while loop while loop ;
	  			//cout <<"size: "<< cZeile.size() << " "<<endl;
	  			if(cZeile.size()!=0)
	  			{
	  				if(cZeile.size()!=max_ncol)
	  				{
	  					string pr_msg="Problem in reading  the  "+change_int_into_string(nLines) +"th  line: \n";
	  					printLIN(pr_msg);
	  					cerr <<"Number of columns in first line : "<< max_ncol << endl;
	  					LIN <<"Number of columns in first line : "<< max_ncol << endl;
	  					//
	  					cerr <<"Number of columns in current line : "<< cZeile.size() << endl;
	  					LIN <<"Number of columns in current line : "<< cZeile.size() << endl;
	  					pr_msg="These two lines should have the same number of columns as in the first line.\n" ;
	  					//cerr_with_delete(vPed,genInfo,pr_msg);
	  					error(pr_msg);
	  				}
	  				//cout <<"max_ncol:"<< max_ncol <<endl;

	  				//now read the genotypes one by one.
	  				//CBSNP* pCBSNP		=new CBSNP;
	  				//cout << "indiv: " << vPed.size()<< ", snp size: "<< genInfo.size()<<endl;
	  				//for(unsigned int i=0; i<genInfo.size();i++)
	  				if(snpAnzhal<=genInfo.size())
	  				{
	  					CBSNP* pCBSNP		=new CBSNP;
	  					pCBSNP->nchr		=cZeile[0];
	  					pCBSNP->rsId		=cZeile[1];
	  					pCBSNP->bp			=atoi(cZeile[2].c_str());
	  					pCBSNP->allele1		=cZeile[3];
	  					pCBSNP->allele2		=cZeile[4];
	  					//now erase the first 5 columns of pedInfo from ZeilenVector
	  					cZeile.erase(cZeile.begin(),cZeile.begin()+5);
	  					//cout << "cZeile size: "<< cZeile.size()<<endl;
	  					CBSNP::no_of_persons=(max_ncol-5)/2;
	  					if(vPed.size()!=CBSNP::no_of_persons)
	  					{
	  							printLIN("genotype information is not given for all samples.\n");
	  							cout <<"No of Samples: "<<vPed.size()<<", No of Genotypes:"<<cZeile.size()/2<<endl;
	  							error("Please check your haps file.\n");
	  					}
	  					//lese_imputegens_file_genInfo( pCBSNP, cZeile,nrow);
	  					CIPED::lese_hapmap_genotype_vec(pCBSNP,cZeile, nLines);
	  					genInfo.push_back(pCBSNP);
	  					++snpAnzhal;
	  				}

	  			}
	  			else
	  			continue;
	  			ncol=0;
	  			cZeile.clear();
	  			cZeile.resize(0);
	  			if(_ifs_eof(pf))break;

	  			//break;
	  		}//end of first while loop
	  		_ifs_close(pf);
	  		nrow	=0;
	  		ccol	=0;
	  		nLines	=0;


	  	}catch(bad_alloc& memoryAllocationException)
	  	{
	   		 cout<<"\n Please inform the programmer that you could not read mlgeno file. Probbaly he will give you a solution. \n ";
	   		 error( "\n bad _alloc()  in reading mach hapamp file \n");
	    	 }

	  	//cout << "test" <<endl; //debug


	//------------------------------------------------------------------------------//
 }
 //shapeit haps file
  void CIPED::lese_shapeit_sample_file(vector<CBPED*>&vPed,const string& fileName)
  {
	  	//cout << "hi from sample file \n";
	  	// here comes the read extra file like structre.
	  	string r_msg = "\n*->Reading SHAPEIT sample file:\""+fileName+"\":\n";
	  	printLIN(r_msg);
	  	string line;
	  	rFile file; file.close(); file.clear();
	  	static int fileLineZaehler =0;
	  	file.open(fileName.c_str(),ios::in);
	  	if(!file){
	  		file.setstate(ios_base::failbit);
	  		file.clear();
	  		error("Your sample file \""+fileName+"\" either does not exit or can't be opened. Have you given the correct path and name of the file?\n");

	  	}
	  	// check the first header and the number of column
	  	do
	  	{
	  		getline(file,line,'\n');
	  		if(line[0]=='#') continue;
	  		if(file.eofbit) break;

	  	}while(line.empty());
	  	string buffer="";
	  	vector<string> header;
	  	stringstream line_parser(line);
	  	if(line_parser.good())
	  	{
	  		while(line_parser>>buffer)
	  			header.push_back(buffer);
	  		++fileLineZaehler;
	  	}
	  	else
	  	{
	  		string msg= "Problem in parsing the first line!\n";
	  		error(msg);

	  	}
	  			// this part is only for only mach hapmap data

	  	if(header.size()<3 || header[0]!="ID_1" || header[1]!="ID_2" || !(header[2]=="missing") )
	  	{//cout << "This part is only for mach hapmap \n";

	  			string msg= "Problem in parsing the first line! The first line should contain at lest three column \"ID_1\", \"ID_2\", \"missing\"  as header.\n";
	  			error(msg);

	  	}
	  	//2nd line is just to given info. This line should be discarded
	  	do
	  	{
	  		getline(file,line,'\n');
	  		if(line[0]=='#') continue;
	  		if(file.eofbit) break;
	  	}while(line.empty());
	  	line.clear();line="";
	  	vector<string> tokens;
	  	//find header size
	  	 const unsigned int NCOLS=header.size();
	  	 //header.resize(0); header.clear();
	  	 while(getline(file,line, '\n')) // until eof hit
	  	 {

		  		//-------------------------------------------------------------------------------//
	  			++fileLineZaehler;
	  			//cout <<"line no: "<< extraMapFileLineZaehler<<endl; //debug
	  			if(line.length()==0)
	  				continue;
	  			if(line[0]=='#')
	  				continue;
	  			//
	  			tokens.resize(0);
	  			tokens.clear();
	  			stringstream line_parser1(line);
	  			//cout <<"\n line_parser1.good(): "<<line_parser1.good() <<endl; // only for test
	  			if(line_parser1.good())
	  			{
	  				while(line_parser1>>buffer)
	  				{
	  					tokens.push_back(buffer);
	  					//	cout << buffer << " "; // only for test
	  				}
	  				if(tokens.size()!=NCOLS)
	  				{
	  					string msg= "Your sample  file \""+fileName+"\"  must contain"+change_int_into_string(NCOLS) +"columns in the line "+change_int_into_string(fileLineZaehler)+\
	  							", however there are "+change_int_into_string(tokens.size())+" columns.\n";
	  					error(msg);
	  				}
	  				else
	  				{
	  					CBPED* pCBPED		=new CBPED;
	  					pCBPED->famId 		=tokens[0];
	  					pCBPED->indId		=tokens[1];
	  					vPed.push_back(pCBPED);
	  					//cout <<" genInfo.size()"<< vPed.size()<<endl; //test
	  					for(unsigned int j=3; j<NCOLS;++j)
	  					{
	  						if(header[j]=="sex")
	  							CBPED::lese_sex_info(pCBPED, tokens[j]);
	  						//pheno type
	  						if(header[j]=="plink_pheno"|| header[j]=="pheno")
	  							CBPED::lese_phenotype_info(pCBPED,tokens[j]);
	  					}


	  				}
	  			}
	 //-------------------------------------------------------------------------------//

	  		}//end of while
  //end
  }
