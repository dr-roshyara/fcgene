/*
 * bimbam.cpp
 *
 *  Created on: May 3, 2012
 *	Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 *     These are the basic parameters
 *

 */
#include"bimbam.h"

CBIMPED::~CBIMPED()
{
	/*	
	for (unsigned int i=0;i<vPed.size();++i)
	{
		delete vPed[i];
	}
		vPed.clear();
		vPed.resize(0);
	*/	
}
//
CBIMSNP::~CBIMSNP()
{
	/*
	for (unsigned int i=0;i<genInfo.size();++i)
	{
		delete genInfo[i];
	}
		genInfo.clear();
		genInfo.resize(0);
	*/	
}

//

//-------------------------------------------------//
void CBIMPED::read_bimbam_data(vector<CBPED*>&pedInfo, vector<CBSNP*>&genInfo, vector<string> cZeile,const string fileName)
{
	printLIN("*->Reading file \""+ fileName+ "\": \n");
	string myString					="";
	string emsg						="";
	bool fatal							=false;
	bool ende							= false;
	unsigned int ccol					=0;
	unsigned int ncol					=0;
	unsigned int nrow					=0;
	unsigned int nLines					=0;
	unsigned int max_ncol				=0;
	bool def_max_ncol_once				=true;
	//unsigned long int	snpAnzahl		=0;
	unsigned long int  nSample			=0;
	bool def_nSample					=true;
	string first_geno					="";
	string second_geno					="";
	string pheno_value					="";
	string sex_value					="";
	string temp_value					="";
	string temp_value1					="";
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
							if(nrow==0|| nrow==1)
							{
								if(ncol==1)
									++nrow;
							}
							else if(def_max_ncol_once &&(ncol!=0))
							{
								max_ncol=ncol;
								def_max_ncol_once=false;
								//cout << "max_ncol: "<<max_ncol <<" " <<"ncol: "<<ncol <<endl; // only test; // only test
							}
							if((ncol==max_ncol)&& (max_ncol!=0))
								++nrow;
								//cout << "nrow : "<< nrow << endl; // debug
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
						if(nrow<4)
						 {
							 if(nrow==1)
							 {
									//save person id ;
								 //cout << cZeile.size()<<endl; // debug
								if(cZeile.size()==1)
								 {
									 if(def_nSample)
										 nSample=atoi(cZeile[0].c_str());
									 def_nSample=false;

								 }
								 else
								 	 {
									 emsg ="The first line of file  \""+fileName+\
											 "\" should contain the only a number representing the number of individuals .\n";
									 error(emsg);

								 	 }
							 }
							 else if(nrow==2)
							 {
								 //save person id ;
								// cout << "nrow: "<< nrow <<endl;
								// cout <<"cZeile.size(): "<< cZeile.size()<<endl; // debug
								if(cZeile.size()==1)
								 {
									//snpAnzahl=atoi(cZeile[0].c_str());
								 }
								 else
								 {
									 emsg ="The second line of file  \""+fileName+\
											 "\" should contain the only a number representing the number of SNPs .\n";
									 error(emsg);

								 }
							 }
							//--
							 else if(nrow==3)
							 {
								 if(cZeile[0]!="IND,")
								 {
									 emsg ="Third line of Bimbam format file \""+fileName+\
											 "\" should start with the word \"IND,\"  or   the word \"id,\" .\n";
									 error(emsg);
								 }
								 cZeile.erase(cZeile.begin(),cZeile.begin()+1);
								 if(cZeile.size()== nSample)
								 {

									 for(unsigned int i=0; i<nSample; ++i)
									 {
										 temp_value=cZeile[i];
										 temp_value1=temp_value.substr(0,temp_value.rfind(","));
										 CBPED* pCBPED = new CBPED;
										 pCBPED->indId=temp_value1;
										 CBPED::lese_sex_info( pCBPED, sex_value);
										 CBPED::lese_phenotype_info(pCBPED,pheno_value);
										 pedInfo.push_back(pCBPED);
									 }
								}
								 else
								 {
									 emsg ="On the   third line in file \""+fileName+\
											 "\" should contain"+ change_int_into_string(nSample)+\
											 "  individual ids but this is the case here.\n";
									 error(emsg);
								 }
							 }

							//-------

		 				 }
						 else
						 {
							 //save genotypes.
							  temp_value			=cZeile[0];
							  const string rsName =temp_value.substr(0,temp_value.rfind(","));
							  cZeile.erase(cZeile.begin(),cZeile.begin()+1);
							  if(cZeile.size()==(nSample))
							  {
								  CBSNP* pCBSNP = new CBSNP;
								  pCBSNP->no_of_persons=nSample;
								  if(pCBSNP->rsId=="")
									  pCBSNP->rsId=rsName;
								  if(pCBSNP->snpName=="")
									 pCBSNP->snpName=rsName;
									 genInfo.push_back(pCBSNP);

									 if(genInfo.size()==(nrow-3))
									 {
										 CBSNP* npCBSNP=genInfo[genInfo.size()-1];
										 for(unsigned int i=0; i<nSample; ++i)
										 {
											 temp_value=cZeile[i];
											 temp_value1=temp_value.substr(0,temp_value.rfind(","));

											 if(temp_value1.size()==2)
											 {
												 first_geno		= temp_value[0];
												 second_geno	=temp_value1[1];
											 }
											 else
											 {
												 emsg = "Each column should contain a genotype of two letters but on the"+\
														 change_int_into_string(nrow)+"th row , this is not the case!";
												 error(emsg);
											 }

											// cout << first_geno <<endl; // debug
											// cout << second_geno <<endl; // debug
											fatal=lese_genotype_info(npCBSNP, first_geno, second_geno,nrow, emsg);
											 if(fatal)
											 {
												 //delete pCBSNP;
												 error(emsg);
											 }

										 }
									  }
									 else
									 {
										 emsg =" problem in saving genotypes line in file \""+fileName+\
												 "\".\n";
										 error(emsg);
									 }

									}
							 else
							 {
								 emsg =" problem in saving genotypes line in file \""+fileName+\
										 "\".\n";
								 error(emsg);
							 }

						  }//end of genotype saving else loop

						 //CBPED* pCBPED= lese_pedfile_pedInfo(cZeile,snpAnzahl, nrow);
						 //pedInfo.push_back(pCBPED);
						 //cout << "nrow: "<< nrow<< " pedInfo.size()"<< pedInfo.size()<<endl;
							 if(max_ncol!=ncol && nrow>2)
							 {
								 string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
								 error(msg);
						 	 }


					 }// while loop
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
				string pmsg="\t**->Total number of successfully read SNPs: ";
				printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
				pmsg="\t**->Total number of successfully read individuals: ";
				//CBPED::display_ped_summary();
				//printLIN(pmsg);
				printWithSpace(pmsg,  change_int_into_string(pedInfo.size()),coutWidth);
				//cout << genInfo.size()<<endl; // debug
			 }catch(bad_alloc& memoryAllocationException)
			 {
				 cout<<"\n Please inform the programmer that you could not read beagle file. Probably he will give you a solution. \n ";
				 error( "\n bad _alloc()  in reading beagle bgl file \n");
			 }

return;}

void CBIMSNP::read_bimbam_pos_data(vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName){
//---------------------------------------------------------//
	printLIN("*->Reading file \""+ fileName+ "\": \n");
	string myString					="";
	string emsg						="";
	//bool fatal							=false;
	bool ende							= false;
	unsigned int ccol					=0;
	unsigned int ncol					=0;
	unsigned int nrow					=0;
	unsigned int nLines					=0;
	unsigned int max_ncol				=0;
	bool def_max_ncol_once				=true;
	//unsigned long int	snpAnzahl		=0;
	//unsigned long int  nSample			=0;
	//bool def_nSample					=true;
	string first_geno					="";
	string second_geno					="";
	string pheno_value					="";
	string sex_value					="";
	string temp_value					="";
	int temp_int						=0;
	try
	{
		IFS pf(fileName.c_str(),"r");
		if(!pf)
		{
			string msg= fileName+ " either does not exits or could not be opend.\n"	;
			error(msg);
		}
		vector<CBSNP*>::iterator  pCBSNP =genInfo.begin();
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
					//cout << "nrow : "<< nrow << endl; // debug
					myString="";
					//cout <<"-----------------\n";//test
					break;
				}
				if(ccol==0)
					continue;
				// this loop is broken only if end  is found.
			}//end of inner while loop
			//------------------------------//
			if(cZeile.size()!=0)
			{
				if(cZeile.size()==2||cZeile.size()==3)
				{
					//save genotypes.
					//vector<CBSNP*>::iterator  pCBSNP =genInfo.begin();
					temp_value			=cZeile[0];
					const string rsName =temp_value.substr(0,temp_value.rfind(","));
					//cZeile.erase(cZeile.begin(),cZeile.begin()+1);
					//--------------------------------------------------//
					if((*pCBSNP)->rsId!="" || (*pCBSNP)->snpName!="")
					{
						if(rsName!=(*pCBSNP)->rsId || rsName!=(*pCBSNP)->snpName)
						{
							emsg= rsName+" does not match with the rsid or snpid of SNP containing on the "+change_int_into_string(nrow)+"th row  of file "+fileName + ".\n" ;
							error(emsg);
						}


					}
					else
					{
						(*pCBSNP)->rsId=rsName;
					if((*pCBSNP)->snpName=="")
						(*pCBSNP)->snpName=rsName;
					}
					// till here
					//--------------------------//
					// base pair position
					temp_value	=cZeile[1];
					if(cZeile.size()==3)
						temp_value	=temp_value.substr(0,temp_value.rfind(","));
					temp_int	=atoi(temp_value.c_str());
					//cout <<temp_int<<endl;
					if((*pCBSNP)->bp==-1 )
					{
						(*pCBSNP)->bp=temp_int;

					}
					else
					{
						emsg= temp_value+" does not match with the base pair position   of SNP containing on the "+change_int_into_string(nrow)+"th row  of file "+fileName + ".\n" ;
						error(emsg);
					}
					// for third column
					if(cZeile.size()==3)
					{

						temp_value	=cZeile[2];
						if(((*pCBSNP)->nchr) =="0")
						{
							(*pCBSNP)->nchr=temp_value;
						}
						else
						{
							if((*pCBSNP)->nchr!=temp_value)
							{
								emsg= temp_value +" does not match with the chromosome number of  SNP containing on the "+change_int_into_string(nrow)+"th row  of file "+fileName + ".\n" ;
								error(emsg);
							}
						}


					}

					++pCBSNP;
					cZeile.clear();
					cZeile.resize(0);

				}//end of genotype saving else loop
				else
				{
					emsg = " SNP location file "+fileName+" should contain either 2 or three columns but this is not the case here in"+change_int_into_string(nrow)+"th row \n";
					error(emsg);
				}


			}// while loop
			else continue;
			ncol=0;
			cZeile.clear();
			cZeile.resize(0);
			if(_ifs_eof(pf))break;
		} // end of outer while loop
		_ifs_close(pf);
		nrow=0; // total rows
		nLines=0; // total lines
	}catch(bad_alloc& memoryAllocationException)
	{
		cout<<"\n Please inform the programmer that you could not read beagle file. Probably he will give you a solution. \n ";
		error( "\n bad _alloc()  in reading beagle bgl file \n");
	}


//---------------------------------------------------------//
return;}
// This file reads the bimbam snpinfo data produced by Bimbam imputation. 
//-------------------------------------------------//

void CBIMSNP::read_bimbam_snpinfo_data(vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName,const float &maf_value)
{
	printLIN("*->Reading file \""+ fileName+ "\": \n");
	string myString					="";
	string emsg						="";
	//bool fatal							=false;
	bool ende							= false;
	unsigned int ccol					=0;
	unsigned int ncol					=0;
	unsigned int nrow					=0;
	unsigned int nLines					=0;
	unsigned int max_ncol				=0;
	bool def_max_ncol_once				=true;
	int	temp_int						=0;
	int	akt_int							=0;
	//unsigned long int  nSample			=0;
	//bool def_nSample					=true;
	string first_geno					="";
	string second_geno					="";
	string pheno_value					="";
	string sex_value					="";
	string temp_value					="";
	try
	{
		IFS pf(fileName.c_str(),"r");
		if(!pf)
		{
			string msg= fileName+ " either does not exits or could not be opend.\n"	;
			error(msg);
		}
		vector<CBSNP*>::iterator  pCBSNP =genInfo.begin();
		while(true)
		{
			if(_ifs_eof(pf))break;
			while(!_ifs_eof(pf))
			{
				ccol=CBSNP::stringLeser(pf, myString, ende);
				ncol+=ccol;
			//	cout << " ncol: "<< ncol << " ccol: "<< ccol << endl; //debug
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
				//cout << "ende "<<boolalpha<< ende << ", _ifs_eof(pf): "<< boolalpha<<(_ifs_eof(pf)) <<endl; // debug
				//cout<< "ende| (_ifs_eof(pf)): "<< (ende| (_ifs_eof(pf)))<<endl;
				if(ende| (_ifs_eof(pf)))
					{
						//ende=false;
						//cout << "nrow: " << nrow << "ncol: "<< ncol << endl;
						++nLines;
						if(def_max_ncol_once &&(ncol!=0))
						{
								max_ncol=ncol;
								def_max_ncol_once=false;
								//cout << "max_ncol: "<<max_ncol <<" " <<"ncol: "<<ncol <<endl; // only test; // only test
						}
						if((ncol==max_ncol)&& (max_ncol!=0))
							++nrow;
						//cout << "nrow : "<< nrow << endl; // debug
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
				if(cZeile.size()!=5 && cZeile.size()!=6)
				{
					myString="There should be 5 or 6 columns in each row of file: "+fileName+". but this is not the case in "+\
							change_int_into_string(nrow)+"th row. \n";
					error(myString);
				}
				else if(nrow==1)
				{

				}
				else
				{
					string rsName= cZeile[0];
					if((*pCBSNP)->rsId!="" || (*pCBSNP)->snpName!="")
					{
						//cout << rsName <<endl;
						//cout << (*pCBSNP)->rsId << " " << (*pCBSNP)->snpName  <<endl;//debug
						if(rsName!=(*pCBSNP)->rsId && rsName!=(*pCBSNP)->snpName)
						{
							emsg= rsName+" does not match with the rsid or snpid of SNP containing on the "+change_int_into_string(nrow)+"th row  of file "+fileName + ".\n" ;
							error(emsg);
						}
					}
					else
					{
						(*pCBSNP)->rsId=rsName;
						if((*pCBSNP)->snpName=="")
							(*pCBSNP)->snpName=rsName;
					}
					//----------------------------------------//
					// here comes the allele info stuff.
					//----------------------------------------//
					string a1	=cZeile[1];
					string a2	=cZeile[2];
					//cout <<" a1 "<<a1 <<", a2 "<< a2 << " \n";
					//cout <<"a1==N "<< (a1=="N") <<endl;
					if(a1=="N")
					{
						a1=a2;
						a2="";
					}
					//cout <<" a1 "<<a1 <<", a2 "<< a2 << " \n";
					//cout <<(*pCBSNP)->allele1 <<" "<<(*pCBSNP)->allele2 <<" \n";
					//first allele
					if((*pCBSNP)->allele1!="" && (*pCBSNP)->allele1!="0")
					{
						//---------------------------------------------------//
						if((*pCBSNP)->allele1 !=a1)
						{
							bool temp_bool1=(*pCBSNP)->allele1 ==a2;
							bool temp_bool2= (*pCBSNP)->allele2 ==a1;
							bool temp_bool3=(*pCBSNP)->allele2 ==a2;
							bool temp_bool4= ((*pCBSNP)->allele2=="");
							if(temp_bool1 & (temp_bool2||temp_bool3 || temp_bool4))
							{
								//cout << "change of allele! allele 1 is allele2 and allele 2 will be allele1 \n"; //ONly for test
								for(unsigned int i=0; i<(*pCBSNP)->geno1.size();++i)
								{
									//homo1.homo2
									if(!((*pCBSNP)->geno1[i])&& !((*pCBSNP)->geno2[i]))
									{
										((*pCBSNP)->geno1[i])	=true;
										((*pCBSNP)->geno2[i])	=true;
									}
									//homo2 .homo1
									else if(((*pCBSNP)->geno1[i])&& ((*pCBSNP)->geno2[i]))
									{
										(*pCBSNP)->geno1[i]=false;
										(*pCBSNP)->geno2[i]=false;
									}
									else if(!((*pCBSNP)->geno1[i])&& ((*pCBSNP)->geno2[i]))
									{
										if((*pCBSNP)->aOrder[i])
											(*pCBSNP)->aOrder[i]=false;
										if(!((*pCBSNP)->aOrder[i]))
											(*pCBSNP)->aOrder[i]=true;
									}
									// change probability genotypes also
									if((*pCBSNP)->pgeno1.size()>0 && (*pCBSNP)->pgeno3.size()>0)
									{
										//cout << no_of_persons <<endl;
										for(unsigned int j=0;j<(*pCBSNP)->pgeno1.size(); ++j)
											swap((*pCBSNP)->pgeno1[j],(*pCBSNP)->pgeno3[j]);
										//swap will exchange prob(AA) and prob(BB).prob(AB) remains the same
									}
								}
								// change the allele also
								(*pCBSNP)->allele1 =a1;
								(*pCBSNP)->allele2 =a2;
								//only for test
								//cout <<"allele1: "<<(*pCBSNP)->allele1 <<endl ;//debug
								//cout <<"allele2: "<<(*pCBSNP)->allele2 <<endl ;//debug
							}
							else
							{
								cout <<"Allele in original files: " << 	(*pCBSNP)->allele1 << " " <<(*pCBSNP)->allele2 <<endl;
								cout <<"Allele in snpinfo files: " << 	a1 << " " <<a2 <<endl;
								string msg= " Allele1 \""+(*pCBSNP)->allele1+\
										"\" in first file  does not match with the allele \""+\
										a1+" given in "+change_int_into_string(nrow)+"th line of  SNP info file\""+fileName+\
										"\" \n.";
								error(msg);
							}
						}

					}
					else
					{
						(*pCBSNP)->allele1 =a1; // this condition never comes normally
						//cout << (*pCBSNP)->allele1 <<"  allele 1: \n" ; // debug
					}
					//second allele
												//allele2
												
					if((*pCBSNP)->allele2!="" && (*pCBSNP)->allele2!="0")
					{
						//----------------------------------------------//
						if((*pCBSNP)->allele2!=a2)
							{

								bool temp_bool1=(*pCBSNP)->allele2 ==a1;
								bool temp_bool2= (*pCBSNP)->allele1 ==a2;
								bool temp_bool3=(*pCBSNP)->allele1 ==a1;
								bool temp_bool4= ((*pCBSNP)->allele1=="");

								if(temp_bool1 & (temp_bool2||temp_bool3 || temp_bool4))
								{
									//cout << "change of allele! allele 1 is allele2 and allele 2 will be allele1 \n"; //ONly for test
									for(unsigned int i=0; i<(*pCBSNP)->geno1.size();++i)
									{
										//homo1->homo2
										if(!((*pCBSNP)->geno1[i])&& !((*pCBSNP)->geno2[i]))
										{
											(*pCBSNP)->geno1[i]=true;
											(*pCBSNP)->geno2[i]=true;
										}
										//homo2 ->homo1
										else if(((*pCBSNP)->geno1[i])&& ((*pCBSNP)->geno2[i]))
										{
											(*pCBSNP)->geno1[i]=false;
											(*pCBSNP)->geno2[i]=false;
										}
										else if(!((*pCBSNP)->geno1[i])&& ((*pCBSNP)->geno2[i]))
										{
											if((*pCBSNP)->aOrder[i])
												(*pCBSNP)->aOrder[i]=false;
											if(!((*pCBSNP)->aOrder[i]))
												(*pCBSNP)->aOrder[i]=true;
										}
										// change probability genotypes also
										if((*pCBSNP)->pgeno1.size()>0 && (*pCBSNP)->pgeno3.size()>0)
										{
											//cout << no_of_persons <<endl;
											for(unsigned int j=0;j<(*pCBSNP)->pgeno1.size(); ++j)
												swap((*pCBSNP)->pgeno1[j],(*pCBSNP)->pgeno3[j]);
													 //swap will exchange prob(AA) and prob(BB).prob(AB) remains the same
										}

									}
									// change the allele also
									(*pCBSNP)->allele1 =a1;
									(*pCBSNP)->allele2 =a2;
									//only for test
									//cout <<"allele1: "<<(*pCBSNP)->allele1 <<endl ;//debug
									//cout <<"allele2: "<<(*pCBSNP)->allele2 <<endl ;//debug


								}

							}
							//------------------------------------------------//
							if((*pCBSNP)->allele2 !=a2 && (*pCBSNP)->allele1 !=a2)
							{
								cout <<"Allele in original files: " << 	(*pCBSNP)->allele1 << " " <<(*pCBSNP)->allele1 <<endl;
								cout <<"Allele in snpinfo files: " << 	a1 << " " <<a2 <<endl;
								string msg= " Allele2 \""+(*pCBSNP)->allele2+\
										"\" in datfile  does not match with the allele \""+\
										a2+" given in "+change_int_into_string(nrow)+"th line of  SNP info file\""+fileName+\
										"\" \n.";
								error(msg);
							}
						}
						else
							(*pCBSNP)->allele2 =a2; // this condition never comes normally
					//----------------------------------------//
					//end of allele info stuff
					//----------------------------------------//
					temp_value			=cZeile[3];
					double temp_maf		=atof(temp_value.c_str());
					//	cout << "af: "<< temp_maf << endl;  // debug
					if((*pCBSNP)->maf<=0.0 )
					{
						
						(*pCBSNP)->maf		=min(temp_maf,(1.0-temp_maf));
						(*pCBSNP)->freq		=max(temp_maf,(1.0-temp_maf));
					
					}else
					{
						if((*pCBSNP)->maf !=temp_maf){
							
							cout <<"Minor allele frequency (maf) in original files: " << 	(*pCBSNP)->maf <<endl;
							cout <<"maf in snpinfo files: " << temp_maf <<endl;
							string msg= " maf  in original file  with the new one given in "+change_int_into_string(nrow)+"th line of  SNP info file\""+fileName+\
										"\" \n.";
								error(msg);
						
						}
					
					}
					// chromosome number
					temp_value			= (*pCBSNP)->nchr ;
					akt_int				= atoi(temp_value.c_str());
					temp_value			=cZeile[4];
					temp_int			= atoi(cZeile[4].c_str());
					//cout << "nChr: "<< temp_int << endl;  // debug 
					if(akt_int <=0 )
					{
						
						(*pCBSNP)->nchr		=temp_value;
					
					}else
					{
						if(akt_int !=temp_int){
							
							cout <<"Chromosome Number (nChr) in original files: " << 	(*pCBSNP)->nchr <<endl;
							cout <<"nChr in snpinfo files: " << temp_int <<endl;
							string msg= " Chromosome number  in original file  with the new one given in "+change_int_into_string(nrow)+"th line of  SNP info file\""+fileName+\
										"\" \n.";
								error(msg);
						
						}
					
					}
					//	
					// base pair position 
					if(cZeile.size()==6)
					{
						temp_value		=cZeile[5];
						// chromosome number
						//akt_int				= (*pCBSNP)->bp ;
						temp_int			= atoi(cZeile[5].c_str());
						akt_int			=(*pCBSNP)->bp;
						//cout <<"bp: "<<akt_int <<" "; //debug
						if(akt_int<=0 )
							(*pCBSNP)->bp	=temp_int;
						else
						{
							if((*pCBSNP)->bp !=temp_int)
							{

								cout <<"Base pair position(bp) in original files: " << 	(*pCBSNP)->nchr <<endl;
								cout <<"Base pair position(bp) in snpinfo files: " << temp_int <<endl;
								string msg= " bp  in original file  with the new one given in "+change_int_into_string(nrow)+"th line of  SNP info file\""+fileName+\
											"\" \n.";
									error(msg);
							
							}

						}
					}
					//decide quality of snps, quality 	
					(*pCBSNP)->quality = (*pCBSNP)->quality &&((*pCBSNP)->maf >=maf_value);
		
				//
				}
				
				//---------------------------------------------------------------//
				 //end of snpinfo saving else loop
				 //---------------------------------------------------------------//
				
				if(max_ncol!=ncol )
				{
					string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
					error(msg);
				}

			}// while loop
			else continue;
					 // assignment of genotypes
					//lese_pedfile_genInfo(genInfo, cZeile,nLines);
					//at the end des loops.
					if(nrow>1)
					{
						++pCBSNP;
					}
					 ncol=0;
					 cZeile.clear();
					 cZeile.resize(0);
					 if(_ifs_eof(pf))break;
		} // end of outer while loop
				_ifs_close(pf);
				nrow=0; // total rows
				nLines=0; // total lines
				//cout << genInfo.size()<<endl; // debug
			 }catch(bad_alloc& memoryAllocationException)
			 {
				 cout<<"\n Please inform the programmer that you could not read beagle file. Probably he will give you a solution. \n ";
				 error( "\n bad _alloc()  in reading beagle bgl file \n");
			 }

return;}



//-------------------------------------------------//
//
void CBIMPED::read_bimbam_gprobs_data(vector<CBPED*>& pedInfo, vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName, const float & thresh)
{

		printLIN("*->Reading file \""+ fileName+ "\": \n");
		string myString					="";
		 string emsg						="";
		 bool ende							= false;
		 unsigned int ccol					=0;
		 unsigned int ncol					=0;
		 unsigned int nrow					=0;
		 unsigned int nLines				=0;
		 unsigned int max_ncol				=0;
		 bool def_max_ncol_once				=true;
		 //unsigned long int	snpAnzahl		=0;
		 unsigned long int  nSample			=0;
		 //bool def_nSample					=true;
		 float first_prob					=0.0;
		 float second_prob					=0.0;
		 float third_prob					=0.0;
		 string pheno_value					="";
		 string sex_value					="";

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
					 //snpAnzahl=genInfo.size();
					 {
						  const string rsName=cZeile[0];
						   string allele_A=cZeile[1];
						   string allele_B=cZeile[2];
						  if(allele_A=="N")
						  {
							  allele_A=allele_B;
						  	  allele_B=allele_A;
						  }
						  cZeile.erase(cZeile.begin(),cZeile.begin()+3);
						  nSample=(cZeile.size())/2;
						  if(cZeile.size()==(2*nSample))
						  {
							  CBSNP* pCBSNP = new CBSNP;
							  if(pCBSNP->rsId=="")
								  pCBSNP->rsId=rsName;
							  if(pCBSNP->snpName=="")
								 pCBSNP->snpName=rsName;
							  if(pCBSNP->allele1=="")
								  pCBSNP->allele1=allele_A;
							  if(pCBSNP->allele2=="")
								  pCBSNP->allele2=allele_B;
							  genInfo.push_back(pCBSNP);
							  if(genInfo.size()==nrow)
								 {

								  CBSNP* pCBSNP=genInfo[genInfo.size()-1];
								  pCBSNP->no_of_persons=nSample;
									 for(unsigned int i=0; i<nSample; ++i)
									 {
										 first_prob  =atof( cZeile[2*i].c_str()); // first is homo1
										 second_prob = atof(cZeile[2*i+1].c_str()); // hetero
										 third_prob = 1.0 - (first_prob+ second_prob); // homo2
										 pCBSNP->pgeno1.push_back(first_prob);
										 pCBSNP->pgeno2.push_back(second_prob);
										 pCBSNP->pgeno3.push_back(third_prob);
									 }
								  }
								 else
								 {
									 emsg =" problem in saving genotypes line in file \""+fileName+\
											 "\".\n";
									 error(emsg);
								 }

						  }
						 else
						 {
							 emsg =" problem in saving genotypes line in file \""+fileName+\
									 "\".\n";
							 error(emsg);
						 }

					  }//end of genotype saving else loop

					 //CBPED* pCBPED= lese_pedfile_pedInfo(cZeile,snpAnzahl, nrow);
					 //pedInfo.push_back(pCBPED);
					 //cout << "nrow: "<< nrow<< " pedInfo.size()"<< pedInfo.size()<<endl;
					  if(max_ncol!=ncol)
					 {
						 string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
						 error(msg);
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
			update_pedInfo_given_gProb(pedInfo,genInfo);
			//string pmsg="\t**->Total number of successfully read SNPs: ";
			//printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
			//pmsg="\t**->Total number of successfully read individuals: ";
			//printLIN(pmsg);
			//printWithSpace(pmsg,  change_int_into_string(pedInfo.size()),coutWidth);
			//cout << genInfo.size()<<endl; // debug
			//changing genoprob into genotypes
			genoProb_to_geno_converter(genInfo, thresh);

		 }catch(bad_alloc& memoryAllocationException)
		 {
			 cout<<"\n Please inform the programmer that you could not read beagle file. Probbaly he will give you a solution. \n ";
			 error( "\n bad _alloc()  in reading beagle bgl file \n");
		 }
return ; }
// next



void CBIMPED::lese_bimbam_bestguess_genInfo(const vector<CBSNP*>&genInfo, const vector<string>& zeilenVec,unsigned const int akt_zeile)
{
	string first	="";
	string second	="";
	string third	="";
	bool fatal=false;
	string EMSG ="\t**->probblem with reading line "+change_int_into_string(akt_zeile)+":\n";
	string errmsg="";
	CBSNP* pCBSNP=genInfo[0];

	if(genInfo.size()==(akt_zeile))
	{
		pCBSNP=genInfo[genInfo.size()-1];
		for(unsigned int i=0;i<zeilenVec.size();++i)
		{
			  third=zeilenVec[i];
			//cout << third << endl; // only for test
			//cout << "/"<<endl;
			if(third.length() !=2) // || third[2]!="/")
			{
				errmsg="Problem in reading  SNP "+ pCBSNP->snpName+" in line"+change_int_into_string(akt_zeile)+". The SNP has genotype [ " + third + " ]. \n";
				error(errmsg);

			}
			// if okay now parse the string
			first=third[0];
			second =third[1];
			//cout << " " << first << "  " <<second <<" "<<endl; // for test only
			if((first!="?"&& second!="?") && ( first!=CBSNP::_miss_gvalue|| second!=CBSNP::_miss_gvalue))
					{
						//cout << "allele : "<< pCBSNP->allele1 << ", "<< pCBSNP->allele2 << " \n";
					if(first==second)
						{
							// add allele
							if(pCBSNP->allele1==""&& pCBSNP->allele2!=first)pCBSNP->allele1=first;
							else if (pCBSNP->allele2==""&& pCBSNP->allele1!=first) pCBSNP->allele2=first;
							// check if allele match with each other
							//cout << "allele 1: "<< pCBSNP->allele1<< " allele 2: "<<pCBSNP->allele2<< "\n";
							// add genotype
							if(first==pCBSNP->allele1&& second==pCBSNP->allele1 )
							{// this is the case where homozygote1 is denoted by 00
								//cout << " reading homo1 \n";
								//cout << first << " "<<second<< pCBSNP->allele1 << "\n";
								pCBSNP->geno1.push_back(false);
								pCBSNP->geno2.push_back(false);
								pCBSNP->aOrder.push_back(true);
							}
							else if(first==pCBSNP->allele2 && second==pCBSNP->allele2)
							{// this is the case where homozygote 2  is denoted by 11
								//cout << " reading homo2\n";// for test only
								//cout << first << " "<<second<< " "<< pCBSNP->allele1<< "\n";// for test only
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
				else if (first==CBSNP::_miss_gvalue|| second==CBSNP::_miss_gvalue)
				{ //In this case we use the code 10 to denote the missing genotypes
						//cout << "reading missing \n"; // for test only
						//cout << first << " "<<second<<"\n"; // for test only
						pCBSNP->geno1.push_back(true);
						pCBSNP->geno2.push_back(false);
						pCBSNP->aOrder.push_back(true);
				}
				else if (first=="?" || second=="?" )
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

		}

				// check if size of geno1 and geno2 are equal they must be equal.
					if(pCBSNP->geno1.size()!= pCBSNP->geno2.size() )
						error(EMSG);


		 //-------------------------------------------------------------#
	}
	else
	{
		string emsg ="problem in saving Genotypes in lese_bimbam_bestguess_genoInf function .\n";
		 error(emsg);
	 }

	// end of first if  loop
		//cout << "It seems like that you have  genotypes not in the form \"A  B\", but in the form \"A/B\". Please let us know if you want to convert \"A/B\" types of genotypes. We will update the program for you to convert such type of genotypes.\n ";
		//cout << "email:nab.raj.roshyara@imise.uni-leipzig.de \n";
		//exit(1);
		// this is the case where snps info is given as  AB

 return ;
}

// final


void CBIMPED::read_bimbam_bestguess_data(vector<CBPED*>& pedInfo, vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName)
{

		printLIN("*->Reading file \""+ fileName+ "\": \n");
		string myString					="";
		 string emsg						="";
		 bool ende							= false;
		 unsigned int ccol					=0;
		 unsigned int ncol					=0;
		 unsigned int nrow					=0;
		 unsigned int nLines				=0;
		 unsigned int max_ncol				=0;
		 bool def_max_ncol_once				=true;
		 //unsigned long int	snpAnzahl		=0;
		 unsigned long int  nSample			=0;
		 bool def_nSample					=true;
		// float first_prob					=0.0;
		 //float second_prob					=0.0;
		 //float third_prob					=0.0;
		 string pheno_value					="";
		 string sex_value					="";

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

						  const string rsName=cZeile[0];
						  cZeile.erase(cZeile.begin(),cZeile.begin()+1);
						  if(def_nSample)
							  nSample=(cZeile.size());
						  def_nSample=false;
						  if(cZeile.size()==nSample)
						  {
							  CBSNP* pCBSNP = new CBSNP;
							  pCBSNP->no_of_persons=nSample;
							  if(pCBSNP->rsId=="")
								  pCBSNP->rsId=rsName;
							  if(pCBSNP->snpName=="")
								 pCBSNP->snpName=rsName;
							  genInfo.push_back(pCBSNP);
							 // cout << "test1"; //debug
							  lese_bimbam_bestguess_genInfo(genInfo, cZeile,nrow);
							 // cout << "test2"; //debug

						  }
						 else
						 {
							 emsg =" problem in saving genotypes line in file \""+fileName+\
									 "\".\n";
							 error(emsg);
						 }

					  if(max_ncol!=ncol)
					 {
						 string msg= "Not enough columns. Only "+ change_int_into_string(ncol) +" columns were recognized in the line of "+change_int_into_string(nLines)+"th row of   file \""+fileName+"\".\n";
						 error(msg);
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
			update_pedInfo_given_gProb(pedInfo,genInfo);
			//string pmsg="\t**->Total number of successfully read SNPs: ";
			//printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
			//pmsg="\t**->Total number of successfully read individuals: ";
			//CBPED::display_ped_summary();
			//printLIN(pmsg);
			//printWithSpace(pmsg,  change_int_into_string(pedInfo.size()),coutWidth);
			//cout << genInfo.size()<<endl; // debug
			//changing genoprob into genotypes

		 }catch(bad_alloc& memoryAllocationException)
		 {
			 cout<<"\n Please inform the programmer that you could not read beagle file. Probbaly he will give you a solution. \n ";
			 error( "\n bad _alloc()  in reading beagle bgl file \n");
		 }
return ; }
























