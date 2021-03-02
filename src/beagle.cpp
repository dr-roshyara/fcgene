/*
 * beagle.cpp
 *
 *  Created on: Jan 6, 2012
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */
#include"beagle.h"
CBEPED::~CBEPED()
{	/*
	for (unsigned int i=0;i<vPed.size();++i)
	{
		delete vPed[i];
	}
		vPed.clear();
		vPed.resize(0);
	*/	
}
CBESNP::~CBESNP()
{
	/*for (unsigned int i=0;i<genInfo.size();++i)
	{
		delete genInfo[i];
	}
		genInfo.clear();
		genInfo.resize(0);
	*/	
}

bool CBESNP::_given_gprobs_header=true;
void CBEPED::read_bgl_gprobs_data(vector<CBPED*>& pedInfo, vector<CBSNP*>& genInfo,  vector<string> cZeile,const string fileName, const float & thresh)
{
	//-----------------------------------------------------------------//
		printLIN("*->Reading file \""+fileName+"\": \n");
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
		 float first_prob					=0.0;
		 float second_prob					=0.0;
		 float third_prob					=0.0;
		 string pheno_value					="";
		 string sex_value					="";
		 bool	go_further					=false;
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
						// cout << nrow<<endl;
						// cout <<"cZeile[0]: "<< cZeile[0]<<endl;
					 if(CBESNP::_given_gprobs_header&& (nrow==1))
					 {
						 //save person id ;
						 // cout << "first: "<<boolalpha << ((cZeile[0]!="Marker") &&cZeile[0]!="marker" && cZeile[0]!="SNP") <<endl;
						 //cout << (cZeile[1]!="alleleA" && cZeile[1]!="A1")<< endl;
						 //cout << (cZeile[1]!="alleleB" && cZeile[1]!="A2")<< endl;
						 if(((cZeile[0]!="Marker")&& cZeile[0]!="marker" && cZeile[0]!="SNP") || (cZeile[1]!="alleleA" && cZeile[1]!="A1") || (cZeile[2]!="alleleB"&& cZeile[2]!="A2"))
						 {
							 emsg =" First line of Beagle format file \""+fileName+\
									 "\" should start with the letter \"marker or Marker\", the second column should contain  the word \"alleleA\" and the thrid  column should contain  the word \"alleleB\" .\n";
							 error(emsg);
						 }
						 else
						 {
							 cZeile.erase(cZeile.begin(),cZeile.begin()+3);
							 if(def_nSample)
							 {
								 nSample=cZeile.size()/3;
								 def_nSample=false;

							 }

							 for(unsigned int i=0; i<nSample; ++i)
							 {
								 if((cZeile[3*i]== cZeile[3*i+1] ) &&(cZeile[3*i]== cZeile[3*i+2] ))
								 {
									 CBPED* pCBPED = new CBPED;
									 pCBPED->indId=cZeile[3*i];
									 pedInfo.push_back(pCBPED);
									 //cout << "pedInfo.size(): "<< pedInfo.size()<<endl; // debug
								}
								else
								{
									 emsg ="On the "+change_int_into_string(3*i)+"th, "+change_int_into_string(3*i+1)+"th and "+change_int_into_string(3*i+2)+"th columns of  first line in file \""+fileName+\
											 "\" should contain the same individual id but this is the case here.\n";
									 error(emsg);
								}
							 }
						 }
					 }
					 else
					 {
						 const string rsName=cZeile[0];
						 const string allele_A=cZeile[1];
						 const string allele_B=cZeile[2];

						 cZeile.erase(cZeile.begin(),cZeile.begin()+3);
						 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
						 // in the case of having no header
						 if(!CBESNP::_given_gprobs_header&&!pedInfo.size()&&(nrow==1))
						 {
							 //CBESNP::_given_gprobs_header=true;
							 if(def_nSample)
							 {
								nSample=cZeile.size()/3;
								 def_nSample=false;

							 }
							 for(unsigned int i=0; i<nSample; ++i)
							 {
								 CBPED* pCBPED = new CBPED;
								 // pCBPED->indId=cZeile[3*i];
								 pedInfo.push_back(pCBPED);

							 }
						 }
						 //xxxxxxxxxxxxxxxxxxxxxxx
						 //cout << cZeile.size() <<", "<< nSample<<endl;
						 if(cZeile.size()==(3*nSample))
						  {
							  CBSNP* pCBSNP = new CBSNP;
							  if(pCBSNP->rsId=="")
								  pCBSNP->rsId=rsName;
							  if(pCBSNP->snpName=="")
								 pCBSNP->snpName=rsName;
							  if((pCBSNP->allele1=="")&&(allele_A!="-"))
								  pCBSNP->allele1=allele_A;
							  if(((pCBSNP->allele2=="")&&(allele_A!=allele_B)&&(allele_B!="-")))
								  pCBSNP->allele2=allele_B;
							  genInfo.push_back(pCBSNP);
							  if(CBESNP::_given_gprobs_header)
								  go_further =(genInfo.size()==(nrow-1));
							  else
								  go_further =(genInfo.size()==nrow);
							  if(go_further) // -1 means first row was indiv info.
								 {

								  	  CBSNP* pCBSNP=genInfo[genInfo.size()-1];
								  	  pCBSNP->no_of_persons=nSample;
								  	  bool temp_summ 	=false;
								  	  bool temp_boole	=false;
								  	  bool _tmp_miss 	=false;
									 for(unsigned int i=0; i<nSample; ++i)
									 {
										// first		=atof(zeilenVec[3*i-3].c_str());
										// second		=atof(zeilenVec[3*i-2].c_str());
										// third 		=atof(zeilenVec[3*i-1].c_str());
										 first_prob  	=atof(cZeile[3*i].c_str());
										 second_prob 	=atof(cZeile[3*i+1].c_str());
										 third_prob 	=atof(cZeile[3*i+2].c_str());
										 temp_summ		=first_prob+second_prob+third_prob;
										 temp_boole		=(first_prob<0.0 ||first_prob>1.0001) || (second_prob<0.0 ||second_prob>1.0001) ||(third_prob<0.0 ||third_prob>1.0001);
										 _tmp_miss  	=(temp_summ==0.0);
										 if((temp_summ!=1 || temp_boole)&& !_tmp_miss)
										 {
											 if(temp_summ!=1.0 && temp_summ<1.1 )
											 {
												 first_prob 	= first_prob/temp_summ;
												 second_prob	= second_prob/temp_summ;
												 third_prob 	= third_prob/temp_summ;
											 }
											 else
											 {
												 printLIN("The probabality distribution of SNP: "+pCBSNP->rsId+" is not correct. \n" );
												 cout <<"Three probabilities: "<< first_prob <<" "<< second_prob <<" "<<third_prob <<" \n";
												 error("Please check your file containing  genotype probabilities\n");
											 }
										 }
										 //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
										 //cout <<first_prob <<" "<< second_prob << " "<< third_prob << " \n";
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
									 "\". There are not enough columns for genotype data.\n";
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
			string pmsg="\t**->Total number of successfully read SNPs: ";
			printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
			pmsg="\t**->Total number of successfully read individuals: ";
			//------------------------------------------------//
			//CBPED::display_ped_summary();
			//printLIN(pmsg);
			printWithSpace(pmsg,  change_int_into_string(pedInfo.size()),coutWidth);
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

void CBEPED::write_beagle_gprobsFile(const vector<CBPED*>&pedInfo, const vector<CBSNP*>genInfo, const string& outFileName)
{
	//cout <<"hi\n";
	ofstream ofs;
	ofs.clear();ofs.close();
	CBPED* pCBPED =pedInfo[0];
	CBSNP* pCBSNP	=genInfo[0];
	string bglFileName=outFileName+".bgl.gprobs";
	ofs.open(bglFileName.c_str());
	if(!ofs)
		error("out file "+outFileName+" does not exit or can not be written out.\n");
	ofs.clear();
	ofs<< "marker" <<" "<< "alleleA"<<" "<<"alleleB";
	for(unsigned int i=0; i<pedInfo.size();++i)
	{
		pCBPED =pedInfo[i];
		if(pCBPED->quality)
		{
			if(pCBPED->indId!="")
				ofs <<" "<<pCBPED->indId<< " "<<pCBPED->indId<< " "<<pCBPED->indId;
			else
				ofs <<" "<<"indv"<<i+1<< " "<<"indv"<<i+1<< " "<<"indv"<<i+1<< " ";
		}
	}
	ofs<<"\n";
		//end of first line
		//ofs<<"A"<<"  "<< "phenotype" << " ";
		//end of second line.
	//writing genotypes
		const unsigned int _SZ =pedInfo.size();
	for(unsigned int i=0;i<genInfo.size();i++)
	{
		pCBSNP	=genInfo[i];
		if(pCBSNP->quality)
		{

			if((_SZ!=pCBSNP->pgeno1.size())||(_SZ!=pCBSNP->pgeno2.size())||(_SZ!=pCBSNP->pgeno3.size()))
			{
				cout<< "pCBSNP->pgeno1.size(): "<<pCBSNP->pgeno1.size()<<", pCBSNP->pgeno2.size(): "<<pCBSNP->pgeno2.size()<<", pCBSNP->pgeno3.size(): "<<pCBSNP->pgeno3.size()<<" \n";
				string msg= " Problem in writing line"+change_int_into_string(i)+"! SNP \""+pCBSNP->snpName+"\" was not saved correctly while reading. Different number of individual alleles\n";
				error(msg);
			}
			if(pCBSNP->rsId!="")
				ofs<< pCBSNP->rsId ;
			else
			ofs<< pCBSNP->snpName ;
			ofs<<" "<<pCBSNP->allele1;
			if(pCBSNP->allele2!="")
			ofs<<" "<<pCBSNP->allele2;
			else
				ofs<<" "<<"-";
			for(unsigned int j=0; j<_SZ;++j)
			{
				pCBPED=pedInfo[j];
				if(pCBPED->quality)
					ofs<<" "<<pCBSNP->pgeno1[j] << " "<<pCBSNP->pgeno2[j] << " "<<pCBSNP->pgeno3[j];
			} //end of for loop j
					ofs << "\n";
		}

	} //end of for i loop
	ofs.clear(); ofs.close();

	string msg="*->Gentoype probability distribution data for Beagle-input has been saved as \""+bglFileName+".\".\n";
		 printLIN(msg);
return;	}
//writing bimbam output geno file
//--------------------------------------------------//



//-------------------------------------------------//
void CBEPED::read_bgl_data(vector<CBPED*>&pedInfo, vector<CBSNP*>&genInfo, vector<string> cZeile,const string fileName)
{
	//-----------------------------------------------------------------//
		printLIN("*->Reading file \""+fileName+"\": \n");
			 string myString					="";
			 string emsg						="";
			 bool fatal							=false;
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
			 string first_geno					="";
			 string second_geno					="";
			 string pheno_value					="";
			 string sex_value					="";
			 //bool  def_nperson_once				true;

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
						// snpAnzahl=genInfo.size();
						 if(nrow<3)
						 {
							 if(nrow==1 ||((cZeile[0]!="I" || cZeile[1]!="id")&& nrow<=1))
							 {
									//save person id ;
								 if(cZeile[0]!="I" || cZeile[1]!="id")
								 {
									 emsg =" First line of Beagle format file \""+fileName+\
											 "\" should start with the letter \"I\" and the second column should contain  the word \"id\" .\n";
									error(emsg);
								 }
								 else
								 {
									 cZeile.erase(cZeile.begin(),cZeile.begin()+2);
									 if(def_nSample) nSample=cZeile.size()/2;
									 def_nSample=false;
									 for(unsigned int i=0; i<nSample; ++i)
									 {
				 						 if(cZeile[2*i]== cZeile[2*i+1])
				 						 {
				 							 CBPED* pCBPED = new CBPED;
										 	 pCBPED->indId=cZeile[2*i];
										     pedInfo.push_back(pCBPED);
										     //cout << "pedInfo.size(): "<< pedInfo.size()<<endl; // debug
				 						 }
				 						 else
				 						 {
											 emsg ="On the "+change_int_into_string(2*i)+"th and "+change_int_into_string(2*i+1)+"th columns of  first line in file \""+fileName+\
													 "\" should contain the same individual id but this is the case here.\n";
											error(emsg);
				 						 }
									 }
								 }
							}
							 else if(nrow==2||( (cZeile[0]!="M"||cZeile[0]!="A" || cZeile[1]!="phenotype") && nrow<=2))
							 {
								//cout << cZeile[0] <<" " << cZeile[1]<<endl;
								//cout << boolalpha <<(cZeile[0]=="A"&& cZeile[1]=="phenotype" )<<endl;
								if(cZeile[0]=="A"&& cZeile[1]=="phenotype")
								 {
									 cZeile.erase(cZeile.begin(),cZeile.begin()+2);
									 if(cZeile.size()!=(2*nSample))
									 {
										 emsg ="The second line of file \""+fileName+\
												 "\" should contain"+change_int_into_string((2*nSample))+ " disease information columns, but this is the case here.\n";
										 error(emsg);
									 }
									 else
									{
										 for(unsigned int i=0; i<nSample; ++i)
										 {
											 CBPED* pCBPED = pedInfo[i];
											 pheno_value = cZeile[2*i];
											 CBPED::lese_phenotype_info(pCBPED,pheno_value);
											 CBPED::lese_sex_info( pCBPED, sex_value);

										 }

									}
								}
								else 
								{ 
									//------------------------------------------------------------------------------------------------------------------------//
									 //save genotypes.
									// cout <<"nth row: " <<nrow << " cZeile: "<< cZeile[0]<< endl; // debug 
									if(cZeile[0]!="M")
										cout << "  There is no M in first column !\n";
									const string rsName=cZeile[1];
									cZeile.erase(cZeile.begin(),cZeile.begin()+2);
									if(cZeile.size()==(2*nSample))
									{
										CBSNP* pCBSNP = new CBSNP;
										pCBSNP->no_of_persons=nSample;
										if(pCBSNP->rsId=="")
											  pCBSNP->rsId=rsName;
										if(pCBSNP->snpName=="")
											pCBSNP->snpName=rsName;
										genInfo.push_back(pCBSNP);
										{
												 CBSNP* npCBSNP=genInfo[genInfo.size()-1];
												for(unsigned int i=0; i<nSample; ++i)
												 {
													 first_geno  =cZeile[2*i];
													 second_geno = cZeile[2*i+1];
													 fatal=lese_genotype_info(npCBSNP, first_geno, second_geno,nrow, emsg);
													 if(fatal)
													 {
														 //delete pCBSNP;
														 error(emsg);
													 }

												}
										}
										

									}
									else
									{
										 emsg ="Problem in saving genotypes line in file \""+fileName+\
												 "\".\n";
										 error(emsg);
									}	 
								}
								//----------------------------------------------------------------------------------------------------------------------//
								
							}
		 				 }
						 else
						 {
							 //save genotypes.
							// cout << "cZeile[0]: "<< cZeile[0]<< endl; // debug 
							 if(cZeile[0]!="M")
							 {
								cout << "  There is no M in first column !\n";
							 }
							  const string rsName=cZeile[1];
							  cZeile.erase(cZeile.begin(),cZeile.begin()+2);
							  //cout << "nth row "<< nrow << ", cZeile.size()"<< cZeile.size() <<endl ; //debug 
							  if(cZeile.size()==(2*nSample))
							  {
								  CBSNP* pCBSNP = new CBSNP;
								  pCBSNP->no_of_persons=nSample;
								  if(pCBSNP->rsId=="")
									  pCBSNP->rsId=rsName;
								  if(pCBSNP->snpName=="")
									 pCBSNP->snpName=rsName;
									 genInfo.push_back(pCBSNP);

									 //if(genInfo.size()==(nrow-2))
									 {
										 CBSNP* npCBSNP=genInfo[genInfo.size()-1];
										 for(unsigned int i=0; i<nSample; ++i)
										 {
											 first_geno  =cZeile[2*i];
											 second_geno = cZeile[2*i+1];
											 fatal=lese_genotype_info(npCBSNP, first_geno, second_geno,nrow, emsg);
											 if(fatal)
											 {
												 //delete pCBSNP;
												 error(emsg);
											 }

										 }
									  }

								}
							 else
							 {
								 cout << "Total genotype elements in "<<nrow<<"th line: " << cZeile.size()<<endl;
								cout << "Total elements which are necessary: "<< (2*nSample) <<endl; 	
								emsg ="Problem in saving genotypes of "+change_int_into_string(nrow)+"th line of  file \""+fileName+\
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
					 // now erase the first 6 columns of pedInfo from ZeilenVector
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
				string pmsg="\t**->Total number of successfully read SNPs: ";
				printWithSpace(pmsg,  change_int_into_string(genInfo.size()),CBPED::coutWidth);
				pmsg="\t**->Total number of successfully read individuals: ";
				printWithSpace(pmsg,  change_int_into_string(pedInfo.size()),coutWidth);
				//cout << genInfo.size()<<endl; // debug
			 }catch(bad_alloc& memoryAllocationException)
			 {
				 cout<<"\n Please inform the programmer that you could not read beagle file. Probably he will give you a solution. \n ";
				 error( "\n bad _alloc()  in reading beagle bgl file \n");
			 }

return;}
// read beagle rsq-info file
void CBESNP::lese_beagle_rsq_score_file(vector<CBSNP*>&genInfo,vector<string>&cZeile, const double& rsq_thresh, const string& fileName) {
//-----------------------------------------------------------------//
	printLIN("*->Reading file \""+fileName+"\": \n");
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
	string emsg							="";
	//long int temp_bp					=0;
	//double  temp_info					=0.0;
	double  temp_thresh					=0.0;
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
				int spaltenZahl=cZeile.size();
				//-----------------------------------------------------------------//
				bool tfValue =(spaltenZahl==2 ); //||(spaltenZahl==7 ) ;
				if(!tfValue)	// Fall mit 2 Zeilen.
				{
					emsg="\t**->Problem on the " + change_int_into_string(nrow)+ "th row.\n";
					printLIN(emsg);
					emsg="There must be 2 columns in each line of file: "+fileName+" but there  are "+change_int_into_string(spaltenZahl)+" columns.\n";
					//delete pCBPED;
					error(emsg);

				}else
				{
					if((*it_CBSNP)->snpName!="" && cZeile[0]!="")
					{
						if((*it_CBSNP)->snpName!=cZeile[0])
						{
							emsg ="on the"+change_int_into_string(nrow)+"th row of  file  \""+fileName+\
									"\" SNP Name"+cZeile[0]+ " does not match with SNP Name "+(*it_CBSNP)->snpName+ " contained in original file.\n";
							error(emsg);
						}
						else
							(*it_CBSNP)->snpName=cZeile[0];

					}
					//for freq1
					//for freq1
					if(cZeile[1]!="")
					{
						if(cZeile[1]=="NaN")
							temp_thresh=0.0;
						else
						temp_thresh	=atof(cZeile[1].c_str());
					}
					if(temp_thresh<=rsq_thresh)
					{
						(*it_CBSNP)->quality=false;
						//cout <<temp_info<< " "<< info_thresh << " "<<( temp_info<=info_thresh)<< endl; // debug
						//cout <<temp_freq<<" "<<(*it_CBSNP)->maf<< " "<< maf_thresh << " "<<( (*it_CBSNP)->maf<=maf_thresh)<< endl; // debug

						//++count_low_info;
					}

					++it_CBSNP;



				}
			}
			else continue;
			ncol=0;
			cZeile.clear();
			cZeile.resize(0);
			if(_ifs_eof(pf))break;
		} // end of outer while loop
		_ifs_close(pf);
		nrow=0;
		nLines=0;
	}catch(bad_alloc& memoryAllocationException)
	{
		cout<<"\n Please inform the programmer that you could not read mlgeno file. Probbaly he will give you a solution. \n ";
		error( "\n bad _alloc()  in reading mach mlgeno file \n");
	}


//------------------------------------------------//
return ;}




