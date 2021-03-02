/*
 * locus.cpp
 *
 *  Created on: 14.11.2011
  *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */
#include "plinksnp.h"
#include "helper.h"

CPSNP::~CPSNP()
{
	/*
	for (unsigned int i=0;i<genInfo.size();++i)
	{
		delete genInfo[i];
		
	}
	genInfo.clear();
	*/	
 	//cout << "\n Good bye from CBSNP destructor.\n";
}
bool CPSNP::_given_three_cols=false;

bool CPSNP::checkSNPinfo(const vector<string> & snp, const long int& countSNPs)
{
	// This function will return true if there are four columns and if not then false:
	bool tfValue=false;
	if((snp.size()!=3) && (snp.size()!=4)){
			error("It is expected to have either 3 or 4 columns in "\
				+ change_int_into_string(countSNPs)+\
				"th line (excluding comments and blank lines) of your map file, but this is not the case here. Please either correct or delete this line. \n");

		}
	if(!_given_three_cols && (snp.size()==3))
		_given_three_cols=true;
				// you check a lot of things for the snp info here
	else if(!_given_three_cols && snp.size()==4)
			tfValue=true;
				// you check a lot of things for the snp info here
return tfValue;}

CBSNP* CPSNP::fillUpLocusInfo(const vector<string>& lines)
{
	CBSNP* pCBSNP=new CBSNP;
	if( _given_three_cols && lines.size()==3)
		{
			pCBSNP->nchr=lines[0];
			pCBSNP->snpName=lines[1];
			if(pCBSNP->rsId=="")
				pCBSNP->rsId=lines[1];
			//pCBSNP->rsId=lines[1];
			pCBSNP->bp=(long int )atoi(lines[2].c_str());
		}
		else if(lines.size()==4)
		{
			pCBSNP->nchr=lines[0];
			pCBSNP->snpName=lines[1];
			if(pCBSNP->rsId=="")
				pCBSNP->rsId=lines[1];
			//pCBSNP->rsId=lines[1];
			pCBSNP->cm_pos=atof(lines[2].c_str());
			pCBSNP->bp=(long int )atoi(lines[3].c_str());

		}
	return pCBSNP;
}

void readMapFile1(const string&  fileName,  vector<CBSNP*>& genInfo)
{
	string line="test";
	rFile ftmp(line );
	const std::string  delim=" ";
	while(!ftmp.eof())
	{

	}
	ftmp.close();
	ftmp.clear();

return ;
}

//----------------------------------------------------------------------------//
void CPSNP::readMapFile(const string&  fileName,  vector<CBSNP*>& genInfo){
	// counting number of SNPs in the mapfile
	static long int countSNPs=1;
	static vector<bool>checkFileColumns(2); // to check number of columns in each row.
	string r_msg="\n*->Reading map file \""+fileName+ "\":\n";
	printLIN(r_msg);
	string line;
	rFile file ;
	//rFile file;
	file.close(); file.clear(); // check if it has been already opend.
	//ios_base::iostate i= file.rdstate();
	//file.open(fileName.c_str());
	// this will check if each row has same number of columns.

	//file(fileName);
	file.open(fileName.c_str(),ios::in);
		if(!file)
		{
			file.setstate(ios_base::failbit);
			file.clear();
			error("your map file "+fileName+ " either does not exit or could not be opend.");
		}

	while(getline(file,line, '\n')) // until eof hit
	{
		//file >> ws; // leave white spaces;
		if (line=="") // leave blank lines
			continue;
		if(line[0]=='#') // leave the comments;
			continue;
		string buffer;
		vector<string> tokens;
		stringstream line_parser(line);
		// check if parsing fails
		// if parsing is good then go further
		if(line_parser.good())
		{
			while(line_parser>>buffer)
			tokens.push_back(buffer);
			//cout << "number of map file columns: " << tokens.size()<<endl;
			bool  mapFileColumns_3or4=checkSNPinfo(tokens, countSNPs);
			//if true then 4 columns if not then three.
			// we have to check if there are 3 or four columns
			// and if there are same columns in every row.
			if(checkFileColumns.empty())
				checkFileColumns.push_back(mapFileColumns_3or4);
			checkFileColumns.push_back(mapFileColumns_3or4);
			if(checkFileColumns[0]!=checkFileColumns[1])
				error(" The " + change_int_into_string(countSNPs)+"th row (excluding comments and blank lines)"+
					"of your map file does not contain the same number of columns as in the previous rows. Please delete or correct this row.\n");
			checkFileColumns.pop_back();
	 		// checking 3 or 4 columns is finished.
			//Now saving each mapfile row as dynamic allocation
			CBSNP* value = fillUpLocusInfo(tokens);
			//cout<<*value<<endl;
			genInfo.push_back(value);
		}

	++countSNPs	;
	}

	// writing what if !file.bad();
	//file.exceptions(file.exceptions()|ios_base::badbit);
	//informing how many snps has successfully been read.
	//display_snp_summary(genInfo);
	//cout <<"\t"<< setw(col_width);
	//printLIN("**->Total number of successfully read SNPs: "+ change_int_into_string(countSNPs-1)+"\n ");

return ;
}


void CPSNP::snpinfoFileLeser(const string& fileName, const vector<CBSNP*> & genInfo)
{
	string r_msg="\n*->Reading SNP info  file \""+fileName+ "\":\n";
	printLIN(r_msg);
	string line;
	rFile file;
	file.close(); file.clear(); 	// check if it has been already opend.

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
		const unsigned int dvalue	=0;
		unsigned int col_nchr		=dvalue; // initialization of columns
		unsigned int col_rsid		=dvalue;
		unsigned int col_snpid		=dvalue;
		unsigned int col_bp			=dvalue;
		unsigned int col_cm_pos		=dvalue;
		unsigned int nCols			=dvalue;
		unsigned int col_allele1	=dvalue;
		unsigned int col_allele2	=dvalue;

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
			//rsid
			bool found	= binary_search(sort_header.begin(), sort_header.end(),"rsid");
			bool found1	= binary_search(sort_header.begin(), sort_header.end(),"snpid");
			if(found|| found1)
			{
				if(found)
				{
					//rsid
					it=find(header.begin(),header.end(),"rsid");
					col_rsid= (int)(distance(header.begin(),it)+1);
				}
				if(found1)
				{
					it=find(header.begin(),header.end(),"snpid");
					col_snpid= (int)(distance(header.begin(),it)+1);
				}

			}
			//else
			//{
			//	string msg="your extra map file \""+fileName+"\"  neither  possesses acolumn containing rsid nor a column containing snpid.\n";
			//	error(msg);
			//}
			found= binary_search(sort_header.begin(), sort_header.end(),"nchr");
			if(found)
			{	// nchr
				it=find(header.begin(),header.end(),"nchr");
				col_nchr= (int)(distance(header.begin(),it)+1);
			}
			found= binary_search(sort_header.begin(), sort_header.end(),"bp");
			if(found)
			{
			// bp
				it=find(header.begin(),header.end(),"bp");
				col_bp= (int)(distance(header.begin(),it)+1);
			}
			found= binary_search(sort_header.begin(), sort_header.end(),"cm_pos");
			if(found)
			{	// cm_pos
				it=find(header.begin(),header.end(),"cm_pos");
				col_cm_pos= (int)(distance(header.begin(),it)+1);
			}
			//allele1
			 found= binary_search(sort_header.begin(), sort_header.end(),"allele1");
						if(found)
						{	// allele1
							it=find(header.begin(),header.end(),"allele1");
							col_allele1= (int)(distance(header.begin(),it)+1);
						}
						//allele2
			 found= binary_search(sort_header.begin(), sort_header.end(),"allele2");
			 if(found)
			 {	// allele2
				 it=find(header.begin(),header.end(),"allele2");
				 col_allele2= (int)(distance(header.begin(),it)+1);
			 }


			//displaying for test
				//cout << "col_rsId: "<< col_rsid << " "<< "col_snpId: "<< col_snpid << " "<< "col_nchr: "<< col_nchr << " "<< "col_bp: "<< col_bp << " " << "col_cm_pos: "<< col_cm_pos << " ";


		}
		else
		{
			string msg= "Problem in parsing the first line!\n";
			error(msg);

		}
		++extraMapFileLineZaehler;
		vector<string> tokens;
		vector<CBSNP*>::const_iterator it_pCBSNP=genInfo.begin();

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

					//rsid
					if(col_rsid!=dvalue)
					{
						if((*it_pCBSNP)->rsId!="")
						{
							if((*it_pCBSNP)->rsId !=tokens[col_rsid-1])
							{
								if((*it_pCBSNP)->rsId !=tokens[col_snpid-1])
								{
									string msg= "rs id  \""+(*it_pCBSNP)->rsId+"\" in datfile  does not match with rs id \""+\
									  tokens[col_rsid-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra map file\""+fileName+\
									  "\" \n.";
									error(msg);
								}
								else
								 (*it_pCBSNP)->rsId =tokens[col_rsid-1];
							}

						}
						else
							(*it_pCBSNP)->rsId =tokens[col_rsid-1];
					}
					//snpid

					if(col_snpid!=dvalue)
					{
						if((*it_pCBSNP)->snpName!="")
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
					if(col_nchr!=dvalue)
					{
						if((*it_pCBSNP)->nchr!="0")
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
					//bp
					if(col_bp!=0)
					{
						if((*it_pCBSNP)->bp!=(-1) )
						{
							if((*it_pCBSNP)->bp !=atoi(tokens[col_bp-1].c_str()))
							{
								string msg= " Base pair position \""+change_int_into_string((*it_pCBSNP)->bp)+
										"\" in datfile  does not match with base pair position \""+\
										tokens[col_bp-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of extra map file\""+fileName+\
										"\" \n.";
								error(msg);
							}
						}
						else
						(*it_pCBSNP)->bp =atoi(tokens[col_bp-1].c_str());
					}
					//cm_pos
					if( col_cm_pos!=dvalue)
					{
						if((*it_pCBSNP)->cm_pos!=-1.0)
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

					//allele1
					if( col_allele1!=dvalue)
					{
						if((*it_pCBSNP)->allele1!="" || (*it_pCBSNP)->allele1!="0")
						{
							if((*it_pCBSNP)->allele1 !=tokens[col_allele1-1])
							{
								string msg= " Allele1 \""+(*it_pCBSNP)->allele1+\
										"\" in datfile  does not match with the allele \""+\
										tokens[col_allele1-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of  SNP info file\""+fileName+\
										"\" \n.";
								error(msg);
							}
						}
						else
							(*it_pCBSNP)->allele1 =tokens[col_allele1-1];
					}
					//allele2
					if( col_allele2!=dvalue)
					{
						if((*it_pCBSNP)->allele2!="" ||(*it_pCBSNP)->allele2!="0" )
						{
							if((*it_pCBSNP)->allele2 !=tokens[col_allele2-1])
							{
								string msg= " Allele2 \""+(*it_pCBSNP)->allele2+\
									"\" in datfile  does not match with the allele \""+\
										tokens[col_allele2-1]+" given in "+change_int_into_string(extraMapFileLineZaehler)+"th line of  SNP info file\""+fileName+\
										"\" \n.";
								error(msg);
							}
						}
						else
							(*it_pCBSNP)->allele2 =tokens[col_allele2-1];
					}

					//vector<CBSNP*>gifo=genInfo;
					//sort(gifo.begin(),gifo.end());
					//bool b= binary_search(gifo.begin()	,gifo.end(),"snp1",CompareString());
					//cout <<"binary_search: "<<endl;
					//cout << b<<endl;
					//cout<< (int)*upper_bound(gifo.begin(),gifo.end(),"snp1",CompareString())<<endl;
					//cout << lower_bound(gifo.begin(),gifo.end(),"snp1",CompareString())<<endl;

				}
				else
				{
					printLIN( "Line parser not good while reading extra snp info file.\n");
					break;
				}
				if(it_pCBSNP==(genInfo.end()-1)) break;
				else
				++it_pCBSNP;
				++extraMapFileLineZaehler;

		}


		// check if parsing fails
		// if parsing is good then go further


		file.close();



return ;}

//end of plink




