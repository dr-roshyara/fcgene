#include"machped.h"
// lese plink formatted covariate file:
 void CBPED::lese_plink_covariate_file(const vector<CBPED*>& pedInfo,const string& fileName)
 {
	 static long int countIndiv=1;
	 static vector<bool>checkFileColumns;
	 string pr_msg="\n*-> Reading  Covariate  file: \""+fileName+ "\":\n";
	 printLIN(pr_msg);
	const int pTotal		=pedInfo.size();
	string line				=""; 
	string buffer			="";
	unsigned int nCols		=0;
	vector<string> tokens;

	rFile file; file.close();file.clear();
	 file.open(fileName.c_str(),ios::in);
	 if(!file)
	 {
		file.setstate(ios_base::failbit);
		file.clear();
		error("Your map file "+fileName +" either does not exit or could not be opend!\n");
	 
	 }
	 //read the fist rwo of the file which should contain
	 // covariate names
	 do
	 {
		getline(file,line,'\n');
		if(line[0]=='#')
			continue;
		if(file.eofbit) break;
	 }while(line.empty());
	 stringstream line_parser(line);
	 if(line_parser.good())
	 {
		while(line_parser>>buffer)
		CBPED::covariate_names.push_back(buffer);
		nCols=covariate_names.size();
		printLIN("\t**->total covariates: "+change_int_into_string((nCols-2))+" \n");
		int tmp_int =0;
		if(_given_cov_names)
			tmp_int =_cov_names_frm_argv.size();
		else 
			tmp_int =nCols-2;
		printLIN("\t**->Number of covariates used for further anyliss: "+change_int_into_string(tmp_int) +"\n")		;
			//cout <<_cov_names_frm_argv.size()<<endl;
		
	 }else
	 {
		pr_msg = "Problem in parsing the first line of the file "+fileName +".\n";
		error(pr_msg);
	 }
	 ++countIndiv;
	 vector<CBPED*>::const_iterator it_pCBPED=pedInfo.begin();
	 while(getline(file,line,'\n'))
	 {
		if(line=="") // this blank line 
			continue;
		if(line[0]=='#') // this is comment line 
			continue;
		stringstream line_parser1(line);

		if(line_parser1.good())
		{
			
			if((countIndiv-1)<=pTotal)
			{
				while(line_parser1>>buffer)
					(*it_pCBPED)->covariate.push_back(buffer);
			}		

			// write eror if individual id and fid does not match with the 
			// indid  and fid of ped file 
			if((*it_pCBPED)->famId=="")
			{
				(*it_pCBPED)->famId=(*it_pCBPED)->covariate[0];
			}
			else if((*it_pCBPED)->covariate[0]!=(*it_pCBPED)->famId)
			{
				
				cout << "Problem reading the follwoing "+change_int_into_string((1+(int)(it_pCBPED-pedInfo.begin())) )+"th row of covariate file: \n";
				for(unsigned int i=0;i<(*it_pCBPED)->covariate.size(); i++)
					cout << (*it_pCBPED)->covariate[i] << " " ;
				cout << "\n original family id: "<< (*it_pCBPED)->famId <<"\n"; 
				pr_msg="The family ID given in covariate file does not match with the family ID   in the ped file.  Please correct it first.\n";
				error(pr_msg);
			}	
			//individual id check 
			if((*it_pCBPED)->indId=="")
			{
				(*it_pCBPED)->indId=(*it_pCBPED)->covariate[1];
			}
			else if((*it_pCBPED)->covariate[1]!=(*it_pCBPED)->indId)
			{
				cout << "Problem reading the follwoing "+change_int_into_string((1+(int)(it_pCBPED-pedInfo.begin())) )+"th row of covariate file: \n";
				for(unsigned int i=0;i<(*it_pCBPED)->covariate.size(); i++)
					cout << (*it_pCBPED)->covariate[i] << " " ;
				cout<< "\n original individual id: "<< (*it_pCBPED)->indId << "\n";
				//cout << "no of pedid: "<< (int)(it_pCBPED-pedInfo.begin())<<endl; //debug 
				pr_msg="The Individual ID given in covariate file does not match with the Inidvidual ID   in the ped file.  Please correct it first.\n";
				error(pr_msg);
			}			
			if((*it_pCBPED)->covariate.size()!=nCols)
			{
				cout << "Problem reading the follwoing row of covariate file: \n";
				cout << (*it_pCBPED)->covariate.size() <<endl;
				for(unsigned int i=0;i<(*it_pCBPED)->covariate.size(); i++)
					cout << (*it_pCBPED)->covariate[i] << " " ;
				cout <<endl;
				pr_msg="The number of columns in header of the file and the number of columns  in the line containing "+change_int_into_string(countIndiv-1)+"th Individual are not same. Please correct it first.\n";
				error(pr_msg);
			}
			
			
		}
		else
		{
			printLIN( "Line parser not good while reading the line containing "+\
			change_int_into_string(countIndiv-1)+"the individual.\n");
			break;
		}
		if(it_pCBPED==(pedInfo.end()-1)) break;
		 //vector<CBSNP*>::const_iterator it_pCBSNP=genInfo.begin();
		if(countIndiv<=pTotal)
			++it_pCBPED;
		++countIndiv;
	}
	
	 
	 file.clear();file.close();
	 if(it_pCBPED!=(pedInfo.end()-1))
	 {
		pr_msg="The covariate file "+fileName+" contains only  "+change_int_into_string(1+it_pCBPED-pedInfo.begin() )+\
		" Individuals in covariate files. However in  the ped file there are "+change_int_into_string(pTotal)+" individuals.\n";
		error(pr_msg);
	 
	 }
	 
 }

//writing snptest sample file
 //individual callrate
 void CBPED::calculate_indiv_callrate(const vector<CBPED*>&pedInfo, const vector<CBSNP*> &genInfo)
 // This function calculate call rate by ignoring excluded snps and individuals or quality of snps and individuals
 {
	 CBPED* pCBPED =pedInfo[0];
	if(!pCBPED->given_indiv_crate)
	 {
	   CBSNP* pCBSNP 	=genInfo[0];
	   unsigned int _nmiss	=0;
	   bool tmp_convert_into_0123 =true;
	   const unsigned int _tmp_good_SNPS =(int) CBSNP::no_of_good_snps(genInfo); //genInfo.size();
	   const unsigned int _INDS =pedInfo.size();

	   for(unsigned int i=0; i<_INDS;++i)
	   {
		   pCBPED =pedInfo[i];
		  // if(pCBPED->quality)
		   {
			 if(tmp_convert_into_0123)
			 {
				 for(unsigned int j=0;j<genInfo.size();++j)
				 {
					 pCBSNP =genInfo[j];
					 CBSNP::convert_snp_genotype_into_0123(pCBSNP);
				 }
				 tmp_convert_into_0123=false;
			 }
			 else{
				 	 for(unsigned int j=0;j<genInfo.size();++j)
					 {
						 pCBSNP =genInfo[j];
						 if(pCBSNP->geno_0123[i]==3) // && pCBSNP->quality)
							 ++_nmiss;
					 }
			 }
				 //cout << _nmiss << " "<< _INDS<< "\n"; //debug
			 pCBPED->ind_CR =(1.0 -(float)(1.0*_nmiss/float(_tmp_good_SNPS)));
			 _nmiss=0;
		  }
	   }
	  //assign true
	   pCBPED->given_indiv_crate=true;
	 }

 }// end of function
//calcualte hwe
 void CBPED::calculate_snp_hwe(const vector<CBSNP*> &genVec,const vector<CBPED*>& vPed)
 {
 	unsigned int j		  	 =0;
 	CBSNP* _pCBSNP			 =genVec[0];
 	unsigned int _tmp_zeros  = 0;
 	unsigned int _tmp_ones   = 0;
 	unsigned int _tmp_twos   = 0;
 	unsigned int i			 =0;
	 for(j=0;j<genVec.size();++j)
	 {

		  _pCBSNP =genVec[j];
		  if(_pCBSNP->quality)
		  {
			  CBSNP::convert_snp_genotype_into_0123(_pCBSNP);
			 	  //Bpar::calculate_hardy =false;
				  //const int ZEROS			 =count( _pCBSNP->geno_0123.begin(),_pCBSNP->geno_0123.end(),0);
				  //const int ONES			 =count( _pCBSNP->geno_0123.begin(),_pCBSNP->geno_0123.end(),1);
				  //const int TWOS 			=count( _pCBSNP->geno_0123.begin(),_pCBSNP->geno_0123.end(),2);
				  //
			  	 // cout << "_pCBSNP->geno0123.size(): "<< _pCBSNP->geno_0123.size()<<endl;
				  for(i=0;i<vPed.size();++i)
				  {
					  if(vPed[i]->quality)
					  {
						  const int _tmp_val =_pCBSNP->geno_0123[i];
						  if(_tmp_val==0)
							  ++_tmp_zeros;
						  else if(_tmp_val==1)
							  ++_tmp_ones;
						  else if(_tmp_val==2)
							  ++_tmp_twos;
						//  cout << _tmp_val << ", ";
					  }

				  }
				  //cout <<endl;
				 // cout << _tmp_zeros << ", "<< _tmp_ones << ", "<<_tmp_twos <<endl;
				  _pCBSNP->pvalue_hwe_exact= CBSNP::SNPHWE (_tmp_ones,_tmp_zeros,_tmp_twos); //(ONES, ZEROS, TWOS);
			      _pCBSNP->given_snp_hwe=true;

		  }
		   _tmp_zeros	=0;
		   _tmp_ones 	=0;
		   _tmp_twos	=0;
		   //cout<<"testing\n";
 	}
 	//Bpar::calculate_hardy=false;
  }
 //
 void CBPED::calculate_maf( const vector<CBSNP*>&genVec, const vector<CBPED*>&vPed)
 {
	 //cout<< count( _pCBSNP->geno_0123.begin(),_pCBSNP->geno_0123.end(),0)<<endl;
	 //determinie total alleles etc
	  unsigned int gpTotal	 =0; // good pedigrees
	  unsigned int gsTotal 	=0; // good snps
	  unsigned int i,j		=0;
	  unsigned int _tmp_ones =0;
	  unsigned int _tmp_twos =0;
	  unsigned int _tmp_threes =0;
	 // cout << vPed.size() <<"," <<genVec.size() <<endl; //debug
	  const CBPED* pCBPED=vPed[0];
	  CBSNP* pCBSNP =genVec[0];
	 // cout << pCBSNP->geno_0123.size()<<endl;
	  for(j=0;j<genVec.size();j++)
	  {
		  pCBSNP =genVec[j];
		//  cout <<"test: "<< j <<endl;
		  CBSNP::convert_snp_genotype_into_0123(pCBSNP);
		 // cout << "pCBSNP->quality: "<<pCBSNP->quality<<endl;
		  if(pCBSNP->quality)
		  {
			  for(i=0;i<vPed.size();i++)
			  {
				  pCBPED=vPed[i];
				 // cout <<"\npCBPED->quality: "<<pCBPED->quality<<endl;
				  if(pCBPED->quality)
				  {
					  if(pCBSNP->geno_0123[i]!=3)
						  ++gpTotal;
					  if(pCBSNP->geno_0123[i]==1)
						  ++_tmp_ones;
					  else if(pCBSNP->geno_0123[i]==2)
						  ++_tmp_twos;
					  else if (pCBSNP->geno_0123[i]==3)
						  ++_tmp_threes ;
					//  cout << pCBSNP->geno_0123[i]<<", ";

				  }
			  }
			 //cout <<endl;
			 //here comes rest for the SNP
			  const int TOTAL_ALLELES		=2*gpTotal;
			  const int ALLELECOUNT 	= _tmp_ones +_tmp_twos*2;
			  //cout << "ALLELECOUNT: "<<ALLELECOUNT<<endl;
			  const float ALLELE_FREQ	=(float) ALLELECOUNT/TOTAL_ALLELES;
			 // cout << "2nd ALLELE_FREQ: "<<ALLELE_FREQ<< endl;
			  pCBSNP->freq=(1.0-ALLELE_FREQ); // freq is frequency of first allele
			  pCBSNP->nchrobs =(TOTAL_ALLELES);
			  if(ALLELE_FREQ<0.5)
			   	{

			   		pCBSNP->maf=ALLELE_FREQ;
			   		pCBSNP->maj_allele=pCBSNP->allele1;
			   		//cout<<"1_:"<<pCBSNP->min_allele<< " "; // minor allele
			   		//if(pCBSNP->allele2!="")
			   			pCBSNP->min_allele=pCBSNP->allele2;
			   	//	else
			   		//	pCBSNP->min_allele="0";
			   		//cout <<"2_:"<<pCBSNP->min_allele<< " "; // minor allele
			   		//cout <<boolalpha <<!(pCBSNP->min_allele=="")<< " "; // minor allele
			   	}else{
			   		pCBSNP->maf=(1.0-ALLELE_FREQ);
					pCBSNP->min_allele=pCBSNP->allele1;
					pCBSNP->maj_allele=pCBSNP->allele2;
					//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),2,99); // value 99 is used just as a temp coding to swap later .
					//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),0,2);
					//replace(pCBSNP->geno_0123.begin(),pCBSNP->geno_0123.end(),99,0);
		   	}
			  //assing all tmp values to zero
			  _tmp_threes =0;
			  _tmp_ones	 =0;
			  _tmp_twos	 =0;
			  gpTotal	 =0;

		  }//end of snp quality
			//cout <<"1_: "<<pCBSNP->maj_allele<< ", "; // minor allele
			//cout <<"2_: "<<pCBSNP->min_allele<< "  \n"; // minor allele
	  } // end of for genvec loop


 }
 //-------------------------------------------//

 void CBPED::schreibe_ind_callrate(const vector<CBSNP*>&genInfo,const vector<CBPED*>&pedInfo, const string & outFileName)
 {
	 	 const CBPED* pCBPED =pedInfo[0];
 		//CBPED::calculate_indiv_callrate(pedInfo, genInfo);
 		string crFileName =outFileName+"_indiv.crate";
 		ofstream ofs;
 		ofs.clear(); ofs.close();
 		ofs.open(crFileName.c_str(),ios::out);
 		if(!ofs)
 			error("out file "+outFileName+" does not exit or can not be written out. Please check the filename and file path.\n");

 		ofs <<"FID"<<" "<< "IID"<<" "<<"CRATE" ;
 		ofs<<"\n";//end of line
 	for(unsigned int i=0; i<pedInfo.size();++i)
 	{
 		pCBPED =pedInfo[i];
 		ofs << pCBPED->famId<<" "<< pCBPED->indId <<" "<< pCBPED->ind_CR;
 		ofs << "\n";//end of line

 	}
 	ofs.clear(); ofs.close();
	printLIN("\t**->File containing individual callrate has been written out and saved as \""+\
			crFileName+"\". \n");
return; }//end of function

void CBPED::write_snptest_sample_file(const vector<CBPED*>&pedInfo, const string& outFileName)
{

	ofstream ofs;
	ofs.clear();ofs.close();
	string sampleFileName=outFileName+".sample";
	ofs.open(sampleFileName.c_str());
	if(!ofs)
		error("out file "+outFileName+" does not exit or can not be written out.\n");
		
	 ofs.clear();
	  vector<int>ncols;
	  
	 CBPED* pCBPED =pedInfo[0];
	 if( _given_cov_types&& _given_cov_names)
	 {
		for(unsigned int i=0;i<_cov_names_frm_argv.size();++i)
		{
			vector<string>::iterator it =find(covariate_names.begin(),covariate_names.end(),_cov_names_frm_argv[i]);
			int mycol 	=(int)(it-covariate_names.begin());
			//cout << mycol << " ";
			ncols.push_back(mycol);
			//cout<< covariate_names[3]<<" ";
			//cout << position_in_vector(covariate_names,"sst")<< " ";
			//ncols.push_back(position_in_vector(covariate_names,_cov_names_frm_argv[i]));
			//cout<< ncols[i]<<" ";
		}
			
	 }
	ofs << "ID_1"<<" "<< "ID_2"<<" "<<"missing"<<" ";
	for(unsigned int i=0;i<ncols.size();++i)
	ofs << covariate_names[ncols[i]]<<" ";
	ofs<<"sex" <<" "<<"status"<<endl; 	
	//write covariate types
	ofs << "0"<<" "<< "0" <<" "<<" "<<"0"<<" ";
	for(unsigned int i=0;i<_cov_types_frm_argv.size();++i)
		ofs <<_cov_types_frm_argv[i]<< " ";
	ofs<<"D" <<" "<<"B"<<endl; 	// this is for sex and disease status
	//writing info form pedigree file
	
	for(unsigned int i=0; i<pedInfo.size();++i)
	{
			pCBPED =pedInfo[i];
		if(pCBPED->quality)
		{
			if(pCBPED->famId!="") ofs<<pCBPED->famId<< " ";
			 else
			{
				string newid="fam"+change_int_into_string(i+1);
			  	ofs << newid << " ";
			}
			if(pCBPED->indId!="") ofs<<pCBPED->indId<< " ";
			else
			{
				string newid="indv"+change_int_into_string(i+1);
			  	ofs << newid << " ";
			}

			//cout <<"indiviudal call rate:" <<pCBPED->ind_CR<<endl; //debug
			ofs<<(1.0-(pCBPED->ind_CR))<< " "; // call rate 
			if(_given_cov_types&& _given_cov_names)
			{
				for(unsigned int i=0;i<ncols.size();++i)
					ofs <<pCBPED->covariate[ncols[i]]<<" ";
	
			}
			//sex 
			if(pCBPED->sex) ofs<<1<< " "; // 1 for male
			else if(!pCBPED->sex && !pedInfo[i]->miss_sex) ofs<<2<<" "; // for female
			 else if(pCBPED->miss_sex) ofs<<"NA"<< " "; // -9 is missing
			//ofs<<pCBPED->sex<<" "<< pCBPED->pheno<<endl; // sex and disease status 
			//status 
			if(pCBPED->pheno)
	  		 	 ofs<<"1"<< " \n";
	  		 else if(!pCBPED->pheno && !pCBPED->miss_pheno)
	  				ofs<<"0"<< " \n";
	  		else if(!pCBPED->pheno && pCBPED->miss_pheno)
	  			ofs<<"NA"<< " \n";
		}// end of if quality loop
	}//end of for loop
	pCBPED =pedInfo[0];
	ofs.clear(); ofs.close(); 
	printLIN( "\t**-> IMPUTE2/SNPTEST *.SAMPLE file has been written out and saved as \""+sampleFileName+"\". \n");
	 
}
// read mach hapmap file
 void CBPED::cerr_with_delete(vector<CBPED*>& vPed, vector<CBSNP*>&genInfo, const string& cerr_message){
	 for(unsigned int i=0;i<vPed.size();++i)
		delete vPed[i];
	vPed.clear();
	 for(unsigned int i=0;i<genInfo.size();++i)
		delete genInfo[i];
	genInfo.clear();
	cerr<<" ERROR: " <<cerr_message;
	LIN<<" Error:"<< cerr_message;
	LIN.close();
	exit(1);	
		
 }
void CMPED::lese_mach_hapmap_file(vector<CBPED*>&vPed,vector<CBSNP*>& genInfo,vector<string>& cZeile,const string& fileName)
{
	//cout << "test"<<endl;//debug 
	printLIN("*->Reading file \""+fileName+"\": \n");
	string myString					="";
	string pheno					="";
	string temp_value1				="";
	string temp_value2				="";
	string mygeno					="";
	bool 	ende 					=false;
	unsigned long int ccol			=0;
	unsigned long int ncol			=0;
	unsigned int 	 nrow 			=0;
	unsigned int nLines				=0;
	unsigned long int max_ncol		=0;
	unsigned long int snpAnzahl		=0;
	bool def_max_ncol_once			=true;
	vector<string> zeile_2nd;
	vector<string>genotype;
	bool define_pCBSNP				=true;	
	
	try
	{
		IFS pf(fileName.c_str(),"r");
		if(!pf){
			string msg =fileName+ " either does not exits or could not be opend.\n";
			error(msg);		
		}
		while(true)
		{
			if(_ifs_eof(pf))
			{
				if(cZeile.size()!=0 && (zeile_2nd.size()==0))
				{
					string pr_msg	="Problem in reading the second line of individual \""+cZeile[0]+"\": in "+change_int_into_string(nrow)+"th row: \n" ;
					printLIN(pr_msg);	
					pr_msg= "In MACH hapmap file, there should be always  two lines representing an individual, but this is not the case here.\n";
					cerr_with_delete(vPed,genInfo,pr_msg);
				}
				else 
				break;
			
			} 
			while(!_ifs_eof(pf))
			{
				ccol	=CBSNP::stringLeser(pf,myString,ende);
				ncol	=ncol+ccol;
				if(ccol!=0)
				{
					if(!myString.empty())
					{
						//cout <<"first modulo: "<<( nrow%2==0)<< endl;//debug 
						if(cZeile.size()<3)
							cZeile.push_back(myString);
						else  
						zeile_2nd.push_back(myString);
					}
						
					else
					{
						string pr_msg= "Error in reading "+change_int_into_string(nLines+1)+"th line of file \""+fileName+"\". \n";
						cerr_with_delete(vPed,genInfo,pr_msg);
					}
					//cout << "test1" << endl; // debug 
					if(ende || _ifs_eof(pf))
					{
						
						++nLines;
						if(def_max_ncol_once && (ccol!=0))
						{
							max_ncol=ncol;
							def_max_ncol_once =false;
							//cout << "max_ncol: "<< max_ncol << endl;//debug 
						}
						if((ncol==max_ncol) && (max_ncol!=0))
							++nrow;
						myString="";
						//cout <<"nrow  : "<< nrow << " at  the end of line and later modulo: "<<( nrow%2==0)<< endl;//debug 		
						break;	
						
					}
					
				
				}
				if(ccol==0)continue;
			}//end of inner while loop while loop ;
			//	cout <<"first size: "<< cZeile.size() << " second size: "<< zeile_2nd.size() << endl;
				//cout <<"first : "<< cZeile[1] << " second size: "<< zeile_2nd[2] << endl;
				//cout <<"nrow: "<< nrow << " Later modulo: "<<( nrow%2==0)<< endl;//debug 	
			if(cZeile.size()!=0&&zeile_2nd.size()!=0) // this condition will read two lines till czeile size ==3 and zeile_2nd size ==3
			{
				if(cZeile.size()==3 && zeile_2nd.size()==3) // if both have size 3 then.
				{
					
					//if(nrow%2==0)
					{
						//cout <<"ind: " <<(nrow)<< endl;
						temp_value1	=cZeile[0];
						temp_value2	=zeile_2nd[0];
						//cout << "first: "<< temp_value1 <<" second: " << temp_value2 <<endl; // debug
						size_t found;
						found=temp_value1.find_last_not_of("_A");
						if(found!=string::npos)
							temp_value1.erase(found+1);
						found=temp_value2.find_last_not_of("_B");
						if(found!=string::npos)
							temp_value2.erase(found+1);
						//cout << temp_value1 << " "<<temp_value2 << endl;
						/*
						 * Normally the follwoing condition is true in mach references but not all
						 * references have this situation, so we just omit the check
						if(temp_value1!=temp_value2)
						{

							string pr_msg="Problem reading individual \""+temp_value1+"\" and \""+temp_value2+"\" from the  "+change_int_into_string(nLines-1) +"th and " +change_int_into_string(nLines)+"th lines: \n";
							printLIN(pr_msg);
							pr_msg="These two lines should have the same individual name with ending \"_A\" and  \"_B\" respectively.\n" ;
							cerr_with_delete(vPed,genInfo,pr_msg);
						}
						*/
						CBPED* pCBPED =new CBPED;
						pCBPED->famId=temp_value1;
						//check next 
						temp_value1	=cZeile[1];
						temp_value2	=zeile_2nd[1];
						found=temp_value1.find_last_not_of("_A");
						if(found!=string::npos)
							temp_value1.erase(found+1);
						found=temp_value2.find_last_not_of("_B");
						if(found!=string::npos)
							temp_value2.erase(found+1);
						//cout << temp_value1 << " "<<temp_value2 << endl;
						/*
						 * same as before this is not necessary to check
						if(temp_value1!=temp_value2)
						{
							delete pCBPED;
							string pr_msg="Problem reading individual \""+temp_value1+"\" and \""+temp_value2+"\" from the  "+change_int_into_string(nLines-1) +"th and " +change_int_into_string(nLines)+"th lines: \n";
							printLIN(pr_msg);
							pr_msg="These two lines should have the same individual name with ending \"_A\" and  \"_B\" respectively.\n" ;
							cerr_with_delete(vPed,genInfo,pr_msg);
						}
						*/
						//cZeile.erase(cZeile.begin(),cZeile.begin()+2);
						//zeile_2nd.erase(zeile_2nd.begin(),zeile_2nd.begin()+2);
						//cout << cZeile.size() << " " <<zeile_2nd.size() <<endl;
						pCBPED->indId=temp_value1;
						lese_phenotype_info(pCBPED, pheno);
						vPed.push_back(pCBPED);
						temp_value1=cZeile[2];
						temp_value2=zeile_2nd[2];
						if(temp_value1.length()!=temp_value2.length())
						{
							cerr<< "size of first haplotype: "<< temp_value1.length() << ", \n";
							cerr<< "size of second haplotype: "<< temp_value2.length() << ". \n";	
							string pr_msg="Problem in reading individual \""+cZeile[0]+"\" from the  "+change_int_into_string(nLines-1) +"th and " +change_int_into_string(nLines)+"th lines: \n";
							printLIN(pr_msg);
							pr_msg="These two lines should have  same number of characters in the third columns.\n" ;
							cerr_with_delete(vPed,genInfo,pr_msg);
						
						}
						else
						{	 if(define_pCBSNP)
							{
								define_pCBSNP= false;
								snpAnzahl=genInfo.size();	
							}
							for(unsigned int i=0;i<temp_value1.length();++i)
							{
								mygeno=temp_value1[i];
								mygeno+="/";
								mygeno+=temp_value2[i];
								//cout <<temp_value1[i]<< ", "<< temp_value2[i]<< ": "<< mygeno << " \n" ;
								genotype.push_back(mygeno);
								
							}
							//cout <<"genotype.size()"<< genotype.size()<<endl;
							if(snpAnzahl!=genotype.size())
							{
								cerr<< "Expected number of SNPs (according as given for first individual): "<< snpAnzahl << ", \n";
								cerr<< "Number of SNPs in current  individual: "<< temp_value1.length() << ". \n";	
								string pr_msg="Problem in reading individual \""+cZeile[0]+"\" from the  "+change_int_into_string(nLines-1) +"th and " +change_int_into_string(nLines)+"th lines: \n";
								printLIN(pr_msg);
								pr_msg="The number of SNPs in these two  lines should be equal to the number of SNPs given for first individual.\n" ;
								cerr_with_delete(vPed,genInfo,pr_msg);
								
							}	
							CBPED::lese_pedfile_genInfo(genInfo, genotype,(int)round(nLines/2+1));
							//cout << "testing"<< " ";
						}
						
					}
				
						//cout <<"vPed size: "<< vPed.size() << endl;

				}
				else
				{
					string prg_msg	=" Problem in reading "+change_int_into_string(nrow)+"th line of file \""+fileName+"\": \n";
					printLIN(prg_msg);	
					prg_msg="Each line of mach HapMap file should contain 3 columns but this is not the case here.\n";
					cerr_with_delete(vPed,genInfo,prg_msg);
				
				}
			}
			else 
			continue;
			ncol=0;
			cZeile.clear();
			cZeile.resize(0);
			zeile_2nd.clear();
			zeile_2nd.resize(0);
			genotype.clear();
			genotype.resize(0);
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



return;}
//
void CBPED::lese_remove_indivs_file(const string& fileName, const vector<CBPED*> &  pedInfo)
{
	printLIN("\t**->Removing individuals mentioned in file \""+fileName+"\": \n");
		string line ;
		string buffer;
		string _tmp_pid="";
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
		vector<string> tmp_vec_pid;
		CBPED* pCBPED ;
		vector<string>::iterator it;
		int ped_index =0;
		bool _found =false;
		if(pedInfo.size()>0)
		{
			for(unsigned int i=0;i<pedInfo.size();++i)
			{
				pCBPED 		=pedInfo[i];
				_tmp_pid 	=pCBPED->indId;
				//cout << _tmp_pid <<" ";
				tmp_vec_pid.push_back(_tmp_pid);
			}
		}
		while(getline(file,line))
		{
			if(line.size()>100)
			{
				file.clear();file.close();
				error("Are you sure about the  file containing the list of individual-ids  mentioned with \"--remove\" command?\n");
			}
			if(line.empty())
				continue;
			if(line[0]=='#')
				continue;
			stringstream line_parser(line);
			if(line_parser.good())
				line_parser>>buffer;
			//cout <<"\nbuffer: "<< buffer <<endl;//debug
			it =find(tmp_vec_pid.begin(),tmp_vec_pid.end(),buffer);
			if( it!=tmp_vec_pid.end())
			{
				_found =true;
				ped_index =(int)(it-tmp_vec_pid.begin());
				//alternative
				pCBPED 		= pedInfo[ped_index];
				pCBPED->quality =false;
				++_tmp_nLine;
			}
			else if(buffer!="")
			{
					++_tmp_nFound;
					if(_tmp_nFound<20)
					 printLIN("\t ... Individual named \""+ buffer+"\" not found in the data \n");
			}
			//cout << boolalpha <<_found<<endl;

			_found=false;
			line="";
		}//end of while getline loop;
		 if(_tmp_nFound>20)
			 printLIN("\tThere are more than 20 Individuals not found in the data");
		printLIN("\t**->Total number of excluded individuals: "+change_int_into_string(_tmp_nLine)+". ");
		//closing file
		file.clear(); file.close();
	return ;
}
///xxxxxxxxxxxxxx
void CBPED::write_plink_dosagesFiles(const vector<CBPED*>&pedInfo,const vector<CBSNP*>genInfo, const string& outFileName){
	//cout << "Hi from write_plink_dosagesFiles(const vector<CBPED*>&pedInfo,const vector<CBSNP*>genInfo, const string& outFileName) \n"; //debug 
	printLIN("*->Writing files in " +CBSNP::_out_format+" format: \n");
	//writing plink map file 
	CBSNP::write_plinkMapFile(genInfo,outFileName);
	//writing plink fam file  and plink dosage file 
	string doseFileName =outFileName+".dose";
	string famFileName =outFileName+".fam";
	bool sex_missing=true; 
	bool pheno_missing=true;
	//writing plink fam file 	
	ofstream ofs; 
	ofs.clear(); ofs.close();
	ofs.open(famFileName.c_str());
	if(!ofs)
		error("out file "+doseFileName+" can not be written out. \n");
	ofs.clear();
	for(unsigned int i=0; i<pedInfo.size(); ++i)
	{
		if(pedInfo[i]->quality)
		{
			if(pedInfo[i]->famId!="")
				ofs<<pedInfo[i]->famId<< " ";
			else if (pedInfo[i]->indId!="")
				ofs<<pedInfo[i]->indId << " ";
				else
					ofs << "famid"<<i+1<< " ";
			if(pedInfo[i]->indId!="")ofs<<pedInfo[i]->indId << " ";
				else ofs<< "indv"<<i+1<<" ";
			if(pedInfo[i]->patId!="")ofs<<pedInfo[i]->patId << " ";
				else ofs<< 0<<" ";
			if(pedInfo[i]->matId!="") ofs <<pedInfo[i]->matId <<" ";
				else ofs<< 0 << " ";
			if(pedInfo[i]->sex) ofs<<"1"<< " "; // 1 is male
				else if(!pedInfo[i]->sex && !pedInfo[i]->miss_sex) ofs<<"2"<<" "; // 2 is female
				else if(pedInfo[i]->miss_sex)
				{
					ofs<<"0"<< " "; // 0 is missing
					if(sex_missing)
					{
						string sout ="\n **->WARNING: Your genotyped data contains individuals with missing sex information.\n";
						printLIN(sout);
						sout ="\t Plink may not  analyse these individiuals. Please make sure that the individuals have correct sex information.\n";
						printLIN(sout);
					}
					sex_missing=false ;

				}
				
			if(pedInfo[i]->pheno)
					 ofs<<"2"<< " ";
				else if(!pedInfo[i]->pheno && !pedInfo[i]->miss_pheno)
						ofs<<"1"<< " ";
				else if(!pedInfo[i]->pheno && pedInfo[i]->miss_pheno)
				{
					ofs<<"-9"<< " "; //-9 here missing
					if(pheno_missing)
					{
						string sout ="\n **->WARNING: Your genotyped data contains individuals with missing phenotype information.\n";
						printLIN(sout);
						sout ="\t Plink may not  analyse these individiuals. Please make sure that the individuals have correct phenotype information.\n";
						printLIN(sout);
					}
					pheno_missing=false ;

				}
			ofs<< "\n";
		}
	}
	ofs.clear(); ofs.close();
	printLIN("\t**->PLINK fam file has been written out and saved as \""+famFileName+\
	  			"\".\n");
	//writing  plink dosage file 
	
	ofs.clear(); ofs.close();
	ofs.open(doseFileName.c_str());
	if(!ofs)
		error("out file "+doseFileName+" can not be written out. \n");
	ofs <<"SNP" << " "<< "A1" <<" "<< "A2"<< " ";
	for(unsigned int i=0; i<pedInfo.size(); ++i)
	{
		if(pedInfo[i]->quality)
		{
			if(pedInfo[i]->indId!="")
			{
				if(pedInfo[i]->famId!="")
					ofs<<pedInfo[i]->famId << " ";
				else
					ofs<<pedInfo[i]->indId << " ";
				ofs<<pedInfo[i]->indId << " "; 
			}
	  		else 
			{
				ofs<< "famid"<<i+1 << " "; //<<"_"<<1<<" "; //two times  prob 1 and prob 2 
				ofs<< "indv"<<i+1<< " "; // <<"_"<<2<<" "; //two times  
				//ofs<<"F1"<< " " <<" I1"<< " ";
			}
		}
	}
	ofs<<"\n";
	//writing genotype probs 
	CBSNP* pCBSNP=genInfo[0];
	CBPED* pCBPED =pedInfo[0];
	for(unsigned int j=0; j<genInfo.size(); j++)
	{
		pCBSNP=genInfo[j];
		if(pCBSNP->quality)
		{
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
				if(pCBSNP->rsId!="") ofs << pCBSNP->rsId<< " ";
				else if(pCBSNP->snpName!="") ofs << pCBSNP->snpName<< " ";
				if(pCBSNP->allele1!="") ofs << pCBSNP->allele1 << " ";
				else if(pCBSNP->allele2!="") ofs << pCBSNP->allele2<< " ";
				else
				{
					ofs << "A" << " ";
				}
				if(pCBSNP->allele2!="") ofs << pCBSNP->allele2;
				else if(pCBSNP->allele1!="") ofs << pCBSNP->allele1;
				 else
					ofs << "B" ;
				for(unsigned int i=0; i<pedInfo.size();++i)
				{
					pCBPED =pedInfo[i];
					if(pCBPED->quality)
					{
							//case1 if genotype probabilities are given
						if((pCBSNP->pgeno1.size()==pedInfo.size()) &&(pCBSNP->pgeno2.size()==pedInfo.size())&&(pCBSNP->pgeno3.size()==pedInfo.size()))
						{
							if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) //if missing
							{
								//ofs <<" "<<Bpar::miss_geno << " "<<Bpar::miss_geno;// << " "<<Bpar::miss_geno;
								// if there are missing , you can not write a dose file with two genotype probability because :
								//	1 	0  	=>AA
								//	0 	1	=>AB
								//	0	0 	=>BB
								//	How to denote missing??
								ofs.clear();
								ofs.close();
								string sout= "problem in writing out "+change_int_into_string(j) +"th SNP:\n";
								printLIN(sout) ;
								sout = "There are missings in original genotype file and expressing missings in plink dose format is not possible.\n";
								error(sout);
							}else
							ofs <<" "<< pCBSNP->pgeno1[i] << " "<< pCBSNP->pgeno2[i] << " "; //<< pCBSNP->pgeno3[i];
						}
						//case 2 if genotype probabilities are not given
						else
						{
							if( !pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if homozygote1  00
								ofs <<" "<< "1" <<" "<< "0"; //<<" "<< "0";
							else if( pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if homozygote2 11
								ofs <<" "<< "0" <<" "<< "0"; //<<" "<< "1";
							else if( !pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if heterozygote
								ofs <<" "<< "0" <<" "<< "1";// <<" "<< "0";
							else if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) //if missing
							{
								//ofs <<" "<<Bpar::miss_geno << " "<<Bpar::miss_geno;// << " "<<Bpar::miss_geno;
								// if there are missing , you can not write a dose file with two genotype probability because :
								//	1 	0  	=>AA
								//	0 	1	=>AB
								//	0	0 	=>BB
								//	How to denote missing??
								ofs.clear();
								ofs.close();
								string sout= "problem in writing out "+change_int_into_string(j) +"th SNP:\n";
								printLIN(sout) ;
								sout = "There are missings in original genotype file and expressing missings in plink dose format is not possible.\n";
								error(sout);
							}
							else
							{
								ofs.clear();
								ofs.close();
								error("problem in writing out "+change_int_into_string(j) +"th SNP.\n");
							}
						}
					}//end of quality if loop
				}//end of pedinfo for loop
				ofs<< "\n";

			}
	
	} //end of if quality loop
	}//end of for j loop
	ofs.clear();
	ofs.close();
	printLIN("\t**->PLINK dose file has been written out and saved as \""+doseFileName+\
	  			"\".\n");

}
//

void CBPED::write_pedlist_file(const vector<CBPED*>& pedInfo, const string& outFileName)
{
	//----------------------------------#--------------------------------------#
	// writing  ped list files.
	//----------------------------------#--------------------------------------#
	ofstream ofs;
	string extraPedFileName=outFileName+"_pedlist.txt";
	ofs.clear();	 ofs.close();
	ofs.open(extraPedFileName.c_str());
	if(!ofs)
		error("out file "+outFileName+" does not exit or can not be written out.\n");
	CBPED* pCBPED =pedInfo[0];
	for(unsigned int i=0; i<pedInfo.size();++i)
	{
		pCBPED =pedInfo[i];
		if(pCBPED->quality)
		{
			//#####
			if(pCBPED->indId!="")ofs<<pCBPED->indId;
			else ofs<< "indv"<<i+1;
			ofs<<"\n"; //end of line
		}
	 }
	ofs.clear();	 ofs.close();
	printLIN("\t**->pedigree list has been written out and saved as \""+extraPedFileName+\
			"\".\n");
	return ;
}


