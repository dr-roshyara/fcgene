#include"genable.h"
#include<iterator>
void CGPED::write_genable_raw_file(const vector<CBPED*>& vPed, const vector<CBSNP*>& genInfo, const string & outFileName,bool bcast)
{
	//this array has been made from the follow r-code copied from GenABEL
	/* alleleID.alleles <- function() {
		a <- list();
		a[[1]] <- c("1","2")
		a[[2]] <- c("A","B")
		alleles <- c("A","T","G","C","-")
		idx <- 3
		for (i in alleles) {
			for (j in alleles) {
				if (i==j) next;
				a[[idx]] <- c(i,j)
				idx <- idx + 1
			}
		}
		a[[idx]] <- c("2","1")
		idx <- idx + 1
		a[[idx]] <- c("B","A")
		idx <- idx + 1
		a[[idx]] <- c("I","D")
		idx <- idx + 1
		a[[idx]] <- c("D","I")
		idx <- idx + 1
		allalleles <- c("1","2","B","I","D","A","T","G","C","-")
		for (jj in allalleles) {
			a[[idx]] <- c(jj,jj)
			idx <- idx + 1
		}
		a
	}
	alleleID.alleles()
*/
	const unsigned int IND_SZ =vPed.size();
	const unsigned int SNP_SZ =genInfo.size();
	string allele_codes[37] ={"12", "AB" ,"AT" ,"AG", "AC", "A-", "TA",	"TG", "TC", "T-" ,"GA", "GT" ,"GC", "G-", "CA", "CT", "CG",
							"C-", "-A", "-T", "-G", "-C", "21", "BA", "ID", "DI", "11", "22", "BB", "II", "DD", "AA", "TT", "GG", "CC" ,"--"
							};
	vector<string> codeset(allele_codes,allele_codes+36);
	//const short unsigned int NCODES =codeset.size(); //this is the sizeof allele_codes
	//cout <<"codeset.size():  "<< codeset.size()<<" "<<codeset[35]<<endl;//debug
	//int verbose = bcast?  1 : 0;
	int offset[4] = {6, 4, 2, 0};
	int nbytes = (int)ceil((double)IND_SZ/4.);
	//int ind;
	vector<unsigned short int> intcoding; //
	vector<unsigned char * > gtype;
	unsigned char* tmp_gtype;
	//define vectors to save SNP info
	vector<string> chrom(SNP_SZ); 					//no of chromosome
	vector<string> snpnm(SNP_SZ); 					// snpnames
	vector<double> genmap(SNP_SZ);
	vector<unsigned long > phymap(SNP_SZ); 			//physical positions
	vector<unsigned short int > strand(SNP_SZ); 	// strand

	for(unsigned int j=0; j<SNP_SZ; ++j)
	{
		 CBSNP*  const pCBSNP =genInfo[j];
		//if(pCBSNP->geno1.size()!=IND_SZ || pCBSNP->geno2.size()!=IND_SZ)
		if(pCBSNP->geno1.size()!=vPed.size() || pCBSNP->geno2.size()!=vPed.size() )
		 {
			 cout <<"geno1.size(): "<<pCBSNP->geno1.size() <<"\n geno2.size():  "<<pCBSNP->geno1.size();
			 //for(unsigned int i=0;i<pCBSNP->geno1.size();++i)
			 //	outfile<< pCBSNP->geno1[i] << " " <<pCBSNP->geno2[i] << " ";
			// outfile.close();
			 error("Problem in writing out "+change_int_into_string(j+1) +"th SNP. geno1 and geno2 has no equal size\n");
		 }
		//write snps name
		if(pCBSNP->rsId!="")
			snpnm[j]=pCBSNP->rsId;
		else
			snpnm[j]=pCBSNP->snpName;
		 chrom[j]=pCBSNP->nchr; // chromosome
		//strand == 0 -> unknown ('u')
		//# strand == 1 -> plus ('+')
		//# strand == 2 -> minus ('-')
		//# strand == 3 -> four columns (name,chr,pos,strand) expected at the beginning
			 strand[j]=pCBSNP->coding_strand;
		 if(pCBSNP->bp!=-1)
			 phymap[j]=pCBSNP->bp;
		 else
			 phymap[j]=j;

		convert_snpCoding_into_genable_coding(pCBSNP,nbytes,codeset,intcoding,offset,gtype);
		//for(unsigned int i=0;i<gtype.size();++i)
		//	cout<< gtype[i]<< " ";
		 // cout <<endl;
		//_add_sz =gtype1.size()+gtype.size();
		//gtype.resize(_add_sz);
		//std::copy( gtype1.begin(), gtype1.end(), gtype.end() );
		//cout <<"gtype.size(): "<< gtype.size()<<" \n";

	}

	// starting here , genable format will be written.
	ofstream outfile ;
	outfile.clear(); outfile.close();
	const string raw_fileName =outFileName+".raw";
	outfile.open(raw_fileName.c_str(),ios::out);
	if(!outfile)
		error("out file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");
	const ios_base::fmtflags hex = ios_base::hex;
	outfile << "#GenABEL raw data version 0.1"; //first line
	outfile << endl;
	for(unsigned int i=0; i<IND_SZ;++i)
		outfile<< vPed[i]->indId << " "; // second line
	outfile << endl;
	//cout << "IND_SZ: "<< IND_SZ<<endl;
	//cout <<"test"<<endl;
	//writing SNP Name
	outfile<<snpnm[0];
	for(unsigned j=1;j<SNP_SZ;++j)
		outfile<<" "<<snpnm[j]; // third line
		outfile<<"\n";
	//writing chromosome name
	outfile<<chrom[0];
	for(unsigned j=1;j<SNP_SZ;++j)
		outfile<<" "<<chrom[j]; // fourth line
	outfile << endl;
	//writign physical maps
	outfile<<phymap[0];
	for(unsigned j=1;j<SNP_SZ;++j)
		outfile<<" "<<phymap[j]; //fifth line
	outfile << endl;
	//writing gentoypes
	outfile.flags(hex);
	for (unsigned long int i=0;i<chrom.size();i++)
	{
		outfile.width(2);
		outfile.fill('0');
		outfile << (unsigned int)intcoding[i] << " "; //sixth line
	}
	outfile << endl;

	for (unsigned long int i=0;i<chrom.size();i++)
	{
		outfile.width(2);
		outfile.fill('0');
		outfile << (unsigned int)strand[i] << " ";
	}
	outfile << endl;

	 //
	 int byte=0;
	 for(unsigned long int i = 0;i<gtype.size();i++)
	{
		tmp_gtype = gtype[i];
		for (byte = 0; byte < nbytes; ++byte) {
			outfile.width(2);
			outfile.fill('0');
			outfile << (unsigned int)tmp_gtype[byte];
			outfile << " ";
		}
		outfile << endl;

		delete [] tmp_gtype;

		//      gtype.pop_front();
	} //while (!gtype.empty());

	outfile.clear();outfile.close();
	printLIN("\t->**GenABEL-formatted raw file has been saved as \""+raw_fileName+"\".\n");
}
void  CGPED::convert_snpCoding_into_genable_coding(CBSNP* const pCBSNP,
						const int nbytes,
						const vector<string>&codeset,
						vector<unsigned short int>&intcoding,
						int*offset,
						vector<unsigned char *>&gtype
							)
{

	if(pCBSNP->change_0123)
		CBSNP::convert_snp_genotype_into_0123(pCBSNP); //convert genotype into 0123
	const int _IND_SZ =pCBSNP->geno1.size();
	int idx=0;
	//vector<unsigned char * > gtype;
	unsigned char* tmp_gtype;
	int byte=0;
	int ccd = -1;
	unsigned int _temp_val1=0;
	unsigned int  _temp_val2=0;
	const unsigned int ncodes =codeset.size();
	//
	string _allele1 =pCBSNP->allele1;  //->maj_allele; // allele 1 is the first allele and allele 2 is the second allele.
	string _allele2 =pCBSNP->allele2;  //->min_allele;
	string tmp_coding="";
	if (_allele1=="" && _allele2=="") tmp_coding="12"; // all genotypes missing
	else if (_allele1!="" && _allele2=="") tmp_coding=_allele1+_allele1;  //sprintf(tmp_chcoding,"%c%c",_allele1,_allele1); // only one allele present
	else if (_allele1=="" && _allele2!="") tmp_coding=_allele2+_allele2;  //sprintf(tmp_chcoding,"%c%c",_allele2,_allele2); // only one allele present
	else if (_allele1!="" && _allele2!="") tmp_coding=_allele1+_allele2;   //sprintf(tmp_chcoding,"%c%c",_allele1,_allele2); // only one allele present

	for (unsigned int i = 0; i < ncodes; i++)
	{
		if (codeset[i].compare(tmp_coding)==0)
		{
			ccd = i + 1;
			intcoding.push_back(ccd);
		}
	}
	//cout << " intcoding.size()"<<intcoding.size() <<endl;
	//intcoding contains indices of types of coding of vec  allele_codes
	if (ccd<0) error("coding "+ tmp_coding+ " for SNP "+pCBSNP->rsId +"not recognized !\n");

	try
	{
		tmp_gtype = new unsigned char [nbytes];
	}
	catch (bad_alloc)
	{
		error("ran out of memory reading SNP "+pCBSNP->rsId+"!\n");
	}
	if (tmp_gtype == NULL)
	{
		error("ran out of memory reading SNP "+pCBSNP->rsId+"!\n");
	}
	///xxxxxxxx
	for (byte = 0; byte < nbytes; ++byte)
	{
		tmp_gtype[byte] = 0;
		for (unsigned int ind = 0; ind < 4; ++ind)
		{
			//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

			if( !pCBSNP->geno1[idx] &&  !pCBSNP->geno2[idx]  )// if homozygote1  00
				{
					//	ofs <<" "<< pCBSNP->allele1 << " "<<pCBSNP->allele1;
						_temp_val1 =1;
						_temp_val2 =1;
				}
				else if( pCBSNP->geno1[idx] &&  pCBSNP->geno2[idx]  ) // if homozygote2 11
				{
					//ofs <<" " << pCBSNP->allele2 << " "<<pCBSNP->allele2;
					_temp_val1 =3;
					_temp_val2 =3;
				}
				else if( !pCBSNP->geno1[idx] &&  pCBSNP->geno2[idx]  ) // if heterozygote
				{
					//ofs << " "<<pCBSNP->allele1 << " "<<pCBSNP->allele2;
					if(pCBSNP->aOrder[idx])
					{
						//ofs <<" "<< pCBSNP->allele1 << " " <<pCBSNP->allele2;
						_temp_val1 =1;
						_temp_val2 =3;
					}
					else
					{
						//ofs<<" " <<pCBSNP->allele2 << " " <<pCBSNP->allele1;
						_temp_val1 =3;
						_temp_val2 =1;
					}
				}
				else if( pCBSNP->geno1[idx] &&  !pCBSNP->geno2[idx]  ) // if missing
					{
						//ofs <<" "<<Bpar::miss_geno << " "<<Bpar::miss_geno;
						_temp_val1 =0;
						_temp_val2 =0;
					}
				else
				{
					error("Problem in writing out "+pCBSNP->rsId+" SNP.\n");
				}

			//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
			//cout << _temp_val1+_temp_val2<< " "; /debux
			switch ( _temp_val1+_temp_val2)
			{
				case 2:
				if(pCBSNP->allele1==pCBSNP->maj_allele)
					//if (ca1 > ca2)
					tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)1 << offset[ind]);
				else
					tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)3 << offset[ind]);
				break;
			case 4:
				tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)2 << offset[ind]);
				break;
			case 6:
			//	if (ca1 > ca2)
			if(pCBSNP->allele1==pCBSNP->maj_allele)
					tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)3 << offset[ind]);
				else
					tmp_gtype[byte] = tmp_gtype[byte] | ((unsigned char)1 << offset[ind]);
				break;
			case 0:
				tmp_gtype[byte] = tmp_gtype[byte] | (0 << offset[ind]); // this does nothing
				break;
			default:
				string  _show_error="Problem with converting format into GenABEL: \n Out of two parts of a genotype of SNP "+pCBSNP->rsId+", half part is missing and only half is genotyped.\nSuch type of genotypes are not accepted.\n";
				error(_show_error);
			}

			++idx  ;

			if (idx >=_IND_SZ ) break;
		//	gnum[idx]+gnum[idx+1]

		}
	}
	//cout << "tmp_gtype: "<<tmp_gtype << endl;
	gtype.push_back(tmp_gtype);
	//for(int i =0; i<gtype.size();++i)
	//	cout <<gtype[i]<< " ";
	//return gtype;
}
//
void CGPED::write_genable_pheno_file(const vector<CBPED*>& vPed, const string & outFileName){

	const CBPED*  pCBPED	 =vPed[0];
	const unsigned int IND_SZ=vPed.size();
	ofstream ofs;
	ofs.clear(); ofs.close();
	const string pheno_fileName=outFileName+"_phe0.dat";
	ofs.open(pheno_fileName.c_str(),ios::out);
	if(!ofs)
		error("output file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");

	ofs<<"id"<<" "<<"sex"<<" "<<"age"<<" "<<"disease";
	ofs<<"\n";
	for(unsigned int i=0;i<IND_SZ;++i){

		pCBPED=vPed[i];
		//cout << pCBPED->famId<<" "<< pCBPED->indId<<" \n";
		/*
		//famid
		if(pCBPED->famId!="")
			ofs<<pCBPED->famId;
		else
			ofs<<"fid_"<<i+1;
		//individual id
		*/
		if(pCBPED->indId!="")
			ofs<<" \""<<pCBPED->indId<<"\"";
		else ofs<<" \""<<"iid_"<<i+1<<"\"";
		/*//patid
		if(pCBPED->patId!="")
		ofs <<" "<<pCBPED->patId;
		else ofs <<" "<<"0";
		//mapid
		if(pCBPED->matId!="")
			ofs <<" "<<pCBPED->mattId;
		else ofs <<" "<<"0";
		*/
		//if male, sex=TRUE;miss_sex=false, if female sex=false miss_sex=false,
		//	if not defined then sex=fasle, miss_sex=true;
		//write sex info
		if(pCBPED->sex&& !pCBPED->miss_sex)
			ofs<<" "<<1; //male
		else if(!pCBPED->sex&& !pCBPED->miss_sex) //female
			ofs<< " "<<0;
		else  ofs<< " "<<"NA"; //missing
		//write age
		if(pCBPED->age!=0.0)
		ofs<< " "<<pCBPED->age;
		else
			ofs<<" "<<"NA";
		//disease
		// if pheno= disease, then pheno=true, miss_state=false,
		//if no diseases, then pheno=false; miss_state=false
		// if undefined, then pheno=false, miss_state=true;
		if(pCBPED->pheno&& !pCBPED->miss_pheno)
			ofs<< " "<<1; //disease
		else if(!pCBPED->pheno&& !pCBPED->miss_pheno)
			ofs<< " "<<0; // no disease
		else
			ofs<<" "<<"NA";//missing
		ofs<<"\n"; //end of line
	}
	ofs.clear();ofs.close();
	printLIN("\t->**A phenotype file used for GenABEL has been saved as \""+pheno_fileName+"\".\n");

}

void CGPED::write_geable_commands(const string & outFileName){
	ofstream ofs;
	ofs.clear(); ofs.close();
	const string commandFileName =outFileName + "_commands.R";
	ofs.open(commandFileName.c_str(),ios::out);
	if(!ofs)
		error("output file "+outFileName+" does not exist or can not be written out.\n Please check if you have given correct file path or name.\n");
	//
	ofs<<"#Command necessary for uploading data into GenABEL: \n";
	ofs<<"#=====================================================\n";
	ofs<<"#load library\n";
	ofs<<"library(GenABEL)\n";
	ofs<< "mydata<- load.gwaa.data(phe = \""<<outFileName+"_phe0.dat\","<<" gen = \"";
	ofs<<outFileName+".raw\",  force =T)\n";
	ofs<<"----------------------------------------------------\n\n";
	//
	ofs<<"#Some examples of data Analysis in GenABEL: \n";
	ofs<<"#==========================================\n";
	ofs<<"#For phenotype data: \n";
	ofs<<"#-------------------#\n";
	ofs<<"#";
	ofs<<"summary(mydata@phdata$age)\n";
	ofs<<"summary(mydata@phdata$sex)  \n";
	ofs<<"summary(mydata@phdata$disease) \n";
	ofs<<"hist(mydata@phdata$sex) \n";
	ofs<<"hist(mydata@phdata$disease)\n";
	ofs<<"#You can use this command if you have more phenotype information.\n";
	ofs<<"#plot(mydata@phdata$qt1 , mydata@phdata$qt2)\n";
	ofs<<"#\n";
	ofs<<"#To investigates the traits in the loaded data:\n";
	ofs<<"descriptives.trait(mydata)\n";
	ofs<<"descriptives.trait(mydata, by.var = mydata@phdata$sex)\n";
	ofs<<"#quality control for phenotypes:\n";
	ofs<<"check.trait(c(\"sex\",\"age\") , mydata@phdata)\n";
	ofs<<"#Simple linear regression: \n";
	ofs<<"lmResults <-lm(mydata@phdata$disease ~ mydata@phdata$sex)\n";
	ofs<<"summary(lmResults)\n";
	ofs<<"#For genotype data: \n";
	ofs<<"#-------------------\n";
	ofs<<"#To generate  human-redable format:\n"
		"print(as.character(mydata@gtdata@coding))\n"
		"as.character(mydata@gtdata@strand)\n"
		"as.character(mydata@gtdata)\n"
		"as.numeric(mydata@gtdata)\n";
	ofs<<"sumgt <- summary(mydata@gtdata)\n";
	ofs<<"#For call rate:\n";
	ofs<<"crate <- sumgt[ ,\"CallRate\"]\n";
	ofs<<"hist(crate)\n";
	ofs<<"#If you would like to produce a summary table, showing how many\n"
	   "#markers had call rate lower than, say, 93%, between 93 and 95%, between 95\n"
		"#and 99% and more than 99%. You can use catable() command for that:\n";
	ofs<<"catable(crate, c(0.93, 0.95, 0.99))\n"
		"catable(crate , c(0.95 , 0.96 , 0.97) , cumulative=TRUE)\n"
	    "#For HWE-extact test: \n"
		"hwp <- sumgt[ ,\"Pexact\"]\n"
		"estlambda(hwp)\n"
		"catable(hwp, c((0.05/mydata@gtdata@nsnps), 0.01, 0.05, 0.1))\n"
		"catable(hwp, c((0.05/mydata@gtdata@nsnps), 0.01, 0.05, 0.1), cum = T)\n"
	 	"#To investigate the minor allele frequency (MAF) distribution,\n"
	 	"#the same logic would apply. \n"
	 	"#Calculating MAF: \n"
	 	"afr<-sumgt[, \"Q.2\"]\n"
	   "maf<-pmin(afr, (1 - afr))\n"
		"par(mfcol = c(2, 1))\n"
		"hist(afr)\n"
		"hist(maf)\n"
		"catable(afr, c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99))\n"
		"catable(maf, c(0, 0.01, 0.05, 0.1, 0.2), cum = T)\n"
		"#Use \"perid.summary\", to produce summary SNP statistics per per-son.\n"
		"perid.summary(mydata[1:10, ])\n"
		"descriptives.marker(mydata)\n"
		"het<-perid.summary(mydata)$Het\n"
		"mean(het)\n"
		"catable(het, c(0.1, 0.25, 0.3, 0.35, 0.5))\n"
		"#To produce genomewide descriptives SNP data\n"
		"descriptives.marker(mydata)\n"
		"descriptives.marker(mydata , ids= (mydata@phdata$sex==1))\n"
		"#Genetic data QC:\n"
		"#First without checking HWE: \n"
		"qc1<-check.marker(mydata, p.level = 0)\n"
		"#or, one can also use the follwoing comand.\n"
		"qc1<-check.marker(mydata, callrate=0.95, perid.call=0.95, maf=0.01,p.level=-1)\n"
		"summary(qc1)\n"
		"print(names(qc1))";
		ofs<<"\n";
		ofs<<"#Second round of Quality control, with HWE:\n";
		ofs<<"data1<-mydata[qc1$idok, qc1$snpok]\n"
		"#data1<-Xfix(data1)\n"
		"qc2<-check.marker(data1)\n";
		ofs<<"print(ls())\n";
	/*#Fast score test for association between a trait
	#and genetic polymorphism (GLM)
	#With no adjustment for binary traits, 1 d.f.,
	#the test is equivalent to the Armitage test
	analysis <- qtscore(dm2 , ge03d2exQC , trait ="binomial")
	#Genome-wide significance for a GWA scan. Analysis function is qtscore
	#Genome-wide significance (permutation)
		qtAnalysis <- emp.qtscore(dm2, ge03d2exQC, trait="binomial")
		plot(analysis)
		qtAnalysis2 <- qtscore(dm2 ~age + sex, ge03d2exQC , trait ="binomial")
		plot(qtAnalysis2)

	#Fast case-control analysis by computing chi-square test
	#from 2x2 (allelic) or 2x3 (genotypic) tables
		ccAnalysis <- ccfast("d2m",data=ge03d2exQC)

#Genome-wide significance for a case-control GWA scan.
#Analysis function is ccfast
#Genome-wide significance (permutation) for ccfast() scan
ccAnalysis2 <- emp.ccfast("d2m",data=ge03d2exQC)
#top 10 results
descriptives.scan(ccAnalysis2)
#######################################################################
*/
	ofs.clear(); ofs.close();
	printLIN("\t->**A phenotype file used for GenABEL has been saved as \""+commandFileName+"\".\n");

}


 void CGSNP::write_probabel_mlinfo_file(const vector<CBSNP*>& genInfo, const string & outFileName){
	 printLIN("\t**->Writing mach-formatted mlinfo file:\n");
	 //------------------------------------------------------------------------------------------------------------//
	 	ofstream ofs;
	 	ofs.clear();ofs.close();
	 	string freqFileName=outFileName+".mlinfo";
	 	ofs.open(freqFileName.c_str(),ios::out);
	 	if(!ofs)
	 		error("out file "+outFileName+" does not exit or can not be written out. Please check the filename and file path.\n");
	 	//ofs <<"CHR"<< " ";
	 	ofs<< "SNP"<<"\t"<< "A11"<<"\t"<< "A12"<<"\t"<<"Freq1"<<"\t"<< "MAF" << "\t"<<"Quality"<<"\t"<<"Rsq\n" ;
	 	//total 7 columns
	 	CBSNP* pCBSNP= genInfo[0];
	 	for(unsigned int i=0; i<genInfo.size();++i)
	 	{
	 		pCBSNP= genInfo[i];
	 		//convert_snp_genotype_into_0123(pCBSNP);
	 		//calculate_maf(pCBSNP);
	 		//cout<<boolalpha<< pCBSNP->quality<<" ";
	 		if(pCBSNP->quality ) // default value of quality is true;
	 		{
	 			if(pCBSNP->rsId!="")
	 				ofs<<pCBSNP->rsId<<"\t";
	 			else
	 				ofs<<pCBSNP->snpName<<"\t";
	 			//allele1
	 			if(pCBSNP->allele1 !="" )
	 				ofs<< pCBSNP->allele1 <<"\t";
	 				else ofs <<"-"<< "\t";
	 			//allele2
	 			if(pCBSNP->allele2 !="" )
	 				ofs<< pCBSNP->allele2 <<"\t";
	 				else ofs <<"-"<< "\t";
	 			//write freq1
	 			ofs<<pCBSNP->freq<<"\t";
	 			ofs<<pCBSNP->maf<<"\t";
	 			//imputation quality value
	 			ofs<<0.0<<"\t";
	 			ofs<<pCBSNP->rsq<<"\t";
	 			ofs<<"\n"; //end of line
	 		}
	 	}
	 	pCBSNP= genInfo[0];
	 	//end of writing
	 	ofs.clear();ofs.close();
	 	printLIN("\t***->File in MLINFO format has been written out and saved as \""+\
	 			freqFileName+"\". \n");
//------------------------------------------------------------------------------------------------------------//
}
void CGPED::write_probabel_mlprob_file(const vector<CBPED*>& vPed, const vector<CBSNP*>& genInfo, const string & outFileName){
	printLIN("\t**->Writing mach-formatted mlprob and mldose files:\n");
	//----------------------------------#--------------------------------------#
	//----------------------------------#--------------------------------------#
	// writing mlprob files.
	//----------------------------------#--------------------------------------#
	const CBPED*  pCBPED	 =vPed[0];
	const unsigned int IND_SZ=vPed.size();
	CBSNP* pCBSNP =genInfo[0];
	//for mlprob
	//-------------------------------------------------------//
	string mlprobFileName=outFileName+".mlprob";
	ofstream ofs;
	ofs.clear();ofs.close();
	ofs.open(mlprobFileName.c_str());
	if(!ofs)
		error("out file "+outFileName+" does not exit or can not be written out.\n");
	ofs.clear();
	//for mldose file
	//-------------------------------------------------------//
	string mldoseFileName=outFileName+".mldose"; //dose file
		ofstream ofs1;
	ofs1.clear();ofs1.close();
	ofs1.open(mldoseFileName.c_str());
	if(!ofs1)
			error("out file "+outFileName+" does not exit or can not be written out.\n");
		ofs1.clear();
	//-------------------------------------------------------//
			//collect bp
	for(unsigned int i=0;i<IND_SZ;i++)
	{
		pCBPED =vPed[i];

		if(pCBPED->quality)
		{

			//----------------------------------------------//
			ofs<<i+1<<"->";
			if(pCBPED->indId!="")
			{
				ofs<<pCBPED->indId<<" ";
			}
			else {
				ofs<<"iid_"<<i+1<<" ";
			}
			//mlprob
			ofs<<"ML_PROB"<<" ";
			// for mldose
			//----------------------------------------------//
			ofs1<<i+1<<"->";
			if(pCBPED->indId!="")
			{
				ofs1<<pCBPED->indId<<" ";
			}
			else {
				ofs1<<"iid_"<<i+1<<" ";
			}
			//mlprob
			ofs1<<"ML_DOSE"<<" ";
			//---------------------------------------------//
			for(unsigned int j=0;j<genInfo.size(); ++j)
			{

				pCBSNP=genInfo[j];
			//cout <<"(pCBSNP->quality): "<< boolalpha << pCBSNP->quality <<endl;
				if(pCBSNP->quality)
				{							// write geno prob if they are given.

					//case1 if genotype probabilities are given
					if((pCBSNP->pgeno1.size()==IND_SZ) &&(pCBSNP->pgeno2.size()==IND_SZ)&&(pCBSNP->pgeno3.size()==IND_SZ))
					{
						ofs <<" "<< pCBSNP->pgeno1[i] << " "<< pCBSNP->pgeno2[i] << " ";
						// for mldose file
						ofs1<<" "<<  (1.0*pCBSNP->pgeno2[i] +2.0*pCBSNP->pgeno1[i])<<" ";
					}
					//case 2 if genotype probabilities are not given
					else
					{
						if( !pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if homozygote1  00
						{
							ofs <<" "<< "1" <<" "<< "0" <<" ";
							// for mldose file
							ofs1<<" "<<  (2.0*1.0 +1.0*0.0)<<" ";
						}
						else if( pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if homozygote2 11
						{
							ofs <<" "<< "0" <<" "<< "0" <<" ";
							// for mldose file
							ofs1<<" "<<  (2.0*0.0 +1.0*0.0)<<" ";
						}
						else if( !pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if heterozygote
							{
								ofs <<" "<< "0" <<" "<< "1" <<" ";
								// for mldose file
								ofs1<<" "<<  (2.0*0.0 +1.0*1.0)<<" ";
							}
						else if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) //if missing
							{
								ofs <<" "<<atoi(CBSNP::_miss_gvalue.c_str()) << " "<<CBSNP::_miss_gvalue << " ";
								// for mldose file
								ofs1<<" "<<-1.0<<" ";
							}
						else
						{
							error("problem in writing out "+change_int_into_string(j) +"th SNP.\n");
						}
						//	ofs <<"TEst\n";
					}
				} //end of if pedind quality

			}//end geninfo loop
			ofs<< "\n"; //end of line
			ofs1<< "\n";//end of line


		} //end of ped quality
	}// end of vped loop
	ofs.clear();	 ofs.close();
	ofs1.clear();	 ofs1.close();
	printLIN( "\t***->mach-formtted mlprob file has been written out and saved as \""+mlprobFileName+"\". \n");
	printLIN( "\t***->mach-formtted mlprob file has been written out and saved as \""+mldoseFileName+"\". \n");

//----------------------------------------------------------------------------------------------------------------------------------//
}







