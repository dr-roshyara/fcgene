#include "vcf.h"
#include <fstream>
#include "bpar.h"
//void write_vcf_format( const CBSNP & pCBSNP, const CBPED & pCBPED , const string & outFile ); 
//VCF::meta_id.push_back("Number");
//meta_id.push_back("Type");
//meta_id.push_back("Description");
VCF::VCF(){
	//create a  default HEADER 
	//enum  meta_id _ids	=ID;
	// We create a defult header for snp data if they are read form snp data and are to be converted SNPdata format 
	pSNPDATA_HEADER =new VCF_HEADER;
	//info 
	//##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
	//AF allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
	pSNPDATA_HEADER->info[0].name ="AF";
	pSNPDATA_HEADER->info[0].number ='A';
	pSNPDATA_HEADER->info[0].type =Float;
	pSNPDATA_HEADER->info[0].description="Allele Frequency"; // for alternative frequency 

	//filter 
	pSNPDATA_HEADER->filter[0].name ="PASS";
	pSNPDATA_HEADER->filter[0].number =1;
	pSNPDATA_HEADER->filter[0].type =String;
	pSNPDATA_HEADER->filter[0].description="Passed variant FILTERs";
	//genotype
	pSNPDATA_HEADER->format[0].name ="GT";
	pSNPDATA_HEADER->format[0].number ="1";
	pSNPDATA_HEADER->format[0].type =String;
	pSNPDATA_HEADER->format[0].description ="Genotype";
	//genotype probabilty
	pSNPDATA_HEADER->format[1].name ="GP";
	pSNPDATA_HEADER->format[1].number ='G';
	//The Number entry is an Integer that describes the number of values that can be included with the INFO field. For example,
	//if the INFO field contains a single number, then this value should be 1; if the INFO field describes a pair of numbers, then this value should be 2 and so on.
	//If the field has one value per alternate allele then this value should be 'A';
	//if the field has one value for each possible genotype (more relevant to the FORMAT tags) then this value should be 'G'.

	pSNPDATA_HEADER->format[1].type =Float;
	pSNPDATA_HEADER->format[1].description ="Genotype probability distribution";
	// please always give allele frequency to the vcf format data 
	//pSNPDATA->info[]
	/*pSNPDATA->filter[ids[0]] ="PASS";
	pSNPDATA->filter[ids[1]] ="PASS";
	pSNPDATA->filter[ids[2]] ="PASS"
	pSNPDATA->filter[ids[3]] ="PASS";
	*/
	
}
VCF::~VCF(){
	//create a  default HEADER
	pSNPDATA_HEADER->info.clear();
	pSNPDATA_HEADER->filter.clear();
	pSNPDATA_HEADER->info.clear();
	delete pSNPDATA_HEADER;
	
	//pSNPDATA->info[]
	//pSNPDATA->filter[]
	
}
void VCF::write_vcf_format(const vector<CBPED*>&pedVec,const vector<CBSNP*>&genVec,const std::string & outFile,const std::vector<std::string>& ids, const bool _uncompress){
	//int i=0;
	//declairing stream 
	std::string oFile =outFile+"_vcf.gz";
	if(_uncompress)
		oFile =outFile+".vcf";
	#ifndef  _HAVE_ZLIB__
	 		oFile =outFile+".vcf";
	#endif //	_HAVE_ZLIB__
	
	printLIN("*->Writing file in VCF format:\n");
	wFile ofs;
	if(!_uncompress)
	 ofstream ofs;
	 
	 ofs.clear(); ofs.close();
		
		
	ofs.open(oFile.c_str(),std::ios::out);
	
	//ofs<<"#test\n";
	//file format 
	ofs<<start_fileformat<<file_format;
	ofs<<"\n"; 		
	// source 
	ofs<<start_source <<source;
	ofs<<"\n";
	// filter  
	map_int_meta::iterator it;

	//##FILTER=<ID=q10,Description="Quality below 10">
	
	std::string _ss;
	//info meta 
	for(map_int_meta::iterator it =pSNPDATA_HEADER->info.begin(); it !=pSNPDATA_HEADER->info.end(); it++)
	{
		
		_ss	=VCF_HEADER::create_meta_info(ids, it->second); 
		ofs<< start_info <<  _ss ; //"<"<<ids[0]<< "="<< (it->second).name; 
		
	}
	
	//filter meta 
	for(map_int_meta::iterator it =pSNPDATA_HEADER->filter.begin(); it !=pSNPDATA_HEADER->filter.end(); it++)
	{
		
		_ss	=VCF_HEADER::create_meta_info(ids, it->second,false); 
		ofs<< start_filter <<  _ss ; //"<"<<ids[0]<< "="<< (it->second).name; 
		
	}
	//format meta 
	// this format is at the moment only two. You can use either genotype
	// or genotype probability
	// or both
	// if we use genotype probability, we can use both
	for(map_int_meta::iterator it =pSNPDATA_HEADER->format.begin(); it !=pSNPDATA_HEADER->format.end(); it++)
	{

			_ss	=VCF_HEADER::create_meta_info(ids, it->second);
			ofs<< start_format <<  _ss ; //"<"<<ids[0]<< "="<< (it->second).name;

	}

	//written should be the following 
	// before calculate Allele frequency 
	// 
	ofs<<"#CHROM"<<"\t"<<"POS" <<"\t" <<"ID"<<"\t"<<"REF"<<"\t"<<"ALT"<<"\t"<<"QUAL"<<"\t"<<"FILTER"<<"\t"<<"INFO"<<"\t"<<"FORMAT" ;
	//now pedVec and genVec 
	//const CBPED* pCBPED	=pedVec[0];
	const CBSNP* pCBSNP	=genVec[0];
	for(unsigned int i=0;i<pedVec.size();++i)
	{
	 ofs <<"\t"<< (pedVec[i])->indId;
	}
	ofs<<"\n";
	// write snp wise now 
	for(unsigned int j=0;j<genVec.size();++j)
	{
		pCBSNP = genVec[j];
		if(pCBSNP->quality)
		{
			ofs<<pCBSNP->nchr <<"\t"<<pCBSNP->bp;
			if(pCBSNP->rsId!="")
				ofs<<"\t"<<pCBSNP->rsId;
			else if(pCBSNP->snpName!="")
				ofs<<"\t"<<pCBSNP->snpName;
			else ofs<<"\t"<<"rsid_"<<j;

			//write first allele
			if(pCBSNP->allele1!="")
				ofs<<"\t"<< pCBSNP->allele1;
			else
				ofs<<"\t"<<".";
			//write second allele
			if(pCBSNP->allele2!="")
				ofs<<"\t"<< pCBSNP->allele2;
			else
				ofs<<"\t"<<".";
			//write QUAL
			ofs <<"\t"<<".";
			//write filter criteria
			 ofs <<"\t"<<"PASS";
			//write INFO
			ofs <<"\t"<<"AF="<<round_nplaces((1.0-pCBSNP->freq),3);
			//write format
			ofs <<"\t"<<"GT:GP";
			//now write genotypes
			//homozygotye
			bool _tmp_bool1=false;
			bool _tmp_bool2=false;
			const CBPED* pCBPED	=pedVec[0];

			for(unsigned int i =0; i<pedVec.size();++i)
			{
				pCBPED	=pedVec[i];

				//--------------------------------------------------------------------//

				 if(pCBPED->quality)
				 {
					_tmp_bool1	=pCBSNP->geno1[i] ;
					_tmp_bool2	=pCBSNP->geno2[i] ;
					//false && false 	=homo1
					//true && true 	=homo2
					// fase & true 	hetero
					//true && false 	=missing
					if(!_tmp_bool1&&!_tmp_bool2)
					{
						ofs<<"\t"<<"0/0";
					}
					else if(_tmp_bool1&&_tmp_bool2)
					{
						ofs<<"\t"<<"1/1";
					}
					else if(!_tmp_bool1&&_tmp_bool2)
					{
						ofs<<"\t"<<"0/1";
					}
					else{
						ofs<<"\t"<<"./.";
					}
					//case1 if genotype probabilities are given
					if((pCBSNP->pgeno1.size()==pedVec.size()) &&(pCBSNP->pgeno2.size()==pedVec.size())&&(pCBSNP->pgeno3.size()==pedVec.size()))
					{
						ofs <<":"<< pCBSNP->pgeno1[i] << ","<< pCBSNP->pgeno2[i] << ","<< pCBSNP->pgeno3[i];
					}
					//case 2 if genotype probabilities are not given
					else
					{
						if( !pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) // if homozygote1  00
							ofs <<":"<< "1" <<","<< "0" <<","<< "0";
						else if( pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if homozygote2 11
							ofs <<":"<< "0" <<","<< "0" <<","<< "1";
						else if( !pCBSNP->geno1[i] &&  pCBSNP->geno2[i]  ) // if heterozygote
							ofs <<": "<< "0" <<","<< "1" <<","<< "0";
						else if( pCBSNP->geno1[i] &&  !pCBSNP->geno2[i]  ) //if missing
							ofs <<":"<<atoi(CBSNP::_miss_gvalue.c_str()) << ","<<CBSNP::_miss_gvalue << ","<<CBSNP::_miss_gvalue;
						else
						{
							error("problem in writing out "+change_int_into_string(j) +"th SNP.\n");
											}
						//	ofs <<"TEst\n";
					}
				} //end of if pedind quality

			}
			ofs<<"\n";
		}
	}
	// write 
    //#CHROM
    //POS
    //ID
    //REF
    //ALT
    //QUAL
    //FILTER
    //INFO

	
	ofs.close();
	printLIN("\t**->VCF format file has been successfully written and saved as \""+oFile+"\".\n");
	

} 

