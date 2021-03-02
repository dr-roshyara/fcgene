#ifndef _VCF_H__
#define _VCF_H__
#include "vcf_header.h"
#include  "helper.h"
#include "gz_stream.h"
#include "machdat.h"
#include "machped.h"
 

class VCF{
public: 

VCF_HEADER HEADER;
VCF_HEADER*  pSNPDATA_HEADER;
VCF();
~VCF();
//void write_vcf_format( const CBSNP & pCBSNP, const CBPED & pCBPED , const string & outFile ); 
 void  write_vcf_format(const vector<CBPED*>&pedVec, const vector<CBSNP*>&genVec, const std::string & outFile , const std::vector<std::string>& ids,const bool _uncompress);
};
#endif //_VCF_H__