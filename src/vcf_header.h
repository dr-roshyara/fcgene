#ifndef _VCF_HEADER__
#define _VCF_HEADER__
#include <ctime> 
#include <iostream>
#include<string>
#include<cstring>
#include<set> 
#include<map>
#include<vector>
#include<algorithm>
#include "helper.h"
#include "gz_stream.h"
/** 	*File meta-information is included after the ##  string and must be key = value pairs.
	*single 'fileformat' field is always required, must be the first line in the file, and details the VCF format version number. 
	*For example, for VCF version 4.1, this line should read:
	*##fileformat=VCFv4.1
	*It is strongly encouraged that information lines describing the 
	*INFO, 
	*FILTER and
	*FORMAT  entries used in the body of the VCF file be included in the meta-information section. 
	*Although they are optional, if these lines are present then they must be completely well-formed.
*/
enum type_name {Integer=0, Float=1, Character=2, String=3, Flag=4};

inline void increment(int & ii){ ++ii;} 

class META_INFO;
typedef 	std::map <int, META_INFO> 				map_int_meta;
typedef		std::map<std::string, META_INFO>		map_string_meta;
typedef 	std::string								string_meta;
//typedef 	std::map <double, META_INFO> 			map_double_meta;
//typedef 	std::map<std::string, META_INFO> 		map_string_meta;
//typedef 	std::map <bool, META_INFO> 				map_bool_meta;
const string_meta		file_format		="VCFv4.1";
const string_meta		source			="fcGENE";
const string_meta 		start_info		="##INFO=";
const string_meta		start_format	="##FILTER=";
const string_meta		start_filter	="##FORMAT=";
const string_meta		start_source	="##source=";
const string_meta		start_fileformat="##fileformat=";

//
class VCF_HEADER{
public: 
	
	string_meta			fileDate;
	string_meta 		vcf_version; 	// vcf version 
	string_meta			source;
	map_int_meta 		format;  	 	//  format can be genotype genotype quality , read depth haplotype quality 
	map_int_meta 		info; 	
	map_int_meta		filter;
	static std::vector<string_meta > meta_id;
//static members
	/* *   Possible Types for INFO fields are: Integer, Float, Flag, Character, and String. 
		*The Number entry is an Integer that describes the number of values that can be included with the INFO field. 
		*For example, if the INFO field contains a single number, then this value should be 1; 
		*if the INFO field describes a pair of numbers, then this value should be 2 and so on. 
		*If the field has one value per alternate allele then this value should be 'A'; 
		*if the field has one value for each possible genotype (more relevant to the FORMAT tags) then this value should be 'G'.  
		*If the number of possible values varies, is unknown, or is unbounded, then this value should be '.'. 
		*The 'Flag' type indicates that the INFO field does not contain a Value entry, and hence the Number should be 0 in this case. 
		*The Description value must be surrounded by double-quotes. 
		*Double-quote character can be escaped with backslash (\") and backslash as \\.
	*/
 public: 
	VCF_HEADER( );
	~VCF_HEADER(){} 
	//string_meta create_info_line(const map_int_meta &info )
	static string_meta  create_meta_info(const std::vector<string_meta>&meta_id, const META_INFO & pMETA_INFO,const bool no_FILTER=true); // no_FILTER assumes no meta information for meta info:filter 
	static string_meta	create_type_string(const enum type_name  typ ); 

	//static string_meta create_info_line(const map_int_meta &info );
	//int add_info_in_header(const std::string & _mystr, const int& _key);
	 
	 //int 	fuege_info_ein(const string_meta & line);
	 //static void screibe_info(const map_int_meta &info); 
	
};

class META_INFO{ // this class describes meta information for each of the header  given after ## in vcf file 
public:
	string_meta			name;
	string_meta 		number; 
	enum type_name 		type;
	string_meta			description;
	int					len; 
	META_INFO(): name(""),number("0"),type(Integer),description(""){  };
	~META_INFO(){ //std::cout<<"bye meta info \n";
	};
public:
static string_meta  create_number_info();
static string_meta  create_type_info();
static string_meta  crate_description_info();

	//##INFO=<ID=ID,Number=number,Type=type,Description=”description”>
};
#endif 
/*
##fileformat=VCFv4.1
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
*/