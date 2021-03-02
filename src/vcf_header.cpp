#include"vcf_header.h"
std::vector<string_meta> VCF_HEADER::meta_id;
VCF_HEADER::VCF_HEADER( ) {
	vcf_version 	="4.1";
	VCF_HEADER::meta_id.resize(0);
	VCF_HEADER::meta_id.push_back("ID");
	VCF_HEADER::meta_id.push_back("Number");
	VCF_HEADER::meta_id.push_back("Type");
	VCF_HEADER::meta_id.push_back("Description");
}
/*
int VCF_HEADER::fuege_info_ein(const string_meta & line){
	std:: cout << "hi fuege:_info_ein: \n"; //debug 
	
return 0;}
*/

string_meta VCF_HEADER::create_meta_info(const std::vector<string_meta>&meta_id, const META_INFO & pMETA_INFO , const bool no_FILTER){
	string_meta _ms;
	_ms		="<";
	_ms		+=meta_id[0]+"="+pMETA_INFO.name+",";
	if(no_FILTER)
	{
		_ms		+=meta_id[1]+"="+pMETA_INFO.number+",";
		_ms		+=meta_id[2]+"="+VCF_HEADER::create_type_string(pMETA_INFO.type)+","; 
	}
	//_ms		+=meta_id[2]+"="+pMETA_INFO.type+",";
	_ms		+=meta_id[3]+"=\""+pMETA_INFO.description+"\""; 
	_ms		+= ">\n";
	
return _ms;}
string_meta VCF_HEADER::create_type_string(const enum type_name  typ){
	string_meta _sm="."; 
	//std::cout<< "type. "<< typ <<std::endl; //debug 
	if(typ==0)
		_sm		="Integer";
	else if(typ==1)	
		_sm		="Float";
	else if(typ==2)
		_sm		="Character";
	else if(typ==3)
			_sm 	="String";
	else if(typ==4)	
	_sm		="Flag";
	/*
	switch(typ)
		{

		case 0:
			_sm		="Integer";
		case 1:
			_sm		="Float";
		case 2:
			_sm		="Character";
		case 3: 
			_sm 	="String";
		case 4: 
			_sm		="Nothing";
	}
	*/
	return _sm;
}	

/*
string_meta VCF_HEADER::create_info_line(const map_int_meta &info )
{
	string_meta ss="";
	ss=ss+ start_info;
	//ss+=
	//NS,Number=1,Type=Integer,Description="Number of Samples With Data">	
return ss;} 

//int VCF_HEADER:: add_info_in_header(const std::string & _mystr, const int& _key){}
string_meta  META_INFO::create_meta_info(){
	string_meta _ss="";
	//ss=	"ID="+name;
return _ss;}
string_meta  META_INFO::create_number_info(){
	string_meta _ss="";
	//ss =""+nu
return _ss;}
string_meta  META_INFO::create_type_info(){
	string_meta _ss="";
	
return _ss;}
string_meta META_INFO::crate_description_info(){
	string_meta _ss="";
	
return _ss;}
*/