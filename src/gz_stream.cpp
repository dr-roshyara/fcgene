#include "gz_stream.h"
#include <stdarg.h>
#include<string>
#include<cstring>
#include<sstream>
bool rFile::checkFileExists( const std::string & f )
{
	rFile inp(f);
  if(inp.fail())
    {
      inp.clear(ios::failbit);
      inp.close();
      string msg = "No any File named [ " + f + " ] exists.";
      error(msg);
    }
  inp.close();
  return true;
}
std::vector<std::string> rFile::split_line_as_vector(const string &sline){
	std::string buffer;
	std::stringstream ss(sline);
	std::vector<std::string> tokens;
	while(ss>>buffer)
		tokens.push_back(buffer);
	return tokens;
}

//#define _HAVE_ZLIB__

#ifdef _HAVE_ZLIB__

//define functions of rFile
std::string rFile::readLine()
{
  std::string sline;
  gzifstream & gg = *this;
  std::getline( gg ,  sline );
  return sline;
}

IFS::IFS(const char * filename, const char * mode)
   {
   // If Zlib cant read file that is larger than 2Gb, then  just use  regular fopen file if gzopen fails and the file name has no name ending with ".gz". 

   _case_gz = true;
   //_file_gz = (void**) gzopen(filename, mode);
   _file_gz = gzopen(filename, mode);
   bool _tfValue=false;
   if (_file_gz == NULL)
      {
      int _das_letzte_char = 0;

      while (filename[_das_letzte_char] != 0) 
		_das_letzte_char++;
		//
		_tfValue =(_das_letzte_char >= 3) ;
		_tfValue=_tfValue &&(	filename[_das_letzte_char - 3] == '.' &&
                           filename[_das_letzte_char - 2] == 'g' &&
                           filename[_das_letzte_char - 1] == 'z');
						   
      if (_tfValue)
			return;

      _case_gz = false;
      _file = fopen(filename, mode);
      }
   }
#else 
std::string rFile::readLine()
{
  std::string sline;
  ifstream & gg = *this;
  std::getline( gg ,  sline );
  return sline;
}

#endif //_HAVE_ZLIB__

/*
int main (int argc , char** argv)
{
	//IFS IFS;
	 const char* filename =argv[1];
	 IFS file;
	 file	= _ifs_open(filename, "rb");
	 if (file == NULL){ printf("Sample Marker Information File [%s] could not be opened\n", filename); exit(1);}
	while (!_ifs_eof(file))
	{
			char myChar;
			myChar= _ifs_getc(file);
			std::cout << myChar << "" ; 
			
	}
	_ifs_close(file);
	//IFS IFS =ifs_open(filename,"rb");
	std::string fs= filename;
	//rFile ofs; 
	rFile::checkFileExists(fs);
	rFile ofs(fs);
	while(!ofs.eof())
	{
		std:: string ss =ofs.readLine();
		//&std::vector<std::string> line = rFile::split_line_as_vector(ss);
		//for(unsigned int i=0; i<line.size(); ++i)
		//	std:: cout << line[i]<< " ";
		//std::cout<< std::endl;
	}	
return 0;}
*/
