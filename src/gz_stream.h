#ifndef _GZSTREAM_H__
#define _GZSTREAM_H__
#include <stdio.h>
#include<iostream>
#include<vector>
#include <fstream>
#include<cstring>
#include<string>
#include "helper.h"

#include <stdio.h>
#include<iostream>
#include<vector>
#include"zfstream.h"

#ifdef  _HAVE_ZLIB__
#include <zlib.h>

//create a class to read files
class rFile : public gzifstream {
private: 
	ifstream _file ;
	gzifstream _gz_file;
public:
     bool _case_gz;
	
	union {
		 //gzifstream _file_gz; 
		 //ifstream   _file;
	};
	rFile (void){}
	rFile (const std:: string & _fileName , ios_base::openmode mode = ios_base::in ) : _gz_file(_fileName.c_str(),mode)
	{
		
		x=11;
		printf("your group is: %i \n", x );	
		//exist_file();
	}
	std:: string readLine();
	std::vector<std::string> split_line_as_vector(const std::string & sline);
	static bool checkFileExists( const std::string & f );
	private: 
	int x ;
};

//writing outfile
class wFile : public gzofstream{
public:
	bool _case_gz;
	wFile (){}
	wFile(const std::string & _fileName, ios_base::openmode mode=ios_base::out):gzofstream(_fileName.c_str(),mode){ }

};

class IFS  {
   public:
      bool _case_gz;
	  int x;
      union
         {
			//A union is a special data type available in C that enables you to store different data types 
			//in the same memory location. You can define a union with many members, but only one member 
			//can contain a value at any given time. Unions provide an efficient way of using the same memory location for multi-purpose.


         gzFile _file_gz;
         FILE * _file;
         };
		 
	
   IFS()
      {
		x =1;
      _case_gz = false;
      _file = NULL;
      //std::cout<<"your number is:" << x <<std::endl;
      }

   IFS(const char * filename, const char * mode);
   // 
   
   operator void * ()
      { return _case_gz ? (void *) _file_gz : (void *) _file; }

   IFS operator = (const IFS & rhs)
      {
      if ((_case_gz = rhs._case_gz) == true)
         _file_gz = rhs._file_gz;
      else
         _file = rhs._file;

      return *this;
      }

   IFS operator = (FILE * rhs)
      {
      _case_gz = false;
      _file = rhs;
      return *this;
      }

   IFS operator = (gzFile & rhs)
      {
      _case_gz = true;
      _file_gz = rhs;
      return *this;
      }

   bool operator == (void * rhs)
      {
      if (rhs != NULL)
         return false;
      return _case_gz ? _file_gz == rhs : _file == rhs;
      }
   };
   //
   

inline IFS _ifs_open(const char * filename, const char * mode)
   { IFS file(filename, mode); return file; }

inline int _ifs_close(IFS & file)
   {
   int result = file._case_gz ? gzclose(file._file_gz) : fclose(file._file);
   file._file_gz = NULL;
   return result;
   }

inline int _ifs_getc(IFS & file)
   { return file._case_gz ? gzgetc(file._file_gz) : fgetc(file._file); }

inline void ifrewind(IFS & file)
   { if (file._case_gz) gzrewind(file._file_gz); else rewind(file._file); }

inline int _ifs_eof(IFS & file)
   { return file._case_gz ? gzeof(file._file_gz) : feof(file._file); }

inline unsigned int _ifs_read(IFS & file, void * buffer, unsigned int size)
   { return file._case_gz ? gzread(file._file_gz, buffer, size) :
                          fread(buffer, 1, size, file._file); }

inline unsigned int ifwrite(IFS & file, void * buffer, unsigned int size)
   { return file._case_gz ? gzwrite(file._file_gz, buffer, size) :
                          fwrite(buffer, 1, size, file._file); }

#else 
#include <stdio.h>
#include<fstream>
class rFile:public std::ifstream{
public:
	rFile (void){ }
	rFile (const std:: string & _fileName , ios_base::openmode mode = ios_base::in ):ifstream(_fileName.c_str(),mode)
	{
		x=22;
		printf("your group is:: %i \n", x );	
		
	}
	std:: string readLine();
	std::vector<std::string> split_line_as_vector(const std::string & sline);
	static bool checkFileExists( const std::string & f );
private: 
int x ;
};

//writing outfile
class wFile:  public std::ofstream
{
public:
	wFile (){}
	wFile(const std::string & _fileName, ios_base::openmode mode=ios_base::out): ofstream(_fileName.c_str(),mode){ }

};

class IFS
   {
   public:
	 int x ;	
	 FILE * _file;

      IFS()
         { _file = NULL;
			x=2 ; 
			//printf("your number is::  %i\n" ,x) ; //debug

		 }
      IFS(const char * filename, const char * mode)
         { _file = fopen64(filename, mode); }
      ~IFS()
         { }

      operator FILE *()
         { return _file; }

      IFS & operator = (FILE * rhs)
         { _file = rhs; return *this; }

      IFS & operator = (const IFS & rhs)
         { _file = rhs._file; return * this; }

      bool operator == (void * rhs)
         {
         if (rhs != NULL)
            return false;
         return _file == rhs;
         }
   };

inline IFS _ifs_open(const char * filename, const char * mode)
   { IFS file(filename, mode); return file; }

inline int _ifs_close(IFS & file)
   {
   int result = fclose(file._file);
   file._file = NULL;
   return result;
   }

inline int _ifs_getc(IFS & file)
   { return fgetc(file._file); }

inline void ifrewind(IFS & file)
   { rewind(file._file); }

inline int _ifs_eof(IFS & file)
   { return feof(file._file); }

inline unsigned int _ifs_read(IFS & file, void * buffer, unsigned int size)
   { return fread(buffer, 1, size, file._file); }

inline unsigned int ifwrite(IFS & file, void * buffer, unsigned int size)
   { return fwrite(buffer, 1, size, file._file); }
   


#endif		//  _HAVE_ZLIB__

#endif 		//GZSTREAM_H__
