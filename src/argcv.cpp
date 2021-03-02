/*
 * argcv.cpp
 *
 *  Created on: 03.11.2011
  *  Created on: Dec 29, 2011
 *  Copyright (C) <2012>  <Nab Raj Roshyara>
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.

 */
#include "argcv.h"
#include<typeinfo>
#include<algorithm>
#include<vector>
using namespace std;
//#include "bpar.h"
extern ofstream LIN;
 string Argcv::msg="";
 //int Argcv::argCount=0;

vector<Bpar*>  Argcv::bparam_vec;
vector<string> Argcv::stcomVec;
vector<string> Argcv::comVec;
vector<string> Argcv::original;
vector< bool > Argcv::parsed;
vector<bool>   Argcv::option;
int 		   Argcv::argCount=0;
int 		   Argcv::bparas_vec_sz=0;

Argcv::Argcv(int arg_count, char** argv)
{
	argCount	=arg_count;
	original.resize(argCount);
	stcomVec.resize(argCount);
	parsed.resize(arg_count,false); // parsed bool
	option.resize(arg_count,false); // this will show the  undefined commands till now
	for( int i=0;i<argCount;++i)
		comVec.push_back(argv[i]);
	 original=comVec;
	 stcomVec = comVec; // sorted command vector
	 stcomVec.erase(stcomVec.begin());
	 sort(stcomVec.begin(),stcomVec.end());
	// cout<< comVec.size()<<endl;
	 //crate a class for basic parameters
	 Bpar * const p0Bpars =new Bpar; // class for basic parameters class . this will save basic parameters.
	 bparam_vec.push_back(p0Bpars);

	  while(compfun("--new-start"))
	 {
		 int _start_idx=scsFind("--new-start");
		 int _end_idx=scsFind("--new-end");
		 if(_start_idx<_end_idx) // if --new-start and --new-end found then create new parameter class .
		 {
			 Bpar* pBpars =new Bpar; // create Bpars
			 //const int _length_idx= _end_idx-_start_idx;
			 vector<string>::iterator it_ori =comVec.begin();
			//ori_merge_commands.assign(it_ori+(_start_idx),it_ori+(_end_idx+1));
			pBpars->original.assign(it_ori+(_start_idx+1),it_ori+(_end_idx)); 	//this will leave --merge-start and --merge-end here
			pBpars->parsed_comms.assign(pBpars->original.size(),false); 		// this will tell if all commands are pased in every parameter class.
			pBpars->nArgs =pBpars->original.size();								//this is the number of parameters in bpar class
			//cout << pBpars->nArgs<<endl;
			 if(pBpars->search_string("--out"))
					 pBpars->find_double_opts(pBpars->output_fileName, "--out");
			//pBpars->stcomVec= pBpars->original; 								// this just saves sorted parameters for every class.
			//sort(pBpars->stcomVec.begin(),pBpars->stcomVec.end() );
			//cout <<"merge commands:\n"; //debug
			//for(unsigned int i =0; i<pBpars->original.size(); ++i)
			//cout<<pBpars->original[i]<<" ";//debug
			//cout <<endl; //debug
			//remove commands in between merge
			comVec.erase(it_ori+(_start_idx),it_ori+(_end_idx+1));
			// debug
		 	//cout << comVec.size()<<endl;
			//for(unsigned int i=0; i<comVec.size();i++)
			//	cout << comVec[i]<< " "; //debug
			// afer erasing the commands in between --new-start and --new-end, we can
			// assign sort the remaining
			stcomVec = comVec; // sorted command vector
			sort(stcomVec.begin(),stcomVec.end());
			vector<bool>::iterator it_parse=parsed.begin()+(_start_idx); 				// this is necessary here
			for(it_parse=(parsed.begin()+(_start_idx));it_parse<(parsed.begin()+(_end_idx+2)); ++it_parse)
				*(it_parse) =true;
			//parsed.assign(it_parse+(_start_idx),it_parse+(_end_idx+1),true); 		//this is necessary here
			//cout<< parsed.size();//debug
			// say that this Bpar class has --new command
			 pBpars->given_new_command=true; // this will tell that given new command
			 bparam_vec.push_back(pBpars);

		 }
		 else
		 {
			 error("Problem in Parsing \"--new-start\" and \"--new-end\" commands. Each \"--new-start\" command should be followed by a \"--new-end\" command.\n");
		 }

	 }
	//now add to origianl
	 Bpar * const _p_0Bpars =bparam_vec[0];
	_p_0Bpars->original=comVec;
	_p_0Bpars->original.erase(_p_0Bpars->original.begin());
	_p_0Bpars->parsed_comms.assign(_p_0Bpars->original.size(),false);
	_p_0Bpars->nArgs =_p_0Bpars->original.size();
	  parsed[0]=true;
	 if( _p_0Bpars->search_string("--out" ))
		 _p_0Bpars->find_double_opts(_p_0Bpars->output_fileName, "--out");
	 Bpar::cpars_sz=bparam_vec.size();
	 //cout << comVec.size() << original.size() <<endl;
	 //comVec.clear();
	 //stcomVec.clear();
	 //comVec.resize(0);
	 //stcomVec.resize(0);
}
// parse arguments for general commands like freq
// bool Argcv::compfun(const string & s )	{bool bv = binary_search(stcomVec.begin(),stcomVec.end(),s); return bv;}
Argcv::~Argcv(){


	if(bparam_vec.size()>0)
	{
		for(vector<Bpar*>::iterator it=bparam_vec.begin(); it!=bparam_vec.end();++it)
			delete *it;
	}
	bparam_vec.clear();
	bparam_vec.resize(0);

}

// parse subArg
vector<string> Argcv::parse_subArg(string& subArg, const string& seperator)
{
	vector<string> myVec;  
	int sLength	=subArg.length();
	size_t found; 
	found 	= subArg.find(seperator);
	if(found==string::npos&&sLength>0)
		myVec.push_back(subArg);
	else
	{
		do
		{
			//cout << subArg.substr(0,found)<<  " "; // debug 
			myVec.push_back(subArg.substr(0,found));
			subArg		=subArg.substr(found+1,sLength);
			sLength		=subArg.length();
			found 		= subArg.find(seperator);
			if((found==string::npos)&&sLength>0)
			{
				//cout<< "last: "<< subArg;
				myVec.push_back(subArg);
			}
		
		}while(	found!=string::npos);
	}
	
	
return myVec;}

// single command string finding
int Argcv::scsFind(const string& s)
{ // single command string finding
	 bool tfValue =compfun(s);// binary_search(stcomVec.begin(),stcomVec.end(),s);
	 int dis=0;
	 if(tfValue){
		vector<string>::iterator it;
		 it = find (comVec.begin(), comVec.end(), s); // finding the string in the argv
		  dis = int(distance(comVec.begin(),it)); // index of the string
		 parsed[dis]=tfValue;
	}
	else{
	 msg = "An argument for "+s+ " is missing.\n";
	 error(msg);
	}
	return dis;

}


void Argcv::dcsFind(string&coms , const string& s)
{
	 bool tfValue =compfun(s);// binary_search(stcomVec.begin(),stcomVec.end(),s);
	//the following code is just to know
	// for(unsigned int i=0; i<stcomVec.size();i++) cout <<" " <<stcomVec[i]<<endl; // test command
	//cout << "tfValue: "<<tfValue<<endl; // delete it in the main program just to check  // test command
	if(tfValue)
	{
		vector<string>::iterator it;
		 it = find (comVec.begin(), comVec.end(), s); // finding the string in the argv
		 int dis = int(distance(comVec.begin(),it)); // index of the string
		 parsed[dis]=tfValue;
		++it;
		dis = int(distance(comVec.begin(),it)); // looking for the next term
		//test
		//cout << comVec[dis].substr(0,2); // this command is only for test
		if(dis<argCount && comVec[dis].substr(0,1)=="-")
		{
					 msg = "Are you sure about the argument after "+ s + "? If yes, we suggest you to rename the argument: " + comVec[dis]+".\n";
					 error(msg); // is giving warning message: control reaches end of non-void function .

		}

		if(dis<argCount && comVec[dis].substr(0,2)!="--")
		{
			 parsed[dis]=tfValue;  // next value is also true. e.g. --ped test.ped  both strings should exit there.
			 coms = comVec[dis];
			 //cout <<"test"<<endl;
		}
		else
		{
			msg = "An argument after "+ s + " is missing.\n";
			error(msg); // is giving warning message: control reaches end of non-void function .
			//return (msg);
		return;
		}
	}
	else{
	 msg = "An argument after "+s+ " is missing.\n";
	 error(msg);
	}
	 //return (msg);
return;
}


 // printing the options used in command line
  //void Argcv::printOptions(const Bpar* pBpar)
  void Argcv::print_given_options()
  {
	  if(argCount>1)
	  {

		  printLIN("*->Given command options:");
		  for ( int i=1;i<argCount;++i)
		  {
			  string s= original[i];
			  if((i==1)||original[i].substr(0,2)=="--")
				  printLIN("\n\t "+original[i]);
			  else
				  printLIN(" "+ original[i]);

		  }
		  printLIN("\n");
		  printLIN("||====================================================================================||\n");
	  }
	return ;
  }

  // printing the options used in command line
   void Argcv::print_new_options(const Bpar* pBpar)
   {

	   if(pBpar->original.size()>1)
 	  {
 		  int cn =0;
 		  cn=count(pBpar->parsed_comms.begin(), pBpar->parsed_comms.end(), false);
 		  // cout<<"counting false in parsed: "<< cn <<endl; // this command is just to test.
 		  // cout<<original[cn]<<endl;
 		  if(cn>0)
 		  {
 			  printLIN("*->Unnecessary  or miss-typed  commands:");
 			  for ( unsigned int i=0; i<pBpar->original.size();++i) // later argCount
 			  {
 				  //cout << !parsed[i];
 				  if(!pBpar->parsed_comms[i])
 				  {
 				  	  string s= pBpar->original[i];
 				  	  if(pBpar->original[i].substr(0,2)=="--")
 					  printLIN("\t"+pBpar->original[i]);
 				  	  else
 					  printLIN("\t"+pBpar->original[i]);

 				  }
 			  }
 			  printLIN("\n");
 			 msg= "\t**->Frequent possible errors:\n"
 			 "\t -necessary argument(s) or command option(s) is(are) missing,\n"
 			 "\t -unnecessary or unacceptable command option at this particular case,\n"
 			 "\t -having a blank space at unnecessary places,(Note that blank space is accepted only before symbol: \"- -\" )\n"
 			 "\t -having no blank space at necessary places,\n"
 			 "\t -too much or too few letters in the command options,\n"
 			 "\t -having not enough  arguments (some commands need two or more arguments),\n"
 			 "\t -having \"-\", instead of \"- - \". (fcGENE  commands are separated with \"- -\" symbols.)\n";
 			 printLIN(msg);
 			// cout <<"test1\n";
 			 msg = "Please remove (or correct) above unnecessary (or miss-typed) commands. If necessary, please add the missing command option(s).\n";
 			  error(msg);

 			  return;
 		  }
 		  else{
 			  printLIN("\t**->Used command options: ");
 			  	  for (unsigned int i=0;i<pBpar->original.size();++i){
 			  		  string s= pBpar->original[i];
 			  		  if(pBpar->original[i].substr(0,2)=="--")
 			  			  printLIN("\n\t "+pBpar->original[i]);
 			  		  else
 			  			  printLIN(" "+ pBpar->original[i]);

 			  	  }
 			  	  printLIN("\n");
 		  }
 	  }
 	  else if( (compfun("--new-start") )||(compfun("--new-end")))

 	  {

 		  msg ="To use \"--new-start\" and \"--new-end\" commands, there should be at least more than one command option in between them.\n";
 		  error(msg);

 	  }
 	  else
 	  {

 		  string msg1= "*->fcGENE is a format converting tool for GENOTYPE DATA. The data may be given in different formats. Interface of this program is inspired by PLINK software.\n"
 				  	   "\t*->To know about more options, please use the command:\n"
 		  	  	  	   "\t**->\"fcgene --help\" for window system.\n"
 				  	   "\t**->\"./fcgene --help\" for linux systems.\n"
 				  	   "\t**->A pdf documentation file can be found in sourceforge site: \"http://sourceforge.net/projects/fcgene/files/\" .\n";
 		  printLIN(msg1);
 		  printLIN("||====================================================================================||");


 	  }
 	return ;
  }

