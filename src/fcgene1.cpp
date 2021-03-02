/*
 * fcgene1.cpp
 *
 *  Created on: 31.12.2013
 *      Author: nroshyar
 */
 #include<iostream>
#include<string>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<algorithm>
#include "fcgene.h"
#include "general.h"
#include "smartpca.h"
#include "matrix_def.h"
#include <cmath>

TimeInfo::TimeInfo(){ // constructor of TimeInfo
		setTime(startTime);
		endTime = 0;
		diffTime=0.0;
		a_hours=0;
		a_minutes=0;
		a_seconds=0;
		displayTime();
}
TimeInfo::~TimeInfo(){setTime(endTime);}//destructor of TimeInfo

void TimeInfo::displayTime(){
		if(endTime==0)
		{
		timeInfo=localtime(&startTime);
		//printf("\n The analysis has started at %s",asctime(timeInfo));
		string msg= "\n*->The analysis has started at "+ string(asctime(timeInfo));
		printLIN(msg);
		//cout << msg;
	}
	else{
		timeInfo=localtime(&endTime);

		string msg= "\n*->The analysis has ended at "+string(asctime(timeInfo));
		printLIN(msg);
		diffTime=difftime(endTime,startTime);
		//
		a_minutes=(int)floor(diffTime/60);
		a_hours=(int)floor(a_minutes/60);
		a_minutes=a_minutes%60;
		a_seconds= ceilf(floor(fmod(diffTime,60)));
		//
		char msg2 [200];
		sprintf(msg2,"\n*->Total time taken for the analysis is: %d hours, %d minutes and %.2f seconds.\n",a_hours, a_minutes,a_seconds);
		printLIN(msg2);
		/*
		char sa_hours[10];
		char sa_mins[10];
		char sa_secs[10];
		itoa(a_hours,sa_hours,10);
		itoa(a_minutes,sa_mins,10);
		itoa(a_seconds,sa_secs,10);
		string msg1= "\n Total time taken for the analysis is: "+string(sa_hours)+" hours " +string(sa_mins)+ " minutes and " + string(sa_secs)+" seconds.\n" ;
		printLIN(msg1);
		*/
		 //msg1= printf(msg," horus: %d",a_hours );  //  Total time taken for the analysis is %d ",a_hours); // "hours"; //, %i minutes and %1f  seconds.",a_hours,a_minutes,a_seconds);
		//printf("\n Total time taken for the analysis is: %i hours, %i minutes and %1f  seconds.",a_hours,a_minutes,a_seconds);
	}
}


void show_data_summary_if_given_merge(const vector<CGENERAL*>& pGEN_VEC){

}


void handel_with_general_commands(CGENERAL* const pGENERAL,Bpar* pBpar)
{
	if(pBpar->code_readType!=""&& pBpar->is_general_command)
		pGENERAL->general_commands(pBpar);

}

int handel_with_merging(vector<CGENERAL*> CGEN_VEC, Bpar*pBpar)
{
	const unsigned int _sz =Bpar::cpars_sz;
	if( (_sz>1))
	{
		if(CGEN_VEC.size()!=_sz)
		{
			for(unsigned int j=0;j<CGEN_VEC.size();++j)
				delete CGEN_VEC[j];
			error("There is internal problem probably in saving merge commands\n ");

		}

		CGENERAL* const pCGENERAL = CGEN_VEC[0];
		const CGENERAL* pCGENERAL_new;
		 unsigned int idx_bparam=1;
		 
		while(idx_bparam<_sz)
		{
			Bpar * _pParams=Argcv::bparam_vec[idx_bparam];
			if(_pParams->merge_data)
			{
				//cout <<"I am a merging function! \n ";
				string _tmp="";
				if (idx_bparam==1)
					_tmp = "first";
				else if (idx_bparam==2)
					_tmp ="second";
				else
					_tmp =change_int_into_string(idx_bparam)+"th";
				printLIN("*->Merging "+_tmp+" data into the base data: \n");
				//cout << CGEN_VEC.size()<< " " << idx_bparam << endl;
				pCGENERAL_new =CGEN_VEC[idx_bparam];
				//merge_second_into_first(pCGENERAL, pCGENERAL_new);
				//cout <<"TEST\n";
				merge_second_data_into_first(pCGENERAL, pCGENERAL_new,pBpar->threshold);
			}
			 ++idx_bparam;
		}
	}

return 0;}

void merge_second_data_into_first(CGENERAL* first, const CGENERAL* second, const float&thresh)
{
	static int pid_count =1;
	vector<string>pids_1(first->pedVec.size());
	vector<string>pids_2(second->pedVec.size());
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//assign pids
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	vector<string>npids_1(first->pedVec.size()); // new pids first
	vector<string>npids_2(second->pedVec.size()); //new pids second
	vector<string>cnpids(first->pedVec.size());	 // common npids in both first and second 
	vector<int>pids_index_1_in_2; //index of first pids in second pid list
	vector<int>pids_index_2_notin_1; // index of second pids which are not in first
	vector<bool>handelled_pids_2(second->pedVec.size(),false); // this is just to check if all pids of second are handelled
	// set differences
	unsigned int i =0;
	for(i =0; i<first->pedVec.size();++i)
	{
		if(first->pedVec[i]->indId!="")
			pids_1[i]=first->pedVec[i]->indId;
			else
			{
				pids_1[i]="pid_1_"+change_int_into_string(pid_count);
				++pid_count;
			}
				
			
	}
	for(i =0; i<second->pedVec.size();++i)
	{
		if(second->pedVec[i]->indId!="")
			pids_2[i]=second->pedVec[i]->indId;
		else 
		{
			pids_2[i]="pid_2_"+change_int_into_string(pid_count);
			++pid_count;
		}	
	}		
	//  copy pids as they are  in original
	vector<string>spids_1	 =pids_1;
	vector<string>spids_2	 =pids_2;
	//make stable sort
	stable_sort(spids_1.begin(),spids_1.end());
	stable_sort(spids_2.begin(),spids_2.end());
	 //pids which are in 2nd but not in first
	// first sort pids and then find intersection , difference in between two data.
	// find those pids which are in second but not in first. 
	vector<string>::iterator it;
	unsigned int _tmp =0;
	it=set_difference(spids_2.begin(),spids_2.end(),spids_1.begin(),spids_1.end(),npids_2.begin());
	npids_2.resize(it-npids_2.begin()); //these are the pids which are in 2nd but not in first
	
	if(npids_2.size()>0)
	{
			pids_index_2_notin_1.resize(npids_2.size());
			for(unsigned int i=0; i<npids_2.size(); ++i)
			{
				//cout << npids_2[i]<<" ";//debug 
				it=find(pids_2.begin(),pids_2.end(),npids_2[i]);
				_tmp =(int)(it-pids_2.begin());
			if((second->pedVec[_tmp])->quality)	
					pids_index_2_notin_1[i]=_tmp;
				handelled_pids_2[_tmp]=true;
			}
	}
	// The indices of pids_index_2_notin_1 may be  sorted according as the sorting names of  pids so  re-arange them according as their position in the second data
	stable_sort(pids_index_2_notin_1.begin(),pids_index_2_notin_1.end());
	//cout <<endl;
	// resize npids_1
	npids_2.clear();
	npids_2.resize(0);
	 pids_index_1_in_2.resize(pids_1.size());
	for(unsigned int i=0; i<pids_1.size();++i)
	 {
			it =find(pids_2.begin(),pids_2.end(),pids_1[i]);
			if(it!=pids_2.end())
			{
				handelled_pids_2[((int)(it-pids_2.begin()))]=true;
				pids_index_1_in_2[i]=(int)(it-pids_2.begin());
			}	
			else
				pids_index_1_in_2[i]=-1; // that means not found
	 }
	// debug 
	//for(unsigned int i=0;i<pids_index_1_in_2.size();++i)
	//	cout <<pids_index_1_in_2[i] <<" ";
	//	cout<<endl;
	//Till now we have determined two vectors namely
	 // 1. pids_index_1_in_2 and
	 //2. pids_index_2_notin_1
	 //cout <<"first pids: "<< pids_index_1_in_2.size() << " and second " << pids_index_2_notin_1.size() <<endl;
	
	 // If there are new  pids in the second then add them in the first
	 //but to add pids we also need to add pids in  snps also.
	 //for this :
	 //1.  determine which SNPS of first are given in second. 
	 //2.  which snps of second are not given in first. 
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//------------------for SNPS-----------------------------------------------------------------
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	//cout <<"test"<<endl;
	//rsid
	vector<string>rsid_1(first->genVec.size());
	vector<string>rsid_2(second->genVec.size());
	//snpname
	vector<string>snpid_1(first->genVec.size());
	vector<string>snpid_2(second->genVec.size());
	vector<string>nrsid(second->genVec.size());
	vector<string>nsnpid(second->genVec.size());
	unsigned int j =0;
	for(j =0; j<first->genVec.size();++j)
	{
		rsid_1[j]	 =first->genVec[j]->rsId;
		snpid_1[j]	=first->genVec[j]->snpName;
	}
	for(j =0; j<second->genVec.size();++j)
	{
		rsid_2[j]	 =second->genVec[j]->rsId;
		snpid_2[j]	=second->genVec[j]->snpName;
	}
	//rsid
	vector<string>srsid_1	=rsid_1;
	vector<string>srsid_2	=rsid_2;
	//snpname
	vector<string>ssnpid_1	=snpid_1;
	vector<string>ssnpid_2	=snpid_2;
	//sorting 
	stable_sort(srsid_1.begin(),srsid_1.end());
	stable_sort(srsid_2.begin(),srsid_2.end());
	//snpname
	stable_sort(ssnpid_1.begin(),ssnpid_1.end());
	stable_sort(ssnpid_2.begin(),ssnpid_2.end());
    //2-1 : find index of first-snps in 2nd->snps
    vector<int> snpindex_1_in_2(rsid_1.size());
    vector<int> snpindex_2_notin_1;
	vector<bool>added_allSNPs_2ndta(second->genVec.size(),false);
   	j								=0;
    string _rsid					="";
	string _snpid					="";
	long int _tmp_int				=-1;
	//bool    _is_ok					=false;
	// the following while loop finds if rsids of first are in rsids of 2nd . if found in 2nd, then it gives the index of that rsid in 2nd  and saves 
	//  in vector  snpidx_1_in_2  otherwise it gives an index -1 -
	// in short fill snpidx_1_in_2 with the follwonig loop 
	while(j<rsid_1.size())
    {
    	_snpid		 =snpid_1[j];
    	_rsid		 = rsid_1[j];
    	if(_rsid!="")
    	{
			it =find(rsid_2.begin(),rsid_2.end(),_rsid);
    		if(it!=rsid_2.end())
    			 _tmp_int =(int)(it-rsid_2.begin());
    	}
		else if(_snpid!="")
		{	
			it	 =find(snpid_2.begin(),snpid_2.end(),_snpid);
    		if(it!=snpid_2.end())
				_tmp_int =(int)(it-snpid_2.begin());
    		
		}
		else
		{
			 error(change_int_into_string(j)+"th SNP has no SNPid and rsid. Probably an internal problem in saving. \n");
		}
		snpindex_1_in_2[j]=_tmp_int;
		if(_tmp_int!=-1&&(!added_allSNPs_2ndta[_tmp_int]))
			added_allSNPs_2ndta[_tmp_int]=true;
		 _tmp_int =-1;
		++j; // increase an increment.	
	}	
	//Now find those SNPs which are in 2nd data but not in the first data.
	// rsid given in 2nd but not in first. 
	it =set_difference(srsid_2.begin(),srsid_2.end(),srsid_1.begin(),srsid_1.end(),nrsid.begin());
	nrsid.resize(it-nrsid.begin());
	// snpname
	it=set_difference(ssnpid_2.begin(),ssnpid_2.end(),ssnpid_1.begin(),ssnpid_1.end(),nsnpid.begin());
	nsnpid.resize(it-nsnpid.begin());
	//cout <<"first pids size():  " <<npids_1.size()<<", npids_2.size(): "<< npids_2.size()<<", "<<"nrsid.size(): "<< nrsid.size() << " nsnpid.size():  "<< nsnpid.size()<<"\n"; // << " nbp.size(): "<< nbp.size()<< " \n";
	j=0;
	while(j<nrsid.size())
	{
		// This while loop finds those snps which are in 2nd but  not in first  and saves the indices in vector snpindex_2_notin_1 ; 
		_rsid	=nrsid[j];
		if(_rsid!="")
		{
			it =find(rsid_2.begin(),rsid_2.end(),_rsid);
			if(it!=rsid_2.end())
			{
				_tmp_int =(int)(it-rsid_2.begin());
				//cout << _tmp_int << ", bool: " << added_allSNPs_2ndta[_tmp_int]<< endl;
				//	cout << added_allSNPs_2ndta.size() << ","<<_tmp_int <<endl;
				if(added_allSNPs_2ndta[_tmp_int])
					error(change_int_into_string(_tmp_int)+"th SNP has already found in data one and it can not be assumed again as SNP which is not in first data.\n");
				added_allSNPs_2ndta[_tmp_int]=true; 	
				if((unsigned int)_tmp_int<snpid_2.size())
				{
					it=find(nsnpid.begin(),nsnpid.end(),snpid_2[_tmp_int]);
					if(it!=nsnpid.end())
						nsnpid.erase(it);
				}	
					
			}
		}
		if((second->genVec[_tmp_int])->quality)
			snpindex_2_notin_1.push_back(_tmp_int);	
		++j;//make increment 
	}
	j=0;
	nrsid.clear();	
	nrsid.resize(0);
	vector<bool>::iterator it_bool =find(added_allSNPs_2ndta.begin(),added_allSNPs_2ndta.end(),false);
	// cout << "(nsnpid.size()>0): "<<(nsnpid.size()>0) << ", "<< (it_bool!=added_allSNPs_2ndta.end() )<<endl;
	if(nsnpid.size()>0 ||it_bool!=added_allSNPs_2ndta.end() )
	{
		// if rsid could not resolve the problem of extra snps		
		error("There is a problem in merging SNPs. Please contact the software provider.\n");
	
	}
	bool _tfvec=false;
	 // Now the following loop takes again  those snps which are in first and 2nd  and  checks if they are really the same  
	 // it checks if they have same genotype format 
	 // it checks if they have same allele order 
	 // it checks if they have same  snp information 
	 for(j=0;j<snpindex_1_in_2.size();++j)
	{
		_tmp_int =snpindex_1_in_2[j];
		if(_tmp_int!=-1)
		{
			CBSNP * fpCBSNP =first->genVec[j];
			CBSNP* spCBSNP =second->genVec[_tmp_int];
			
			// check if they are same snps 
			_tfvec=CGENERAL::are_they_same_snp(fpCBSNP,spCBSNP);
			//cout << boolalpha<< _tfvec <<endl; 
			if(_tfvec) // if they are same snps 
			{
				//cout << boolalpha << _tfvec << " \n" ;
				//_tfvec =false;
				//checks if they have same allele order if not change the allele order 
				CGENERAL::check_make_same_alleleOrder(fpCBSNP,spCBSNP);
				//check and update snp information if necessary 
				CGENERAL::check_n_update_snpInfo(fpCBSNP,spCBSNP);
				//cout << _tfvec << " \n";
				// it checks if they have same format . if not make the same format 
				CGENERAL::check_n_make_same_genoFormat(fpCBSNP, spCBSNP, thresh);
			}	
		}
		
	}
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	// end for SNPs--------------------------------------------------------------------------------
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	// this follwoing loop adds the genotypes in first data for  those  individudals which are not in first and for those snps which are in  first   or and second  	
		//check how many pids are there 
	// cout <<"pids_index_2_notin_1: \n";
	// for(unsigned int i=0; i<pids_index_2_notin_1.size()>0;++i) //debug 
	//cout << pids_index_2_notin_1[i]<< " " ; //debug 
 	//cout << endl; // deubg 
		
	if(pids_index_2_notin_1.size()>0)
	{
		 // create new CBPED
		 // first->pedVec.resize(pedVec.size()+index_2_notin_1.size());
		 CBSNP * _npCBSNP=first->genVec[0];
		 int _tmp_int =-1;
		 unsigned int _tmp_pid_indx =0;
		for(unsigned int i=0;  i<pids_index_2_notin_1.size(); ++i)
		{
			 _tmp_pid_indx =pids_index_2_notin_1[i];
			//first add new person to pedVec of first data.
			 CBPED * npCBPED =new CBPED;
			 *npCBPED =*second->pedVec[_tmp_pid_indx];
			 // update sex info 
			 bool update_sex =true; 
			 bool update_pheno=true;
			 CBPED::add_sex_n_phenotype_info(npCBPED, update_sex,update_pheno);
			 first->pedVec.push_back(npCBPED);
			 //second , add all genotypes to the genVec of first data.
			 // for this we do the following.
			// cout <<"test1"<<endl; 
			 for(unsigned int _j=0;_j<first->genVec.size();++_j)
			{
				 _npCBSNP =first->genVec[_j];
				 _tmp_int =snpindex_1_in_2[_j];  //snpindex_1_in_2 has eitehr index -1 or some other value . if -1 then this snp is not in 2nd 
				if(_tmp_int!=-1)
				{
				 	 const CBSNP* _second_psnp =second->genVec[_tmp_int];
					//cout << _second_psnp->rsId << ", first: " << _npCBSNP->rsId <<" \n" ;
					//cout <<_npCBSNP->geno1.size()<<endl;
					//cout <<_npCBSNP->pgeno1.size()<<endl;
					//cout <<_npCBSNP->pgeno2.size()<<endl;
					//cout <<_npCBSNP->pgeno3.size()<<endl;
					//cout << _tmp_pid_indx<<endl; 
					//cout << _second_psnp->pgeno1.size()<<endl;
					if(_npCBSNP->geno1.size()>0 && _tmp_pid_indx<_npCBSNP->geno1.size()) // this checks if there are hard calls . if given in first then it must be in second also . we have already checked  and updated this 
					{
						_npCBSNP->geno1.push_back(_second_psnp->geno1[_tmp_pid_indx]);
						_npCBSNP->geno2.push_back(_second_psnp->geno2[_tmp_pid_indx]);
						
					}
					else if(_npCBSNP->geno1.size()<_tmp_pid_indx){
						error("Problem in merging: "+change_int_into_string(_tmp_pid_indx)+"th SNP of second data\n");
					}	
					//		cout <<"testing\n";
					if(_npCBSNP->pgeno1.size()>0) // this checks if there are probabilities of genotypes in first .  if given in first then it must be in second also . we have already checked  and updated this 
					{
						_npCBSNP->pgeno1.push_back(_second_psnp->pgeno1[_tmp_pid_indx]);
						_npCBSNP->pgeno2.push_back(_second_psnp->pgeno2[_tmp_pid_indx]);
						_npCBSNP->pgeno3.push_back(_second_psnp->pgeno3[_tmp_pid_indx]);
					}	
					_npCBSNP->aOrder.push_back(_second_psnp->aOrder[_tmp_pid_indx]);
					
				}	
				else // this is the case where snps   of first not given in second.
				{
						//set missings 
					if(_npCBSNP->geno1.size()>0) // this checks if there are hard calls . if given in first then it must be in second also . we have already checked  and updated this 
					{
						_npCBSNP->geno1.push_back(true);
						_npCBSNP->geno2.push_back(false);
						_npCBSNP->aOrder.push_back(true);
					}	
					if(_npCBSNP->pgeno1.size()>0) // this checks if there are probabilities of genotypes in first .  if given in first then it must be in second also . we have already checked  and updated this 
					{
						_npCBSNP->pgeno1.push_back(0);
						_npCBSNP->pgeno2.push_back(0);
						_npCBSNP->pgeno3.push_back(0);
						_npCBSNP->aOrder.push_back(true);
					}	
				}
				 
			}
		}
		
		
		
		}
		
	//cout << "first->pedVec.size():"<<first->pedVec.size() << endl; 
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	////now work in those SNPS which are in 2nd data but not given in first data 

		//cout <<"test2"<<endl;

	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	if(snpindex_2_notin_1.size()>0)
	{
		//CBPED * _pCBPED =first->pedVec[0];
		int _tmp_int =-1;
		int _tmp_snp_indx =0;
		for(unsigned int j=0; j<snpindex_2_notin_1.size();++j)
		{
			 _tmp_snp_indx =snpindex_2_notin_1[j];
			//first add new person to pedVec of first data.
			CBSNP * npCBSNP 	=new CBSNP;
			CBSNP * _tmp_pCBSNP =second->genVec[_tmp_snp_indx];
			*npCBSNP 			=*_tmp_pCBSNP;
			//first make the same format of the SNPs  which is given in first 
			bool _is_first_geno =((first->genVec[0])->geno1.size()>0);		
			bool _is_first_pgeno =((first->genVec[0])->pgeno1.size()>0);		
			bool _tmp_cond1 = _is_first_geno && (_tmp_pCBSNP->pgeno1.size()>0) &&(_tmp_pCBSNP->geno1.size()==0);
			bool _tmp_cond2 =_is_first_pgeno && (_tmp_pCBSNP->geno1.size())	&& (_tmp_pCBSNP->pgeno1.size()==0);
			if(_tmp_cond1)
				CBSNP::convert_genoprob_into_hardcalls(_tmp_pCBSNP, thresh);
			if(_tmp_cond2)
				CBSNP::convert_hardcalls_into_genoprob(_tmp_pCBSNP);
			 // fix total individuals 
			 unsigned const int	 _tot_indiv    =first->pedVec.size();
			 unsigned const int  _first_indivs =pids_index_1_in_2.size();
			// cout <<"first pids: "<< pids_index_1_in_2.size() << " and second " << pids_index_2_notin_1.size() <<endl;
			// If there are new  pids in the second then add them in the first
				if(_tot_indiv!=( pids_index_1_in_2.size() +pids_index_2_notin_1.size()))
				{ 
					delete npCBSNP;
					error("problem");
				}
				else 
				{
					//crate for order 
					vector<bool> naOrder(_tot_indiv);
					bool update_order_once=true;
					if(_is_first_geno)
					{
						// create  geno vectors here. 	
						vector<bool> ngeno1(_tot_indiv)	;
						vector<bool> ngeno2(_tot_indiv)	;
						for(unsigned int i=0;i<pids_index_1_in_2.size(); ++i)
						{
							_tmp_int  =pids_index_1_in_2[i];
							if(_tmp_int!=-1)
							{
								ngeno1[i]  = _tmp_pCBSNP->geno1[_tmp_int];
								ngeno2[i]  = _tmp_pCBSNP->geno2[_tmp_int];
								naOrder[i] = _tmp_pCBSNP->aOrder[_tmp_int];
							}
							else
							{
								// individual not found in second . set missing 
								ngeno1[i]  = true;
								ngeno2[i]  = false;
								naOrder[i] =true;
							}	
						}
						// now update pids which are not in first but only in second
						
						for(unsigned int i=0;i<pids_index_2_notin_1.size(); ++i)
						{
							_tmp_int  =pids_index_2_notin_1[i];
							ngeno1[_first_indivs+i]   = _tmp_pCBSNP->geno1[_tmp_int];
							ngeno2[_first_indivs+i]   = _tmp_pCBSNP->geno2[_tmp_int];
							naOrder[_first_indivs+i]  = _tmp_pCBSNP->aOrder[_tmp_int];
						}
						//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
						update_order_once=false;
						npCBSNP->geno1 =ngeno1;
						npCBSNP->geno2 =ngeno2;
						
					}	
					if(_is_first_pgeno )
					{
						//crates pgeno vector here	
						vector<double> npgeno1(_tot_indiv)	;
						vector<double> npgeno2(_tot_indiv)	;
						vector<double> npgeno3(_tot_indiv)	;
						for(unsigned int i=0;i<pids_index_1_in_2.size(); ++i)
						{
							_tmp_int  =pids_index_1_in_2[i];
							if(_tmp_int!=-1)
							{
								npgeno1[i]  = _tmp_pCBSNP->pgeno1[_tmp_int];
								npgeno2[i]  = _tmp_pCBSNP->pgeno2[_tmp_int];
								npgeno3[i]  = _tmp_pCBSNP->pgeno3[_tmp_int];
								if(update_order_once)
								 naOrder[i] = _tmp_pCBSNP->aOrder[_tmp_int];
							}
							else
							{
								// individual not found in second . set missing 
								npgeno1[i]  = 0.0;
								npgeno2[i]  = 0.0;
								npgeno3[i]  = 0.0;
								if(update_order_once)
									naOrder[i] =true;
							}	
						}
						//update the individuals which are only in second 
						for(unsigned int i=0;i<pids_index_2_notin_1.size(); ++i)
						{
							_tmp_int  =pids_index_2_notin_1[i];
							npgeno1[_first_indivs+i]  = _tmp_pCBSNP->pgeno1[_tmp_int];
							npgeno2[_first_indivs+i]  = _tmp_pCBSNP->pgeno2[_tmp_int];
							npgeno3[_first_indivs+i]  = _tmp_pCBSNP->pgeno3[_tmp_int];
							if(update_order_once)
							 naOrder[_first_indivs+i] = _tmp_pCBSNP->aOrder[_tmp_int];
						}
						npCBSNP->pgeno1 =npgeno1;
						npCBSNP->pgeno2 =npgeno2;
						npCBSNP->pgeno3 =npgeno3;
						
					}	
				
					npCBSNP->aOrder =naOrder;	
				}
					//create aorder here. 
					
				first->genVec.push_back(npCBSNP);	
		}
				
			// update sex info 
			// CBPED::add_sex_n_phenotype_info(npCBPED, update_sex,update_pheno);
		
			//second , add all genotypes to the genVec of first data.
			// for this we do the following.
	}
	
	//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
return; }
 
// splitting
void handel_with_splitting(const CGENERAL* const pCGENERAL, const Bpar* const pBpar)
{
	  // splitting option
	       //snp split
	if(pBpar->ssplitt|| pBpar->isplitt||pBpar->bpsplitt)
	{
		printLIN("*->Splitting the genotype data:\n ");
		//first check the size of split vec
			const int _sSize 				=pBpar->v_ssplitt.size();
			const int _iSize 				=pBpar->v_isplitt.size();
			unsigned const int _genSize 	=pCGENERAL->genVec.size();	//snp upper index
			unsigned const int _pedSize 	=pCGENERAL->pedVec.size();	//snp upper index
			vector<CGENERAL*>	 SPLITT_CGEN;
			vector<Bpar*> 		 SPLITT_Bpar;
			unsigned int 		_sl_indx	=0;							//snp lower index
			unsigned int 		_su_indx	=_genSize;	//snp upper index
			unsigned int 		_il_indx	=0;							//indiv lower index
			unsigned int 		_iu_indx	=_pedSize;	//indiv upper index
			bool				_check		=false;
			//this loop creates the new parameter class and new General class
			unsigned const int _nSize =pBpar->v_split_suffix.size(); // this contains ending name of splitted data
			//Therefore we must create the same number of new cgeneral classes and basic parameter classes which represent
			//each of the new split classes.
			vector<int> _cur_snp_indices;
			for(unsigned int j=0;j<_nSize; ++j)	//everything is done inside now
			{
				CGENERAL* _tpCGENERAL =new CGENERAL; // creating new Class
				_cur_snp_indices.clear();
				_cur_snp_indices.resize(0);
				// do something for  _tpCGENERAL
				//check if the given upper and lower limits with ssplit, bpsplit and isplit are correct and
				//meaningful. Meaningful means they should lie in between the snp and indiviudal indices
				//first for snpwise split
				//----------------------------------------------------------------------//
				if(pBpar->v_ssplitt.size()>0)
				{
					if(pBpar->v_ssplitt.size()==2)//this is a case where only one split task is given.
					{
						if(pBpar->ssplitt) // if snp splitting is given as indices of snps
						{
							bool _b1 =pBpar->v_ssplitt[0]<1 || pBpar->v_ssplitt[0]>_genSize;
							_b1 	 =_b1|| (pBpar->v_ssplitt[1]<1 || pBpar->v_ssplitt[1]>_genSize);
							if(_b1){
								string _msg ="indices given with \"--ssplitt\" must lie in between 1 and total number of SNPs: "+change_int_into_string(_genSize)+ "\n";
								error(_msg);
							}
						}
						_sl_indx =pBpar->v_ssplitt[0];
						_su_indx =pBpar->v_ssplitt[1];

					}
					else{ //this is a case where multiple split commands are given.

							if(pBpar->ssplitt) // if ssplit , then check if the indices lie in between 1 and size of genoVec
							{
								bool _b1 =pBpar->v_ssplitt[(2*j)]<1 || pBpar->v_ssplitt[(2*j)]>_genSize;
								_b1 	 =_b1|| (pBpar->v_ssplitt[(2*j+1)]<1 || pBpar->v_ssplitt[(2*j+1)]>_genSize);
								if(_b1){
								string _msg ="indices given with \"--ssplitt\" must lie in between 1 and total number of SNPs: "+change_int_into_string(_genSize)+ "\n";
								error(_msg);
								}
							}

							_sl_indx =pBpar->v_ssplitt[(2*j)];
							_su_indx =pBpar->v_ssplitt[(2*j+1)];

					}
					//now if _sl_indx and _su_index are used for bp then we should find a vector containing bps that lie in between
					//otherwise we should find a vector which contains the indices that lie in between
					//If bp are given
					if(pBpar->bpsplitt)
					{
						unsigned int _x_ind =0;
						bool _bool_tmp	=false;
						for(unsigned int j=0; j<pCGENERAL->genVec.size();++j)
						{
							_x_ind =pCGENERAL->genVec[j]->bp;
							//cout <<j<< ", bp: "<< pCGENERAL->genVec[j]->bp<<", " <<_x_ind<< ",  "<<_sl_indx<<", "<<_su_indx<<", "<<boolalpha <<(_x_ind>=_sl_indx&& _x_ind<=_su_indx)<<endl;
							if(_x_ind>=_sl_indx&& _x_ind<=_su_indx)
							{
								_cur_snp_indices.push_back(j);
								 _bool_tmp=true;
							}

						}
						if(!_bool_tmp)
						{
							//this is the case where no bp found in the range
							string _msg =" There is no base pair position in the range \""+change_int_into_string(_sl_indx)+"-"+change_int_into_string(_su_indx)+"\", you have given with command \"--bpsplit\".\n"
									"Please give the correct range of base pair position for splitting genotype data.\n";
							error(_msg);
						}
					}else{
						//then its  snp indx
						_cur_snp_indices.resize((_su_indx-_sl_indx +1));
						for(unsigned int j=0;j<_cur_snp_indices.size();++j)
							_cur_snp_indices[j]= (_sl_indx -1 +j);
					}

				}
				//----------------------------------------------------------------------//

				//indices of iindiviudals
				if(pBpar->v_isplitt.size()>0)
				{
					if(pBpar->v_isplitt.size()==2)
					{
						bool _b1 =pBpar->v_isplitt[0]<1 || pBpar->v_isplitt[0]>_pedSize;
							_b1 =_b1|| (pBpar->v_isplitt[1]<1 || pBpar->v_isplitt[1]>_pedSize);
							//cout <<boolalpha << _b1 <<endl;
						if(_b1){
							string _msg ="indices given with --issplit must lie in between 1 and total number of individuals: "+change_int_into_string(_pedSize)+ "\n";
							error(_msg);
						}else{
							_il_indx =(pBpar->v_isplitt[0]-1);
							_iu_indx =pBpar->v_isplitt[1];
						}
					}
					else{
							bool _b1 =pBpar->v_isplitt[(2*j)]<1 || pBpar->v_isplitt[(2*j)]>_pedSize;
								_b1 =_b1|| (pBpar->v_isplitt[(2*j)+1]<1 || pBpar->v_isplitt[(2*j)+1]>_pedSize);
							if(_b1)
							{
								string _msg ="indices given with --issplit must lie in between 1 and total number of individuals:"+change_int_into_string(_pedSize)+"\n";
								error(_msg);
							}else
							{
								_il_indx =pBpar->v_isplitt[(2*j)];
								_iu_indx =pBpar->v_isplitt[(2*j+1)];
							}
					}
				}
				// In above . we do not need to check more coz everything has been checked already.
				//---------------------------------------------------------------------------------//
				//cout <<_sl_indx<<", "<<_su_indx<<", "<<_il_indx<<", "<<_iu_indx<<endl;
				// create ped Vec
				for(unsigned int _ii=(_il_indx);_ii<_iu_indx;++_ii)
				{
					// take only those pCBSNP whose index are given
					CBPED* npCBPED =new CBPED;
					 *npCBPED = *(pCGENERAL->pedVec[_ii]);
					 _tpCGENERAL->pedVec.push_back(npCBPED);
				}
				unsigned int _tmp_indx =0;
				for(unsigned int _si=0;_si<_cur_snp_indices.size();++_si) //
				{
					_tmp_indx =_cur_snp_indices[_si];

					// her we can update snps but
					 CBSNP* _p_lhs  =new CBSNP;
					 *_p_lhs		=*(pCGENERAL->genVec[_tmp_indx]);
					 _p_lhs->geno1.clear();
					 _p_lhs->geno2.clear();
					 _p_lhs->aOrder.clear();
					 //
					 _p_lhs->pgeno1.clear();
					 _p_lhs->pgeno2.clear();
					 _p_lhs->pgeno3.clear();
					 //
					 _p_lhs->geno_0123.clear();
					 _p_lhs->geno_dose.clear();
					 // All vectors related to saving indiviudals
					 // information should be cleared and resized to zero first
					 // later only the necessary individuals are assigned.
					 //its important to check if this is possible ;
					for(unsigned int _ii=(_il_indx);_ii<_iu_indx;++_ii)
					{


						//selected individual loop
						// here indiviudal can be updated
						//1. add geno1 and geno2
						_p_lhs->geno1.push_back(pCGENERAL->genVec[_tmp_indx]->geno1[_ii] );
						_p_lhs->geno2.push_back(pCGENERAL->genVec[_tmp_indx]->geno2[_ii] );
						_p_lhs->aOrder.push_back(pCGENERAL->genVec[_tmp_indx]->aOrder[_ii] );
						//copy   genotypeprobabilites

						 if(pCGENERAL->given_pgeno)
						{
							_p_lhs->pgeno1.push_back(pCGENERAL->genVec[_tmp_indx]->pgeno1[_ii] );
							_p_lhs->pgeno2.push_back(pCGENERAL->genVec[_tmp_indx]->pgeno2[_ii] );
							_p_lhs->pgeno3.push_back(pCGENERAL->genVec[_tmp_indx]->pgeno3[_ii] );

						}

						 // copy 012
						if(!pCGENERAL->genVec[_tmp_indx]->change_0123){
							_p_lhs->geno_0123.push_back(pCGENERAL->genVec[_tmp_indx]->geno_0123[_ii] );

						}
						//copy geno dose
						if(!pCGENERAL->genVec[_tmp_indx]->change_dose)
						{
							_p_lhs->geno_dose.push_back(pCGENERAL->genVec[_tmp_indx]->geno_dose[_ii] );

						}

						// cout <<"["<<_tmp_indx<<","<<_ii<<"]"<<"\t";
						 // pedVec can be
					}
					//set maf, hwe etc to zero because new snp has different number of indiviudals
					_p_lhs->given_both_crate	=false;
					_p_lhs->given_indiv_crate	=false;
					_p_lhs->given_snp_crate		=false;
					_p_lhs->given_snp_hwe		=false;
					//check if the no of individuals and SNP vector related to individuals (e.g. geno1 geno2) have same 
					//number of elements 
					//	check =_p_lhs->geno1.size()==
						
					_tpCGENERAL->genVec.push_back(_p_lhs);
					//cout << endl;
					//cout<<"pCGENERAL->genVec[_tmp_indx]->geno1.size(): "<<pCGENERAL->genVec[_tmp_indx]->geno1.size()<<endl;

				}

				Bpar* pSplittBpar =new Bpar;
				*pSplittBpar =*pBpar;
				if(!pSplittBpar->change_format)
				{
					pSplittBpar->change_format =true;
					if(pSplittBpar->is_code_plink)
					{
						if(pSplittBpar->code_readType=="plink")
							pSplittBpar->oFormat	="plink";
						if(pSplittBpar->code_readType=="bplink")
							pSplittBpar->oFormat	="plink-bed";
						else if(pSplittBpar->code_readType=="plink-dosage")
							pSplittBpar->oFormat	="plink-dosage";
						else if (pSplittBpar->code_readType=="plink-rawA")	
							pSplittBpar->oFormat	="plink-recodeA";
						else if(pSplittBpar->code_readType=="plink-rawAD")	
							pSplittBpar->oFormat		="plink-recodeAD" ;
						
					}
					if(pSplittBpar->is_code_mach)
					{					
						pSplittBpar->oFormat	="mach";
					}
					else if(pSplittBpar->is_code_minimac)
					{		
						pSplittBpar->oFormat	="minimac";
					}
					else if(pSplittBpar->is_code_impute)
					{					
						pSplittBpar->oFormat ="impute";
					}
					else if(pSplittBpar->is_code_snptest)
					{				
						pSplittBpar->oFormat ="snptest";
					}
					else if(pSplittBpar->is_code_beagle)
					{	
						pSplittBpar->oFormat ="beagle";
					}
					else if(pSplittBpar->is_code_bimbam)
					{					
					}					
					else if(pSplittBpar->is_code_rformat)
					{
						if(pSplittBpar->code_readType=="r")
							pSplittBpar->oFormat="r";
						
					}
				
						
						
				}
				pSplittBpar->output_fileName+=pBpar->v_split_suffix[j];
				printLIN("||-----------------------------------------------------------------||\n");
				printLIN("*->Working on "+ change_int_into_string(j+1)+"th Splitted data:\n");
				handel_with_general_commands(_tpCGENERAL,pSplittBpar);
				//printLIN("||-----------------------------------------------------------------||\n");
	
					
				// check if all CBPED and CBSNP are are included into new or not;
				//cout<<"	_tpCGENERAL->pedVec.size(): "<<	_tpCGENERAL->pedVec.size()<<endl;
				//cout<<"	_tpCGENERAL->genVec.size(): "<<	_tpCGENERAL->genVec.size()<<endl;
	
				SPLITT_CGEN.push_back(_tpCGENERAL);
				//pSplittBpar->output_fileName+=pBpar->v_split_suffix[j];
				//cout << pSplittBpar->output_fileName<<endl;
				SPLITT_Bpar.push_back(pSplittBpar);
			}

			/*


			for(unsigned int i=0;i<pBpar->v_ssplitt.size();++i)
			{
				cout <<"splitting "<< pBpar->v_ssplitt[i]<<" ";
	       	}
	       	cout <<endl;
	       	*/
	}
}

