/*
 * 	gebnable.h
 *
 *  Created on: 09.01.2013
 *  Author: Nab Raj Roshyara
 *  Copyright (C) <2012>  <Nab Raj Roshyara> and GenABEL community
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU
 *  General PublicLicense as published by the Free Software Foundation, either version 3 of the License,
 *  or  any later version.This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the  GNU General Public License for more details.You should have received a copy of the GNU General Public License  along with this program.
 *  If not, see <http://www.gnu.org/licenses/>.
 */
// GenABLE
/**
 * 	Some codes of given below are used from GenABEL package of R. If format of Genable package changes, then this codes may not be useful unless the changes in the following codes accordingly.
 * 	R  is commonly referred to as \Base R plus Recommended Packages" and is released in both source code and binary executable forms under the Free Software Foundation's GNU Public License
 *	(hereafter referred to as the GPL)
 *  For more information on GenABEL, please visit http://www.genabel.org/packages/GenABEL
 *  Acknowledgments:
 *	The GenABEL project community
 *	Above all, we appreciate continuous contribution from the GenABEL users and developers [ People ]. The R and the R-Forge
 *	The GenABEL project is in large part made possible by continuous efforts of the R and R-Forge teams [ R | R-Forge ]
 *	Current support:  The Laboratory of Recombination and Segregation Analysis (Prof. Tatiana I. Axenovich and Prof. Pavel M. Borodin),
 *	Institute of Cytology and Genetics SD RAS, Novosibirsk, Russia [ the Laboratory | ICG SD RAS ]
 *  the BioUML/ISB/DotE (Dr. Fedor Kolpakov), Novosibirsk, Russia [ BioUML ]
 *  Former support:
 *   The Genetic Epidemiology Unit (Prof. Cornelia M. van Duijn), Department of Epidemiology, Erasmus MC Rotterdam, #
 *   The Netherlands [ Genetic Epidemiology Unit | Department of Epidemiology | Erasmus MC Rotterdam ]
 *   the Netherlands Bioinformatics Centre (Dr. Marc van Driel, Dr. Morris Swertz) [ NBIC ]
 *	Declaimer
 *	At this page, we acknowledge these who supported the GenABEL project as the whole, e.g. through sponsoring salary of GenABEL-core team members,
 *	 direct financial or labor support of developments, critical to the project as a whole, etc. For acknowledgments of specific
 *	sub-projects and packages from the GenABEL suite, see Packages pages [ Packages ]. It is up to specific package maintainers to provide this information.
 *
 */

#ifndef _GENABLEFORMAT_H__
#define _GENABLEFORMAT_H__
#include "helper.h"
#include "machdat.h"
#include "machped.h"
// for locus information
class CGSNP:public CBSNP{
public:
	//probabel
	static void write_probabel_mlinfo_file(const vector<CBSNP*>& genInfo, const string & outFileName);


};

// for pedigree information
class CGPED:public CBPED
{
 public:

	static void write_genable_raw_file(const vector<CBPED*>& vPed, const vector<CBSNP*>& genInfo, const string & outFileName,bool bcast=true);
	static void write_probabel_mlprob_file(const vector<CBPED*>& vPed, const vector<CBSNP*>& genInfo, const string & outFileName);

	static void convert_snpCoding_into_genable_coding(CBSNP* const pCBSNP,const int nbyte,
									const vector<string>&codeset,
									vector<unsigned short int>&intcoding,
									int *offset,
									vector<unsigned char *>& gtype
									);
	static void write_genable_pheno_file(const vector<CBPED*>& vPed,  const string & outFileName);
	static void write_geable_commands(const string & outFileName);
//probabel
	//static void write_probabel_mlinfo_file(const vector<CBSNP*>& genInfo, const string & outFileName,bool bcast=true);
	//static void write_genable_pheno_file(const vector<CBPED*>& vPed,  const string & outFileName);

};

#endif /*_GENABLEFORMAT_H__ */
