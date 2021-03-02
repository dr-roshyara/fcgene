/*

 * matrix_def.h
 *
 *  Created on: 31.12.2013
 *      Author: nroshyar
 */
#ifndef _MATRIX_H__
#define _MATRIX_H__
#include <vector>

	template<class T>
	class Matrix{
	public:
		std::vector<std::vector< T > > data; // This is a data matrix for nrow and ncol.
		Matrix(int nrow, int ncol, T val=0)
		{
			if(nrow>0)
			{
				data.resize(nrow); // resize the vector with nrow
				if(ncol>0)
				{
					for(int i=0;i<nrow;++i )
					{
					data[i].resize(ncol, val); // make ncols with default values;
					}
				}
			}
		}


		int nRows(){return (int )data.size();}
		int nCols(){ return (data.size()==0)?0: (int)data[0].size();}
	};

#endif //_MATRIX_H__


