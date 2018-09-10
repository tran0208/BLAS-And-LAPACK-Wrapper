
#ifdef  _DEBUG
#include <cstdio>
#else
#define NDEBUG
#endif


#ifndef _LapackMatrix_
#define _LapackMatrix_
//
// Prototypes for the only two LAPACK BLAS routines used by this class
//
extern "C" void dgemm_(char* TRANSA,char* TRANSB,long* M, long*N ,long* K,double* ALPHA,
                       double* A,long* LDA,double* B, long* LDB,double* BETA,double* C,long* LDC);

extern "C" void dgemv_(char* TRANS, long* M, long* N, double* alpha, double* Aptr,
                       long* LDA, double* Xptr, long* INCX, double* BETA, double* Yptr, long* INCY);

//
//

#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
/**
 *  A minimal matrix class that facilitates the invocation of LAPACK
 *  routines. This class
 *
 *  (+) stores data by columns :  FORTRAN convention
 *  (+) indexing starts at 0   :  C/C++   convention
 *
 *
 * Calling Fortran versions directly using appropriately modified prototypes.
 *
 * Data type mapping used:
 *
 * C++  char   ==  Fortran character
 * C++  int    ==  Fortran LOGICAL
 * C++  long   ==  Fortran INTEGER
 * C++  double ==  Fortran DOUBLE PRECISION
 *
 * Linking to the Fortran routines using -llapack -lblas
 *
 * April 15, 2016
 *
 * Updated July 27, 2018 (CRA)
 *
 *
 */
/*
#############################################################################
#
# Copyright  2016-2018 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/

namespace SCC
{
class LapackMatrix
{
public:

	LapackMatrix()
	{
	dataPtr        = nullptr;
	externDataFlag = false;
	rows           = 0;
	cols           = 0;
	}

	LapackMatrix(long rows, long cols)
	{
	dataPtr        = nullptr;
	externDataFlag = false;
	initialize(rows,cols);
	}

	// Constructing an instance with externally
	// defined data. Deleting or re-initializing
	// this instance will not delete the external data.
	//
	// double* must point to a single double array
	// of size rows*cols. The data is assumeed to
	// be stored in column major order (Fortran convention)
	//

	LapackMatrix(long rows, long cols, double* dataPtr)
	{
	initialize(rows,cols,dataPtr);
	}

	LapackMatrix(const LapackMatrix& M)
	{
	dataPtr        = nullptr;
	externDataFlag = false;

	if(M.dataPtr == nullptr)
	{rows = 0; cols = 0; return;}

	initialize(M);
	}


	LapackMatrix(LapackMatrix&& V)
    {
      dataPtr        = V.dataPtr;
      externDataFlag = V.externDataFlag;
      rows           = V.rows;
      cols           = V.cols;
      V.dataPtr      = nullptr;
      V.rows         = 0;
      V.cols         = 0;
    }


	~LapackMatrix()
	{
		if((dataPtr != nullptr)&&(not externDataFlag)) delete [] dataPtr;
	}

	//
	// Initialize member functions with local memory
	// allocation. Re-use existing allocation if possible
	//
	void initialize()
	{
		if((dataPtr != nullptr)&&(not externDataFlag))
		{
			delete [] dataPtr;
		}
		dataPtr          = nullptr;
		externDataFlag   = false;
		rows             = 0;
		cols             = 0;
	}

    //
    // This initialize always creates (or uses) a
    // local memory allocation.
    //
	void initialize(long rows, long cols)
	{
	    // Re-use existing allocation if possible

		if((dataPtr != nullptr)&&(not externDataFlag))
		{
			if( (this->rows*this->cols) != rows*cols)
			{
				delete [] dataPtr;
				dataPtr = new double[rows*cols];
			}
		}

		if((dataPtr == nullptr)||(externDataFlag))
		{
		dataPtr        = new double[rows*cols];
		externDataFlag = false;
		}

		this->rows       = rows;
		this->cols       = cols;
		for(long i = 0; i < rows*cols; i++) {dataPtr[i] =0.0;}
	}

    //
    // This initialize always creates (or uses) a
    // local memory allocation.
    //

	void initialize(const LapackMatrix& M)
	{
		if((dataPtr != nullptr)&&(not externDataFlag))
		{
			if( (this->rows*this->cols) != M.rows*M.cols)
			{
				delete [] dataPtr;
				dataPtr = new double[M.rows*M.cols];
			}
		}

		if((dataPtr == nullptr)||(externDataFlag))
		{
		dataPtr        = new double[M.rows*M.cols];
		externDataFlag = false;
		}

		this->rows = M.rows;
		this->cols = M.cols;
		for(long i = 0; i < rows*cols; i++)
		{
			dataPtr[i] = M.dataPtr[i];
		}
	}


	//
	// Initialize with externally defined data.
	//
	// The data is assumed to be stored in column major
	// order (Fortran convention)
	//
	// Deleting or re-initializing this instance
	// will not delete the data.
	//
	// double* must point to a single double array
	// of size rows*cols. There is no error checking
	// on the size of the double* array.
	//

	void initialize(long rows, long cols, double* dataPtr)
	{
		externDataFlag = true;
		this->rows     = rows;
		this->cols     = cols;
		this->dataPtr  = dataPtr;
	}

#ifdef _DEBUG
    double&  operator()(long i1, long i2)
    {
    assert(boundsCheck(i1, 0, rows-1,1));
    assert(boundsCheck(i2, 0, cols-1,2));
    return *(dataPtr +  i1 + i2*rows);
    };

    const double&  operator()(long i1, long i2) const
    {
    assert(boundsCheck(i1, 0, rows-1,1));
    assert(boundsCheck(i2, 0, cols-1,2));
    return *(dataPtr +   i1  + i2*rows);
    };
#else

    /*!
    Returns a reference to the element with index (i1,i2) - indexing
    starting at (0,0).
    */
    inline double&  operator()(long i1, long i2)
    {
    return *(dataPtr +  i1 + i2*rows);
    };

    /*!
    Returns a reference to the element with index (i1,i2) - indexing
    starting at (0,0).
     */
    inline const double&  operator()(long i1, long i2) const
    {
    return *(dataPtr +   i1  + i2*rows);
    };


#endif


    inline void operator=(const LapackMatrix& B)
	{
    	if(dataPtr == nullptr)
    	{
    		rows    = B.rows;
    		cols    = B.cols;
    		dataPtr = new double[rows*cols];
    	}

        assert(sizeCheck(this->rows,B.rows));
    	assert(sizeCheck(this->cols,B.cols));
    	for(long i = 0; i < rows*cols; i++)
    	{
    		dataPtr[i] = B.dataPtr[i];
    	}
	}


    inline void operator+=(const  LapackMatrix& B)
    {
    	assert(sizeCheck(this->rows,B.rows));
    	assert(sizeCheck(this->cols,B.cols));
    	for(long i = 0; i < rows*cols; i++)
    	{
    		dataPtr[i] += B.dataPtr[i];
    	}
    }

    LapackMatrix operator+(const LapackMatrix& B)
    {
    	assert(sizeCheck(this->rows,B.rows));
    	assert(sizeCheck(this->cols,B.cols));

    	LapackMatrix     C(*this);
    	for(long i = 0; i < rows*cols; i++)
    	{
    		C.dataPtr[i] = dataPtr[i] +  B.dataPtr[i];
    	}
        return C;
    }

    LapackMatrix operator-(const LapackMatrix& B)
    {
    	assert(sizeCheck(this->rows,B.rows));
    	assert(sizeCheck(this->cols,B.cols));

    	LapackMatrix     C(*this);
    	for(long i = 0; i < rows*cols; i++)
    	{
    		C.dataPtr[i] = dataPtr[i] -  B.dataPtr[i];
    	}
    	return C;
    }


    inline void operator-=(const  LapackMatrix& D)
    {
      assert(sizeCheck(this->rows,D.rows));
      assert(sizeCheck(this->cols,D.cols));
      for(long i = 0; i < rows*cols; i++)
      {
    		dataPtr[i] -= D.dataPtr[i];
      }
    }

    inline void operator*=(const double alpha)
    {
    for(long i = 0; i < rows*cols; i++)
    {
    		dataPtr[i] *= alpha;
    }
    }

    LapackMatrix operator*(const double alpha)
    {
    LapackMatrix R(*this);
    R *= alpha;
    return R;
    }

    friend LapackMatrix operator*(const double alpha, const LapackMatrix& B)
    {
    LapackMatrix R(B);
    R *= alpha;
    return R;
    }


    void setToValue(double val)
    {
      for(long i = 0; i < rows*cols; i++)
      {
    		dataPtr[i] = val;
      }
    }

    void setToIdentity()
    {
    setToValue(0.0);
    long dimMin = (rows < cols) ? rows : cols;

    for(long k = 0; k < dimMin; k++)
    {
    	this->operator()(k,k) = 1.0;
    }

    }

    void setDiagonal(vector<double>& diag)
    {
    	long kMax = (rows < cols) ? rows : cols;
    	kMax      = ((long)diag.size() < kMax) ? diag.size() : kMax;
    	for(long k = 0; k < kMax; k++)
    	{
    		this->operator()(k,k) = diag[k];
    	}
    }

    //  C := alpha*op( A )*op( B ) + beta*C,

LapackMatrix operator*(const LapackMatrix& B) const
{
    assert(sizeCheck(this->cols,B.rows));

    LapackMatrix C(this->rows,B.cols);

    char TRANSA = 'N';
    char TRANSB = 'N';

    long M       = this->rows;
    long N       = B.cols;
    long K       = this->cols;
    double ALPHA = 1.0;
    double BETA  = 0.0;
    double*Aptr  = dataPtr;
    double*Bptr  = B.dataPtr;
    double*Cptr  = C.dataPtr;
    long LDA     = this->rows;
    long LDB     = B.rows;
    long LDC     = C.rows;

    dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA, Aptr,&LDA,Bptr,&LDB,&BETA,Cptr,&LDC);
    return C;
}


vector<double> operator*(vector<double>& x)
{
	vector<double> y(rows,0.0);

    char TRANS     = 'N';
    double ALPHA   = 1.0;
    double BETA    = 0.0;
    long INCX      = 1;
    long INCY      = 1;

    dgemv_(&TRANS,&rows,&cols,&ALPHA,dataPtr,&rows,&x[0],&INCX,&BETA,&y[0],&INCY);
	return y;
}

vector<double> applyTranspose(vector<double>& x)
{
	vector<double> y(cols,0.0);

    char TRANS     = 'T';
    double ALPHA   = 1.0;
    double BETA    = 0.0;
    long INCX      = 1;
    long INCY      = 1;

    dgemv_(&TRANS,&rows,&cols,&ALPHA,dataPtr,&rows,&x[0],&INCX,&BETA,&y[0],&INCY);
	return y;
}


LapackMatrix transpose() const
{
	LapackMatrix R(cols,rows);
	for(long i = 0; i < rows; i++)
	{
		for(long j = 0; j < cols; j++)
		{
			R(j,i) = this->operator()(i,j);
		}
	}
	return R;
}

double normFrobenius() const
{
	double valSum = 0.0;
    for(long i = 0; i < rows*cols; i++)
    {
    		valSum += dataPtr[i]*dataPtr[i];
    }
    return sqrt(valSum);
}

double elementMaxAbs() const
{
    if(rows*cols == 0) return 0.0;
	double val = abs(dataPtr[0]);
    for(long i = 0; i < rows*cols; i++)
    {
    		val = (val > abs(dataPtr[i])) ? val : abs(dataPtr[i]);
    }
    return val;
}

vector<double> getColumn(long colIndex)
{
	vector<double> r(rows,0.0);
	for(long i = 0; i < rows; i++)
	{
	r[i]=operator()(i,colIndex);
    }
	return r;
}

void insertColumn(vector<double> colVals, long colIndex)
{
	for(long i = 0; i < rows; i++)
	{
	operator()(i,colIndex) = colVals[i];
    }
}

void swapColumns(long colA, long colB)
{
	double val;

    for(long i = 0; i < rows; i++)
	{
    val = operator()(i,colA);
	operator()(i,colA) = operator()(i,colB);
	operator()(i,colB) = val;
    }
}

LapackMatrix getRowSlice(long rowStartIndex, long rowEndIndex)
{
	LapackMatrix M((rowEndIndex-rowStartIndex)+1,this->cols);

	for(long i = rowStartIndex; i <= rowEndIndex; i++)
	{
	for(long j = 0; j < this->cols; j++)
	{
	M(i-rowStartIndex,j) = this->operator()(i,j);
	}}

	return M;
}

LapackMatrix getColSlice(long colStartIndex, long colEndIndex)
{
	LapackMatrix M(this->rows, (colEndIndex-colStartIndex)+1);

	for(long i = 0; i < this->rows; i++)
	{
	for(long j = colStartIndex; j <= colEndIndex; j++)
	{
	M(i,j-colStartIndex) = this->operator()(i,j);
	}}

	return M;
}




/*!  Outputs the matrix values to the screen with the (0,0) element in the upper left corner  */

friend ostream& operator<<(ostream& outStream, const LapackMatrix&  V)
{
        long i; long j;

        for(i = 0;  i < V.rows; i++)
        {
        for(j = 0; j <  V.cols; j++)
        {
          outStream <<   std::scientific << setprecision(3) <<  std::right << setw(10) << V(i,j) << " ";
        }
        outStream << endl;
        }
        return outStream;
}


long getRowDimension() const
{return rows;}

long getColDimension() const
{return cols;}

double* getDataPointer() const
{return dataPtr;}

bool getExternalDataFlag() const
{
	return externDataFlag;
}


//###################################################################
//                      Bounds Checking
//###################################################################
//



#ifdef _DEBUG
        bool boundsCheck(long i, long begin, long end,int coordinate) const
        {
        if((i < begin)||(i  > end))
        {
        cerr << "LapackMatrix index " << coordinate << " out of bounds " << endl;
        cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "] " << endl;
        return false;
        }
        return true;
        }
#else
        bool boundsCheck(long, long, long,int) const {return true;}
#endif

#ifdef _DEBUG
    bool sizeCheck(long size1, long size2)
    {
    if(size1 != size2)
    {
    cerr << "LapackMatrix sizes are incompatible : " << size1 << " != " << size2 << " ";
    return false;
    }
    return true;
    }

    bool sizeCheck(long size1, long size2) const
    {
    if(size1 != size2)
    {
    cerr << "LapackMatrix sizes are incompatible : " << size1 << " != " << size2 << " ";
    return false;
    }
    return true;
    }
#else
    bool sizeCheck(long, long) {return true;}
    bool sizeCheck(long, long) const{return true;}
#endif



    double*       dataPtr;
	long             rows;
	long             cols;

	bool   externDataFlag;

};


}
#endif
