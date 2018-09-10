#include "LapackMatrix.h"
#include "Lapack.h"
//
// SCC::LapackMatrixRoutines
//
// A miscellaneous collection of classes whose functionality is
// based upon LAPACK routines. These routine are meant
// to be used with instances of LapackMatrix.
//
// Author          : Chris Anderson
// Version         : July, 27, 2018
/*
#############################################################################
#
# Copyright  2015-2018 Chris Anderson
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


#include <iostream>
#include <vector>
#include <cassert>
#include <complex>
#include <cstdlib>
#include <limits>
#include <cstring>
#include <stdexcept>
using namespace std;

#ifndef  _LAPACKMATRIXROUTINES_
#define  _LAPACKMATRIXROUTINES_

namespace SCC
{

// C++  int    ==  Fortran LOGICAL
// C++  long   ==  Fortran INTEGER
// C++  double ==  Fortran DOUBLE PRECISION

/*
	 DGELSY computes the minimum-norm solution to a real linear least
     squares problem:
         minimize || A * X - B ||
     using a complete orthogonal factorization of A.  A is an M-by-N
     matrix which may be rank-deficient.

     Several right hand side vectors b and solution vectors x can be
     handled in a single call; they are stored as the columns of the
     M-by-NRHS right hand side matrix B and the N-by-NRHS solution
     matrix X.

     The routine first computes a QR factorization with column pivoting:
         A * P = Q * [ R11 R12 ]
                     [  0  R22 ]
     with R11 defined as the largest leading submatrix whose estimated
     condition number is less than 1/RCOND.  The order of R11, RANK,
     is the effective rank of A.

     Then, R22 is considered to be negligible, and R12 is annihilated
     by orthogonal transformations from the right, arriving at the
     complete orthogonal factorization:
        A * P = Q * [ T11 0 ] * Z
                    [  0  0 ]
     The minimum-norm solution is then
        X = P * Z**T [ inv(T11)*Q1**T*B ]
                     [        0         ]
     where Q1 consists of the first RANK columns of Q.

     This routine is basically identical to the original xGELSX except
     three differences:
       o The call to the subroutine xGEQPF has been substituted by the
         the call to the subroutine xGEQP3. This subroutine is a Blas-3
         version of the QR factorization with column pivoting.
       o Matrix B (the right hand side) is updated with Blas-3.
       o The permutation of matrix B (the right hand side) is faster and
         more simple.
 */

class DGELSY
{
public:

	DGELSY()
	{
		initialize();
	}

    void initialize()
	{
	this->A.initialize();
	this->X.clear();
    this->WORK.clear();
    this->JPVT.clear();
    RCOND = 10.0*numLimits.epsilon();
    RANK  = 0;
	}


	void initialize(const DGELSY& dgelsy)
	{
	this->A    = dgelsy.A;
	this->X    = dgelsy.X;
    this->JPVT = dgelsy.JPVT;
    RCOND      = dgelsy.RCOND;
    RANK       = dgelsy.RANK;
	}

    long getRank()
    {
    return RANK;
    }


	vector<double> qrSolve(const vector<double> B, const LapackMatrix& A, double rcondCutoff = -1.0)
	{
	    if(rcondCutoff < 0) {RCOND = 10.0*numLimits.epsilon();}
	    else                {RCOND = rcondCutoff;}

		// Copy A, since A is overwritten or destroyed by
	    // the routine.

		if(not equalMatrixDimensions(this->A,A))
		{
		this->A.initialize(A);
		}
		else
		{
		this->A = A;
		}


		long M  = this->A.getRowDimension();
		long N  = this->A.getColDimension();

		// Capture B and increase size if N > M

        X = B;
        if(M < N) {X.resize(N,0.0);}

		JPVT.resize(N);

		long NRHS = 1;
		long LDA  = M;
		long LDB  = X.size();


	   // Query to obtain optimal work array size

		long LWORK = -1;
		INFO       =  0;

		double WORKtmp;

		dgelsy_(&M, &N, &NRHS, this->A.getDataPointer(), &LDA,&X[0],&LDB,
		&JPVT[0], &RCOND, &RANK, &WORKtmp, &LWORK,&INFO);


	    LWORK = (long)(WORKtmp + 100);
		WORK.resize(LWORK);


		// Second call to create qr solution

	    INFO       = 0;

		dgelsy_(&M, &N, &NRHS, this->A.getDataPointer(), &LDA,&X[0],&LDB,
		&JPVT[0], &RCOND, &RANK, &WORK[0], &LWORK,&INFO);

		if(INFO != 0)
        {
        cerr << "DGELSY  Failed : INFO = " << INFO  << endl;
        exit(1);
        }

        // Set X to be the right dimension
        X.resize(N);

        return X;
	}


    //
    // A convenience solve interface. It is assumed that Bptr points to contiguous
    // data of size of the number of rows of A. No bounds checking is performed.
    //

	vector<double> qrSolve(double* Bptr, const LapackMatrix& Ain, double rcondCutoff = -1.0)
	{
		long M     = Ain.getRowDimension();
		vector<double>                 B(M);
		std::memcpy(B.data(),Bptr,M*sizeof(double));
		return qrSolve(B,Ain,rcondCutoff);
	}

	bool equalMatrixDimensions(const LapackMatrix& A, const LapackMatrix& B)
	{
		if((A.rows != B.rows)||(A.cols != B.cols)) return false;
		return true;
	}


	LapackMatrix           A;
	vector<double>         X;
    vector<long>        JPVT;
    double             RCOND;
    long                RANK;

	long                INFO;
	vector<double>      WORK;

	numeric_limits<double>   numLimits;
};


/*
*> DGESVD computes the singular value decomposition (SVD) of a real
*> M-by-N matrix A, optionally computing the left and/or right singular
*> vectors. The SVD is written
*>
*>      A = U * SIGMA * transpose(V)
*>
*> where SIGMA is an M-by-N matrix which is zero except for its
*> min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
*> V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
*> are the singular values of A; they are real and non-negative, and
*> are returned in descending order.  The first min(m,n) columns of
*> U and V are the left and right singular vectors of A.
*>
*> Note that the routine returns V**T, not V.
*>
*> SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
*> WORK, LWORK, INFO )
*>
*/
class DGESVD
{
public:

	DGESVD()
	{
		initialize();
	}


	DGESVD(const DGESVD& dgesvd)
	{
		initialize(dgesvd);
	}

	void initialize()
	{
	this->A.initialize();
    this->U.initialize();

    this->singularValues.clear();
    this->WORK.clear();

    this->VT.initialize();
    this->svdDim = 0;
	}

	void initialize(const DGESVD& dgesvd)
	{
	this->A = dgesvd.A;
    this->U = dgesvd.U;

    this->singularValues = dgesvd.singularValues;
    this->VT             = dgesvd.VT;
    this->svdDim         = dgesvd.svdDim;
	}


	bool equalMatrixDimensions(const LapackMatrix& A, const LapackMatrix& B)
	{
		if((A.rows != B.rows)||(A.cols != B.cols)) return false;
		return true;
	}

	void computeSVD(const LapackMatrix& A)
	{
		// Copy A, since A is overwritten or destroyed by
	    // the routine.


		if(not equalMatrixDimensions(this->A,A))
		{
		this->A.initialize(A);
		}
		else
		{
		this->A = A;
		}


		long M  = this->A.getRowDimension();
		long N  = this->A.getColDimension();

		// S dimension (min(M,N))

		long minMN = (M < N) ? M : N;

	    if((U.rows != M) || (U.cols != M)) {U.initialize(M,M);}

		singularValues.resize(minMN,0.0);

		if((VT.rows != N) || (VT.cols != N)) {VT.initialize(N,N);}


		long  LDA = M;
		long  LDU = M;
		long LDVT = N;

		JOBU  = 'A';  //  All M columns of U are returned in array U:
	    JOBVT = 'A';  //  All N rows of V**T are returned in the array VT;

	   // Query to obtain optimal work array size

		long LWORK = -1;
		INFO       = 0;

		double WORKtmp;
	    dgesvd_(&JOBU, &JOBVT, &M, &N, this->A.dataPtr, &LDA, &singularValues[0], U.dataPtr, &LDU, VT.dataPtr, &LDVT,
               &WORKtmp, &LWORK, &INFO);

	    LWORK = (long)(WORKtmp + 100);
		WORK.resize(LWORK);

		double* WORKptr =    &WORK[0];

		// Second call to create svd

	    INFO       = 0;

	    dgesvd_(&JOBU, &JOBVT, &M, &N, this->A.dataPtr, &LDA, &singularValues[0], U.dataPtr, &LDU, VT.dataPtr, &LDVT,
               WORKptr, &LWORK, &INFO);

		if(INFO != 0)
        {
        cerr << "DGESVD  Failed : INFO = " << INFO  << endl;
        exit(1);
        }
	}

	void computeThinSVD(const LapackMatrix& A)
	{
		// Copy A, since A is overwritten or destroyed by
	    // the routine.


		if(not equalMatrixDimensions(this->A,A))
		{
		this->A.initialize(A);
		}
		else
		{
		this->A = A;
		}


		long M  = this->A.getRowDimension();
		long N  = this->A.getColDimension();

		// S dimension (min(M,N))

		long minMN = (M < N) ? M : N;

	    if((U.rows != M) || (U.cols != minMN)) {U.initialize(M,minMN);}

		singularValues.resize(minMN,0.0);

		if((VT.rows != minMN) || (VT.cols != N)) {VT.initialize(minMN,N);}


		long  LDA = M;
		long  LDU = M;
		long LDVT = minMN;

		JOBU  = 'S';  //  minMN columns of U are returned in array U:
	    JOBVT = 'S';  //  minMN rows of V**T are returned in the array VT;

	   // Query to obtain optimal work array size

		long LWORK = -1;
		INFO       = 0;

		double WORKtmp;
	    dgesvd_(&JOBU, &JOBVT, &M, &N, this->A.dataPtr, &LDA, &singularValues[0], U.dataPtr, &LDU, VT.dataPtr, &LDVT,
               &WORKtmp, &LWORK, &INFO);

	    LWORK = (long)(WORKtmp + 100);
		WORK.resize(LWORK);

		double* WORKptr =    &WORK[0];

		// Second call to create svd

	    INFO       = 0;

	    dgesvd_(&JOBU, &JOBVT, &M, &N, this->A.dataPtr, &LDA, &singularValues[0], U.dataPtr, &LDU, VT.dataPtr, &LDVT,
               WORKptr, &LWORK, &INFO);

		if(INFO != 0)
        {
        cerr << "THIN DGESVD  Failed : INFO = " << INFO  << endl;
        exit(1);
        }
	}

	long getSVDdim(){return svdDim;}

	vector<double> applyPseudoInverse(vector<double>& b, double svdCutoff = -1.0)
	{

		vector<double> x(A.cols,0.0);

		svdDim = 0;

		if(svdCutoff < 0){svdDim = singularValues.size();}
		else
		{
	    for(long i = 0; i < (long)singularValues.size(); i++)
	    {
	    	if(singularValues[i] > svdCutoff){svdDim++;}
	    }
		}

		if(svdDim == 0) return x;

		//
		// Construct pseudo-inverse using components of SVD
		//

		vector<double> xStar;

        /*
		LapackMatrix   Ustar = U.getColSlice(0,svdDim-1);
		xStar = Ustar.applyTranspose(b);
        */

		char TRANS     = 'T';
		long Mstar     = U.rows;
		long Nstar     = svdDim;
		double ALPHA   = 1.0;
		double BETA    = 0.0;
		long LDA       = Mstar;
		long INCX      = 1;
		long INCY      = 1;

		xStar.resize(svdDim,0.0);

		dgemv_(&TRANS,&Mstar,&Nstar,&ALPHA,U.dataPtr,&LDA,&b[0],&INCX,&BETA,&xStar[0],&INCY);

		for(long i = 0; i < svdDim; i++)
		{
			xStar[i] /= singularValues[i];
		}

		LapackMatrix Vstar;

	    /*
	    Vstar = VT.getRowSlice(0,svdDim-1);
	    x = Vstar.applyTranspose(xStar);
        */

		if(svdDim == VT.rows)
		{
	    TRANS     = 'T';
		Mstar     = VT.rows;
		Nstar     = VT.cols;
		ALPHA   = 1.0;
		BETA    = 0.0;
		LDA       = Mstar;
		INCX      = 1;
		INCY      = 1;
		dgemv_(&TRANS,&Mstar,&Nstar,&ALPHA,VT.dataPtr,&LDA,&xStar[0],&INCX,&BETA,&x[0],&INCY);
		}
		else // Extract columns of Vstar as first svdDim rows of VT
		{
	    Vstar.initialize(VT.cols,svdDim);
	    for(long i = 0; i < VT.cols; i++)
	    {
	    for(long j = 0; j < svdDim;  j++)
	    {
	    Vstar(i,j) = VT(j,i);
	    }}

	    TRANS     = 'N';
		Mstar     = Vstar.rows;
		Nstar     = Vstar.cols;
		ALPHA   = 1.0;
		BETA    = 0.0;
		LDA       = Mstar;
		INCX      = 1;
		INCY      = 1;
		dgemv_(&TRANS,&Mstar,&Nstar,&ALPHA,Vstar.dataPtr,&LDA,&xStar[0],&INCX,&BETA,&x[0],&INCX);
		}

		return x;
	}

    LapackMatrix     A;
    LapackMatrix     U;

    vector<double> singularValues;

    LapackMatrix           VT;

	char                JOBU;
	char               JOBVT;
	long                INFO;
	vector<double>      WORK;

	long               svdDim;


};

/*
   Upon construction or initialization with an input LapackMatrix instance A
   this class provides an interface to various methods for computing 
   the solution of A * X = B 

 */


// Calling Fortran versions directly using appropriately modified prototypes.
//

// Data type mapping used:
//
// C++  int    ==  Fortran LOGICAL
// C++  long   ==  Fortran INTEGER
// C++  double ==  Fortran DOUBLE PRECISION
//
// Linking to the Fortran routines  using -llapack -lblas
//
// Oct. 15, 2016
//


/*
*       SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*
*
*       extern "C" void dsyev_(char* JOBZ,char* UPLO, long*N, double* Aptr, long* LDA, double* Wptr,
*		double* WORKptr, long* LWORK, long* INFO);
*
*
*/


class DSYEV
{
public:

	DSYEV()
	{
		JOBZ = 'N';
		UPLO = 'U';            // Using upper triangular part of A
	}

	void computeEigenvalues(const LapackMatrix& A, vector<double>& eigenValues)
	{
		assert(A.sizeCheck(A.rows,A.cols));

		U.initialize(A);
		JOBZ = 'N';
		UPLO = 'U';            // Using upper triangular part of A

	    long N       = A.rows;
		double* Uptr = U.dataPtr;

		long LDA = N;

		eigenValues.resize(N);
		double*Wptr = &eigenValues[0];

		long LWORK = -1;

		double WORKtmp;

		long INFO = 0;

		// First call to get optimal workspace

		dsyev_(&JOBZ,&UPLO,&N, Uptr, &LDA, Wptr, &WORKtmp, &LWORK, &INFO);

		LWORK = (long)(WORKtmp + 100);
		vector<double>    WORK(LWORK);
		double* WORKptr =    &WORK[0];

		// Second call to create eigensystem

	    dsyev_(&JOBZ,&UPLO,&N, Uptr, &LDA, Wptr, WORKptr, &LWORK, &INFO);

		if(INFO != 0)
        {
        cerr << "dsyev  Failed : INFO = " << INFO  << endl;
        exit(1);
        }

	}

	void computeEigensystem(const LapackMatrix& A, vector<double>& eigenValues, vector < vector < double> >& eigenVectors)
	{
		assert(A.sizeCheck(A.rows,A.cols));
		computeEigensystem(A,eigenValues, U);

		// Pack eigenvectors into return argument

		vector <double> eigVector(A.rows);

		eigenVectors.resize(A.cols,eigVector);

		for(long j = 0; j < A.cols; j++)
		{
			for(long i = 0; i < A.rows; i++)
			{
				eigenVectors[j][i] = U(i,j);
			}
		}

	}

	void computeEigensystem(const LapackMatrix& A, vector<double>& eigenValues, LapackMatrix& eigenVectors)
	{
		assert(A.sizeCheck(A.rows,A.cols));

        eigenVectors.initialize(A);
        JOBZ = 'V';
        UPLO = 'U';            // Using upper triangular part of A

        long N       = A.rows;
		double* Uptr = eigenVectors.dataPtr;

		long LDA = N;

		eigenValues.resize(N);
		double*Wptr = &eigenValues[0];

		long LWORK = -1;

		double WORKtmp;

		long INFO = 0;

		// First call to get optimal workspace

		dsyev_(&JOBZ,&UPLO,&N, Uptr, &LDA, Wptr, &WORKtmp, &LWORK, &INFO);

		LWORK = (long)(WORKtmp + 100);
		vector<double>    WORK(LWORK);
		double* WORKptr =    &WORK[0];

		// Second call to create eigensystem

	    dsyev_(&JOBZ,&UPLO,&N, Uptr, &LDA, Wptr, WORKptr, &LWORK, &INFO);

		if(INFO != 0)
        {
        cerr << "dsyev  Failed : INFO = " << INFO  << endl;
        exit(1);
        }
	}

	LapackMatrix    U;
	vector<double>  eigValues;
	char            JOBZ;
	char            UPLO;

};


/*
*  SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
*                          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
*                          WORK, IWORK, INFO )



extern "C" void dgesvx_(char* FACT, char* TRANS, long* N, long* NRHS, double* Aptr, long* LDA, double* AFptr, long* LDAF, long* IPIVptr,
		                char* EQED, double* Rptr, double* Cptr, double* Bptr, long* LDB, double* Xptr, long* LDX,  double* RCOND,
						double* FERR, double* BERR, double* WORKptr,
		                long* IWORKptr, long* INFO);

*/

class DGESVX
{
public:

	DGESVX()
	{}


    void applyInverse(const LapackMatrix& A,vector <double >& b)
	{
			applyInverse(A,&b[0]);
	}

    void applyInverse(const LapackMatrix& A,LapackMatrix& b)
	{
    		applyInverse(A,b.dataPtr,b.cols);
	}

	void applyInverse(const LapackMatrix& A, double* b, long NRHS = 1)
	{
		assert(A.sizeCheck(A.rows,A.cols));

		char FACT  = 'E'; // Equilibrate, then factor
		char TRANS = 'N'; // No transpose
		long N     = A.rows;

		// Allocate temporaries

		this->A.initialize(A);
		this->AF.initialize(N,N);

		double* Aptr  =  A.dataPtr;
		double* AFptr = AF.dataPtr;

		long LDA   = N;
		long LDAF  = N;

		vector <long >   IPIV(N);
		long* IPIVptr = &IPIV[0];

		char  EQED;

		vector<double>   R(N);
		double* Rptr  = &R[0];

		vector<double>    C(N);
		double* Cptr  =  &C[0];

		vector<double>   B(N*NRHS);
		double* Bptr  =      &B[0];
	    long LDB      =          N;


		// b will be overwritten with the solution
	    // so no need to declare X separately

		double* Xptr = b;
		long LDX     = N;

		FERR.resize(NRHS);
		BERR.resize(NRHS);

		vector<double>   WORK(4*N);
		double* WORKptr = &WORK[0];

		vector<long>       IWORK(N);
		long* IWORKptr  = &IWORK[0];

		long   INFO = 0;


		// Assign right hand side to B

		for(long i = 0; i < N*NRHS; i++)
		{
			Bptr[i] = b[i];
		}

		dgesvx_(&FACT, &TRANS, &N, &NRHS, Aptr, &LDA, AFptr, &LDAF, IPIVptr,
		        &EQED, Rptr, Cptr, Bptr,&LDB, Xptr, &LDX, &RCOND,
				&FERR[0], &BERR[0], WORKptr, IWORKptr, &INFO);


		if(INFO != 0)
        {
        cerr << "dgesvx  Failed : INFO = " << INFO  << endl;
        exit(1);
        }
	}

/*
*  RCOND is DOUBLE PRECISION
*  The estimate of the reciprocal condition number of the matrix
*  A after equilibration (if done).  If RCOND is less than the
*  machine precision (in particular, if RCOND = 0), the matrix
*  is singular to working precision.  This condition is
*  indicated by a return code of INFO > 0.
*/

	double getReciprocalConditionNumber()
	{
		return RCOND;
	}

/*
*  FERR is DOUBLE PRECISION array, dimension (NRHS)
*  The estimated forward error bound for each solution vector
*  X(j) (the j-th column of the solution matrix X).
*  If XTRUE is the true solution corresponding to X(j), FERR(j)
*  is an estimated upper bound for the magnitude of the largest
*  element in (X(j) - XTRUE) divided by the magnitude of the
*  largest element in X(j).  The estimate is as reliable as
*  the estimate for RCOND, and is almost always a slight
*  overestimate of the true error.
*/

	double getSolutionErrorEst()
	{
		return FERR[0];
	}

	vector<double> getMultipleSolutionErrorEst()
	{
		return FERR;
	}

/*
*  BERR is DOUBLE PRECISION array, dimension (NRHS)
*  The componentwise relative backward error of each solution
*  vector X(j) (i.e., the smallest relative change in
*   any element of A or B that makes X(j) an exact solution).
*/
	double getSolutionBackwardErrorEst()
	{
		return BERR[0];
	}

	vector<double> getMultipleSolutionBackwardErrorEst()
	{
		return BERR;
	}



    double         RCOND;
	vector<double>  FERR;
	vector<double>  BERR;

	LapackMatrix       A;
	LapackMatrix      AF;

};



/*
*    DPOSV computes the solution to a real system of linear equations
*       A * X = B,
*    where A is an N-by-N symmetric positive definite matrix and X and B
*    are N-by-NRHS matrices.
*
*    The Cholesky decomposition is used to factor A as
*       A = U**T* U,  if UPLO = 'U', or
*       A = L * L**T,  if UPLO = 'L',
*    where U is an upper triangular matrix and L is a lower triangular
*    matrix.  The factored form of A is then used to solve the system of
*    equations A * X = B.
*
*    SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )

*/

extern "C" void dposv_(char* UPLO, long* N, long* NRHS, double* Aptr, long* LDA, double* Bptr, long* LDB, long* INFO);


class DPOSV
{
public:

	DPOSV()
	{}


    void applyInverse(const LapackMatrix& A,vector <double >& b)
	{
			applyInverse(A,&b[0]);
	}

    void applyInverse(const LapackMatrix& A,LapackMatrix& b)
	{
    		applyInverse(A,b.dataPtr,b.cols);
	}

	void applyInverse(const LapackMatrix& A, double* b, long NRHS = 1)
	{
		assert(A.sizeCheck(A.rows,A.cols));

		char UPLO  = 'U'; // A = U**T* U
		long N     = A.rows;

		double* Aptr  =  A.dataPtr;

		long LDA    = N;
		double* Bptr = b;
	    long LDB    = N;
		long   INFO = 0;

		dposv_(&UPLO, &N, &NRHS,Aptr,&LDA,Bptr,&LDB,&INFO);

		if(INFO != 0)
        {
        cerr << "dposv Failed : INFO = " << INFO  << endl;
        exit(1);
        }
	}
};

//
// This class creates the solution of the normal equations with the singular value
// parameter cut-off value svdCutoff. Since the singular values of the normal
// equations are square roots of the eigenvalues of the normal equations, the
// accuracy of the solution is limited to problems where it is sufficient to
// construct a solution using singular vector components associated with
// singular values where (singular values)^2 > accuracy of the
// eigensystem routine. Typically this limits the application to problems where
// svdCutoff is approximately > 1.0e-07. The returned solution can only be
// expected to have single precision accuracy for the same reasons.
//
// In the implementation below, the solution components are only accumulated
// when the singular values are greater than the max(svdCutoff,1.0e-07).
//
// When M < N, so the solution is underdetermined, the minimum L2
// norm solution is returned and only M approximate singular
// values are computed.
//
class NORMALEQ
{
public:

	NORMALEQ(){};

	long getSVDdim(){return svdDim;}


	// Normal equation solution of
	//
	// A*x = b
	//
	// with svdCutoff.
	//
	// Only components with singular values > svdCutoff or 1.0e-07 are accumulated.
	//
	// When A is M x N and M < N, then only M singular values are evaluated.
	//
	//
	vector<double> computeNormalEquationSolution(vector<double>& b, LapackMatrix& AB, double svdCutoff)
	{
		// Form right hand side

		bBar = AB.applyTranspose(b);

		// Form normal equations

		long M = AB.getRowDimension();
		long N = AB.getColDimension();

		vector<double> x(N,0.0);

		ABnormal.initialize(N,N);

		// Call dgemm to form normal equations. This is done using a direct call to dgemm so
		// that forming A^T isn't required.

	    char TRANSA  = 'T';
	    char TRANSB  = 'N';
	    long  Mstar  = N;
	    long  Nstar  = N;
	    long Kstar   = M;
	    double ALPHA = 1.0;
	    double BETA  = 0.0;
	    long LDA     = M;
	    long LDB     = M;
	    long LDC     = N;
	    long INCX    = 1;
	    long INCY    = 1;

		dgemm_(&TRANSA,&TRANSB,&Mstar,&Nstar,&Kstar,&ALPHA, AB.dataPtr,&LDA,AB.dataPtr, &LDB,&BETA,ABnormal.dataPtr,&LDC);

		// Diagonalize the normal equations

	    eigenValues.resize(N,0.0);

	    if(M >= N)
	    {
	    singularValues.resize(N,0.0);
	    }
	    else
	    {
	    singularValues.resize(M,0.0);
	    }

	    eigenVectors.initialize(N,N);
		diagonalizer.computeEigensystem(ABnormal, eigenValues, eigenVectors);


	    // The eigenvalues from the diagonalization routine
        // are returned smallest to largest, so pack the singular values in reverse order
        // to agree with the SVD convention

		long      stopIndex = 0;
        if(M < N) stopIndex = N-M;

        long svIndex             = 0;
        long firstComponentIndex = N-1;
        long componentCount      = 0;
        for(long i= N-1; i >= stopIndex; i--)
        {
        	singularValues[svIndex] = sqrt(abs(eigenValues[i]));
        	if((singularValues[svIndex] > svdCutoff) && (singularValues[svIndex] > 1.0e-7))
        	{
        	  firstComponentIndex = i;
        	  componentCount++;
        	}
        	svIndex++;
        }

        /*
		cout << std::setprecision(15) << endl;
        for(long i = 0; i < eigenValues.size(); i++)
        {
        	cout << sqrt(abs(eigenValues[i])) <<  " " << abs(eigenValues[i]) << endl;
        }
        */


        // Project solution onto selected subspace

        // A direct call to dgemv to create
        //
        // bStar = eigenVectors.applyTranspose(bBar);
        //
        // in order to avoid accumulating the values known
        // a-priori to be set to zero
        //

        bStar.resize(N,0.0);

        char TRANS = 'T';
        Mstar      =  N;
        Nstar      =  componentCount;
        ALPHA      = 1.0;
        BETA       = 0.0;
        INCX       = 1;
        INCY       = 1;
        LDA        = N;

        double* bPtr   = &bStar[firstComponentIndex];
        double* eigPtr = eigenVectors.dataPtr + firstComponentIndex*N;

        dgemv_(&TRANS,&Mstar,&Nstar,&ALPHA,eigPtr,&LDA,&bBar[0],&INCX,&BETA,bPtr,&INCY);

        // Solve he diagonal system for the components being kept

        for(long i= N-1; i >= firstComponentIndex; i--)
        {
            bStar[i] /= abs(eigenValues[i]);
        }

        // Zero out entries in case M < N

        for(long i = firstComponentIndex-1; i >= 0; i--)
        {
        	bStar[i] = 0;
        }

        //
        // Evaluate the solution. Using direct call to dgemv to avoid
        // accumulating columns with zero weights.
        //
        // x = eigenVectors*bStar;

        x.resize(M,0.0);

        TRANS      = 'N';
        Mstar      =  N;
        Nstar      =  componentCount;
        ALPHA      = 1.0;
        BETA       = 0.0;
        INCX       = 1;
        INCY       = 1;
        LDA        = N;

        bPtr   = &bStar[firstComponentIndex];
        eigPtr = eigenVectors.dataPtr + firstComponentIndex*N;

        dgemv_(&TRANS,&Mstar,&Nstar,&ALPHA,eigPtr,&LDA,bPtr,&INCX,&BETA,&x[0],&INCY);

        svdDim = componentCount;

        return x;

	}

	vector<double>          bStar;
	vector<double>           bBar;
	LapackMatrix          ABnormal;

	DSYEV             diagonalizer;
	vector<double>    eigenValues;
	LapackMatrix      eigenVectors;

	vector<double> singularValues;
		long               svdDim;

};

//
// Class QRutility can be used to create a solution to
//
//               A x = b
//
// where A is a full rank M x N matrix with M >= N.
// If M > N then a least squares approximation to the
// solution is determined.
//
// The solution process is decomposed into two
// steps; a create QR factors step followed by
// a create solution step. The QR factorization
// is retained internally so that, for a given
// system and multiple right hand sides, the QR
// factorization need only be computed once.
//
// The procedure used is essentially that of
// the LAPACK routine DGELS performed without scaling.
// That routine, as well as this class, creates
// a QR factorization without pivoting, so no
// rank determination is performed. For the solution
// of a system with undetermined properties, one
// might consider using SCC::DGELSY or another class
// that uses Lapack routine DGELSY to construct
// QR solutions.
//
class QRutility
{
    public :

    QRutility()
    {
    initialize();
    }

    QRutility(const QRutility& Q)
    {
    initialize(Q);
    }

    void initialize()
    {
    QRfactors.initialize();
    TAU.clear();
    WORK.clear();
    bTemp.clear();
    dormqrWorkSize = -1;
    }

    void initialize(const QRutility& Q)
    {
    QRfactors.initialize(Q.QRfactors);
    TAU            = Q.TAU;
    WORK           = Q.WORK;
    bTemp          = Q.bTemp;
    dormqrWorkSize = Q.dormqrWorkSize;
    }

    vector<double> createQRsolution(vector<double>& b)
    {
    try {if(QRfactors.getRowDimension() != (long)b.size()) {throw std::runtime_error("createQRsolution");}}
    catch (std::runtime_error& e)
    {
     cerr << "Runtime exception in QRutility member function " <<  e.what() << endl;
     cerr << "createSolution called be fore createQRfactors " << endl;
     exit(1);
    }

    // Capture right hand side

    bTemp = b;

    char SIDE    = 'L';
    char TRANS   = 'T';
    long M       = (long)b.size();
    long NRHS    = 1;
    long K       = M;
    long LDA     = QRfactors.getRowDimension();
    double* Aptr = QRfactors.getDataPointer();
    long LDC     = M;

    long LWORK  = -1;
    long INFO   =  0;

    double WORKDIM;

    // Obtain the optimal work size if it has not already
    // been determined

    if(dormqrWorkSize < 0)
    {
    dormqr_(&SIDE, &TRANS, &M, &NRHS, &K,Aptr,& LDA, &TAU[0], &bTemp[0],
    &LDC, &WORKDIM, &LWORK, &INFO);

    dormqrWorkSize = (long)WORKDIM+1;
    }

    LWORK = dormqrWorkSize;
    WORK.resize(LWORK);

    dormqr_(&SIDE, &TRANS, &M, &NRHS, &K,Aptr,& LDA, &TAU[0], &bTemp[0],
    &LDC, &WORK[0], &LWORK, &INFO);

    try {if(INFO != 0) {throw std::runtime_error("createQRsolution");}}
    catch (std::runtime_error& e)
    {
     cerr << "Runtime exception in QRutility member function " <<  e.what() << endl;
     cerr << "DORMQR Failed : INFO = " << INFO  << endl;
     exit(1);
    }

    // Note: only using QRfactors.cols elements of bTemp.

	// Backsolve upper trignular system to obtain a solution

	char UPLO    = 'U';
	TRANS        = 'N';
	char DIAG    = 'N';
	M            = QRfactors.getColDimension();
	long N       = M;
	NRHS         = 1;
	LDA          = QRfactors.getRowDimension();
	long LDB     = M;

	dtrtrs_(&UPLO, &TRANS, &DIAG, &N, &NRHS, Aptr, &LDA,&bTemp[0],&LDB,&INFO);

    try {if(INFO != 0) {throw std::runtime_error("createQRsolution");}}
    catch (std::runtime_error& e)
    {
     cerr << "Runtime exception in QRutility member function " <<  e.what() << endl;
     cerr << "DTRTRS detected singular matrix : INFO = " << INFO  << endl;
     exit(1);
    }

	return bTemp;
    }


    //
    // A convenience solve interface. It is assumed that Bptr points to contiguous
    // data of size of the number of rows of A. No bounds checking is performed.
    //

	vector<double> createQRsolution(double* Bptr)
	{
		long M  = QRfactors.getRowDimension();
		vector<double>                   B(M);
		std::memcpy(B.data(),Bptr,M*sizeof(double));
		return createQRsolution(B);
	}
//
//  This member function creates and stores internally the QR factorization
//  of the M x N input matrix A. It is assumed that M >= N and A is of
//  full rank.
//
	void createQRfactors(const SCC::LapackMatrix& A)
	{
    long M   = A.getRowDimension();
    long N   = A.getColDimension();

    try
    {
    if( M < N ) throw std::runtime_error("createQRfactors");
    }
    catch (std::runtime_error& e)
    {
     cerr << "Runtime exception in QRutility member function " <<  e.what() << endl;
     cerr << "Input matrix rows < cols " << '\n';
     cerr << "rows : " << M << " cols : " << N << endl;
     exit(1);
    }

	// Capture system

	QRfactors.initialize(A);

    // Create QR factors

    long LDA = M;

    TAU.resize(N);

    long LWORK  = -1;
    long INFO   =  0;

    double WORKDIM;

    // First call to obtain the optimal work size

    dgeqrf_(&M, &N, A.getDataPointer(), &LDA, &TAU[0],&WORKDIM, &LWORK, &INFO);

    LWORK = (long)WORKDIM+1;
    WORK.resize(LWORK);

    // Create QR factors

    dgeqrf_(&M, &N, QRfactors.getDataPointer(), &LDA, &TAU[0],&WORK[0], &LWORK, &INFO);

    try {if(INFO != 0) {throw std::runtime_error("createQRfactors");}}
    catch (std::runtime_error& e)
    {
     cerr << "Runtime exception in QRutility member function " <<  e.what() << endl;
     cerr << "DGEQRF Failed : INFO = " << INFO  << endl;
     exit(1);
    }

    // Reset work sizes to -1  to trigger a re-computation
    // of work storage requirements on first use
    // of createQRsolution

    dormqrWorkSize = -1;
	}

    SCC::LapackMatrix QRfactors; // QR factor component as returned by dgeqrf
    vector<double>          TAU; // QR factor component as returned by dgeqrf
    vector<double>         WORK;
    vector<double>        bTemp;

    long         dormqrWorkSize;
};


} // End of SCC Namespace declaration

#endif
