#include "Matrix.h" 
#include "Lapack.h"
#include "Blas.h" 
#include <stdio.h>

#pragma region "Constructor"
	Matrix::Matrix(void)
	{
		this->ni = 0;
		this->mj = 0;
	} 


	//Ain matrix in read in row by row
	Matrix::Matrix( double Ain[], int nd )
	{
		this->isTridiagonal = false;
        this->isSquareMatrix = true;
		this->ni = nd;
		this->mj = nd; 
		//store the matrix in transpose form, ie store by column (convention used by lapack)
		At = new double[nd*nd];
		int k = 0;
		for (int i = 0; i<ni; i++)
		{
			for (int j = 0; j<ni; j++)
			{
				At[i + ni * j] = Ain[k];
				k = k + 1;
			}
		}
	}


	Matrix::Matrix(double Ain[], int nd, int md, bool isTranspose)
	{
		this->isTridiagonal = false;
		this->ni = nd;
		this->mj = md;
		if (nd == md)
		{
			this->isSquareMatrix = true;
		}
		else
		{
			this->isSquareMatrix = false;
		}
		//x = new double[md];
		At = new double[nd*md];
		if (isTranspose)
		{
			int k = 0;
			for (int i = 0; i<ni; i++)
			{
				for (int j = 0; j<mj; j++)
				{
					At[i + ni * j] = Ain[k];
					k = k + 1;
				}
			}
		}
		else
		{
			for (int j = 0; j<(nd*md); j++)  //the ith index of each row
			{
				At[j] = Ain[j];
			}
		}
	}


	Matrix::Matrix(double Ain[], int nd, int md, bool isTranspose, int flag)
    {
		this->isTridiagonal = false;
        this->info = flag;
		this->ni = nd;
		this->mj = md;
		if (nd == md)
		{
			this->isSquareMatrix = true;
		}
		else
		{
			this->isSquareMatrix = false;
		}
		//x = new double[md];
		At = new double[nd*md];
		if (isTranspose)
		{   
			int k = 0;
			for (int i = 0; i<ni; i++)
			{
				for (int j = 0; j<mj; j++)
				{
					At[i + ni * j] = Ain[k];
					k = k + 1;
				}
			} 
		}
		else
		{   
			for (int j = 0; j<(nd*md); j++)  //the ith index of each row
			{
				At[j] = Ain[j];
			}
		}  
    }


	//Ain matrix in read in row by row
	Matrix::Matrix(double Ain[], int nd , int md)
	{  
		this->isTridiagonal = false;
        if (nd == md)
        { 
			this->isSquareMatrix = true;
        }
        else
        { 
			this->isSquareMatrix = false;
        }
		this->ni = nd;
		this->mj = md; 
		//store the matrix in transpose form, ie store by column
		At = new double[nd*md];
		//x = new double[md]; 
		int k = 0;
		for (int i = 0; i<ni; i++)
		{
			for (int j = 0; j<mj; j++)
			{
				At[i + ni * j] = Ain[k];
				k = k + 1;
			}
		}  
	}
 

	Matrix::Matrix(double ain[], double bin[], double cin[], int nd) 
	{  
        this->isTridiagonal = true;
		this->isSquareMatrix = true;
		this->ni = nd;
		this->mj = nd;
		alpha = new double[nd-1];
        beta = new double[nd]; 
		gamma = new double[nd-1];
		for (int i = 0; i<(nd-1); i++)
		{
			alpha[i] = ain[i];
			beta[i] = bin[i];
			gamma[i] = cin[i];
		}
		beta[nd-1] = bin[nd-1];
    }
     

	Matrix::Matrix(double Ain[], double din[], double Hin[], double fin[], int nd , int md, int pd)
	{
		this->isTridiagonal = false;
		this->ni = nd;
		this->mj = md;
		this->pk = pd;
		//x = new double[md];
		At = new double[nd*md];
		d = new double[nd];
        Ht = new double[pd*md];
        f = new double[pd]; 
		int k = 0;
		for (int i = 0; i<ni; i++)
		{
			for (int j = 0; j<mj; j++)
			{
				At[i + ni * j] = Ain[k];
				k = k + 1;
			}
		} 
		k = 0;
		for (int i = 0; i<pk; i++)
		{
			for (int j = 0; j<mj; j++)
			{
				Ht[i + pk * j] = Hin[k];
				k = k + 1;
			}
		}  
		for (int j = 0; j<nd; j++)
        {
			d[j] = din[j];
        }
		for (int j = 0; j<pd; j++)
		{
			f[j] = fin[j];
		} 
	}


	Matrix::Matrix(double Ain[], double Bin[], double din[], int nd, int md, int pd)
	{
		this->isTridiagonal = false;
		this->ni = nd;
		this->mj = md;
		this->pk = pd;
		At = new double[nd*md];
		Bt = new double[nd*pd];
		d = new double[nd];
		//x = new double[md];
		//y = new double[pd];
		for (int j = 0; j<nd; j++)
		{
			d[j] = din[j];
		}
		int k = 0;
		for (int i = 0; i<ni; i++)
		{
			for (int j = 0; j<mj; j++)
			{
				At[i + ni * j] = Ain[k];
				k = k + 1;
			}
		} 
	    k = 0;
		for (int i = 0; i<nd; i++)
		{
			for (int j = 0; j<pd; j++)
			{
				Bt[i + nd * j] = Bin[k];
				k = k + 1;
			}
		} 
    }
#pragma endregion


#pragma region "Methods"    

	int Matrix::solve(double xinout[])
    {
		int lwork = -1;
		double wkopt;
		double* work; 
		dgglse_(&ni, &mj, &pk, At, &ni, Ht, &pk, d, f, xinout, &wkopt, &lwork, &info);
		lwork = (int)wkopt;
		work = new double[lwork];
		dgglse_(&ni, &mj, &pk, At, &ni, Ht, &pk, d, f, xinout, work, &lwork, &info);
        return info; 
    }
    
 
	int Matrix::solve(double xinout[], double yinout[])
	{ 
		int lwork = -1;
		double wkopt;
		double* work;
		dggglm_(&ni, &mj, &pk, At, &ni, Bt, &ni, d, xinout, yinout, &wkopt, &lwork, &info);
		lwork = (int)wkopt;
		work = new double[lwork];
		dggglm_(&ni, &mj, &pk, At, &ni, Bt, &ni, d, xinout, yinout, work, &lwork, &info);
		return info;
    }


	void Matrix::printMatrix(double input[], int nd, int md)
	{
		int k = 0;
		for (int i = 0; i<nd; i++)
		{
			for (int j = 0; j<md; j++)
			{
				printf(" %lf  ", input[k]);
				k = k + 1;
			}
			printf("\n");
		}
	}
 

	void Matrix::printMatrix()

	{
		printf("\n");
		int k = 0;
		for (int i = 0; i<ni; i++)
		{
			for (int j = 0; j<mj; j++)
			{
				printf("  %lf   ", At[i + ni * j]);
				k = k + 1;
			}
			printf("\n\n");
		} 
	}
#pragma endregion


#pragma region "Operators"    
	Matrix Matrix::operator + (Matrix operand)
    {   
		int n = ni * mj;
        double* out = new double[n];
        for( int i = 0; i < n ; i++ )
        {
            out[i] = operand.At[i];
        } 
		double one = 1;
		int inc = 1; 
		daxpy_(&n, &one, At, &inc, out, &inc);
        return Matrix(out,ni,mj, false);
    }


	Matrix Matrix::operator - (Matrix operand)
	{
		int n = ni * mj;
		double* P = new double[n];
		double* Q = new double[n];
		for (int i = 0; i < n; i++)
		{
			P[i] = operand.At[i];
			Q[i] = At[i];
		}
		double one_1 = -1;
		int inc = 1;
		daxpy_(&n, &one_1,P, &inc, Q, &inc);
		return Matrix(Q, ni, mj, false);
	}


	Matrix Matrix::operator * (Matrix operand)
	{
        int n = ni * operand.mj;
		double* Q = new double[n]; 
		char trans[] = { "N" };
		double one = 1.0;
		double zero = 0.0;
		dgemm_(trans, trans, &ni, &operand.mj, &mj, &one, At, &ni, operand.At, &mj, &zero, Q, &ni); 
		return Matrix(Q, ni, operand.mj,false);

	} 


	Matrix Matrix::operator * (double operand[])
	{
		double* out = new double[ni];
		char trans[] = { "N" };
		double one = 1.0;
		double zero = 0.0;
		dgemm_(trans, trans, &ni, &mj, &mj, &one, At, &ni, operand, &mj, &zero, out, &ni);
		return Matrix(out, ni, 1, false);
	}


	//scaling: y = s*y
	Matrix Matrix::operator * (double s)
    { 
		int n = ni * mj;
		int incx = 1;
		double* inputAt = new double[ni*mj];
		for (int j = 0; j<(ni*mj); j++)
		{
			inputAt[j] = At[j];
		}
		dscal_(&n, &s, inputAt, &incx);
		return Matrix(inputAt,  ni, mj, false);
    }

    //scaling: y = s*y
	Matrix operator * (double s, Matrix operand)
	{
		int n = operand.ni * operand.mj;
		int incx = 1;
		double* inputAt = new double[operand.ni*operand.mj];
		for (int j = 0; j<(operand.ni*operand.mj); j++)
		{
			inputAt[j] = operand.At[j];
		}
		dscal_(&n, &s, inputAt, &incx);
		return Matrix(inputAt, operand.ni, operand.mj, false);
	}


	Matrix Matrix::operator | (double r[]) 
    {
        if (isSquareMatrix)
        { 
            if (isTridiagonal)
            {
				int bl = 1;
				double* aa = new double[ni - 1];
				double* bb = new double[ni];
				double* cc = new double[ni - 1];
				double* dd = new double[ni];
				for (int i = 0; i<(ni - 1); i++)
				{
					aa[i] = alpha[i];
					bb[i] = beta[i];
					cc[i] = gamma[i];
					dd[i] = r[i];
				}
				bb[ni - 1] = beta[ni - 1];
				dd[ni - 1] = r[ni - 1];
				dgtsv_(&ni, &bl, aa, bb, cc, dd, &ni, &info);
				return Matrix(dd, mj, bl, false, info); 
            }
            else
            {
				int nrh = 1;
				int* ip = new int[ni];
				//the subroutine dgesv changes the matrix input, need to create a copy of the input matrix 
				double* inputAt = new double[ni*mj];
				for (int j = 0; j<(ni*mj); j++)
				{
					inputAt[j] = At[j];
				}
				double* out = new double[ni];
				for (int j = 0; j<ni; j++)
				{
					out[j] = r[j];
				}
				dgesv_(&ni, &nrh, inputAt, &ni, ip, out, &ni, &info);
				return Matrix(out, mj, nrh, false, info); 
            } 
        }
        else
        { 
			int nrhs = 1; 
			int lwork = -1;
			double wkopt;
			double* work;
			char trans[] = { "N" };
			double* inputAt = new double[ni*mj];
			for (int j = 0; j<(ni*mj); j++)
			{
				inputAt[j] = At[j];
			}
			double* out = new double[ni];
			for (int j = 0; j<ni; j++)
			{
				out[j] = r[j];
			}
			dgels_(trans, &ni, &mj, &nrhs, inputAt, &ni, out , &ni, &wkopt, &lwork, &info);
			lwork = (int)wkopt;
			work = new double[lwork];
			dgels_(trans, &ni, &mj, &nrhs, inputAt, &ni, out, &ni, work, &lwork, &info);
			return Matrix(out, mj, nrhs, false,info);
        }
    }


	Matrix Matrix::operator | (Matrix operand)
	{
        if (isSquareMatrix)
        {
			if (isTridiagonal)
			{
				int bl = operand.mj;
				double* aa = new double[ni - 1];
				double* bb = new double[ni];
				double* cc = new double[ni - 1];
				double* dd = new double[operand.ni*operand.mj];
				for (int i = 0; i<(ni - 1); i++)
				{
					aa[i] = alpha[i];
					bb[i] = beta[i];
					cc[i] = gamma[i];
				}
				bb[ni - 1] = beta[ni - 1];
				for (int i = 0; i<(operand.ni*operand.mj); i++)
				{
					dd[i] = operand.At[i];
				}
				dgtsv_(&ni, &bl, aa, bb, cc, dd, &ni, &info);
				return Matrix(dd, mj, operand.mj, false, info); 
			}
			else
			{
				int k = 0;
				double* out = new double[ni*ni];
				int* ip = new int[ni];
				//the subroutine dgesv changes the matrix input, need to create a copy of the input matrix 
				double* inputAt = new double[ni*mj];
				for (int j = 0; j<(ni*mj); j++)
				{
					inputAt[j] = At[j];
				}
				double* inputBt = new double[operand.ni*operand.mj];
				for (int j = 0; j<(operand.ni*operand.mj); j++)
				{
					inputBt[j] = operand.At[j];
				}
				dgesv_(&ni, &(operand.mj), inputAt, &ni, ip, inputBt, &operand.ni, &info);
				return Matrix(inputBt, mj, operand.mj); 
			} 
        }
        else
        {
			int nrhs = operand.mj; //number of column on the right hand side 
			int lwork = -1;
			double wkopt;
			double* work;
			char trans[] = { "N" };
			double* inputAt = new double[ni*mj];
			for (int j = 0; j<(ni*mj); j++)
			{
				inputAt[j] = At[j];
			}
			double* inputBt = new double[operand.ni*operand.mj];
			for (int j = 0; j<(operand.ni*operand.mj); j++)
			{
				inputBt[j] = operand.At[j];
			} 
			dgels_(trans, &ni, &mj, &nrhs, inputAt, &ni, inputBt, &operand.ni, &wkopt, &lwork, &info);
			lwork = (int)wkopt;
			work = new double[lwork];
			dgels_(trans, &ni, &mj, &nrhs, inputAt, &ni, inputBt, &operand.ni, work, &lwork, &info);
			return Matrix(inputBt, mj, operand.mj,false,info);
        }
	}


	Matrix Matrix::operator ^ (int one_1)
    {
		if (one_1 == -1)
		{
			int k = 0;
			double* Bout = new double[ni*ni];
			for (int i = 0; i<ni; i++)
			{
				for (int j = 0; j<ni; j++)
				{
					if (i == j)
					{
						Bout[k] = 1.0;
					}
					else
					{
						Bout[k] = 0.0;
					}
					k++;
				}
			}
			int* ip = new int[ni];
			//the subroutine dgesv changes the matrix input, need to create a copy of the input matrix 
			double* inputAt = new double[ni*ni];
			for (int j = 0; j<(ni*ni); j++)
			{
				inputAt[j] = At[j];
			}
			dgesv_(&ni, &ni, inputAt, &ni, ip, Bout, &ni, &info);
			return Matrix(Bout, ni, ni, false,info);
        }
        else
        {
			double* Bout = new double[ni*ni];
            int k = 0;
			for (int i = 0; i<ni; i++)
			{
				for (int j = 0; j<ni; j++)
				{
					if (i == j)
					{
						Bout[k] = 1.0;
					}
					else
					{
						Bout[k] = 0.0;
					}
					k++;
				}
			}
			return Matrix(Bout, ni, ni, false);
        } 
    }
#pragma endregion


