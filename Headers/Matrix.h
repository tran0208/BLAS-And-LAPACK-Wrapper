#pragma once
  

class Matrix
{
	private:
		//parameters for tridiagonal matrix : subdiagonal
		double* alpha;
		//parameters for tridiagonal matrix : main digonal
		double* beta;
		//parameters for tridiagonal matrix : super diagonal
		double* gamma;
 

        double* Ht;
		double* Bt;
		double* f;
		double* d;

        bool isSquareMatrix = false;
		bool isTridiagonal = false;
  
	public:  
		//This matrix store entries in column format
		double* At;
        //Consider a matrix stored in the format of row, this matrix has ni number of rows
		int ni;
		//Consider a matrix stored in the format of row, this matrix has mj number of columns
		int mj;

		int pk;

		int info;

		//Constructor
		Matrix();
		//~Matrix();
		Matrix(double[], int );
		Matrix(double[],int,int);
		Matrix(double[], double[], double[], int);
		Matrix(double[], double[], double[], double[], int, int, int);
		Matrix(double[], double[], double[], int, int, int);
		Matrix(double[], int, int,bool,int);
		Matrix(double[], int, int, bool);

		int solve(double[]);
		int solve(double[],double[]);
 

		static void printMatrix(double input[], int nd, int md);
	    void printMatrix(void);
 

		Matrix operator + (Matrix);
		Matrix operator - (Matrix);
		Matrix operator * (Matrix);
		Matrix operator * (double[]);
	    Matrix operator * (double);
		friend Matrix operator * (double, Matrix);
		Matrix operator | (double[]);
		Matrix operator | (Matrix);
		Matrix operator ^ (int);


};






