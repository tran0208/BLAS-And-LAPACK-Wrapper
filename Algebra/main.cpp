#include <stdio.h>
#include <conio.h>
#include <Lapack.h>
#include "Matrix.h"

int main(void)
{ 

	printf("\n"); 
	double A1[] = { 5 ,  7 ,  6 , 5 ,
                    7 , 10 ,  8 , 7 ,
                    6 ,  8 , 10 , 9 ,
                    5 ,  7 ,  9 , 10 };
	double b1[] = { 23 , 32 , 33 , 31 };
    Matrix P1(A1,4);   
    Matrix example1 = (5*P1+P1*15) | b1;

	if (example1.info == 0)
	{
		printf("X1 = \n"); 
        example1.printMatrix();
	}
	else
	{
		fprintf(stderr, "Info = %d\n", example1.info);
	}
//Matlab code:
//A1 = [5, 7, 6, 5; 7, 10, 8, 7; 6, 8, 10, 9; 5, 7, 9, 10];
//b1 = [23, 32, 33, 31]';
//example1 = inv(5*A1 + A1*15)*b1
//		example1 =
//		0.05000
//		0.05000
//		0.05000
//		0.05000

	printf("\n\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\n\n");


	double A2[] = { 0 , 1 , 2 ,
		                3 , 4 , 5 , 
                        6 , 7 , 0  };
	double B2[] = { 1 , 2 , 3 }; 
	Matrix P2(A2, 3);
    Matrix example2 = (P2*P2*P2) | B2;  
	if (example2.info == 0)
	{
		printf("X2 = \n");
		example2.printMatrix();
	}
	else
	{
		fprintf(stderr, "Info = %d\n", example2.info );
	} 
	printf("\n\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\n\n");
	//Matlab code:
	//Alpha1 = [0, 1, 2; 3, 4, 5; 6, 7, 0];
	//Beta1 = [1; 2; 3]; 
	//example2 = inv(Alpha1 * Alpha1*Alpha1)*Beta1 
	//   example2 = 
	//	-3.08796296296297
	//	2.69444444444445
	//	- 0.569444444444446




	double a[] = { -1, -1, -1 };
	double b[] = { 0, 0, 0, 0};
	double c[] = { 1, 1, 1 };
	double d[] = { 1, 2, 3, 4 };
	Matrix Q3(d, 4,1); 
	Matrix P3(a,b,c,4); 
	Matrix example3 = P3 | d;
	Matrix example4 = P3 | Q3;

	if (example3.info == 0 && example3.info == 0)
	{
		printf("X3 = \n");
		example3.printMatrix();

		printf("X4 = \n");
		example4.printMatrix();
	}
	else
	{
		fprintf(stderr, "Info = %d\n", example3.info);
		fprintf(stderr, "Info = %d\n", example4.info);
	}


	printf("\n\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\n\n");

	double A4[30] = {  -74 ,  80 ,  18 , -11 , -4 ,
                        14 , -69 ,  21 ,  28 , 0 , 
                        66 , -72 ,  -5 ,   7 , 1 ,
                       -12 ,  66 , -30 , -23 , 3 ,
                         3 ,   8 ,  -7 ,  -4 , 1 ,
                         4 , -12 ,   4 ,   4 , 0};
	double b4[6] = { 51 , -61 , -56 , 69 , 10 , -12 };
	Matrix P4(A4,6,5); 
	Matrix example5 = P4 | b4;
	Matrix Q4(b4, 6, 1);
	Matrix example6 = P4 | Q4;

	printf("\n X5 = \n ");
	example5.printMatrix();

	printf("\n X6 = \n ");
	example6.printMatrix();

	printf("\n\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\n\n");
 
	double *x5 = new double[4];
	double *y5 = new double[3];
	double A5[] = { 1 ,  2 ,  1 , 4  ,
                   -1 ,  1 ,  1 , 1  , 
                   -1 , -2 , -1 , 1  , 
                   -1 ,  2 , -1 , -1 , 
                    1 ,  1 ,  1 , 2} ;
	double B5[] = { 1 ,  2 , 2 ,
                   -1 ,  1 , -2 , 
                    3 ,  1 , 6 ,
                    2 , -2 , 4 ,
                    1 , -1 , 2 };
	double d5[] = { 7.99 , 0.98 , -2.98 ,3.04 , 4.02 };
    Matrix P5(A5, B5, d5, 5, 4, 3); 
	int info = P5.solve(x5,y5);
	if (info == 0)
	{
		printf("X7 = [ %f , %f , %f , %f ] \n", x5[0], x5[1], x5[2], x5[3]);
		printf("Y7 = [ %f , %f , %f  ] \n", y5[0], y5[1], y5[2] );
	}
	else
	{
		fprintf(stderr, "dggglm_ fails %d\n", info);
	} 

	printf("\n\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\n\n");
 
	double A6[] = { 1 ,  1 , 1 , 
                    1 ,  3 , 1 , 
                    1 , -1 , 1 , 
                    1 ,  1 , 1};
	double H6[] = { 1 , 1 , 1 ,
		            1 , 1 , -1 };
	double D6[] = { 1 , 
                    2 , 
                    3 , 
                    4 };
	double F6[] = { 7 ,
		            4 };
	double *X6 = new double[3];
	Matrix P6(A6, D6, H6, F6, 4, 3, 2);
    P6.solve(X6);
	printf("\n X8 = \n ");
	Matrix::printMatrix(X6,3,1);

	printf("\n\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\n\n");
   
	double G[] = { 1 , 1 , 1 ,
		           1 , 1 , 1 ,
		           1 , 1 , 1 ,
		           1 , 1 , 1 ,
		           1 , 1 , 1 };
    Matrix P7(G,5,3);
    Matrix W1 = P7*20;
	printf("\n X9 = \n ");
    W1.printMatrix();
	printf("-------------------------------------------------------------------"); 
	Matrix W2 = 40*P7;
	printf("\n X10 = \n ");
	W2.printMatrix(); 
    
    
	printf("\n\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\n\n");

	printf("Press any key to exit\n");
	_getch(); 
	return 0;


}
