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

	double A5[30] = {  -74 ,  80 ,  18 , -11 , -4 ,
                        14 , -69 ,  21 ,  28 , 0 , 
                        66 , -72 ,  -5 ,   7 , 1 ,
                       -12 ,  66 , -30 , -23 , 3 ,
                         3 ,   8 ,  -7 ,  -4 , 1 ,
                         4 , -12 ,   4 ,   4 , 0};
	double b5[6] = { 51 , -61 , -56 , 69 , 10 , -12 };
	Matrix P5(A5,6,5); 
	Matrix example5 = P5 | b5;
	Matrix Q5(b5, 6, 1);
	Matrix example6 = P5 | Q5;

	printf("\n X5 = \n ");
	example5.printMatrix();

	printf("\n X6 = \n ");
	example6.printMatrix();

	printf("\n\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\n\n");
 
	double *x7 = new double[4];
	double *y7 = new double[3];
	double A7[] = { 1 ,  2 ,  1 , 4  ,
                   -1 ,  1 ,  1 , 1  , 
                   -1 , -2 , -1 , 1  , 
                   -1 ,  2 , -1 , -1 , 
                    1 ,  1 ,  1 , 2} ;
	double B7[] = { 1 ,  2 , 2 ,
                   -1 ,  1 , -2 , 
                    3 ,  1 , 6 ,
                    2 , -2 , 4 ,
                    1 , -1 , 2 };
	double d7[] = { 7.99 , 0.98 , -2.98 ,3.04 , 4.02 };
    Matrix P7(A7, B7, d7, 5, 4, 3); 
	int info7 = P7.solve(x7,y7);
	if (info7 == 0)
	{
		printf("X7 = [ %f , %f , %f , %f ] \n", x7[0], x7[1], x7[2], x7[3]);
		printf("Y7 = [ %f , %f , %f  ] \n", y7[0], y7[1], y7[2] );
	}
	else
	{
		fprintf(stderr, "dggglm_ fails %d\n", info7);
	} 

	printf("\n\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\n\n");
 
	double A8[] = { 1 ,  1 , 1 , 
                    1 ,  3 , 1 , 
                    1 , -1 , 1 , 
                    1 ,  1 , 1};
	double H8[] = { 1 , 1 , 1 ,
		            1 , 1 , -1 };
	double D8[] = { 1 , 
                    2 , 
                    3 , 
                    4 };
	double F8[] = { 7 ,
		            4 };
	double *x8 = new double[3];
	Matrix P8(A8, D8, H8, F8, 4, 3, 2);
    int info8 = P8.solve(x8);
	printf("\n X8 = \n ");
	Matrix::printMatrix(x8,3,1);

	printf("\n\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	printf("\n\n");
   

	printf("Press any key to exit\n");
	_getch(); 
	return 0;


}
