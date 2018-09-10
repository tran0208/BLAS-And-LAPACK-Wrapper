#pragma once
  
//x := a*x
extern "C" void dscal_(int* n, double* a, double*x, int* incx);
 
//y = y + a*x 
extern "C" void daxpy_(int* n, double* a, double* x, int* incx, double* y, int* incy);    

// matrix matrix multiply
//C := alpha * A*B + beta * C 
extern "C" void dgemm_(char* transa, char* transb, int* n, int* m, const int* p, double* alpha, double* A, const int* lda, double* B, int* ldb, double* beta, double* C, const int* ldc);
  
//matrix vector multiply
//y := alpha*A*x + beta*y
extern "C" void dgemv_(char* trans, int* n, int* m, double* alpha, double* A, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
