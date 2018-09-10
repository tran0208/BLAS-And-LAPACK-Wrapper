//Solving the system A*x = d, 
//where A is a n*n general Matrix (in and out)
//"d" is the right hand side and is the output
extern "C" void dgesv_(const int *n, const int *nrhs, double *A, const int *lda, int *ipiv, double *d, const int *ldb, int* flag);

//Solving the system A*x = b, where A is a tridiagonal matrix, 
//with "a" being the subdiagonal, (in and out)
//"b" being the main diagonal, (in and out) 
//"c" being the super diagonal, (in and out)
//"d" is the right hand side and the output will be stored in this matrix, , (in and out)
extern "C" void dgtsv_(int *n, int *NRHS, double *a, double *b, double *c, double *d, int *LDBp, int *flag);

//Solving the system A*x = d, 
//"A" is a n*m general matrix, (in and out)
//"d" is a n*1 matrix,Note: db also store the output
extern "C" void dgels_(char* trans, int* n, int* m, int* nrhs, double* A, int* lda, double* d, int* ldb, double* work, int* lwork, int* flag);
 
//Solving the system min |A*x = d|_{2} with constraints Hx = f, 
//"A" is a n*m general matrix,  (in and out)
//"H" is a p*m, (in and out)
//"d" is n*1 matrix 
//"f" is a p*1 matrix
//"x" is a m*1 matrix
extern "C" void dgglse_(int* n, int* m, int* p, double* A, int* lda, double* H, int* ldb, double* d, double* f, double* x, double* work, int* lwork, int* flag); 
 
//min |y|_{2} s.t x,y with constraints d = A*x + B*y
//"A" is a n*m general matrix, (in and out)
//"B" is a n*p matrix (in and out)
//"d" is n*1 matrix
//"x" is a m*1 matrix
//"y" is a p*1 matrix
extern "C" void dggglm_(int* n, int* m, int* p, double* A, int* lda, double* B, int* ldb, double* d, double* x, double* y, double* work, int* lwork, int* flag); 
