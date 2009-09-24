#include "DenseVector.h"
///////////////////////////////////////////////////////////////////////////////
//* performs a matrix-vector multiply with optional transposes BLAS version
void MultMv(const Matrix<double> &A, const Vector<double> &v, DenseVector<double> &c,
            const bool At, double a, double b)
{
  static char t[2] = {'N','T'};
  char *ta=t+At;
  int sA[2] = {A.nRows(), A.nCols()};  // sizes of A
  int sV[2] = {v.size(), 1};           // sizes of v
  
  GCK(A, v, sA[!At]!=sV[0], "MultAB<double>: matrix-vector multiply");
  if (c.size() != sA[At])
  {
    c.resize(sA[At]); // set size of C to final size
    c.zero();
  }
  // get pointers to the matrix sizes needed by BLAS
  int *M = sA+At;  // # of rows in op[A]  (op[A] = A' if At='T' else A)
  int *N = sV+1;   // # of cols in op[B]
  int *K = sA+!At; // # of cols in op[A] or # of rows in op[B]

  double *pa=A.get_ptr(), *pv=v.get_ptr(), *pc=c.get_ptr();

#ifdef COL_STORAGE
  dgemm_(ta, t, M, N, K, &a, pa, sA, pv, sV, &b, pc, M); 
#else
  dgemm_(t, ta, N, M, K, &a, pv, sV+1, pa, sA+1, &b, pc, N); 
#endif
}
