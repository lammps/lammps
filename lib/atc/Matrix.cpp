#include "DenseMatrix.h"
#include "Solver.h"

//-----------------------------------------------------------------------------
//* performs a matrix-matrix multiply with optional transposes BLAS version
// C = b*C + a*A*B
//-----------------------------------------------------------------------------
void MultAB(const MATRIX &A, const MATRIX &B, DENS_MAT &C, 
            const bool At, const bool Bt, double a, double b)
{ 
  static char t[2] = {'N','T'};
  char *ta=t+At, *tb=t+Bt;
  int sA[2] = {A.nRows(), A.nCols()};  // sizes of A
  int sB[2] = {B.nRows(), B.nCols()};  // sizes of B
  
  GCK(A, B, sA[!At]!=sB[Bt], "MultAB<double>: matrix-matrix multiply");
  if (!C.is_size(sA[At],sB[!Bt]))
  {
    C.resize(sA[At],sB[!Bt]);
    if (b != 0.0) C.zero();
  }
  // get pointers to the matrix sizes needed by BLAS
  int *M = sA+At;  // # of rows in op[A]  (op[A] = A' if At='T' else A)
  int *N = sB+!Bt; // # of cols in op[B]
  int *K = sA+!At; // # of cols in op[A] or # of rows in op[B]

  double *pa=A.get_ptr(), *pb=B.get_ptr(), *pc=C.get_ptr();

#ifdef COL_STORAGE
  dgemm_(ta, tb, M, N, K, &a, pa, sA, pb, sB, &b, pc, M); 
#else
  dgemm_(tb, ta, N, M, K, &a, pb, sB+1, pa, sA+1, &b, pc, N); 
#endif
  
}

//-----------------------------------------------------------------------------
//* returns the inverse of a double precision matrix
//-----------------------------------------------------------------------------
DenseMatrix<double> inv(const MATRIX& A)
{
  SQCK(A, "DenseMatrix::inv(), matrix not square"); // check matrix is square
  DENS_MAT invA(A);  // Make copy of A to invert

  // setup for call to BLAS
  int m, info, lwork=-1;
  m = invA.nRows();

  int *ipiv = new int[m<<1]; // need 2m storage
  int *iwork=ipiv+m;

  dgetrf_(&m,&m,invA.get_ptr(),&m,ipiv,&info); // compute LU factorization
  GCK(A,A,info<0,"DenseMatrix::inv() dgetrf error: Argument had bad value.");
  GCK(A,A,info>0,"DenseMatrix::inv() dgetrf error: Matrix not invertable.");

  // LU factorization succeeded
  // Compute 1-norm of original matrix for use with dgecon
  char norm = '1'; // Use 1-norm
  double rcond, anorm, *workc = new double[4*m];
  anorm = dlange_(&norm,&m,&m,A.get_ptr(),&m,workc);

  // Condition number estimation (warn if bad)
  dgecon_(&norm,&m,invA.get_ptr(),&m,&anorm,&rcond,workc,iwork,&info);
  GCK(A,A,info<0, "DenseMatrix::inv(): dgecon error: Argument had bad value.");
  GCK(A,A,rcond<1e-14,"DenseMatrix::inv(): Matrix nearly singular, RCOND<e-14");

  // Now determine optimal work size for computation of matrix inverse
  double work_dummy[2] = {0.0,0.0};
  dgetri_(&m, invA.get_ptr(), &m, ipiv, work_dummy, &lwork, &info);
  GCK(A,A,info<0,"DenseMatrix::inv() dgetri error: Argument had bad value.");
  GCK(A,A,info>0,"DenseMatrix::inv() dgetri error: Matrix not invertable.");
  
  // Work size query succeded
  lwork = (int)work_dummy[0];  
  double *work = new double[lwork];    // Allocate vector of appropriate size

  // Compute and store matrix inverse
  dgetri_(&m,invA.get_ptr(),&m,ipiv,work,&lwork,&info);
  GCK(A,A,info<0,"DenseMatrix::inv() dgetri error: Argument had bad value.");
  GCK(A,A,info>0,"DenseMatrix::inv() dgetri error: Matrix not invertable.");

  // Clean-up
  delete [] ipiv;
  delete [] workc;
  delete [] work;
  return invA;
}
//-----------------------------------------------------------------------------
//* computes the determinant of a square matrix by LU decomposition (if n>3)
//-----------------------------------------------------------------------------
double det(const MATRIX& A)
{
  static const double sign[2] = {1.0, -1.0};
  SQCK(A, "Matrix::det(), matrix not square"); // check matrix is square
  int m = A.nRows();
  switch (m)   // explicit determinant for small matrix sizes
  {
    case 1: return A(0,0);
    case 2: return A(0,0)*A(1,1)-A(0,1)*A(1,0);
    case 3:
      return A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1))
           + A(0,1)*(A(1,2)*A(2,0)-A(1,0)*A(2,2))
           + A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
    default: break; 
  }
  // First compute LU factorization
  int info, *ipiv = new int[m];
  double det = 1.0;
  DENS_MAT PLUA(A);
  dgetrf_(&m,&m,PLUA.get_ptr(),&m,ipiv,&info);

  GCK(A,A,info>0,"Matrix::det() dgetrf error: Bad argument value.");
  if (!info)  // matrix is non-singular
  {
    // Compute det(A) = det(P)*det(L)*det(U) = +/-1 * det(U)
    int i, OddNumPivots;

    det = PLUA(0,0);
    OddNumPivots = ipiv[0]!=(1);
    for(i=1; i<m; i++) 
    { 
      det *= PLUA(i,i); 
      OddNumPivots += (ipiv[i]!=(i+1)); // # pivots even/odd
    }
    det *= sign[OddNumPivots&1];  
  }
  delete [] ipiv;  // Clean-up
  return det;
}
//-----------------------------------------------------------------------------
//* Returns the maximum eigenvalue of a matrix.
//-----------------------------------------------------------------------------
double max_eigenvalue(const Matrix<double>& A)
{
  GCK(A,A,!is_size(3,3), "max_eigenvalue only implimented for 3x3");
  const double c0 = det(A);
  const double c1 = A(1,0)*A(0,1) + A(2,0)*A(0,2) + A(1,2)*A(2,1) 
                  - A(0,0)*A(1,1) - A(0,0)*A(2,2) - A(1,1)*A(2,2);
  const double c2 = trace(A);
  double c[4] = {c0, c1, c2, -1.0}, x[3];
  int num_roots = ATC::solve_cubic(c, x);
  double max_root = 0.0;
  for (int i=0; i<num_roots; i++) max_root = std::max(x[i], max_root);
  return max_root;
}

