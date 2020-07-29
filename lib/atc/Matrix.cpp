#include "Matrix.h"
#include "DenseMatrix.h"

#include "PolynomialSolver.h"

namespace ATC_matrix {



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

  double *pa=A.ptr(), *pb=B.ptr(), *pc=C.ptr();

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
  
  dgetrf_(&m,&m,invA.ptr(),&m,ipiv,&info); // compute LU factorization
  GCK(A,A,info<0,"DenseMatrix::inv() dgetrf error: Argument had bad value.");
  GCK(A,A,info>0,"DenseMatrix::inv() dgetrf error: Matrix not invertible.");
  if (info > 0) {
    delete [] ipiv;
    invA = 0; 
    return invA;
  }

  // LU factorization succeeded
  // Compute 1-norm of original matrix for use with dgecon
  char norm = '1'; // Use 1-norm
  double rcond, anorm, *workc = new double[4*m];
  anorm = dlange_(&norm,&m,&m,A.ptr(),&m,workc);

  // Condition number estimation (warn if bad)
  dgecon_(&norm,&m,invA.ptr(),&m,&anorm,&rcond,workc,iwork,&info);
  GCK(A,A,info<0, "DenseMatrix::inv(): dgecon error: Argument had bad value.");
  GCK(A,A,rcond<1e-14,"DenseMatrix::inv(): Matrix nearly singular, RCOND<e-14");

  // Now determine optimal work size for computation of matrix inverse
  double work_dummy[2] = {0.0,0.0};
  dgetri_(&m, invA.ptr(), &m, ipiv, work_dummy, &lwork, &info);
  GCK(A,A,info<0,"DenseMatrix::inv() dgetri error: Argument had bad value.");
  GCHK(info>0,"DenseMatrix::inv() dgetri error: Matrix not invertible.");
  
  // Work size query succeeded
  lwork = (int)work_dummy[0];  
  double *work = new double[lwork];    // Allocate vector of appropriate size

  // Compute and store matrix inverse
  dgetri_(&m,invA.ptr(),&m,ipiv,work,&lwork,&info);
  GCK(A,A,info<0,"DenseMatrix::inv() dgetri error: Argument had bad value.");
  GCHK(info>0,"DenseMatrix::inv() dgetri error: Matrix not invertible.");

  // Clean-up
  delete [] ipiv;
  delete [] workc;
  delete [] work;
  return invA;
}
//-----------------------------------------------------------------------------
//* returns all eigenvalues & e-vectors of a pair of double precision matrices
//-----------------------------------------------------------------------------

DenseMatrix<double>  eigensystem(const MATRIX& AA, const MATRIX & BB, 
  DenseMatrix<double> & eVals, bool normalize)
{
  DENS_MAT A(AA);  // Make copy of A 
  DENS_MAT B(BB);  
  int m = A.nRows(); // size
  eVals.resize(m,1); // eigenvectors 
  //A.print("A");
  //B.print("B");
  SQCK(A, "DenseMatrix::eigensystem(), matrix not square"); // check matrix is square
  SQCK(B, "DenseMatrix::eigensystem(), matrix not square"); // check matrix is square
  SSCK(A, B, "DenseMatrix::eigensystem(), not same size");// check same size

  // workspace 
  int lwork=-1; //1+(NB+6+2*NMAX)*NMAX) 
  double tmp[1];
  double *work = tmp;
  int liwork = -1; // 3+5*NMAX
  int  itmp[1];
  int *iwork = itmp;

  // Solve the generalized symmetric eigenvalue problem
  // A*x = lambda*B*x (ITYPE = 1)
  // only accesses upper triangle
  char vectors[] = "Vectors", upper[] = "Upper";
  int type = 1, info;
  // query optimal sizes
  dsygvd_(&type,vectors,upper,&m,A.ptr(),&m,B.ptr(),&m,
    eVals.ptr(),work,&lwork,iwork,&liwork,&info);
  // returns optimal sizes  LWOPT = WORK(1), LIWOPT = IWORK(1)
  lwork = int(work[0]);
  liwork = iwork[0];
  work = new double[lwork];
  iwork = new int[liwork]; 
  dsygvd_(&type,vectors,upper,&m,A.ptr(),&m,B.ptr(),&m,
    eVals.ptr(),work,&lwork,iwork,&liwork,&info);
  GCK(A,B,info!=0,"DenseMatrix::eigensystem(), error");
  //eVals.print("e-values");
  //(A.transpose()).print("e-vectors");

  // normalize
  if (normalize) {
    for (int j = 0; j < A.nCols(); j++) {
      double scale = 0.0;
      for (int i = 0; i < A.nRows(); i++) {
        scale += A(i,j)*A(i,j);
      }
      scale = 1.0/sqrt(scale);
      for (int i = 0; i < A.nRows(); i++) {
        A(i,j) *= scale;
      }
    }
    //(A.transpose()).print("normalized e-vectors");
  }  
  delete [] work;
  delete [] iwork;

  //return A.transpose();
  return A; // column storage
}

//-----------------------------------------------------------------------------
//* returns (1-norm) condition number
//-----------------------------------------------------------------------------
double condition_number(const MATRIX& AA)
{
  DenseMatrix<double> eVals, I;
  I.identity(AA.nRows());
  eigensystem(AA, I, eVals); 
// [1] William W. Hager, "Condition Estimates," SIAM J. Sci. Stat. Comput. 5, 1984, 311-316, 1984.
// [2] Nicholas J. Higham and Françoise Tisseur, "A Block Algorithm for Matrix 1-Norm Estimation with an Application to 1-Norm Pseudospectra, "SIAM J. Matrix Anal. Appl., Vol. 21, 1185-1201, 2000.
  double max = eVals.maxabs();
  double min = eVals.minabs();
  return max/min;
}
//-----------------------------------------------------------------------------
//* returns polar decomposition of a square double precision matrix via SVD
//-----------------------------------------------------------------------------
DenseMatrix<double>  polar_decomposition(const MATRIX& AA, 
  DenseMatrix<double> & rotation,
  DenseMatrix<double> & stretch,
  bool leftRotation)
{
  DENS_MAT A(AA);  // Make copy of A 
  SQCK(A, "DenseMatrix::polar_decomposition(), matrix not square"); 
  int m = A.nRows(); // size
  DENS_MAT D(m,1); 
  DENS_MAT U(m,m), VT(m,m); // left and right SVD rotations

  // workspace 
  int lwork=-1; 
  double tmp[1];
  double *work = tmp;
  
  // calculate singular value decomposition A = U D V^T 
  char type[] = "A"; // all columns are returned
  int info;
  // query optimal sizes
  dgesvd_(type,type,&m,&m,A.ptr(),&m,D.ptr(),
    U.ptr(),&m,VT.ptr(),&m,
    work,&lwork,&info); // simple: svd, div&conq: sdd
  lwork = int(work[0]); // returns optimal size  LWOPT = WORK(1)
  work = new double[lwork];
  // compute SVD
  dgesvd_(type,type,&m,&m,A.ptr(),&m,D.ptr(),
    U.ptr(),&m,VT.ptr(),&m,
    work,&lwork,&info);
  
  //GCK(A,B,info!=0,"DenseMatrix::polar_decomposition(), error");
  GCK(A,A,info!=0,"DenseMatrix::polar_decomposition(), error");
  delete [] work;


  //rotation.resize(m,m);
  rotation = U*VT;
  // A = R' U' = (U V^T) (V D V^T)
  stretch.resize(m,m);
//if (leftRotation) { stretch = (VT.transpose())*D*VT; }
  if (leftRotation) { 
    DENS_MAT V = VT.transpose();
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
        stretch(i,j) = V(i,j)*D(j,0);
      }
    }
    stretch = stretch*VT;
  }
  // A = V' R' = (U D U^T) (U V^T)
//else              { stretch = U*D*(U.transpose()); }
  else { 
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
        stretch(i,j) = U(i,j)*D(j,0);
      }
    }
    stretch = stretch*(U.transpose());
  }
  return D; 
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
  dgetrf_(&m,&m,PLUA.ptr(),&m,ipiv,&info);

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

  GCK(A,A,!A.is_size(3,3), "max_eigenvalue only implemented for 3x3");
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

} // end namescape

