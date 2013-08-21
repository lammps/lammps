#ifndef MATRIXDEF_H 
#define MATRIXDEF_H

/******************************************************************************
* Common definitions for Matrix and Vector classes
* This header file contains macros and inline functions needed for the matrix 
* classes.  All error checking should be defined here as a macro so that it is
* neatly disabled when ATC_PRINT_DEBUGGING is not defined
******************************************************************************/

/******************************************************************************
* Headers and namespaces used by Matrix and Vector classes
******************************************************************************/
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <cstring>
#include <string>
#include <iomanip>
#include <cmath>
#include "Utility.h"

namespace ATC_matrix {

/******************************************************************************
* Typedefs used by Matrix and Vector classes
******************************************************************************/
//* @typedef INDEX
//* @brief indexing type (default: unsigned) for matrix classes
// can switch typedef back to unsigned to be more precise, but will cause warnings everywhere
//typedef unsigned INDEX;
typedef int INDEX;
//* @typedef CLONE_TYPE
//* @brief dimension of matrix to clone
enum CLONE_TYPE { CLONE_ROW=0, CLONE_COL=1, CLONE_DIAG=2 };
//* @struct TRIPLET
//* @brief Triplet output entity
template <typename T> 
struct TRIPLET { TRIPLET<T>(int _i=0, int _j=0, T _v=T(0)) : i(_i), j(_j), v(_v) {}
                 INDEX i, j; T v; };

/******************************************************************************
* Definitions for row/column major storage
******************************************************************************/
#define COL_STORAGE  /* <--- comment out this line for row-major storage*/
#ifdef COL_STORAGE
#define DATA(i,j) _data[(i)+_nRows*(j)]
#else
#define  ROW_STORAGE
#define DATA(i,j) _data[(i)*_nCols+(j)]
#endif

/******************************************************************************
* error checking macros
*  MICK:  checks if index (i,j) is in range                       MATRIX ONLY
*  VICK:  checks if index (i) is in range                         VECTOR ONLY
*  MICM:  checks if index (i,j) is in range, displays message     MATRIX ONLY
*  VICM:  checks if index (i) is in range, displays message       VECTOR ONLY
*  SQCK:  checks if matrix is square, displays message            MATRIX ONLY
*  SSCK:  checks if a has the same size as b                      VECTOR/MATRIX
*  GCK:   generic two object check, checks if c is true           VECTOR/MATRIX
*  GCHK:  generic check, checks if c is true                      ANYTHING
******************************************************************************/
#define ERROR_FOR_BACKTRACE /**/
#define MICK(i,j)     /**/
#define VICK(i)       /**/
#define MICM(i,j,m)   /**/
#define VICM(i,m)     /**/
#define SQCK(a,m)     /**/
#define SICK(a,b,m)   /**/
#define SSCK(a,b,m)   /**/
#define GCK(a,b,c,m)  /**/
#define GCHK(c,m)     /**/

// the following two convert __LINE__ to a string
#define STRING2(x) #x
#define STRING(x) STRING2(x)
// prints file and line number for error messages
#define ERROR(x) __FILE__":"STRING(__LINE__)" "x

/******************************************************************************
* BLAS and LAPACK definitions
******************************************************************************/
#ifdef MKL
#include "mkl.h"
#define dgemv_  dgemv
#define dgemm_  dgemm
#define dgetrf_ dgetrf
#define dgetri_ dgetri
#define dgecon_ dgecon
#define dlange_ dlange
#define dsygvd_ dsygvd
#define dgesvd_ dgesvd
#define dgesdd_ dgesdd
#else
extern "C"
{
extern void dgemv_(char*,int*,int*,double*,const double*,int*,const double*,int *,double*,double*,int*);
extern void dgemm_(char*,char*,int*,int*,int*,double*,const double*,int*,const double*,int*,double*,double*,int*);
extern void dgetrf_(int*,int*,double*,int*,int*,int*);
extern void dgetri_(int*,double*,int*,int*,double*,int*,int*);
extern void dgecon_(char*,int*,double*,int*,double*,double*,double*,int*,int*);
extern double dlange_(char*,int*,int*,const double*,int*,double*);
extern double dsygvd_(int*,char*,char*,int*,double*,int*,double*,int*,double*,double*,int*,int*,int*,int*);
extern double dgesvd_(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);
extern double dgesdd_(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);
};
#endif

// forward declarations of matrix and vector classes
template<typename T> class Matrix;
template<typename T> class DenseMatrix;
template<typename T> class ParDenseMatrix;
template<typename T> class SparseMatrix;
template<typename T> class ParSparseMatrix;
template<typename T> class SparseVector;
template<typename T> class DiagonalMatrix;
template<typename T> class ParDiagonalMatrix;
template<typename T> class Vector;
template<typename T> class DenseVector;
template<typename T> class CloneVector;
template<typename T> class WrapMatrix;
template<typename T> class WrapVector;

//* forward declaration of operations
//@{
template<class T> DenseVector<T>  operator*(const Matrix<T> &M, const SparseVector<T> &v);
template<class T> DenseVector<T>  operator*(const SparseVector<T> &v, const Matrix<T> &M);
template<class T> SparseVector<T> operator*(const SparseMatrix<T> &M, const SparseVector<T> &v);
template<class T> SparseVector<T> operator*(const SparseVector<T> &v, const SparseMatrix<T> &M);
template<class T> DenseVector<T>  operator*(const SparseMatrix<T> &A, const Vector<T>& x);
template<class T> DenseVector<T>  operator*(const Vector<T> &A, const SparseMatrix<T>& x);
template<class T> DenseMatrix<T>  operator*(const SparseMatrix<T> &A, const Matrix<T>& D);
template<class T> SparseMatrix<T> operator*(const SparseMatrix<T> &A, const DiagonalMatrix<T>& D);
template<class T> SparseMatrix<T> operator*(const SparseMatrix<T> &A, const SparseMatrix<T> &B);
template<class T> T dot(const SparseVector<T> &a, const SparseVector<T> &b);
//@}

template<class T> CloneVector<T> column(Matrix<T> &c, INDEX i) {
  return CloneVector<T>(c, CLONE_COL, i); 
}
template<class T> CloneVector<T> row(Matrix<T> &c, INDEX i) {
  return CloneVector<T>(c, CLONE_ROW, i);
}
template<class T> CloneVector<T> diagonal(Matrix<T> &c) {
  return CloneVector<T>(c, CLONE_DIAG); 
}
template<class T> const CloneVector<T> column(const Matrix<T> &c, INDEX i) {
  return CloneVector<T>(c, CLONE_COL, i); 
}
template<class T> const CloneVector<T> row(const Matrix<T> &c, INDEX i) {
  return CloneVector<T>(c, CLONE_ROW, i); 
}
template<class T> const CloneVector<T> diagonal(const Matrix<T> &c) { 
  return CloneVector<T>(c, CLONE_DIAG); 
}
template<class T> const SparseMatrix<T> *sparse_cast(const Matrix<T> *m);
template<class T> const DiagonalMatrix<T> *diag_cast(const Matrix<T> *m);
template<class T> void copy_sparse_to_matrix(const SparseMatrix<T> *s, Matrix<T> &m);
template<typename T> DenseMatrix<T> operator*(const DiagonalMatrix<T>& A, const Matrix<T> &B);
template<typename T> DenseMatrix<T> operator*(const Matrix<T> &B, const DiagonalMatrix<T>& A);

// double precision shortcuts
typedef Matrix<double>         MATRIX;          // matrix of double
typedef Vector<double>         VECTOR;          // vector of double
typedef DenseMatrix<double>    DENS_MAT;        // dense matrix of double type
typedef ParDenseMatrix<double> PAR_DENS_MAT;    // parallel dense matrix of doubles
typedef CloneVector<double>    CLON_VEC;        // cloned vector of double type
typedef DenseVector<double>    DENS_VEC;        // dense vector of double type
typedef DiagonalMatrix<double> DIAG_MAT;        // diagonal matrix of double type
typedef ParDiagonalMatrix<double> PAR_DIAG_MAT; // diagonal matrix of double type
typedef SparseMatrix<double>   SPAR_MAT;        // sparse matrix of double type
typedef ParSparseMatrix<double> PAR_SPAR_MAT;   // parallel sparse matrix of doubles
typedef SparseVector<double>   SPAR_VEC;        // sparse matrix of double type
typedef std::vector<DenseMatrix<double> > DENS_MAT_VEC;
typedef std::vector<SparseMatrix<double> * > SPAR_MAT_VEC;

// int containers
typedef DenseMatrix<int>       INT_ARRAY; // to become vector<int> or Array2D
//typedef SparseMatrix<int>      SPAR_INT_ARRAY; // to become ?
typedef DenseVector<int>       INT_VECTOR; // to become vector<int> or Array

// forward declaration of error messages
template<typename T> void ierror(const Matrix<T> &a, const char *FILE, int LINE, INDEX i, INDEX j=0);
template<typename T> void ierror(const Matrix<T> &a, INDEX i, INDEX j, const std::string m);
template<typename T> void merror(const Matrix<T> &a, const Matrix<T> &b, const std::string m);
inline void gerror(const std::string m) { std::cout<<"Error: "<<m<<"\n"; ERROR_FOR_BACKTRACE ; exit(EXIT_FAILURE); }

} // end namespace
#endif
