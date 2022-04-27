#ifndef PARSPARSEMATRIX_H
#define PARSPARSEMATRIX_H

#include "mpi.h"
#include "MPI_Wrappers.h"
#include "SparseMatrix.h"
#include "DiagonalMatrix.h"
#include <algorithm>

namespace ATC_matrix {

  /**
   *  @class  ParSparseMatrix
   *  @brief  Parallelized version of SparseMatrix class.
   *
   *  ParSparseMatrix<double>::MultMv is used in LinearSolver, which is then
   *  used in NonLinearSolver, PoissonSolver, and SchrodingerSolver. These
   *  parallelized solvers are used in the following locations:
   *
   *  - LinearSolver
   *    - ExtrinsicModelDriftDiffusion.cpp (lines 511 and 522)
   *    - AtomicRegulator.cpp (line 926)
   *    - TransferLibrary.cpp (lines 72 and 260)
   *  - PoissonSolver
   *    - ExtrinsicModelDriftDiffusion.cpp (line 232)
   *  - SchrodingerSolver
   *    - ExtrinsicModelDriftDiffusion.cpp (line 251)
   *  - SliceSchrodingerSolver
   *    - ExtrinsicModelDriftDiffusion.cpp (line 246)
   */

  template <typename T>
  class ParSparseMatrix : public SparseMatrix<T>
  {
    public:
      ParSparseMatrix(MPI_Comm comm, INDEX rows = 0, INDEX cols = 0)
        : SparseMatrix<T>(rows, cols), _comm(comm){}

      ParSparseMatrix(MPI_Comm comm, const SparseMatrix<T> &c)
        : SparseMatrix<T>(c), _comm(comm){}

      ParSparseMatrix(MPI_Comm comm, INDEX* rows, INDEX* cols, T* vals,
                      INDEX size, INDEX nRows, INDEX nCols, INDEX nRowsCRS)
        : SparseMatrix<T>(rows, cols, vals, size, nRows, nCols, nRowsCRS)
         ,_comm(comm){}

      ParSparseMatrix(MPI_Comm comm)
        : SparseMatrix<T>(), _comm(comm){}

    virtual void operator=(const SparseMatrix<T> &source)
    {
      copy(source);
    }

    template<typename U>
    friend void ParMultAB(MPI_Comm comm, const SparseMatrix<U>& A,
                          const Matrix<U>& B, DenseMatrix<U>& C);

    private:
        MPI_Comm _comm;
  };

  template <>
  class ParSparseMatrix<double> : public SparseMatrix<double>
  {
  public:
    // All the same constructors as for SparseMatrix
    ParSparseMatrix(MPI_Comm comm, INDEX rows = 0, INDEX cols=0);
    ParSparseMatrix(MPI_Comm comm, const SparseMatrix<double> &c);
    ParSparseMatrix(MPI_Comm comm, INDEX* rows, INDEX* cols, double* vals, INDEX size,
        INDEX nRows, INDEX nCols, INDEX nRowsCRS);

    // Parallel sparse matrix multiplication functions
    void MultMv(const Vector<double>& v, DenseVector<double>& c) const;
    DenseVector<double> transMat(const Vector<double>& v) const;
    void MultAB(const Matrix<double>& B, DenseMatrix<double>& C) const;
    DenseMatrix<double> transMat(const DenseMatrix<double>& B) const;
    DenseMatrix<double> transMat(const SparseMatrix<double>& B) const;

    virtual void operator=(const SparseMatrix<double> &source);

    template<typename U>
    friend void ParMultAB(MPI_Comm comm, const SparseMatrix<U>& A, const Matrix<U>& B, DenseMatrix<U>& C);

  private:
    void partition(ParSparseMatrix<double>& A_local) const;
    void finalize();
    MPI_Comm _comm;
  };

  // The SparseMatrix versions of these functions will call the correct
  //   MultMv/MultAB:
  // DenseVector<double> operator*(const ParSparseMatrix<double> &A, const Vector<double> &v);
  // DenseVector<double> operator*(const Vector<double> &v, const ParSparseMatrix<double> &A);
  // DenseMatrix<double> operator*(const ParSparseMatrix<double> &A, const Matrix<double> &B);


  template<typename T>
  void ParMultAB(MPI_Comm comm, const SparseMatrix<T>& A, const Matrix<T>& B, DenseMatrix<T>& C)
  {
    SparseMatrix<T>::compress(A);

    INDEX M = A.nRows(), N = B.nCols();
    if (!C.is_size(M, N))
    {
      C.resize(M, N);
      C.zero();
    }

    // Temporarily put fields into a ParSparseMatrix for distributed multiplication
    ParSparseMatrix<T> Ap(comm);
    Ap._nRows       = A._nRows;
    Ap._nCols       = A._nCols;
    Ap._size        = A._size;
    Ap._nRowsCRS    = A._nRowsCRS;
    Ap._val         = A._val;
    Ap._ja          = A._ja;
    Ap._ia          = A._ia;
    Ap.hasTemplate_ = A.hasTemplate_;

    // MultAB calls compress(), but we hope that does nothing because we just
    // compressed A. If it did something, it might mess up other members
    // (e.g. _tri).
    Ap.MultAB(B, C);

    // We're not changing the matrix's values, so we can justify calling A const.
    SparseMatrix<T> &Avar = const_cast<SparseMatrix<T> &>(A);
    Avar._nRows       = Ap._nRows;
    Avar._nCols       = Ap._nCols;
    Avar._size        = Ap._size;
    Avar._nRowsCRS    = Ap._nRowsCRS;
    Avar._val         = Ap._val;
    Avar._ja          = Ap._ja;
    Avar._ia          = Ap._ia;
    Avar.hasTemplate_ = Ap.hasTemplate_;

    // Avoid catastrophe
    Ap._val = nullptr;
    Ap._ja  = nullptr;
    Ap._ia  = nullptr;
  }


  /*SparseMatrix<double> operator*(const ParSparseMatrix<double> &A, const SparseMatrix<double> &B);

  SparseMatrix<double> operator*(const ParSparseMatrix<double> &A, const DiagonalMatrix<double> &B);
  */
} // end namespace
#endif

