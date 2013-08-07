#ifndef PARDENSEMATRIX_H
#define PARDENSEMATRIX_H

#include "MatrixDef.h"
#include "DenseMatrix.h"
#include "DenseVector.h"
#include "MPI_Wrappers.h"

#include "ATC_Error.h"
using ATC::ATC_Error;

#include <algorithm>
#include <sstream>

using namespace MPI_Wrappers;

namespace ATC_matrix {

  /**
   *  @class  ParDenseMatrix 
   *  @brief  Parallelized version of DenseMatrix class.
   */

  template <typename T>
  class ParDenseMatrix : public DenseMatrix<T> {
  public:
    MPI_Comm _comm;

    ParDenseMatrix(MPI_Comm comm, INDEX rows=0, INDEX cols=0, bool z=1)
      : DenseMatrix<T>(rows, cols, z), _comm(comm) {}
    ParDenseMatrix(MPI_Comm comm, const DenseMatrix<T>& c)
      : DenseMatrix<T>(c), _comm(comm) {}
    ParDenseMatrix(MPI_Comm comm, const SparseMatrix<T>& c)
      : DenseMatrix<T>(c), _comm(comm) {}
    ParDenseMatrix(MPI_Comm comm, const Matrix<T>& c)
      : DenseMatrix<T>(c), _comm(comm) {}

    //////////////////////////////////////////////////////////////////////////////
    //* performs a matrix-vector multiply
    void ParMultMv(const Vector<T> &v,
             DenseVector<T> &c, const bool At, T a, T b)
    {
      // We can't generically support parallel multiplication because the data
      // types must be specified when using MPI
      MultMv(*this, v, c, At, a, b);
    }

  };

  template<>
  class ParDenseMatrix<double> : public DenseMatrix<double> {
  public:
    MPI_Comm _comm;

    ParDenseMatrix(MPI_Comm comm, INDEX rows=0, INDEX cols=0, bool z=1)
      : DenseMatrix<double>(rows, cols, z), _comm(comm) {}
    ParDenseMatrix(MPI_Comm comm, const DenseMatrix<double>& c)
      : DenseMatrix<double>(c), _comm(comm) {}
    ParDenseMatrix(MPI_Comm comm, const SparseMatrix<double>& c)
      : DenseMatrix<double>(c), _comm(comm) {}
    ParDenseMatrix(MPI_Comm comm, const Matrix<double>& c)
      : DenseMatrix<double>(c), _comm(comm) {}


    void ParMultMv(const Vector<double> &v, DenseVector<double> &c,
        const bool At, double a, double b) const
    {
      // We don't support parallel vec-Mat multiplication yet
      if (At) {
        MultMv(*this, v, c, At, a, b);
        return;
      }

      const INDEX nRows = this->nRows();
      const INDEX nCols = this->nCols();

      if (c.size() != nRows) {
        c.resize(nRows);             // set size of C
        c.zero();                // do not add result to C
      } else c *= b;

      // Determine how many rows will be handled on each processor
      int nProcs = MPI_Wrappers::size(_comm);
      int myRank = MPI_Wrappers::rank(_comm);



      int *majorCounts = new int[nProcs];
      int *offsets = new int[nProcs];

#ifdef COL_STORAGE // Column-major storage
      int nMajor = nCols;
      int nMinor = nRows;
      int ParDenseMatrix::*majorField = &ParDenseMatrix::_nCols;
      int ParDenseMatrix::*minorField = &ParDenseMatrix::_nRows;
#else // Row-major storage
      int nMajor = nRows;
      int nMinor = nCols;
      int ParDenseMatrix::*majorField = &ParDenseMatrix::_nRows;
      int ParDenseMatrix::*minorField = &ParDenseMatrix::_nCols;
#endif

      for (int i = 0; i < nProcs; i++) {
        // If we have an uneven row-or-col/processors number, or too few rows
        // or cols, some processors will need to receive fewer rows/cols.
        offsets[i] = (i * nMajor) / nProcs;
        majorCounts[i] = (((i + 1) * nMajor) / nProcs) - offsets[i];
      }

      int myNMajor = majorCounts[myRank];
      int myMajorOffset = offsets[myRank];

      // Take data from an offset version of A
      ParDenseMatrix<double> A_local(_comm);
      A_local._data = this->_data + myMajorOffset * nMinor;
      A_local.*majorField = myNMajor;
      A_local.*minorField = nMinor;

#ifdef COL_STORAGE // Column-major storage

      // When splitting by columns, we split the vector as well, and sum the
      // results.

      DenseVector<double> v_local(myNMajor);
      for (int i = 0; i < myNMajor; i++)
        v_local(i) = v(myMajorOffset + i);

      // Store results in a local vector
      DenseVector<double> c_local = A_local * v_local;

      // Sum all vectors onto each processor
      allsum(_comm, c_local.ptr(), c.ptr(), c_local.size());

#else // Row-major storage

      // When splitting by rows, we use the whole vector and concatenate the
      // results.

      // Store results in a small local vector
      DenseVector<double> c_local(myNMajor);
      for (int i = 0; i < myNMajor; i++)
        c_local(i) = c(myMajorOffset + i);

      MultMv(A_local, v, c_local, At, a, b);

      // Gather the results onto each processor
      allgatherv(_comm, c_local.ptr(), c_local.size(), c.ptr(),
                 majorCounts, offsets);

#endif

      // Clear out the local matrix's pointer so we don't double-free
      A_local._data = NULL;

      delete [] majorCounts;
      delete [] offsets;

    }

  };

  // Operator for dense Matrix - dense vector product
  template<typename T>
  DenseVector<T> operator*(const ParDenseMatrix<T> &A, const Vector<T> &b)
  {
    DenseVector<T> c;
    A.ParMultMv(b, c, 0, 1.0, 0.0);
    return c;
  }


} // end namespace
#endif

