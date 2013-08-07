#ifndef PARDIAGONALMATRIX_H
#define PARDIAGONALMATRIX_H

#include "mpi.h"
#include "MPI_Wrappers.h"
#include "MatrixLibrary.h"
#include "DiagonalMatrix.h"
#include "DenseMatrix.h"

using namespace MPI_Wrappers;

namespace ATC_matrix {

  /**
   *  @class  ParDiagonalMatrix 
   *  @brief  Parallelized version of DiagonalMatrix class.
   */

  template <typename T>
  class ParDiagonalMatrix : public DiagonalMatrix<T> {
  public:
    ParDiagonalMatrix(MPI_Comm comm, INDEX rows=0, bool z=0)
      : DiagonalMatrix<T>(rows, z), _comm(comm) {}
    ParDiagonalMatrix(MPI_Comm comm, const DiagonalMatrix<T>& c)
      : DiagonalMatrix<T>(c), _comm(comm) {}
    ParDiagonalMatrix(MPI_Comm comm, const Vector<T>& v)
      : DiagonalMatrix<T>(v), _comm(comm) {}

    void MultAB(const Matrix<T> &B, DenseMatrix<T> &C) const;

    private:
        MPI_Comm _comm;
  };

  template<>
  void ParDiagonalMatrix<double>::MultAB(const Matrix<double> &B, DenseMatrix<double> &C) const;

} // end namespace
#endif

