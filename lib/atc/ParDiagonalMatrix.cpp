#include "ParDiagonalMatrix.h"

using MPI_Wrappers::allgatherv;

namespace ATC_matrix {

  // template<>
  // void ParDiagonalMatrix<double>::MultAB(const Matrix<double> &B, DenseMatrix<double> &C) const
  // {
  //   //SparseMatrix<T>::compress(*this);
  //   GCK(*this, B, this->nCols()!=B.nRows(), "ParDiagonalMatrix * Matrix");

  //   const INDEX nRows = this->nRows();
  //   const INDEX nCols = this->nCols();

  //   // Determine which rows will be handled on this processor
  //   int nProcs = MPI_Wrappers::size(_comm);
  //   int myRank = MPI_Wrappers::rank(_comm);

  //   INDEX startIndex = (myRank * nRows) / nProcs;
  //   INDEX endIndex = ((myRank + 1) * nRows) / nProcs;

  //   // Calculate the scaled rows associated with this processor
  //   for (INDEX i = startIndex; i < endIndex; i++) {
  //     double value = (*this)[i];
  //     for (INDEX j = 0; j < nCols; j++)
  //       C(i, j) = value * B(i, j);
  //   }

  //   // Collect results on all processors

  //   //       consider sending only owned rows from each processor
  //   allsum(_comm, MPI_IN_PLACE, C.ptr(), C.size());
  // }

  template<>
  void ParDiagonalMatrix<double>::MultAB(const Matrix<double> &B, DenseMatrix<double> &C) const
  {
    //SparseMatrix<T>::compress(*this);
    GCK(*this, B, this->nCols()!=B.nRows(), "ParDiagonalMatrix * Matrix");

    const INDEX nRows = this->nRows();
    const INDEX nCols = this->nCols();

    int nProcs = MPI_Wrappers::size(_comm);
    int myRank = MPI_Wrappers::rank(_comm);

#ifdef COL_STORAGE // Column-major storage
    int nMajor = nCols;
    int nMinor = nRows;
#else // Row-major storage
    int nMajor = nRows;
    int nMinor = nCols;
#endif

    int *majorCounts = new int[nProcs];
    int *majorOffsets = new int[nProcs];

    // Determine which rows/columns will be handled on this processor
    for (int i = 0; i < nProcs; i++) {
      majorOffsets[i] = (i * nMajor) / nProcs;
      majorCounts[i] = (((i + 1) * nMajor) / nProcs) - majorOffsets[i];
    }

    INDEX myNMajor = majorCounts[myRank];
    INDEX myMajorOffset = majorOffsets[myRank];

    // Calculate the scaled values associated with this processor, in row chunks

#ifdef COL_STORAGE // Column-major storage
    for (INDEX i = 0; i < nRows; i++) {
      double value = (*this)[i];
      for (INDEX j = myMajorOffset; j < myMajorOffset + myNMajor; j++)
        C(i, j) = value * B(i, j);
    }
#else // Row-major storage
    for (INDEX i = myMajorOffset; i < myMajorOffset + myNMajor; i++) {
      double value = (*this)[i];
      for (INDEX j = 0; j < nCols; j++)
        C(i, j) = value * B(i, j);
    }
#endif

    for (int i = 0; i < nProcs; i++) {
      majorCounts[i] *= nMinor;
      majorOffsets[i] *= nMinor;
    }
    // Collect results on all processors
    allgatherv(_comm, C.ptr() + myMajorOffset * nMinor, myNMajor * nMinor,
               C.ptr(), majorCounts, majorOffsets);

    delete[] majorCounts;
    delete[] majorOffsets;
  }


} // end namespace
