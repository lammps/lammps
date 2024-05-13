#ifndef MATRIXLIBRARY_H
#define MATRIXLIBRARY_H

#include "DenseMatrix.h"
#include "ParDenseMatrix.h"
#include "DenseVector.h"
#include "CloneVector.h"
#include "DiagonalMatrix.h"
#include "ParDiagonalMatrix.h"
#include "SparseMatrix.h"
#include "ParSparseMatrix.h"
#include "SparseVector.h"


namespace ATC_matrix {
template<typename T>
const SparseMatrix<T> *sparse_cast(const Matrix<T> *m)
{
  return dynamic_cast<const SparseMatrix<T>*>(m);
}

template<typename T>
void copy_sparse_to_matrix(const SparseMatrix<T> *s, Matrix<T> &m)
{
  m.zero();
  TRIPLET<T> triplet;
  for (INDEX i=0; i<s->size(); i++)
  {
    triplet = s->triplet(i);
    m(triplet.i, triplet.j) = triplet.v;
  }
}

} // end namespace

using namespace ATC_matrix;

#endif
