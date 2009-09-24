#ifndef MATRIXLIBRARY_H
#define MATRIXLIBRARY_H

#include "DenseMatrix.h"
#include "DenseVector.h"
#include "CloneVector.h"
#include "DiagonalMatrix.h"
#include "SparseMatrix.h"
#include "SparseVector.h"

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
    triplet = s->get_triplet(i);
    m(triplet.i, triplet.j) = triplet.v; 
  }
}

#endif
