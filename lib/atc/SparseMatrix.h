#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
#include <exception>
#include "MatrixLibrary.h"
#include <algorithm>

namespace ATC_matrix {

/**
 * @struct TRI_COORD
 * @brief Triplet SparseMatrix entry
 */
template <typename T> 
struct TRI_COORD
{
  TRI_COORD<T>(INDEX row=0, INDEX col=0);
  TRI_COORD<T>(INDEX row, INDEX col, T val, bool add_to=0);
  INDEX i, j;
  T v;
  bool add;
};

template<typename T>
void ParMultAB(MPI_Comm comm, const SparseMatrix<T>& A, const Matrix<T>& B, DenseMatrix<T>& C);

/**
 * @class SparseMatrix
 * @brief Stores data in triplet format or CRS format 
 */
template<typename T>
class SparseMatrix : public Matrix<T>
{
  //* SparseMatrix-Vector multiplication  (S * v)
  friend DenseVector<T>  operator*<T>(const SparseMatrix<T> &A, const Vector<T>& x);
  //* SparseMatrix-DenseMatrix multiplication (S * F)
  friend DenseMatrix<T>  operator*<T>(const SparseMatrix<T> &A, const Matrix<T>& D);
  //* SparseMatrix-DiagonalMatrix multiplication (S * D)
  friend SparseMatrix<T> operator*<T>(const SparseMatrix<T> &A, const DiagonalMatrix<T>& D);
  //* SparseMatrix-SparseMatrix multiplication (S * S)
  friend SparseMatrix<T> operator*<T>(const SparseMatrix<T> &A, const SparseMatrix<T> &B);
  //* computes the product of a SparseMatrix transpose with a SparseVector (M'*v).
  friend SparseVector<T> operator*<T>(const SparseMatrix<T> &M, const SparseVector<T> &v);
  //* computes the product of a SparseMatrix transpose with a SparseVector (M'*v).
  friend SparseVector<T> operator*<T>(const SparseVector<T> &v, const SparseMatrix<T> &M);

  template<typename U>
  friend void ParMultAB(MPI_Comm comm, const SparseMatrix<U>& A, const Matrix<U>& B, DenseMatrix<U>& C);

public:
  SparseMatrix(INDEX rows=0, INDEX cols=0);
  SparseMatrix(const SparseMatrix<T>& c);
  SparseMatrix(const DenseMatrix<T>& c);
  SparseMatrix(INDEX* rows, INDEX* cols, T* vals, INDEX size, 
               INDEX nRows, INDEX nCols, INDEX nRowsCRS);
  virtual ~SparseMatrix() { _delete(); }
    
  //*  General index by value (requires a binary search on the row)
  T  operator()(INDEX i, INDEX j)  const;
  //*  General index by reference (requires a binary search on the row)
  T& operator()(INDEX i, INDEX j);
  //* General flat index by value operator (by nth nonzero)
  T  operator[](INDEX i) const;
  //* General flat index by reference operator (by nth nonzero)
  T& operator[](INDEX i);
  //* sets a value to index i,j
  void set(INDEX i, INDEX j, T v);
  //* adds a value to index i,j
  void add(INDEX i, INDEX j, T v);
  //* return a triplet value of the ith nonzero
  TRIPLET<T> triplet(INDEX i) const;

  //* full reset - completely wipes out all SparseMatrix data
  void reset(INDEX rows=0, INDEX cols=0, bool zero=true);
  //* only changes the bounds of the matrix, no deletion
  void resize(INDEX rows=0, INDEX cols=0, bool zero=true);
  //* reset - from DenseMatrix - this will be SLOW
  void reset(const DenseMatrix<T>& D, double TOL=-1.0); 
  //* copy data
  void copy(const T * ptr, INDEX rows=0, INDEX cols=0);
  
  void dense_copy(DenseMatrix<T>& D) const;
  DenseMatrix<T>  dense_copy(void) const;

  //* returns true if the matrix has no nonzero elements
  bool  empty()          const;
  //* returns the user-specified number of rows
  INDEX nRows()          const;
  INDEX nRowsCRS()          const;
  //* returns the user-specified number of cols
  INDEX nCols()          const; 
  //* returns the number of non-zero elements
  INDEX size()           const;
  //* returns the number of non-zeros in a row
  INDEX RowSize(INDEX r) const;
  //* returns a pointer to the CRS list of rows
  inline INDEX* rows() const;
  //* returns a pointer to the CRS list of cols
  inline INDEX* cols() const;
  //* returns a pointer to the nonzero data
  inline T* ptr()    const;
  //* checks if the index i,j falls in the user-specified range
  bool in_range(INDEX i, INDEX j) const;
  //* check if the total matrix has a value set for an index pair
  bool has_entry(INDEX i, INDEX j) const;
  //* check if the uncompressed part of the matrix has a value set for an index pair
  bool has_entry_uncompressed(INDEX i, INDEX j) const;
  //* check if the compressed part matrix has a value set for an index pair
  bool has_entry_compressed(INDEX i, INDEX j) const;
  //* check if the matrix has been compressed at least once
  bool has_template(void) const;

/*
 *  \section assignment operators
 */
  //* copies SparseMatrix R to this
  SparseMatrix<T>& operator=(const SparseMatrix &R);
  //* sets all nonzero values to a constant
  SparseMatrix<T>& operator=(const T v);
  //* scales all nonzero values by a constant
  SparseMatrix<T>& operator*=(const T &a);
  //* calls operator*= of base class
  SparseMatrix<T>& operator*=(const SparseMatrix<T> &a);
  // Adds two matrices together.
  SparseMatrix<T>& operator+=(const SparseMatrix & R);

/*
 *  \section Multiplication operations
 */

  //-----------------------------------------------------------------------------
  // multiply sparse matrix by a vector
  //-----------------------------------------------------------------------------
  virtual void MultMv(const Vector<T>& v, DenseVector<T>& c) const
  {
    compress(*this);
    GCK(*this, v, this->nCols() != v.size(), "SparseMatrix * Vector")

    // resize c if necessary
    if (c.size() != this->nRows()) {
      c.resize(this->nRows());
      c.zero();
    }

    INDEX i, j;
    for (i = 0; i < this->_nRowsCRS; i++)
      for (j = this->_ia[i]; j < this->_ia[i + 1]; j++)
        c(i) += this->_val[j] * v(this->_ja[j]);
  }

  //-----------------------------------------------------------------------------
  // multiply sparse matrix by dense matrix
  //-----------------------------------------------------------------------------
  virtual void MultAB(const Matrix<T>& B, DenseMatrix<T>& C) const
  {
    GCK(*this, B, this->nCols() != B.nRows(), "SparseMatrix * DenseMatrix")

    const INDEX J = B.nCols();

    INDEX i, ik, j;
    for (i = 0; i < this->_nRowsCRS; i++)
      for (ik = this->_ia[i]; ik < this->_ia[i + 1]; ik++)
        for (j = 0; j < J; j++)
          C(i, j) += this->_val[ik] * B(this->_ja[ik], j);  // C(i,j) = S(i,k) * B(k, j)
  }

  //-----------------------------------------------------------------------------
  // Multiplies this SparseMatrix transposed times a vector
  //-----------------------------------------------------------------------------
  virtual DenseVector<T> transMat(const Vector<T> &x) const
  {
    compress(*this);
    DenseVector<T> y(nCols(), true);
    GCK(*this, x, nRows()!=x.size(),"operator *: Sparse matrix incompatible with Vector.")

    INDEX i, ij;
    for(i=0; i<_nRowsCRS; i++)
      for(ij=_ia[i]; ij<_ia[i+1]; ij++)
        y(_ja[ij]) += _val[ij]*x(i);
    return y;
  }

  //-----------------------------------------------------------------------------
  // Matrix Transpose/DenseMatrix multiply
  //-----------------------------------------------------------------------------
  virtual DenseMatrix<T> transMat(const DenseMatrix<T> &D) const
  {
    compress(*this);
    GCK(*this, D, nRows()!=D.nRows(),"transMat: Sparse matrix incompatible with DenseMatrix.")
    DenseMatrix<T> C(nCols(), D.nCols(), true);  // initialized to zero
    INDEX j, k, ki;

    for (k=0; k<_nRowsCRS; k++)
      for (ki=_ia[k]; ki<_ia[k+1]; ki++)
        for (j=0; j<D.nCols(); j++)
          C(_ja[ki], j) += _val[ki]*D(k,j);     // C(i,j) = S(k,i) * D(k, j)

    return C;
  }

  //-----------------------------------------------------------------------------
  // Matrix Transpose/SparseMatrix multiply - IS THIS REALLY NEEDED??
  //-----------------------------------------------------------------------------
  virtual DenseMatrix<T> transMat(const SparseMatrix<T> &D) const
  {
    compress(*this);
    GCK(*this, D, nRows()!=D.nRows(),"transMat: Sparse matrix incompatible with DenseMatrix.")
    DenseMatrix<T> C(nCols(), D.nCols(), true); // initialized to zero

    INDEX k, ki, kj;
    for (k=0; k<_nRowsCRS; k++)
      for (kj=D._ia[k]; kj<D._ia[k+1]; kj++)
        for (ki=_ia[k]; ki<_ia[k+1]; ki++)
          C(_ja[ki], D._ja[kj]) += _val[ki]*D._val[kj]; // C(i,j) = S(k,i)*D(k,j)

    return C;
  }

  SparseMatrix<T> transpose() const;
  SparseMatrix<T>& row_scale(const Vector<T> &v);
  SparseMatrix<T>& col_scale(const Vector<T> &v);
  DenseVector<T> col_sum()                               const;
  DenseVector<INDEX> column_count()                      const;
  DiagonalMatrix<T> diag()                           const;
  DiagonalMatrix<T> row_sum_lump()                       const;
  void row(INDEX i, DenseVector<T>& row, DenseVector<INDEX>& indx) const;
  void weighted_least_squares(const SparseMatrix<T> &N, const DiagonalMatrix<T> &D);
  void set_all_elements_to(const T &v);

  T row_max(INDEX row) const;
  T row_min(INDEX row) const;

/*
 *  \section I/O functions
 */
  //* outputs this SparseMatrix to a formatted string
  std::string to_string() const;
  using Matrix<T>::matlab;
  //* writes a command to recreate this matrix in matlab to a stream
  void matlab(std::ostream &o, const std::string &name="S")  const;
  //* prints a row histogram for each row
  void print_row_histogram(const std::string &name, INDEX nbins = 10) const;
  //* prints a histogram of the values in a row
  void print_row_histogram(INDEX row, INDEX nbins) const;
  //* prints the current triplets
  void print_triplets() const;
  //! Writes the matrix to a binary file (after a compress).
  void binary_write(std::fstream& f) const;
  //! Reads a SparseMatrix from a binary file.  (wipes out any original data)
  void binary_read(std::fstream& f);
  //* Dump templated type to disk; operation not safe for all types
  void write_restart(FILE *f) const; 

/*
 *  \section Utility functions
 */
  //* converts all triplets and merges with CRS
  void compress();          
  //* converts T to CRS
  static void compress(const SparseMatrix<T> &C); 
  //* sorts and returns the # of unique triplets
  INDEX CountUniqueTriplets();       

private:
  //* creates a CRS structure
  void _create(INDEX size, INDEX nrows);
  //* clears all memory and nulls references
  void _delete();   
  //* copies all data from another SparseMatrix
  void _copy(const SparseMatrix<T> &C);
  //* general sparse matrix assignment
  void _set_equal(const Matrix<T> &r);
  //* returns the first row with a nonzero in it (from the CRS structure only)
  int _first_nonzero_row_crs() const;

/*
 *  \section CRS storage variables
 */
protected:
  T * _val;                    // matrix non-zeros
  INDEX *_ia, *_ja;            // ptrs to rows, column indexes
  INDEX _size, _nRowsCRS;      // # of non-zeros, rows
  bool hasTemplate_;  

  void copy(const SparseMatrix<T> &C);

  //* new (unsorted triplet values - won't intersect CRS values)
  mutable std::vector<TRI_COORD<T> > _tri;
/*
 *  \section User specified variables
 */  
  INDEX _nRows, _nCols;
  static T _zero;
};

} // end namespace

#include "SparseMatrix-inl.h"
#endif

