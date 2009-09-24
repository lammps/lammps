#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
#include <exception>
#include "MatrixLibrary.h"
#include <algorithm>

/**
 * @struct TRI_COORD
 * @brief Triplet SparseMatrix entry
 */
template <typename T> 
struct TRI_COORD
{
  TRI_COORD<T>(unsigned row=0, unsigned col=0);
  TRI_COORD<T>(unsigned row, unsigned col, T val, bool add_to=0);
  unsigned i, j;
  T v;
  bool add;
};

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
  //* computes the product of a SparseMatrix tranpose with a SparseVector (M'*v).
  friend SparseVector<T> operator*<T>(const SparseMatrix<T> &M, const SparseVector<T> &v);
  //* computes the product of a SparseMatrix tranpose with a SparseVector (M'*v).
  friend SparseVector<T> operator*<T>(const SparseVector<T> &v, const SparseMatrix<T> &M);
public:
  SparseMatrix(INDEX rows=0, INDEX cols=0);
  SparseMatrix(const SparseMatrix<T>& c);
  SparseMatrix(const DenseMatrix<T>& c);
  ~SparseMatrix();
    
  //*  General index by value (requires a binary search on the row)
  T  operator()(INDEX i, INDEX j)  const;
  //*  General index by reference (requires a binary search on the row)
  T& operator()(INDEX i, INDEX j);
  //* General flat index by value operator (by nth nonzero)
  T  operator[](INDEX i) const;
  //* General flat index by reference operator (by nth nonzero)
  T& operator[](INDEX i);
  //* adds a value to index i,j
  void set(INDEX i, INDEX j, T v);
  //* sets a value to index i,j
  void add(INDEX i, INDEX j, T v);
  //* return a triplet value of the ith nonzero
  TRIPLET<T> get_triplet(INDEX i) const;

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
  //* returns the user-specified number of cols
  INDEX nCols()          const; 
  //* returns the number of non-zero elements
  INDEX size()           const;
  //* returns the number of non-zeros in a row
  INDEX RowSize(INDEX r) const;
  //* returns a pointer to the nonzero data
  inline T* get_ptr()    const;
  //* checks if the index i,j falls in the user-specified range
  bool in_range(INDEX i, INDEX j) const;

/*
 *  \section assignment operators
 */
  //* copies SparseMatrix R to this
  SparseMatrix<T>& operator=(const SparseMatrix &R);
  //* sets all nonzero values to a constant
  SparseMatrix<T>& operator=(const T v);
  //* scales all nonzero values by a constant
  SparseMatrix<T>& operator*=(const T &a);

/*
 *  \section Multiplication operations
 */

  //* S' * F
  DenseMatrix<T> transMat(const DenseMatrix<T> &D) const;
  //* S' * S
  DenseMatrix<T> transMat(const SparseMatrix<T> &S) const;
  //* S' * v
  DenseVector<T> transMat(const Vector<T> &x) const;

  SparseMatrix<T> transpose() const;
  SparseMatrix<T>& row_scale(const Vector<T> &v);
  SparseMatrix<T>& col_scale(const Vector<T> &v);
  DenseVector<T> col_sum()                               const;
  DenseVector<INDEX> column_count()                      const;
  DiagonalMatrix<T> get_diag()                           const;
  void get_row(INDEX i, DenseVector<T>& row, DenseVector<INDEX>& indx) const;
  void WeightedLeastSquares(const SparseMatrix<T> &N, const DiagonalMatrix<T> &D);

  T get_row_max(INDEX row) const;
  T get_row_min(INDEX row) const;

/*
 *  \section I/O functions
 */
  //* outputs this SparseMatrix to a formatted string
  string tostring() const;
  using Matrix<T>::matlab;
  //* writes a command to recreate this matrix in matlab to a stream
  void matlab(ostream &o, const string &name="S")  const;
  //* prints a row histogram for each row
  void print_row_histogram(const string &name, INDEX nbins = 10) const;
  //* prints a histogram of the values in a row
  void print_row_histogram(INDEX row, INDEX nbins) const;
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
  //* copys all data from another SparseMatrix
  void _copy(const SparseMatrix<T> &C);
  //* general sparse matrix assignment
  void _set_equal(const Matrix<T> &r);
  //* returns the first row with a nonzero in it (from the CRS structure only)
  int _first_nonzero_row_crs() const;

/*
 *  \section CRS storage variables
 */
  T * _val;                    // matrix non-zeros
  INDEX *_ia, *_ja;            // ptrs to rows, column indexes
  INDEX _size, _nRowsCRS;      // # of non-zeros, rows

  //* new (unsorted triplet values - won't intersect CRS values)
  mutable vector<TRI_COORD<T> > _tri;
/*
 *  \section User specified variables
 */  
  INDEX _nRows, _nCols;
  static T _zero;
};

#include "SparseMatrix-inl.h"
#endif

