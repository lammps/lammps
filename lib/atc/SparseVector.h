#ifndef SPARSEVECTOR_H
#define SPARSEVECTOR_H
#include "MatrixLibrary.h"

namespace ATC_matrix {

// No C++ templated typedefs, so use a define, gets cleaned up at end,
// so don't use outside of this class
#define STORE typename std::map<INDEX, T>

template<class T> class SparseVector;
template<class T> T dot(const SparseVector<T> &a, const SparseVector<T> &b);

  /**
   *  @class  SparseVector 
   *  @brief  Class for vectors that contain a majority of zero elements and provides relevant operations 
   */

template<class T>
class SparseVector : public Vector<T> {
  //* Multiplies a Matrix by a SparseVector (M*v) and returns a DenseVector.
  friend DenseVector<T> operator*<T>(const Matrix<T> &M, const SparseVector<T> &v);
  //* Multiplies a SparseVector by a Matrix (M'*v) and returns a DenseVector.
  friend DenseVector<T> operator*<T>(const SparseVector<T> &v, const Matrix<T> &M);
  //* Computes the dot product between two SparseVectors of equal length.
#ifdef __INTEL_COMPILER
  // for use on Intel compilers
  template<class T> friend T dot(const SparseVector<T> &a, const SparseVector<T> &b);
#else
  // for use with gcc
  friend T dot<T>(const SparseVector<T> &a, const SparseVector<T> &b);
#endif
  //* computes the product of a SparseMatrix transpose with a SparseVector (M'*v).
  friend SparseVector<T> operator*<T>(const SparseMatrix<T> &M, const SparseVector<T> &v);
  //* computes the product of a SparseMatrix transpose with a SparseVector (M'*v).
  friend SparseVector<T> operator*<T>(const SparseVector<T> &v, const SparseMatrix<T> &M);
public:
  //* Constructor - sets length of vector (NOT # of nonzeros).
  SparseVector(INDEX length=0);
  //* Copies another SparseVector 
  SparseVector(const SparseVector<T> &c); 
  //* Copies a general Vector (avoid if possible, its going to be slow).
  SparseVector(const Vector<T> &c); 

  //* Overrides output to string function to list only nonzeros and indices.
  std::string to_string() const;
  //* Indexing operators (w/ const overloading).
  //@{
  T  operator()(INDEX i, INDEX j=0) const;
  
  T& operator()(INDEX i, INDEX j=0);
  T  operator[](INDEX i) const;
  T& operator[](INDEX i);
  //* Returns a pair (index, value) for a nonzero in the vector.
  std::pair<INDEX, T> pair(INDEX i) const;
  //@}
  //* assignment operators
  //@{
  SparseVector<T>& operator=(const SparseVector<T> &c);
  SparseVector<T>& operator=(Vector<T> &c);
  //@}
  //* Return the number of rows in the Vector.
  INDEX nRows() const;
  //* Returns the number of columns - always 1.
  INDEX nCols() const { return 1; }
  //* Change # of Vector rows Vector and optionally keeps nonzeros (ignores nCols).
  void resize(INDEX nRows, INDEX nCols=1, bool copy=0);
  //* Return the number of nonzeros in the Vector.
  INDEX size() const;
  //* Changes size of Vector rows and optionally removes nonzeros.
  void reset (INDEX nRows, INDEX nCols=1, bool zero=0);
  //* zeros out all elements while preserving sparcity pattern
  void zero();
  //* TODO implement copy (or maybe not necessary)
  void copy(const T* ptr, INDEX nRows, INDEX nCols=1);
  //* Writes a restart file (TODO implement this if needed/wanted).
  void write_restart(FILE *F) const;
  //* Adds SparseVector x, scaled by s to this one.  Can be different sparcity.
  void add_scaled(SparseVector<T>& x, const T& s);

  // output to matlab (is this needed)
//  using Matrix<T>::matlab;
  //* Writes a matlab string to a stream that creates this object with a name.
  void matlab(std::ostream &o, const std::string &s="v") const;

protected:
  //* Banned operators
  //@{
  SparseVector(const Matrix<T> &c);
  SparseVector<T>& operator=(Matrix<T> &c);
  T* ptr() const {return NULL; } 
  //@}
  
  STORE data_;  //*> sparse data structure
  INDEX length_;             //*> number of rows
};

} // end namespace

#include "SparseVector-inl.h"
#undef STORE 
#endif

