#ifndef SPARSEVECTOR_H
#define SPARSEVECTOR_H
#include "MatrixLibrary.h"

// No C++ templated typedefs, so use a define, gets cleaned up at end,
// so don't use outside of this class
#define STORE typename std::map<unsigned, T>

/** SparseVector class - an implimentation of a vector that contains a
 **    majority of zero element, and provides the relevant operations 
 **    specified by the base class Vector.         
 **/
template<class T>
class SparseVector : public Vector<T> {
  //* Multiplies a Matrix by a SparseVector (M*v) and returns a DenseVector.
  friend DenseVector<T> operator*<T>(const Matrix<T> &M, const SparseVector<T> &v);
  //* Multiplies a SparseVector by a Matrix (M'*v) and returns a DenseVector.
  friend DenseVector<T> operator*<T>(const SparseVector<T> &v, const Matrix<T> &M);
  //* Computes the dot product between two SparseVectors of equal length.
  friend T dot<T>(const SparseVector<T> &a, const SparseVector<T> &b);
  //* computes the product of a SparseMatrix tranpose with a SparseVector (M'*v).
  friend SparseVector<T> operator*<T>(const SparseMatrix<T> &M, const SparseVector<T> &v);
  //* computes the product of a SparseMatrix tranpose with a SparseVector (M'*v).
  friend SparseVector<T> operator*<T>(const SparseVector<T> &v, const SparseMatrix<T> &M);
public:
  //* Constructor - sets length of vector (NOT # of nonzeros).
  SparseVector(unsigned length=0);
  //* Copies another SparseVector 
  SparseVector(const SparseVector<T> &c); 
  //* Copies a general Vector (avoid if possible, its going to be slow).
  SparseVector(const Vector<T> &c); 

  //* Overrides output to string function to list only nonzeros and indices.
  std::string tostring() const;
  //* Indexing operators (w/ const overloading).
  //@{
  T  operator()(unsigned i, unsigned j=0) const;
  // NOTE reading a non-const SparseVector will call this and add zero entries.
  T& operator()(unsigned i, unsigned j=0);
  T  operator[](unsigned i) const;
  T& operator[](unsigned i);
  //@}
  //* assignment operators
  //@{
  SparseVector<T>& operator=(SparseVector<T> &c);
  SparseVector<T>& operator=(Vector<T> &c);
  //@}
  //* Return the number of rows in the Vector.
  unsigned nRows() const;
  //* Returns the number of columns - always 1.
  unsigned nCols() const { return 1; }
  //* Change # of Vector rows Vector and optionally keeps nonzeros (ignores nCols).
  void resize(unsigned nRows, unsigned nCols=1, bool copy=0);
  //* Return the number of nonzeros in the Vector.
  unsigned size() const;
  //* Changes size of Vector rows and optionally removes nonzeros.
  void reset (unsigned nRows, unsigned nCols=1, bool zero=0);
  //* zeros out all elements while preserving sparcity pattern
  void zero();
  //* TODO impliment copy (or maybe not necessary)
  void copy(const T* ptr, unsigned nRows, unsigned nCols=1);
  //* Writes a restart file (TODO impliment this if needed/wanted).
  void write_restart(FILE *F) const;

  // output to matlab (is this needed)
//  using Matrix<T>::matlab;
  //* Writes a matlab string to a stream that creates this object with a name.
  void matlab(ostream &o, const string &s="v") const;

protected:
  //* Banned operators
  //@{
  SparseVector(const Matrix<T> &c);
  SparseVector<T>& operator=(Matrix<T> &c);
  T* get_ptr() const {return NULL; } 
  //@}
  
  STORE data_;  //*> sparse data structure
  unsigned length_;             //*> number of rows
};

#include "SparseVector-inl.h"
#undef STORE 
#endif

