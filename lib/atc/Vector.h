#ifndef VECTOR_H
#define VECTOR_H

#include "Matrix.h"

namespace ATC_matrix {

///////////////////////////////////////////////////////////////////////////////
// forward declarations ///////////////////////////////////////////////////////

//* Matrix-vector product
//template<typename T>
//void MultMv(const Matrix<T> &A, const Vector<T> &v, DenseVector<T> &c, 
//            const bool At=0, T a=1, T b=0);

/******************************************************************************
* abstract class Vector
******************************************************************************/


template<typename T>
class Vector : public Matrix<T>
{
public:
  Vector() {}
  Vector(const Vector<T> &c); // do not implement!
  virtual ~Vector() {}

  std::string to_string() const;

  // pure virtual functions
  virtual T  operator()(INDEX i, INDEX j=0) const=0;
  virtual T& operator()(INDEX i, INDEX j=0)      =0;
  virtual T  operator[](INDEX i)            const=0;
  virtual T& operator[](INDEX i)                 =0;
  virtual INDEX nRows()                     const=0;
  virtual T* ptr()                      const=0;
  virtual void resize(INDEX nRows, INDEX nCols=1, bool copy=0)=0;
  virtual void  reset(INDEX nRows, INDEX nCols=1, bool zero=0)=0;
  virtual void copy(const T * ptr, INDEX nRows, INDEX nCols=1)=0;
  void write_restart(FILE *f)               const; // will be virtual

  
  // output to matlab
  using Matrix<T>::matlab;
  void matlab(std::ostream &o, const std::string &s="v") const;

  using Matrix<T>::operator=;
  INDEX nCols()                   const;
  bool in_range(INDEX i)          const;
  bool same_size(const Vector &m) const;
  static bool same_size(const Vector &a, const Vector &b);

 protected:
  void _set_equal(const Matrix<T> &r);
  //* don't allow this
  Vector& operator=(const Vector &r);
};

///////////////////////////////////////////////////////////////////////////////
//* performs a matrix-vector multiply with default naive implementation
template<typename T>
void MultMv(const Matrix<T> &A, const Vector<T> &v, DenseVector<T> &c, 
            const bool At, T /* a */, T b)
{
  const INDEX sA[2] = {A.nRows(), A.nCols()};  // m is sA[At] k is sA[!At]
  const INDEX M=sA[At], K=sA[!At];
  GCK(A, v, v.size()!=K, "MultAb<T>: matrix-vector multiply");
  if (c.size() != M)
  {
    c.resize(M);             // set size of C
    c.zero();                // do not add result to C
  }
  else c *= b;
  for (INDEX p=0; p<M; p++)
    for (INDEX r=0; r<K; r++)
       c[p] += A(p*!At+r*At, p*At+r*!At) * v[r];
}
///////////////////////////////////////////////////////////////////////////////
//* Operator for Matrix-vector product
template<typename T>
DenseVector<T> operator*(const Matrix<T> &A, const Vector<T> &b)
{
  DenseVector<T> c;
  MultMv(A, b, c, 0, 1.0, 0.0);
  return c;
}
///////////////////////////////////////////////////////////////////////////////
//* Operator for Vector-matrix product
template<typename T>
DenseVector<T> operator*(const Vector<T> &a, const Matrix<T> &B)
{
  DenseVector<T> c;
  MultMv(B, a, c, 1, 1.0, 0.0);
  return c;
}
///////////////////////////////////////////////////////////////////////////////
//* Multiply a vector by a scalar 
template<typename T>
DenseVector<T> operator*(const Vector<T> &v, const T s)
{
  DenseVector<T> r(v);
  r*=s; 
  return r;
}
///////////////////////////////////////////////////////////////////////////////
//* Multiply a vector by a scalar - communitive
template<typename T>
DenseVector<T> operator*(const T s, const Vector<T> &v)
{
  DenseVector<T> r(v);
  r*=s; 
  return r;
}
///////////////////////////////////////////////////////////////////////////////
//* inverse scaling operator - must always create memory
template<typename T>
DenseVector<T> operator/(const Vector<T> &v, const T s)
{
  DenseVector<T> r(v);
  r*=(1.0/s); // for integer types this may be worthless
  return r;
}
///////////////////////////////////////////////////////////////////////////////
//* Operator for Vector-Vector sum
template<typename T>
DenseVector<T> operator+(const Vector<T> &a, const Vector<T> &b)
{
  DenseVector<T> c(a);
  c+=b;
  return c;
}
///////////////////////////////////////////////////////////////////////////////
//* Operator for Vector-Vector subtraction
template<typename T>
DenseVector<T> operator-(const Vector<T> &a, const Vector<T> &b)
{
  DenseVector<T> c(a);
  c-=b;
  return c;
}

///////////////////////////////////////////////////////////////////////////////
// Template definitions ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//* output operator
template<typename T>
std::string Vector<T>::to_string() const
{
  std::string s;
  int sz = this->size(); 
  for (INDEX i = 0; i < sz; i++) 
    s += std::string(i?"\t":"") + ATC_Utility::to_string((*this)[i],myPrecision);
  return s;
}
///////////////////////////////////////////////////////////////////////////////
//* Writes a matlab script defining the vector to the stream
template<typename T>
void Vector<T>::matlab(std::ostream &o, const std::string &s) const
{
  o << s <<"=zeros(" << this->size() << ",1);\n";
  int sz = this->size(); 
  for (INDEX i = 0; i < sz; i++) 
    o << s << "("<<i+1<<") = " << (*this)[i] << ";\n";
}
///////////////////////////////////////////////////////////////////////////////
//* writes the vector data to a file
template <typename T>
void Vector<T>::write_restart(FILE *f)                                    const
{
  INDEX size = this->size();
  fwrite(&size, sizeof(INDEX),1,f);
  if (size) fwrite(this->ptr(), sizeof(T), this->size(), f);
}
///////////////////////////////////////////////////////////////////////////////
//* returns the number of columns; always 1
template<typename T>
inline INDEX Vector<T>::nCols()                                           const
{
  return 1; 
}
///////////////////////////////////////////////////////////////////////////////
//* returns true if INDEX i is within the range of the vector
template<typename T>
bool Vector<T>::in_range(INDEX i)                                         const  
{
  return i<this->size(); 
}
///////////////////////////////////////////////////////////////////////////////
//* returns true if m has the same number of elements this vector
template<typename T>
bool Vector<T>::same_size(const Vector &m)                                const 
{
  return this->size() == m.size(); 
}
///////////////////////////////////////////////////////////////////////////////
//* returns true if a and b have the same number of elements
template<typename T>
inline bool Vector<T>::same_size(const Vector &a, const Vector &b)
{
  return a.same_size(b); 
}
//----------------------------------------------------------------------------
// general matrix assignment (for densely packed matrices)
//----------------------------------------------------------------------------
template<typename T>
void Vector<T>::_set_equal(const Matrix<T> &r)
{
  this->resize(r.nRows(), r.nCols());
  const Matrix<T> *pr = &r;
#ifdef OBSOLETE
  if (const SparseMatrix<T> *ps = dynamic_cast<const SparseMatrix<T>*>(pr))//sparse_cast(pr)) 
    copy_sparse_to_matrix(ps, *this);
  
  else if (dynamic_cast<const DiagonalMatrix<T>*>(pr))//diag_cast(pr))  // r is Diagonal?
  {
    this->zero();
    for (INDEX i=0; i<r.size(); i++) (*this)(i,i) = r[i];
  }
  else memcpy(this->ptr(), r.ptr(), r.size()*sizeof(T));
#else
  const Vector<T>         *pv = dynamic_cast<const Vector<T>*>         (pr);
  if (pv)  this->copy(pv->ptr(),pv->nRows());
  else
  {
    std::cout <<"Error in general vector assignment\n";
    exit(1);
  }
#endif
}
} // end namespace


#endif
