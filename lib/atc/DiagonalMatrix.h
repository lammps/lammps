#ifndef DIAGONALMATRIX_H
#define DIAGONALMATRIX_H

#include "MatrixDef.h"

template<typename T>
class DiagonalMatrix : public Matrix<T>
{
 public: 
  explicit DiagonalMatrix(INDEX nRows=0, bool zero=0);
  DiagonalMatrix(const DiagonalMatrix<T>& c);
  DiagonalMatrix(const Vector<T>& v);
  virtual ~DiagonalMatrix();
 
  //* resizes the matrix, ignores nCols, optionally zeros 
  void  reset(INDEX rows, INDEX cols=0, bool zero=true);
  //* resizes the matrix, ignores nCols, optionally copies what fits
  void resize(INDEX rows, INDEX cols=0, bool copy=false);
  //* resets based on full copy of vector v
  void reset(const Vector<T>& v);
  //* resets based on full copy of a DiagonalMatrix
  void reset(const DiagonalMatrix<T>& v);
  //* resizes the matrix, ignores nCols, optionally copies what fits
  void copy(const T * ptr, INDEX rows, INDEX cols=0);

  //* resets based on a "shallow" copy of a vector
  void shallowreset(const Vector<T> &v);
  //* resets based on a "shallow" copy of a DiagonalMatrix
  void shallowreset(const DiagonalMatrix<T> &c);

  T& operator()(INDEX i, INDEX j);
  T  operator()(INDEX i, INDEX j) const;
  T& operator[](INDEX i);
  T  operator[](INDEX i)          const;
  INDEX nRows()                   const;
  INDEX nCols()                   const;
  T* get_ptr()                    const;
  void write_restart(FILE *f)     const;

  // Dump matrix contents to screen (not defined for all datatypes)
  string tostring() const { return _data->tostring(); }

  using Matrix<T>::matlab;
  void matlab(ostream &o, const string &s="D") const;

  // overloaded operators
  DiagonalMatrix<T>& operator=(const T s);
  DiagonalMatrix<T>& operator=(const DiagonalMatrix<T> &C);
  //DiagonalMatrix<T>& operator=(const Vector<T> &C);

  INDEX size() const { return _data->size(); }
  //* computes the inverse of this matrix
  DiagonalMatrix<T>& inv_this(); 
  //* returns a copy of the inverse of this matrix
  DiagonalMatrix<T>  inv() const;

protected:
  void _set_equal(const Matrix<T> &r);
  DiagonalMatrix& operator=(const Vector<T> &c) {}
  DiagonalMatrix& operator=(const Matrix<T> &c) {}

private: 
  void _delete();
  Vector<T> *_data;
}; 

//-----------------------------------------------------------------------------
// DiagonalMatrix-DiagonalMatrix multiplication
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T> operator*(const DiagonalMatrix<T>& A, const DiagonalMatrix<T>& B) 
{
  SSCK(A, B, "DiagonalMatrix-DiagonalMatrix multiplication");
  DiagonalMatrix<T> R(A); 
  for (INDEX i=0; i<R.nRows(); i++) R[i] *= B[i];
  return R;
}
//-----------------------------------------------------------------------------
// DiagonalMatrix-matrix multiplication
//-----------------------------------------------------------------------------
template<typename T>
DenseMatrix<T> operator*(const DiagonalMatrix<T>& A, const Matrix<T> &B)
{
  GCK(A, B, A.nCols()!=B.nRows(), "DiagonalMatrix-Matrix multiplication");
  DenseMatrix<T> R(B);  // makes a copy of r to return
  for (INDEX i=0; i<R.nRows(); i++)
    for (INDEX j=0; j<R.nCols(); j++)
      R(i,j) *= A[i];
  return R;
}
//-----------------------------------------------------------------------------
// matrix-DiagonalMatrix multiplication
//-----------------------------------------------------------------------------
template<typename T>
DenseMatrix<T> operator*(const Matrix<T> &B, const DiagonalMatrix<T>& A)
{
  GCK(B, A, B.nCols()!=A.nRows(), "Matrix-DiagonalMatrix multiplication");
  DenseMatrix<T> R(B);  // makes a copy of r to return
  for (INDEX j=0; j<R.nCols(); j++)
    for (INDEX i=0; i<R.nRows(); i++)
      R(i,j) *= A[j];
  return R;
}
//-----------------------------------------------------------------------------
// DiagonalMatrix-vector multiplication
//-----------------------------------------------------------------------------
template<typename T>
DenseVector<T> operator*(const DiagonalMatrix<T>& A, const Vector<T> &b)
{
  GCK(A, b, A.nCols()!=b.size(), "DiagonalMatrix-Vector multiplication");
  DenseVector<T> r(b);  // makes a copy of r to return
  for (INDEX i=0; i<r.size(); i++)
    r[i] *= A[i];
  return r;
}
//-----------------------------------------------------------------------------
// vector-DiagonalMatrix multiplication
//-----------------------------------------------------------------------------
template<typename T>
DenseVector<T> operator*(const Vector<T> &b, const DiagonalMatrix<T>& A)
{
  GCK(b, A, b.size()!=A.nRows(), "Matrix-DiagonalMatrix multiplication");
  DenseVector<T> r(b);  // makes a copy of r to return
  for (INDEX i=0; i<r.size(); i++)
    r[i] *= A[i];
  return r;
}
//-----------------------------------------------------------------------------
// DiagonalMatrix-SparseMatrix multiplication
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T> operator*(const DiagonalMatrix<T> &A, const SparseMatrix<T>& B) 
{
  GCK(A, B, A.nCols()!=B.nRows() ,"DiagonalMatrix-SparseMatrix multiplication");
  SparseMatrix<T> R(B);
  CloneVector<T> d(A);
  R.row_scale(d);
  return R;
}
//-----------------------------------------------------------------------------
// DiagonalMatrix-scalar multiplication
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T> operator*(DiagonalMatrix<T> &A, const T s)
{
  DiagonalMatrix<T> R(A);
  R *= s;
  return R;
}
//-----------------------------------------------------------------------------
// Commute with DiagonalMatrix * double
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T> operator*(const T s, const DiagonalMatrix<T>& A) 
{
  DiagonalMatrix<T> R(A);
  R *= s;
  return R;
}
//-----------------------------------------------------------------------------
// DiagonalMatrix addition
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T> operator+(const DiagonalMatrix<T> &A, const DiagonalMatrix<T> &B)
{
  DiagonalMatrix<T> R(A);
  R+=B;
  return R;
}
//-----------------------------------------------------------------------------
// DiagonalMatrix subtraction
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T> operator-(const DiagonalMatrix<T> &A, const DiagonalMatrix<T> &B)
{
  DiagonalMatrix<T> R(A);
  return R-=B;
}
//-----------------------------------------------------------------------------
// template member definitions
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Default constructor - optionally zeros the matrix
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T>::DiagonalMatrix(INDEX rows, bool zero)
 : _data(NULL)
{
  reset(rows, zero);
}
//-----------------------------------------------------------------------------
// copy constructor - makes a full copy
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T>::DiagonalMatrix(const DiagonalMatrix<T>& c)
 : _data(NULL)
{
  reset(c);
}
//-----------------------------------------------------------------------------
// copy constructor from vector
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T>::DiagonalMatrix(const Vector<T>& v)
 : _data(NULL)
{
  reset(v);
}
//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T>::~DiagonalMatrix() 
{
  _delete();
} 
//-----------------------------------------------------------------------------
// deletes the data stored by this matrix
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::_delete()
{
  if (_data) delete _data;
}
//-----------------------------------------------------------------------------
// resizes the matrix, ignores nCols, optionally zeros 
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::reset(INDEX rows, INDEX cols, bool zero)
{
  _delete();
  _data = new DenseVector<T>(rows, zero);
}
//-----------------------------------------------------------------------------
// resizes the matrix, ignores nCols, optionally copies what fits
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::resize(INDEX rows, INDEX cols, bool copy)
{
  _data->resize(rows, copy);  
}
//-----------------------------------------------------------------------------
// changes the diagonal of the matrix to a vector v (makes a copy)
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::reset(const Vector<T>& v) 
{
  if (&v == _data) return; // check for self-reset
  _delete();
  _data = new DenseVector<T>(v);
}
//-----------------------------------------------------------------------------
// copys from another DiagonalMatrix
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::reset(const DiagonalMatrix<T>& c) 
{
  reset(*(c._data));
}
//-----------------------------------------------------------------------------
// resizes the matrix and copies data
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::copy(const T * ptr, INDEX rows, INDEX cols)
{
  if (_data) _data->reset(rows, false);
  else _data = new DenseVector<T>(rows, false);
  memcpy(_data, ptr, this->size()*sizeof(T));
}
//-----------------------------------------------------------------------------
// shallow reset from another DiagonalMatrix
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::shallowreset(const DiagonalMatrix<T> &c)
{
  _delete();
  _data = new CloneVector<T>(*(c._data));
}
//-----------------------------------------------------------------------------
// shallow reset from Vector
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::shallowreset(const Vector<T> &v)
{
  _delete();
  _data = new CloneVector<T>(v);
}
//-----------------------------------------------------------------------------
// reference indexing operator - must throw an error if i!=j
//-----------------------------------------------------------------------------
template<typename T>
T& DiagonalMatrix<T>::operator()(INDEX i, INDEX j)
{
  GCK(*this,*this,i!=j,"DiagonalMatrix: tried to index off diagonal");
  return VIDX(i);
}
//-----------------------------------------------------------------------------
// value indexing operator - returns 0 if i!=j
//-----------------------------------------------------------------------------
template<typename T>
T DiagonalMatrix<T>::operator()(INDEX i, INDEX j)                         const 
{
  return (i==j) ? (*_data)(i) : (T)0; 
}
//-----------------------------------------------------------------------------
// flat reference indexing operator
//-----------------------------------------------------------------------------
template<typename T>
T& DiagonalMatrix<T>::operator[](INDEX i)                
{
  return (*_data)(i); 
}
//-----------------------------------------------------------------------------
// flat value indexing operator
//-----------------------------------------------------------------------------
template<typename T>
T DiagonalMatrix<T>::operator[](INDEX i)                                  const          
{ 
  return (*_data)(i); 
}
//-----------------------------------------------------------------------------
// returns the number of rows
//-----------------------------------------------------------------------------
template<typename T>
INDEX DiagonalMatrix<T>::nRows()                                          const          
{ 
  return _data->size(); 
}
//-----------------------------------------------------------------------------
// returns the number of columns (same as nCols())
//-----------------------------------------------------------------------------
template<typename T>
INDEX DiagonalMatrix<T>::nCols()                                          const  
{
  return _data->size(); 
}
//-----------------------------------------------------------------------------
// returns a pointer to the diagonal values, dangerous!
//-----------------------------------------------------------------------------
template<typename T>
T* DiagonalMatrix<T>::get_ptr()                                           const
{
  return _data->get_ptr(); 
} 
//-----------------------------------------------------------------------------
// writes the diagonal to a binary data restart file
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::write_restart(FILE *f)                            const 
{
  _data->write_restart(f);
}
//-----------------------------------------------------------------------------
// sets the diagonal to a constant
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(const T v) 
{
  this->set_all_elements_to(v);
  return *this;
}
//-----------------------------------------------------------------------------
// assignment operator with another diagonal matrix
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(const DiagonalMatrix<T>& C) 
{
  reset(C);
  return *this;
}
//-----------------------------------------------------------------------------
// writes a matlab command to duplicate this sparse matrix
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::matlab(ostream &o, const string &s)               const 
{
  _data->matlab(o, s);
  o << s <<"=diag("<<s<<",0);\n";
}
//-----------------------------------------------------------------------------
// inverts this matrix, returns a reference
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T>& DiagonalMatrix<T>::inv_this()
{
  for(INDEX i=0; i<nRows(); i++) 
  {
     if (VIDX(i)!=T(0)) VIDX(i) = 1.0/VIDX(i);
     else 
     {
        cout<<"DiagonalMatrix::inv(): ("<<i<<","<<i<<")=0\n";
     }
  }  
  // Error check info
  const double min_max = _data->minabs() / _data->maxabs();
  if (min_max > 1e-14)   return *this;
  cout << "DiagonalMatrix::inv_this(): Warning: Matrix is badly scaled.";
  cout << "  RCOND = "<<min_max<<"\n";
  return *this;
}
//-----------------------------------------------------------------------------
// returns the inverse of this matrix
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T> DiagonalMatrix<T>::inv() const
{
  DiagonalMatrix<T> invA(*this); // Make copy of A to invert

  for(INDEX i=0; i<invA.nRows(); i++) 
  {
     if (VIDX(i)!=T(0)) invA[i]=1.0/VIDX(i);
     else 
     {
        cout<<"DiagonalMatrix::inv(): ("<<i<<","<<i<<")=0\n";
     }
  }
  // Error check info
  const double min_max = _data->minabs() / _data->maxabs();
  if (min_max > 1e-14)   return invA;
  cout << "DiagonalMatrix::inv(): Warning: Matrix is badly scaled.";
  cout << "  RCOND = "<<min_max<<"\n";
  return invA;
}
//-----------------------------------------------------------------------------
// computes a matrix inverse
//-----------------------------------------------------------------------------
inline DiagonalMatrix<double> inv(const DiagonalMatrix<double>& A)
{
  return A.inv();
}
//-----------------------------------------------------------------------------
// general diagonalmatrix assigment 
//-----------------------------------------------------------------------------
template<typename T>
void DiagonalMatrix<T>::_set_equal(const Matrix<T> &r)
{
  this->resize(r.nRows(), r.nCols());
  const Matrix<T> *pr = &r;

  const SparseMatrix<T>   *ps = dynamic_cast<const SparseMatrix<T>*>   (pr);
  const DiagonalMatrix<T> *pd = dynamic_cast<const DiagonalMatrix<T>*> (pr);
  const Vector<T>         *pv = dynamic_cast<const Vector<T>*>         (pr);

  if (ps)       this->reset(ps->get_diag());
  else if (pd)  this->reset(*pd);
  else if (pv)  this->reset(*pv);
  else
  {
    cout <<"Error in general sparse matrix assignment\n";
  }
}
//-----------------------------------------------------------------------------
// casts a generic matrix pointer into a DiagonalMatrix pointer - null if fail 
//-----------------------------------------------------------------------------
template<typename T>
const DiagonalMatrix<T> *diag_cast(const Matrix<T> *m)
{
  return dynamic_cast<const DiagonalMatrix<T>*>(m);
}

#endif
