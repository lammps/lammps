#ifndef DENSEVECTOR_H
#define DENSEVECTOR_H

#include "Vector.h"

namespace ATC_matrix {

template<typename T>

  /**
   *  @class  DenseVector
   *  @brief  Class for storing data in a "dense" vector form
   */

class DenseVector : public Vector<T>
{
public:
  explicit DenseVector(INDEX n=0, bool z=1)          { _create(n,z); }
  DenseVector(const DenseVector<T> &c) : Vector<T>(), _data(nullptr) { _copy(c); }
  DenseVector(const Vector<T> &c)      : Vector<T>(), _data(nullptr) { _copy(c); }
  DenseVector(const T * ptr, INDEX nrows) : Vector<T>(), _data(nullptr) { copy(ptr,nrows); }
  virtual ~DenseVector()               { _delete();    }

  //* resizes the Vector, ignores nCols, optionally copys what fits
  void resize(INDEX rows, INDEX cols=1, bool copy=false);
  //* resizes the Vector, ignores nCols, optionally zeros it out
  void reset (INDEX rows, INDEX cols=1, bool zero=true);
  //* resizes the Vector and copies data, ignores nCols
  void copy(const T * ptr, INDEX rows, INDEX cols=1);

  // overloaded inline virtual functions
  T  operator[](INDEX i) const { VICK(i) return _data[i]; }
  T& operator[](INDEX i)       { VICK(i) return _data[i]; }
  T  operator()(INDEX i, INDEX /* j */) const { VICK(i) return _data[i]; }
  T& operator()(INDEX i, INDEX /* j */)       { VICK(i) return _data[i]; }
  T  operator()(INDEX i) const { VICK(i) return _data[i]; }
  T& operator()(INDEX i)       { VICK(i) return _data[i]; }
  void set_all_elements_to(const T &v)    {
                                            int sz = this->size();
                                            for (INDEX i = 0; i < sz; i++) _data[i] = v;
                                          }
  INDEX nRows()    const { return _size; }

  T* ptr() const     { return _data; }

  DenseVector<T>& operator=(const T &v);
  DenseVector<T>& operator=(const Vector<T> &c);
  DenseVector<T>& operator=(const DenseVector<T> &c);

  void write_restart(FILE *f) const;

private:
  void _delete();
  void _create(INDEX n, bool zero=0);
  void _copy(const Vector<T> &c);

  T *_data;
  INDEX _size;
};

///////////////////////////////////////////////////////////////////////////////
// Template definitions ///////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// resizes the matrix and optionally copies over what still fits, ignores cols
//-----------------------------------------------------------------------------
template <typename T>
  void DenseVector<T>::resize(INDEX rows, INDEX /* cols */, bool copy)
{
  if (_size==rows) return;  // if is correct size, done
  if (!copy)
  {
     _delete();
     _create(rows);
     return;
  }
  DenseVector<T> temp(*this);
  _delete();
  _create(rows);
  int sz = this->size();
  for (INDEX i = 0; i < sz; i++)
    _data[i] = i<temp.size() ? temp[i] : T(0.0);
  return;
}
///////////////////////////////////////////////////////////////////////////////
//* resizes the matrix and optionally zeros it out
template <typename T>
void DenseVector<T>::reset(INDEX rows, INDEX /* cols */, bool zero)
{
  if (_size!=rows)
  {
     _delete();
     _create(rows);
  }
  if (zero) this->zero();
}
///////////////////////////////////////////////////////////////////////////////
//* resizes the matrix and optionally zeros it out
template <typename T>
void DenseVector<T>::copy(const T * ptr, INDEX rows, INDEX /* cols */)
{
  resize(rows, 1, false);
  memcpy(_data, ptr, this->size()*sizeof(T));
}
///////////////////////////////////////////////////////////////////////////////
//* writes the matrix data to a file
template <typename T>
void DenseVector<T>::write_restart(FILE *f)                               const
{
  fwrite(&_size, sizeof(INDEX),1,f);
  if(_size) fwrite(_data, sizeof(T), _size, f);
}
///////////////////////////////////////////////////////////////////////////////
//* clears allocated memory
template <typename T>
inline void DenseVector<T>::_delete()
{
  if (_data) delete [] _data;
  _size = 0;
}
///////////////////////////////////////////////////////////////////////////////
//* allocates memory for an rows by cols DenseMatrix
template <typename T>
inline void DenseVector<T>::_create(INDEX n, bool zero)
{
  _size=n;
  _data = _size ? new T [_size] : nullptr ;
  if (zero) this->zero();
}
///////////////////////////////////////////////////////////////////////////////
//* creates a deep memory copy from a general matrix
template <typename T>
inline void DenseVector<T>::_copy(const Vector<T> &c)
{
  if (!_data || _size!=c.size())
  {
    _delete();
    _create(c.size(), false);
  }
  else _size = c.size();
  memcpy(_data, c.ptr(), _size*sizeof(T));
}
///////////////////////////////////////////////////////////////////////////////
//* assigns v to all values in the vector
template <typename T>
DenseVector<T>& DenseVector<T>::operator=(const T &v)
{
  int sz = this->size();
  for (INDEX i = 0; i < sz; i++) (*this)[i] = v;
  return *this;
}
///////////////////////////////////////////////////////////////////////////////
//* copys c with a deep copy
template <typename T>
DenseVector<T>& DenseVector<T>::operator=(const Vector<T> &c)
{
  _copy(c);
  return *this;
}
///////////////////////////////////////////////////////////////////////////////
//* copys c with a deep copy
template <typename T>
DenseVector<T>& DenseVector<T>::operator=(const DenseVector<T> &c)
{
  _copy(c);
  return *this;
}

} // end namespace

#endif
