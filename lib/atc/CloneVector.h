#ifndef CLONEVECTOR_H
#define CLONEVECTOR_H

#include "Vector.h"

namespace ATC_matrix {

  /**
   *  @class  CloneVector 
   *  @brief  Class for creating objects that wrap matrix data for manipulation through vector operations
   */

template<typename T>
class CloneVector : public Vector<T>
{
public:
  CloneVector(); // do not implement
  CloneVector(const Vector<T> &c);
  CloneVector(const Matrix<T> &c, int dim, INDEX idx=0);
  CloneVector(const DiagonalMatrix<T> &c, INDEX idx=0);

  // overloaded virtual functions
  T& operator[](INDEX i);
  T  operator[](INDEX i)              const;
  T  operator()(INDEX i, INDEX j=0)   const;
  T& operator()(INDEX i, INDEX j=0);
  INDEX nRows()                       const;

  CloneVector<T>& operator=(const T &v);
  CloneVector<T>& operator=(const CloneVector<T> &C);
  CloneVector<T>& operator=(const Matrix<T> &C);

  virtual bool memory_contiguous()    const;
  T* ptr()             const;
  void resize(INDEX nRows, INDEX nCols=0, bool copy=false);
  void  reset(INDEX nRows, INDEX nCols=0, bool zero=true);
  void copy(const T * ptr, INDEX nRows, INDEX nCols=0);

private:
  void _resize(INDEX nRows, INDEX nCols, bool copy, bool zero);
 
  Vector<T> * const _baseV; // ptr to a base vector
  Matrix<T> * const _baseM; // ptr to a base matrix
  int _clone_type;          // what to clone (see enum CLONE_TYPE)
  INDEX _idx;               // index of matrix dimension to clone
};
///////////////////////////////////////////////////////////////////////////////
// Template definitions ///////////////////////////////////////////////////////
//-----------------------------------------------------------------------------
// Construct from another vector
//-----------------------------------------------------------------------------
template<typename T>
CloneVector<T>::CloneVector(const Vector<T> &c)
 : Vector<T>(), _baseV(const_cast<Vector<T>*>(&c)), _baseM(NULL) 
{}
//-----------------------------------------------------------------------------
// Construct from a matrix, the const_cast isn't pretty
/* CloneVector(const Matrix<T> &c, int dim, INDEX idx)
/   attaches to a slice of a matrix
/   Arguments: c = pointer to the matrix
/              dim = type of slice CLONE_ROW, CLONE_COL, CLONE_DIAG
/              idx = index of row or column (no effect on diag currently)
*/
//-----------------------------------------------------------------------------
template<typename T>
CloneVector<T>::CloneVector(const Matrix<T> &c, int dim, INDEX idx) 
 : Vector<T>(), _baseV(NULL), _baseM(const_cast<Matrix<T>*>(&c))
 , _clone_type(dim), _idx(idx)
{}
//-----------------------------------------------------------------------------
// Construct from a DiagonalMatrix
//-----------------------------------------------------------------------------
template<typename T>
CloneVector<T>::CloneVector(const DiagonalMatrix<T> &c, INDEX idx) 
 : Vector<T>(), _baseV(NULL), _baseM(const_cast<DiagonalMatrix<T>*>(&c))
 , _clone_type(CLONE_DIAG), _idx(0)
{}
//-----------------------------------------------------------------------------
// value (const) indexing operator
//-----------------------------------------------------------------------------
template<typename T>
T CloneVector<T>::operator()(INDEX i, INDEX j) const 
{
  return (*this)[i];
}
//-----------------------------------------------------------------------------
// reference index operator
//-----------------------------------------------------------------------------
template<typename T>
T& CloneVector<T>::operator()(INDEX i, INDEX j)
{
  return (*this)[i];
}
//-----------------------------------------------------------------------------
// Indexes the cloned vector either from another vector or a matrix
//-----------------------------------------------------------------------------
template<typename T>
T CloneVector<T>::operator[](INDEX i)                                     const
{
  if (_baseV)                         return (*_baseV)(i);
  if      (_clone_type == CLONE_ROW)  return (*_baseM)(_idx, i);
  else if (_clone_type == CLONE_COL)  return (*_baseM)(i,_idx);
  else if (_clone_type == CLONE_DIAG) return (*_baseM)(i,i);
  return 0;
}
//-----------------------------------------------------------------------------
// Indexes the cloned vector either from another vector or a matrix
//-----------------------------------------------------------------------------
template<typename T>
T& CloneVector<T>::operator[](INDEX i)
{
  if (_baseV)                    return (*_baseV)(i);
  if (_clone_type == CLONE_ROW)  return (*_baseM)(_idx, i);
  if (_clone_type == CLONE_COL)  return (*_baseM)(i,_idx);
  if (_clone_type == CLONE_DIAG) return (*_baseM)(i,i);
  return (*_baseV)(i);
}
//-----------------------------------------------------------------------------
// Returns the size of the base vector or of the row/col of the base matrix
//-----------------------------------------------------------------------------
template<typename T>
INDEX CloneVector<T>::nRows()                                             const
{
  using std::min;
  if (_baseV)                    return _baseV->size();
  if (_clone_type == CLONE_ROW)  return _baseM->nCols();
  if (_clone_type == CLONE_COL)  return _baseM->nRows();
  if (_clone_type == CLONE_DIAG) return min(_baseM->nRows(), _baseM->nCols());
  return 0;
}
//-----------------------------------------------------------------------------
// assigns all elements to a constant
//-----------------------------------------------------------------------------
template<typename T>
CloneVector<T>& CloneVector<T>::operator=(const T &v)
{
  this->set_all_elements_to(v);   
  return *this;
}
//-----------------------------------------------------------------------------
// assigns all elements to the corresponding elements in C
//-----------------------------------------------------------------------------
template<typename T>
CloneVector<T>& CloneVector<T>::operator=(const CloneVector<T> &C)
{
  GCK(*this, C, this->size()!=C.size(), "Error in CloneVector:operator=");
  int sz = this->size(); 
  for (INDEX i = 0; i < sz; i++) (*this)[i] = C[i];
  return *this;
}
//-----------------------------------------------------------------------------
// assigns all elements to the corresponding elements in C
//-----------------------------------------------------------------------------
template<typename T>
CloneVector<T>& CloneVector<T>::operator=(const Matrix<T> &C)
{
  GCK(*this, C, this->size()!=C.size(), "Error in CloneVector:operator=");
  int sz = this->size(); 
  for (INDEX i = 0; i < sz; i++) (*this)[i] = C[i];
  return *this;
}
//-----------------------------------------------------------------------------
// returns true only if its guaranteed memory is contiguous
//-----------------------------------------------------------------------------
template<typename T>
bool CloneVector<T>::memory_contiguous()  const
{
  // drill down through clone of clones
  if (_baseV) return _baseV->memory_contiguous();
  // could be okay if DiagonalMatrix, but can't guarantee this
  if (_clone_type == CLONE_DIAG) return false; 
#ifdef ROW_STORAGE
  return _clone_type == CLONE_ROW;
#else
  return _clone_type == CLONE_COL;
#endif
}
//-----------------------------------------------------------------------------
// Returns a pointer to the data unless the data is a column of a matrix
//-----------------------------------------------------------------------------
template<typename T>
T* CloneVector<T>::ptr()                                              const
{
  if (_baseV) return _baseV->ptr();
#ifdef ROW_STORAGE
  if (_clone_type == CLONE_ROW)  return  _baseM->ptr() + this->size()*_idx;
  if (_clone_type == CLONE_COL)  return _baseM->ptr() + this->size();
  if (_clone_type == CLONE_DIAG) return _baseM->ptr();
#else
  if (_clone_type == CLONE_COL)  return _baseM->ptr() + this->size()*_idx;
  if (_clone_type == CLONE_ROW)  return _baseM->ptr() + this->size();
  if (_clone_type == CLONE_DIAG) return _baseM->ptr();
#endif
  return 0;
}
//-----------------------------------------------------------------------------
// general resize function, can handle parents that are matrices or vectors
//-----------------------------------------------------------------------------
template<typename T>
void CloneVector<T>::_resize(INDEX nRows, INDEX nCols, bool copy, bool zero)
{
  if (_baseV) 
  {
    if (copy) _baseV->resize(nRows, nCols, copy);
    else      _baseV->reset (nRows, nCols, zero);
    return;
  }
  // parent is a matrix, need to decide what the Vector is cloning
  switch (_clone_type)
  {
    case CLONE_ROW:  // now the leading dimension is rows
      nCols = nCols ? nCols : _baseM->nCols();
      break;
    case CLONE_COL:  // now the leading dimension is columns
      nCols = nCols ? nCols : _baseM->nRows();
      ATC_Utility::swap(nRows, nCols);
      break;
    case CLONE_DIAG: // lets just hope you knew what you were doing
      break;
    default:
      return;
  }
  if (zero)  _baseM->reset(nRows, nCols, zero);  // zero overrides copy
  else      _baseM->resize(nRows, nCols, copy);
}
//-----------------------------------------------------------------------------
// resizes the matrix and optionally copies what fits
//-----------------------------------------------------------------------------
template<typename T>
void CloneVector<T>::resize(INDEX nRows, INDEX nCols, bool copy)
{
  _resize(nRows, nCols, copy, false);
}
//-----------------------------------------------------------------------------
// resizes the matrix and optionally zeros it out 
//-----------------------------------------------------------------------------
template<typename T>
void CloneVector<T>::reset(INDEX nRows, INDEX nCols, bool zero)
{
  _resize(nRows, nCols, false, zero);
}
//-----------------------------------------------------------------------------
// resizes the matrix and copies data
//-----------------------------------------------------------------------------
template<typename T>
void CloneVector<T>::copy(const T * ptr, INDEX nRows, INDEX nCols)
{
  _resize(nRows, nCols, false, false);
  memcpy(this->ptr(), ptr, this->size()*sizeof(T));
}

} // end namespace
#endif
