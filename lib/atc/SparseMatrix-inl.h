#ifndef SPARSEMATRIX_INL_H
#define SPARSEMATRIX_INL_H

#include "mpi.h"
#include "DenseVector.h"

namespace ATC_matrix {

template <typename T> 
TRI_COORD<T>::TRI_COORD(INDEX row, INDEX col) : i(row), j(col) {}
template <typename T> 
TRI_COORD<T>::TRI_COORD(INDEX row, INDEX col, T val, bool add_to) 
  : i(row), j(col), v(val), add(add_to) {}
  
//-----------------------------------------------------------------------------
// default constructor - creates an empty sparsematrix with specified size
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T>::SparseMatrix(INDEX rows, INDEX cols) 
 : _val(NULL), _ia(NULL), _ja(NULL), _size(0), _nRowsCRS(0), hasTemplate_(false),
   _nRows(rows),_nCols(cols) {}
//-----------------------------------------------------------------------------
// copy constructor
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix<T>& C)
 : _val(NULL), _ia(NULL), _ja(NULL), hasTemplate_(false)
{
  _copy(C); 
}
//-----------------------------------------------------------------------------
// copy constructor - converts from DenseMatrix
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T>::SparseMatrix(const DenseMatrix<T>& C) 
: _val(NULL), _ia(NULL), _ja(NULL), hasTemplate_(false)
{
  reset(C); 
}

//-----------------------------------------------------------------------------
// constructor - creates a sparse matrix given an array of row indeces,
// an array of col indeces, and an array of nonzero values.
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T>::SparseMatrix(INDEX* rows, INDEX* cols, T* vals, 
                              INDEX size, INDEX nRows, INDEX nCols, INDEX nRowsCRS)
 : hasTemplate_(true) 
{
  _val = vals;
  _ia = rows;
  _ja = cols;
  _size = size;
  _nRows = nRows;
  _nCols = nCols;
  _nRowsCRS = nRowsCRS;
}

//-----------------------------------------------------------------------------
// assigns internal storage for CRS
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::_create(INDEX size, INDEX nrows)
{
  _size     = size;
  _nRowsCRS = nrows;
  // assign memory to hold matrix
  try
  {
    _val = (_size*nrows) ? new T     [_size]        : NULL;
    _ia  = (_size*nrows) ? new INDEX [_nRowsCRS+1]  : NULL;
    _ja  = (_size*nrows) ? new INDEX [_size]        : NULL;
  }
  catch (std::exception &e)
  {
    cout << "Could not allocate SparseMatrix of "<< _size << " nonzeros.\n";
    ERROR_FOR_BACKTRACE
      exit(EXIT_FAILURE);
  }
  if (!_ia) return;  
  // automatically handle the ends of rowpointer
  *_ia = 0;               // first non-zero is the zero index
  _ia[_nRowsCRS] = _size; // last row pointer is the size
}
//-----------------------------------------------------------------------------
// cleans up internal storage, but retains nRows & nCols
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::_delete()
{
  
  vector<TRI_COORD<T> >().swap(_tri); // completely deletes _tri
  if (_val) delete [] _val;
  if (_ia)  delete [] _ia;
  if (_ja)  delete [] _ja;
  _size = _nRowsCRS = 0;
  _val = NULL;
  _ia = _ja = NULL;
}
//-----------------------------------------------------------------------------
// full memory copy of C into this
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::_copy(const SparseMatrix<T> &C) 
{
  compress(C);
  _delete();
  _create(C.size(), C._nRowsCRS);
  if (_size) {
    std::copy(C._val, C._val+_size, _val);
    std::copy(C._ja,  C._ja+_size,  _ja);
  }
  if (_nRowsCRS)  {
    std::copy(C._ia,  C._ia+_nRowsCRS+1, _ia);
  }
  _nCols = C._nCols;
  _nRows = C._nRows;
  if (_nCols > 0 && _nRows > 0) hasTemplate_ = true; // needs if since map seems to call the copy instead of the default constructor
} 
// this version is accessible to derived classes
template<typename T>
void SparseMatrix<T>::copy(const SparseMatrix<T> &C)
{
  _copy(C);
}
//----------------------------------------------------------------------------
// general sparse matrix assignment
//----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::_set_equal(const Matrix<T> &r)
{
  this->resize(r.nRows(), r.nCols());
  const Matrix<T> *ptr_r = &r;

  const SparseMatrix<T> *s_ptr = dynamic_cast<const SparseMatrix<T>*>(ptr_r);
  if (s_ptr) this->reset(*s_ptr);
  else if (dynamic_cast<const DiagonalMatrix<T>*>(ptr_r))
    for (INDEX i=0; i<r.size(); i++) set(i,i,r[i]);
  else if (dynamic_cast<const DenseMatrix<T>*>(ptr_r))  this->reset(r);
  else
  { 
    cout <<"Error in general sparse matrix assignment\n";
    exit(1);
  }
} 
// General flat index by value operator (by nth nonzero)
template <typename T> inline T SparseMatrix<T>::operator[](INDEX i) const 
{
  VICK(i); return _val[i]; 
}

// General flat index by reference operator (by nth nonzero)
template <typename T> inline T& SparseMatrix<T>::operator[](INDEX i)
{
  VICK(i); return _val[i]; 
}

template<typename T>
T SparseMatrix<T>::_zero = T(0);

//-----------------------------------------------------------------------------
// triplet comparison operator returns true if x < y
//-----------------------------------------------------------------------------
template <typename T>
bool triplet_comparision(const TRI_COORD<T> &x, const TRI_COORD<T> &y)
{ 
  const bool row_less  = (x.i) <  (y.i);
  const bool row_equal = (x.i) == (y.i);
  const bool col_less  = (x.j) <  (y.j);
  return (row_less || (row_equal && col_less));
}
//-----------------------------------------------------------------------------
// triplet comparison operator returns true if x == y
//-----------------------------------------------------------------------------
template <typename T>
bool triplets_equal(const TRI_COORD<T> &x, const TRI_COORD<T> &y)
{ 
  return x.i==y.i && x.j==y.j;
}
//-----------------------------------------------------------------------------
// multiply sparse matrix by a vector
//-----------------------------------------------------------------------------
template<typename T>
DenseVector<T> operator*(const SparseMatrix<T> &A, const Vector<T>& x)
{
  DenseVector<T> y(A.nRows(), true);

  A.MultMv(x, y);

  return y;
}
//-----------------------------------------------------------------------------
// multiply a vector by a sparse matrix
//-----------------------------------------------------------------------------
template<typename T>
DenseVector<T> operator*(const Vector<T>& x, const SparseMatrix<T> &A)
{
  return A.transMat(x);
}
//-----------------------------------------------------------------------------
// multiply sparse matrix by dense matrix
//-----------------------------------------------------------------------------
template<typename T>
DenseMatrix<T> operator*(const SparseMatrix<T> &A, const Matrix<T>& D)
{
  DenseMatrix<T> C(A.nRows(), D.nCols(), true);  // initialized to zero

  A.MultAB(D, C);

  return C;
}
//-----------------------------------------------------------------------------
// multiply sparse matrix by a diagonal matrix - scales each column
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T> operator*(const SparseMatrix<T> &A, const DiagonalMatrix<T>& D)
{
  GCK(A, D, A.nCols()!=D.nRows(),"SparseMatrix * DiagonalMatrix")
  SparseMatrix<T> C(A);  // C has same sparcity as A
 
  // C(i,j) = A(i,k) * D(k, j) * j==k
  INDEX i, ij;
  for (i=0; i<A._nRowsCRS; i++) 
    for (ij=A._ia[i]; ij<A._ia[i+1]; ij++) 
      C[ij] = A._val[ij]*D(A._ja[ij],A._ja[ij]); 

  return C;
}
//-----------------------------------------------------------------------------
// multiplies two sparse matrices - assumes their output is sparse
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T> operator*(const SparseMatrix<T> &A, const SparseMatrix<T> &B)
{
  SparseMatrix<T> At(A.transpose());
  SparseMatrix<T>::compress(B);

  GCK(A, B, A.nCols()!=B.nRows(), "SparseMatrix * SparseMatrix");

  SparseMatrix<T> C(A.nRows(), B.nCols());
  if (At.empty() || B.empty()) return C;
  
  INDEX k, ki, kj;
  INDEX K = std::min(At._nRowsCRS, B._nRowsCRS);
  for (k=0; k<K; k++)   // loop over rows of A or B (smallest)
    for (ki=At._ia[k]; ki<At._ia[k+1]; ki++) // loop over row nonzeros of A
      for (kj=B._ia[k]; kj<B._ia[k+1]; kj++) // loop over row nonzeros of B
        C.add(At._ja[ki], B._ja[kj], At[ki]*B[kj]); // C(i,j) = At(k,i)*B(k, j)

  C.compress();
  return C;
}
//-----------------------------------------------------------------------------
// returns the first row number with a nonzero entry or -1 if no rows
//-----------------------------------------------------------------------------
template<typename T>
int SparseMatrix<T>::_first_nonzero_row_crs()                             const
{
  if (!_nRowsCRS) return -1;
  INDEX r;
  for (r=0; r<_nRowsCRS; r++)
    if (_ia[r+1]>0) return r;
  return -1;
}
//-----------------------------------------------------------------------------
// converts T to CRS
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::compress(const SparseMatrix<T> &C)
{
  const_cast<SparseMatrix<T>*>(&C)->compress();
}
//-----------------------------------------------------------------------------
// merges all the _tri triples with CRS storage
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::compress()
{
  if (_tri.empty()) return;

  // Sort and find the number of unique triplets.
  // Triplet values will all be not present in existing CRS structure.
  const INDEX nUnique = CountUniqueTriplets();  

  // Max number of rows in new CRS structure.
  const INDEX nRows = std::max((INDEX)_tri.back().i+1, _nRowsCRS);
  
  // make a new CRS structure
  INDEX *ia = new INDEX [nRows+1];
  INDEX *ja = new INDEX [nUnique];
  T    *val = new     T [nUnique];

  // Set first and last row ptr to 0 and nnz respectively.
  // Set all else to a flagvalue MAX_UNSIGNED (~0).
  ia[0] = 0;
  INDEX i;
  for (i=1; i<nRows; i++) ia[i]=~0;  // ~0 is max(INDEX)
  ia[nRows] = nUnique;
  
  INDEX crs_pt, crs_row;
  unsigned tri_ct; // must be unsigned to interface with std::vector without warnings

  // Get the first CRS and triplet coordinates (if they exist).
  TRI_COORD<T> nextCRS, nextTRI(_tri[0]), next;
  int first_row = _first_nonzero_row_crs();
  if (first_row != -1)  nextCRS = TRI_COORD<T>(first_row, _ja[0], _val[0]);

  // merge sorted triplets into a new CRS structure
  crs_pt = crs_row = tri_ct = 0;  // initialize counters
  for (i=0; i<nUnique; i++)
  {
    // is the next non-zero in the new triplet vector
    if (tri_ct < _tri.size() 
        && (triplet_comparision(nextTRI, nextCRS) || crs_pt>=_size)) {
      next = nextTRI;
      // advance the triplet counter, and skip voided TRIPLET entries
      do tri_ct++;
      while ( tri_ct<_tri.size() && _tri[tri_ct].j == ~0 );

      // if not at the end of the vector, set the next triplet
      if (tri_ct<_tri.size()) nextTRI = _tri[tri_ct];
    }
    // is the next nonzero in the old CRS data
    else if (crs_pt < _size) {
      next = nextCRS;
      // Advance the CRS counter, don't set next if we are at the end.
      if (++crs_pt < _size) {  
        // advance to the row corresponding to this value
        while (crs_pt >= _ia[crs_row+1]) {
          crs_row++;
        }
        nextCRS = TRI_COORD<T>(crs_row, _ja[crs_pt], _val[crs_pt]);
      }
    }
    else cout << "SparseMatrix - Error in compressing CRS\n";

    // Add next to the new CRS structure.
    // Is this a new row (is j>0 and is ja[j] == 0)?
    if (ia[next.i]==~0) ia[next.i] = i;
    ja[i]  = next.j;
    val[i] = next.v;
  }
  // sweep backwards through row pointers and check for skipped rows
  for (i=nRows-1; i>0; i--) ia[i] = (ia[i]==~0) ? ia[i+1] : ia[i];

  _delete();
  _val      = val;
  _ia       = ia;
  _ja       = ja;
  _size     = nUnique;
  _nRowsCRS = nRows;
  hasTemplate_=true;
}
//-----------------------------------------------------------------------------
// Sorts the triplets, condenses duplicates, and returns the # of unique values
//-----------------------------------------------------------------------------
template<typename T>
INDEX SparseMatrix<T>::CountUniqueTriplets()
{
  if (_tri.empty()) return _size;
  std::sort(_tri.begin(), _tri.end(), triplet_comparision<T>);
  INDEX nUnique=1 + _size;
  
  typename vector<TRI_COORD<T> >::reverse_iterator t;  
  // Loop backwards over all new triplets.
  for (t = _tri.rbegin(); t+1!=_tri.rend();  ++t) {
    // If this triplet is the same as the preceding one.
    if (triplets_equal(*(t+1), *t))  {
      if (t->add) (t+1)->v += t->v;   // Add to previous
      else        (t+1)->v  = t->v;   // Replace previous -- DOES THIS WORK?
      t->j = ~0;                      // Void this entry's column pointer
    }
   else nUnique++;
  }
  return nUnique;
}
//-----------------------------------------------------------------------------
// Checks if a value has been set
//-----------------------------------------------------------------------------
template<typename T>
bool SparseMatrix<T>::has_entry(INDEX i, INDEX j) const
{
  if (has_entry_compressed(i,j)) return true;
  if (has_entry_uncompressed(i,j)) return true;
  return false;
}

template<typename T>
bool SparseMatrix<T>::has_entry_uncompressed(INDEX i, INDEX j) const
{
  for (unsigned k=0; k<_tri.size() ; k++) {  
    if (_tri[k].i == i && _tri[k].j == j) return true;
  }
  return false;
}

template<typename T>
bool SparseMatrix<T>::has_entry_compressed(INDEX i, INDEX j) const
{
  if (_size == 0) return false;
  if (i >= _nRowsCRS) return false;
  if (_ia[i] < _ia[i+1]) {
    return -1 < ATC_Utility::search_sorted(_ja, j, _ia[i], _ia[i+1]);
  }
  return false;
}
//-----------------------------------------------------------------------------
// check if the matrix has been compressed at least once
//-----------------------------------------------------------------------------
template<typename T>
bool SparseMatrix<T>::has_template(void) const
{
  return hasTemplate_;
}
//-----------------------------------------------------------------------------
// Index by copy operator - return zero if not found
//-----------------------------------------------------------------------------
template<typename T>
T SparseMatrix<T>::operator()(INDEX i, INDEX j) const 
{
  MICK(i,j);  // Matrix Index ChecKing
  compress(*this);
  if (i>=_nRowsCRS || _ia[i+1]==_ia[i]) return 0.0;
  INDEX f = std::lower_bound(_ja+_ia[i], _ja+_ia[i+1]-1, j) - _ja;
  if (f>=_ia[i] && f<_ia[i+1] && _ja[f] == j) return _val[f];
  return 0.0;
}
//-----------------------------------------------------------------------------
// Index by reference operator - add to _tri if not found
//-----------------------------------------------------------------------------
template<typename T>
T& SparseMatrix<T>::operator()(INDEX i, INDEX j)
{
  MICK(i,j);  // Matrix Index ChecKing
  compress(*this);
  if (i < _nRowsCRS && _ia[i+1]>_ia[i]) {
    INDEX f = std::lower_bound(_ja+_ia[i], _ja+_ia[i+1]-1, j) - _ja;
    if (f>=_ia[i] && f<_ia[i+1] && _ja[f] == j) return _val[f];
  }
  // NEVER use index operator as LHS to modify values not already in the
  // sparcity pattern - the crude check below will only catch this on the 
  // second infraction.
  if (_zero != T(0)) cout << "Use add or set for SparseMatrix\n";
  return _zero;
}
//-----------------------------------------------------------------------------
// Sets (i,j) to value
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::set(INDEX i, INDEX j, T v)
{
  MICK(i,j);  // Matrix Index ChecKing
  if (i < _nRowsCRS)
  {
    const int loc = ATC_Utility::search_sorted(_ja, j, _ia[i], _ia[i+1]);
    if (loc >=0 )
    {
      _val[loc] = v;
      return;
    }
  }
  _tri.push_back(TRI_COORD<T>(i,j,v,false));
}
//-----------------------------------------------------------------------------
// Adds (i,j) to value
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::add(INDEX i, INDEX j, T v)
{
  MICK(i,j);  // Matrix Index ChecKing
  if (i < _nRowsCRS)
  {
    const int loc = ATC_Utility::search_sorted(_ja, j, _ia[i], _ia[i+1]);
    if (loc >=0 )
    {
      _val[loc] += v;
      return;
    }
  }
  _tri.push_back(TRI_COORD<T>(i,j,v,true));
}
//-----------------------------------------------------------------------------
// returns a triplet value of the ith nonzero
//-----------------------------------------------------------------------------
template<typename T>
TRIPLET<T> SparseMatrix<T>::triplet(INDEX i) const
{
  compress(*this);
  if (i >= _ia[_nRowsCRS]) {
    gerror("ERROR: tried indexing triplet of sparse matrix beyond range");
  }
    
  INDEX row(std::lower_bound(_ia, _ia+_nRowsCRS, i)-_ia);
  row -= _ia[row] != i;
  return TRIPLET<T>(row, _ja[i], _val[i]);
}
//-----------------------------------------------------------------------------
// full reset - completely wipes out all SparseMatrix data, zero is ignored
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::reset(INDEX rows, INDEX cols, bool zero)
{
  _delete();
  _nRows = rows;
  _nCols = cols;
}
//-----------------------------------------------------------------------------
// resize - changes the _nRows and _nCols without changing anything else if
//          the matrix is being enlarged, other wise wipes it
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::resize(INDEX rows, INDEX cols, bool copy)
{
//if (copy) throw; 
  if (_nRowsCRS>rows) {
    _delete();
  }
  if (copy) 
  _nRows = rows;
  _nCols = cols;  // a check on this would be expensive
}
//-----------------------------------------------------------------------------
// get sparsity from DenseMatrix, if TOL < 0, then only zero values are added
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::reset(const DenseMatrix<T>& D, double TOL) 
{
  _delete(); // clears all values
  // if TOL is specified then TOL = TOL^2 * max(abs(D))^2
  if (TOL > 0.0)
  {
    TOL *= D.maxabs();
    TOL *= TOL;
  }
  _nRows = D.nRows();
  _nCols = D.nCols();
  for (INDEX i=0; i<D.nRows(); i++)
    for (INDEX j=0; j<D.nCols(); j++)
      if (D(i,j)*D(i,j) >= TOL)  // if TOL wasn't specified then TOL < 0
        set(i, j, D(i,j));

  compress();
}
//-----------------------------------------------------------------------------
// copy - dangerous: ignores rows & columns
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::copy(const T * ptr, INDEX rows, INDEX cols)
{
  cout << "SparseMatrix<T>::copy() has no effect.\n";
  throw;
}
//-----------------------------------------------------------------------------
// dense_copy - copy to dense matrix
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::dense_copy(DenseMatrix <T> & D ) const
{
  SparseMatrix<T>::compress(*this);
  D.reset(nRows(),nCols());
  for (INDEX i=0; i<_nRowsCRS; i++) 
    for (INDEX j=_ia[i]; j<_ia[i+1]; j++)  
      D(i, _ja[j]) = _val[j];
}
template<typename T>
DenseMatrix <T> SparseMatrix<T>::dense_copy(void) const
{
  DenseMatrix<T> D;
  dense_copy(D);
  return D;
}
//-----------------------------------------------------------------------------
// returns true if the matrix has no non-zero elements
//-----------------------------------------------------------------------------
template<typename T>
bool SparseMatrix<T>::empty()                                             const 
{
  return _size==0 && _tri.empty(); 
}
//-----------------------------------------------------------------------------
// returns the number of rows specified by the user 
//-----------------------------------------------------------------------------
template<typename T>
inline INDEX SparseMatrix<T>::nRows()                                     const 
{
  return _nRows;
}
//-----------------------------------------------------------------------------
// returns ??????????????????????
//-----------------------------------------------------------------------------
template<typename T>
inline INDEX SparseMatrix<T>::nRowsCRS()                                     const 
{
  return _nRowsCRS;
}
//-----------------------------------------------------------------------------
// returns the number of columns specified by the user 
//-----------------------------------------------------------------------------
template<typename T>
inline INDEX SparseMatrix<T>::nCols()                                     const
{
  return _nCols;
}
//-----------------------------------------------------------------------------
// returns the number of non-zeros in the matrix
//-----------------------------------------------------------------------------
template<typename T>
INDEX SparseMatrix<T>::size()                                             const 
{
  compress(*this);
  return _size; 
}
//-----------------------------------------------------------------------------
// returns the number of nonzero elements in a row
//-----------------------------------------------------------------------------
template<typename T>
INDEX SparseMatrix<T>::RowSize(INDEX r)    const 
{
  compress(*this);
  GCHK(r>=_nRows, "Rowsize: invalid row");
  if (r >= _nRowsCRS) return 0;
  return _ia[r+1]-_ia[r]; 
}
//-----------------------------------------------------------------------------
// returns a pointer to the data, causes a compress
//-----------------------------------------------------------------------------
template<typename T>
T* SparseMatrix<T>::ptr()                                             const
{
  compress(*this);
  return _val; 
}
template<typename T>
INDEX* SparseMatrix<T>::rows()                                             const
{
  compress(*this);
  return _ia; 
}
template<typename T>
INDEX* SparseMatrix<T>::cols()                                             const
{
  compress(*this);
  return _ja; 
}
//-----------------------------------------------------------------------------
// returns true if (i,j) falls in the user specified range 
//-----------------------------------------------------------------------------
template<typename T>
bool SparseMatrix<T>::in_range(INDEX i, INDEX j)                          const
{
  return i < nRows() && j < nCols();
}
//-----------------------------------------------------------------------------
// assigns this sparsematrix from another one - full memory copy
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator=(const SparseMatrix<T> &C) 
{
  _delete();
  _copy(C);
  return *this;
}
//-----------------------------------------------------------------------------
// assigns existing sparsematrix to a value, preserving structure
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator=(const T v)
{
  this->set_all_elements_to(v);
  return *this;
}
//-----------------------------------------------------------------------------
// scales this sparse matrix by a constant
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::set_all_elements_to(const T &a) 
{
  compress(*this);
  for (INDEX i=0; i<size(); i++)  _val[i] = a;
}
//-----------------------------------------------------------------------------
// scales this sparse matrix by a constant
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator*=(const T &a) 
{
  compress(*this);
  for (INDEX i=0; i<size(); i++)  _val[i] *= a;
  return *this;
}

template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator*=(const SparseMatrix<T> &a)
{
  compress(*this);
  Matrix<T>::operator*=(a);
  return *this;
}

//-----------------------------------------------------------------------------
// Adds two sparse matrices together. 
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::operator+=(const SparseMatrix & R) 
{

  compress(R);

  int *Ria  = R.rows();
  int *Rja  = R.cols();
  T *Rval = R.ptr();

  int nRowsCRS = R.nRowsCRS();
    
  int rowR, colR;
  T valR;
  for (rowR = 0; rowR < nRowsCRS; ++rowR) { 

    for (int j = Ria[rowR]; j < Ria[rowR+1]; ++j) {
      colR = Rja[j];
      valR = Rval[j];
      
      // Because we simply want to add the value, we call add and let compress 
      // take care of the rest--we don't have to worry about extant entries.
      add(rowR, colR, valR);
    }
  }
  return *this;
}

//-----------------------------------------------------------------------------
// Return matrix transpose 
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T> SparseMatrix<T>::transpose() const 
{
  compress(*this);
  SparseMatrix<T> At(nCols(), nRows());  

  for (INDEX i=0; i<_nRowsCRS; i++)
    for (INDEX ij=_ia[i]; ij<_ia[i+1]; ij++)
      At.set(_ja[ij], i, _val[ij]);
  compress(At);
  return At;
}

//-----------------------------------------------------------------------------
// multiplies each row by the corresponding element in Vector scale
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::row_scale(const Vector<T> &v)
{
  compress(*this);
  INDEX i,ij;
  GCK(*this, v, v.size()!=nRows(), "Incompatible Vector length in row_scale.");
  for(i=0; i<_nRowsCRS; i++) 
    for(ij=_ia[i]; ij<_ia[i+1]; ij++) _val[ij] *= v[i];
  return *this;
}
//-----------------------------------------------------------------------------
// multiples each column by the corresponding element in Vector scale
//-----------------------------------------------------------------------------
template<typename T>
SparseMatrix<T>& SparseMatrix<T>::col_scale(const Vector<T> &v)
{
  compress(*this);
  INDEX i,ij;
  GCK(*this, v, v.size()!=nCols(), "Incompatible Vector length in col_scale.");
  for(i=0; i<_nRowsCRS; i++) 
    for(ij=_ia[i]; ij<_ia[i+1]; ij++) _val[ij] *= v[_ja[ij]];
  return *this;
}
//-----------------------------------------------------------------------------
// Returns a vector of the sums of each column
//-----------------------------------------------------------------------------
template<typename T>
DenseVector<T> SparseMatrix<T>::col_sum() const
{
  compress(*this);
  INDEX i,ij;
  GCHK(!nRows(), "SparseMatrix::Matrix not initialized in col_sum.")
  DenseVector<T> csum(nCols());
  for(i=0; i<_nRowsCRS; i++) 
    for(ij=_ia[i]; ij<_ia[i+1]; ij++)  csum(_ja[ij]) += _val[ij];
  return(csum);
}
//-----------------------------------------------------------------------------
// Returns a vector with the number of nonzeros in each column
//-----------------------------------------------------------------------------
template<typename T>
DenseVector<INDEX> SparseMatrix<T>::column_count() const
{
  compress(*this);
  INDEX i,j;
  DenseVector<INDEX> counts(nCols());
  
  for (i=0; i<_nRowsCRS; i++) 
    for(j=_ia[i]; j<_ia[i+1]; j++)  counts(_ja[j])++;
  return(counts);
}
//-----------------------------------------------------------------------------
// Writes a the nonzeros of a row to a vector
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::row(INDEX i, DenseVector<T>& row, DenseVector<INDEX>& indx) const 
{
  compress(*this);
  GCHK(i>=nRows(), "get_row() - invalid row number");
  if (i >= _nRowsCRS) {
    row.resize(0);
    indx.resize(0);
    return;
  }
  row.resize(RowSize(i));
  indx.resize(row.size());
  INDEX idx=0, ij;
  for(ij=_ia[i]; ij<_ia[i+1]; ij++) 
  {
    row(idx)    = _val[ij];
    indx(idx++) = _ja[ij];
  }
}
//-----------------------------------------------------------------------------
// Computes the product of N'DN
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::
weighted_least_squares(const SparseMatrix<T> &N, const DiagonalMatrix<T> &D) 
{
  compress(N);
  GCK(N,D,N.nRows()!=D.nRows(),"SparseMatrix::WeightedLeastSquares()");
  INDEX k, ki, kj;

  resize(N.nCols(), N.nCols()); // set size of this matrix
  for (k=0; k<_size; k++) _val[k] = 0.0;
  // compute R(i,j) = N(k,i) D(k,q) N(i,j) = N(k,i)*D(k,k)*N(k,j) (sum on k)
  for (k=0; k<N._nRowsCRS; k++)
    for (ki=N._ia[k]; ki<N._ia[k+1]; ki++)
      for (kj=N._ia[k]; kj<N._ia[k+1]; kj++)
        add(N._ja[ki],N._ja[kj], D[k]*N[kj]*N[ki]);
  compress();
}
//-----------------------------------------------------------------------------
// Return a diagonal matrix containing the diagonal entries of this matrix 
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T> SparseMatrix<T>::diag() const 
{
  compress(*this);
  DiagonalMatrix<T> D(nRows(), true); // initialized to zero
  INDEX i, ij;
  for (i=0; i<_nRowsCRS; i++)
  { 
    for(ij=_ia[i]; ij<_ia[i+1]; ij++) 
    {
      if (_ja[ij]>=i)  // have we reached or passed the diagonal?
      {
        if (_ja[ij]==i) D[i]=_val[ij];  // this this the diagonal?
        break;  // D[i] is already zero if there is no diagonal
      }
    }
  }
  return D;
}
//-----------------------------------------------------------------------------
// Return a diagonal matrix containing row-sum lumped entries of the matrix
//-----------------------------------------------------------------------------
template<typename T>
DiagonalMatrix<T> SparseMatrix<T>::row_sum_lump() const 
{
  compress(*this);
  DiagonalMatrix<T> D(nRows(), true); // initialized to zero
  INDEX i, ij;
  for (i=0; i<_nRowsCRS; i++)
  { 
    for(ij=_ia[i]; ij<_ia[i+1]; ij++) 
    {
      D(i,i) += _val[ij];
    }
  }
  return D;
}
//-----------------------------------------------------------------------------
// output function - builds a string with each nonzero triplet value
//-----------------------------------------------------------------------------
template<typename T>
string SparseMatrix<T>::to_string() const 
{
  compress(*this);
  string out;  
  INDEX i, ij;
  for(i=0; i<_nRowsCRS; i++)
  {
    for(ij=_ia[i]; ij<_ia[i+1]; ij++) 
    {
      if (ij) out += "\n";               // append newline if not first nonzero
      out += "(" + ATC_Utility::to_string(i) + ", ";          // append "(i,"
      out +=       ATC_Utility::to_string(_ja[ij])  + ") = "; // append  "j) = "
      out +=       ATC_Utility::to_string(_val[ij]);          // append "value"
    }
  }
  return out;  // return the completed string
}
//-----------------------------------------------------------------------------
// returns the maximum value in the row
//-----------------------------------------------------------------------------
template<typename T>
T SparseMatrix<T>::row_max(INDEX row) const 
{
  compress(*this);
  if (!RowSize(row)) return (T)0; // if there are no nonzeros in the row
  INDEX ij;
  T max = _val[_ia[row]];
  for(ij=_ia[row]+1; ij<_ia[row+1]; ij++) max = std::max(max,_val[ij]);
  return max;
}
//-----------------------------------------------------------------------------
// returns the minimum value in the row
//-----------------------------------------------------------------------------
template<typename T>
T SparseMatrix<T>::row_min(INDEX row) const 
{
  compress(*this);
  if (!RowSize(row)) return (T)0; // if there are no nonzeros in the row
  INDEX ij;
  T min = _val[_ia[row]];
  for(ij=_ia[row]+1; ij<_ia[row+1]; ij++) min = std::min(min,_val[ij]);
  return min;
}
//-----------------------------------------------------------------------------
// prints a histogram of the values of a row to the screen
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::print_row_histogram(const string &name, INDEX nbins) const 
{
  compress(*this);
  cout << "Begin histogram " << name << "\n";
  cout << "#  rows: " << _nRows << " columns: " << _nCols 
       << " size:  " << _size << "\n";
  for(INDEX i=0; i<_nRows; i++) 
  {
    print_row_histogram(i, nbins);
    cout << "\n";
  }
  cout << "End histogram " << name << "\n";
}
//-----------------------------------------------------------------------------
// prints a histogram of the values of a row to the screen
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::print_row_histogram(INDEX row, INDEX nbins) const 
{
  compress(*this);
  if (!nbins) nbins++;
  vector<INDEX> counts(nbins, 0);
  const T min = row_min(row);
  const T max = row_max(row);
  const T range = max-min;
  const double bin_size = range/double(nbins);
  if (range<=0.0) counts[nbins-1]=RowSize(row);
  else 
  {
    for(INDEX ij=_ia[row]; ij<_ia[row+1]; ij++) 
    {
      INDEX bin = INDEX((_val[ij]-min)/bin_size);
      counts[bin-(bin==nbins)]++;
    }
  }
  cout<<showbase<<scientific;
  cout<<"# Histogram: row "<<row<<" min "<<min<<" max "<<max<<" cnt " <<RowSize(row)<<"\n";
  T bin_start = min;
  for(INDEX i=0; i<nbins; i++) 
  {
    cout << "(" << bin_start << ",";
    bin_start += bin_size;
    cout << bin_start << ") " << counts[i] << "\n";
  }
}
//-----------------------------------------------------------------------------
// prints the triplets the screen
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::print_triplets() const 
{
  typename vector<TRI_COORD<T> >::const_iterator t;
  string out;
  out += "==================BEGIN TRIPLETS=======================\n";
  // Loop backwards over all new triplets.
  for (t = _tri.begin(); t!=_tri.end();  ++t) {
    out += "(" + ATC_Utility::to_string(t->i) + ", ";          // append "(i,"
    out +=       ATC_Utility::to_string(t->j)  + ") = "; // append  "j) = "
    out +=       ATC_Utility::to_string(t->v);          // append "value"
    out += "\n";
  }
  out += "===================END TRIPLETS========================\n";
  cout << out;
}
//-----------------------------------------------------------------------------
// Outputs a string to a sparse Matlab type
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::matlab(ostream &o, const string &s) const 
{
  compress(*this);
  INDEX i, ij;
  o << s <<" = sparse(" << nRows() << "," << nCols() << ");\n";
  o << showbase << scientific;
  for(i=0; i<_nRowsCRS; i++) 
    for(ij=_ia[i]; ij<_ia[i+1]; ij++)
      o<<s<<"("<<i+1<<","<<_ja[ij]+1<<")="<<_val[ij]<<";\n";
}

//-----------------------------------------------------------------------------
// Writes the matrix to a binary file (after a compress).
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::binary_write(std::fstream& f) const
{
  compress(*this);  
  f.write((char*)&_size,     sizeof(INDEX));  // writes number of nonzeros  
  f.write((char*)&_nRowsCRS, sizeof(INDEX));  // writes number of rows in crs
  f.write((char*)&_nRows,    sizeof(INDEX));  // write matrix rows
  f.write((char*)&_nCols,    sizeof(INDEX));  // write number of columns
  if (!_size) return;
  f.write((char*)_val, sizeof(T)    *_size); 
  f.write((char*)_ja,  sizeof(INDEX)*_size); 
  f.write((char*)_ia,  sizeof(INDEX)*(_nRowsCRS+1));  
}
//-----------------------------------------------------------------------------
// Reads a SparseMatrix from a binary file.  (wipes out any original data)
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::binary_read(std::fstream& f)
{
  _delete();
  f.read((char*)&_size,     sizeof(INDEX));
  f.read((char*)&_nRowsCRS, sizeof(INDEX));
  f.read((char*)&_nRows,    sizeof(INDEX));
  f.read((char*)&_nCols,    sizeof(INDEX));
  if (!_size) return;
  _create(_size,_nRowsCRS);
  f.read((char*)_val,      sizeof(T)*_size);
  f.read((char*)_ja,       sizeof(INDEX)*_size);
  f.read((char*)_ia,       sizeof(INDEX)*(_nRowsCRS+1)); 
}

//-----------------------------------------------------------------------------
// Writes the sparse matrix to a file in a binary format
//-----------------------------------------------------------------------------
template<typename T>
void SparseMatrix<T>::write_restart(FILE *f) const
{
  compress(*this);
  fwrite(&_size,     sizeof(INDEX), 1 ,f);  // write number of nonzeros
  fwrite(&_nRowsCRS, sizeof(INDEX), 1 ,f);  // write number of rows
  fwrite(&_nRows,    sizeof(INDEX), 1 ,f);  // write number of columns
  fwrite(&_nCols,    sizeof(INDEX), 1 ,f);  // write number of columns
  if (!_size) return;
  fwrite(_val, sizeof(T), _size ,f);
  fwrite(_ja,  sizeof(T), _size ,f);
  fwrite(_ia,  sizeof(INDEX), _nRowsCRS+1 ,f);
}

} // end namespace
#endif

