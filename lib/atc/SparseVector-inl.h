// SparseVector-inl.h: provides templated functions for SparseVector in
// separate header


namespace ATC_matrix {

template<class T>
SparseVector<T> sparse_rand(INDEX n, INDEX fill, int seed=1234)
{
  srand(seed);
  const double rmax_inv = 1.0/double(RAND_MAX);
  SparseVector<T> r(n);
  while (r.size()<fill) 
    r(std::rand()%r.nRows())=double(std::rand()*rmax_inv);
  return r;
}

// Multiplies a Matrix by a SparseVector (M*v) and returns a DenseVector.
template<class T>
DenseVector<T> operator*(const Matrix<T> &M, const SparseVector<T> &v)
{
  DenseVector<T> y(M.nRows());
  STORE::const_iterator it=v.data_.begin();
  for (; it!=v.data_.end(); it++) {
    const INDEX j  = it->first;
    const T& vj = it->second;
     for (INDEX i=0; i<M.nRows(); i++) y(i)+=M(i,j)*vj;
  }
  return y;
}

// Multiplies a SparseVector by a Matrix (M'*v) and returns a DenseVector.
template<class T>
DenseVector<T> operator*(const SparseVector<T> &v, const Matrix<T> &M)
{
  DenseVector<T> y(M.nCols());
  STORE::const_iterator it=v.data_.begin();
  for (; it!=v.data_.end(); it++) {
    const INDEX i  = it->first;
    const T& vi = it->second;
     for (INDEX j=0; j<M.nCols(); j++) y(j)+=vi*M(i,j);
  }
  return y;
}

// Computes the dot product between two SparseVectors of equal length.
template<class T>
T dot(const SparseVector<T> &a, const SparseVector<T> &b)
{
  T v = 0.0;
  for (STORE::const_iterator ai=a.data_.begin(); ai!=a.data_.end(); ai++) {
    STORE::const_iterator bi=b.data_.find(ai->first);
    if (bi == b.data_.end()) continue;
    v += ai->second * bi->second;
  } 
  return v;
}
// Computes the product of a SparseMatrix tranpose with a SparseVector (M'*v).
template<class T>
SparseVector<T> operator*(const SparseMatrix<T> &M, const SparseVector<T> &v)
{
  SparseVector<T> y(M.nRows());
  for (INDEX i=0; i<M.nRows(); i++) {
    double yi=0.0;
    for (INDEX ij=M._ia[i]; ij<M._ia[i+1]; ij++) {
      const INDEX j = M._ja[ij];
      STORE::const_iterator it = v._data.find(j);
      if (it == v._data.end()) continue;
      yi += M._v[ij] * it->second;
    }
    if (yi!=0.0) y(i)+=yi;
  }
  return y; 
}

// computes the product of a SparseMatrix tranpose with a SparseVector (M'*v).
template<class T>
SparseVector<T> operator*(const SparseVector<T> &v, const SparseMatrix<T> &M)
{
  SparseVector<T> y(M.nCols());
  for (INDEX i=0; i<M.nRows(); i++) {
    STORE::const_iterator it = v._data.find(i);
    if (it == v._data.end()) continue;
    for (INDEX ij=M._ia[i]; ij<M._ia[i+1]; ij++) 
      y(M._ja[ij]) += it->second * M._v[ij];
  }
  return y;
}

// General constructor - sets the vector length.
template<class T>
SparseVector<T>::SparseVector(INDEX n) : length_(n){}

// Outputs the vector to string
template<class T>
std::string SparseVector<T>::to_string() const
{
  if (size() > nRows()/2) return Vector<T>::to_string();
  STORE::const_iterator it=data_.begin();
  std::string str;
  using ATC_Utility::to_string;
  for (; it!=data_.end(); it++)
    str += to_string(it->first)+": "+to_string(it->second)+'\n';
  return str;
}

// Indexes the ith component of the vector or returns zero if not found. 
template<class T>
T SparseVector<T>::operator()(INDEX i, INDEX j) const
{
  STORE::const_iterator it = data_.find(i);
  if (it == data_.end()) return 0.0;
  return it->second;
}

// Indexes the ith component of the vector or returns zero if not found. 
template<class T>
T& SparseVector<T>::operator()(INDEX i, INDEX j)
{
  return data_[i]; 
}

// Indexes the ith component of the vector or returns zero if not found. 
template<class T> T SparseVector<T>::operator[](INDEX i) const
{
  return (*this)(i);
}

// Indexes the ith component of the vector or returns zero if not found. 
template<class T> T& SparseVector<T>::operator[](INDEX i)
{
  return (*this)[i];
}

// Returns a pair (index, value) for a nonzero in the vector.
template<class T> 
std::pair<INDEX, T> SparseVector<T>::pair(INDEX i) const
{
  STORE::const_iterator it=data_.begin() + i;
  return std::pair<INDEX, T>(it->first, it->second);
}

//* Adds SparseVector x, scaled by s to this one.  Can be different sparcity.
template<class T> 
void SparseVector<T>::add_scaled(SparseVector<T>& x, const T& s)
{
  STORE::const_iterator it;
  for (it=x.data_.begin(); it!=x.data_.end(); it++)  {
   data_[it->first] += it->second * s;
  }
}

// Returns the number of nonzeros in the sparse vector.
template<class T> inline INDEX SparseVector<T>::size() const
{ return data_.size(); }

// Returns the number of nonzeros in the sparse vector.
template<class T> inline INDEX SparseVector<T>::nRows() const
{ return length_; }

// Copy constructor for sparse vector.
template<typename T>
SparseVector<T>::SparseVector(const SparseVector<T> &c)
{
  length_ = c.length_;
  data_ = c.data_;
}

// operator equal for sparse vector.
template<class T>
SparseVector<T>& SparseVector<T>::operator=(const SparseVector<T> &c)
{
  length_ = c.length_;
  data_ = c.data_;
  return *this;
}

// Changes the size of the SparseVector
template<class T>
void SparseVector<T>::resize(INDEX nRows, INDEX nCols, bool copy)
{
  length_ = nRows;
  STORE::iterator it;
  for (it=data_.begin(); it!=data_.end(); it++) {
    if (it->second >= length_) data_.erase(it);
    else if (!copy) it->second=T(0);
  }
}
// same as resize, but zero rather than copy
template<class T>
void SparseVector<T>::reset(INDEX nRows, INDEX nCols, bool zero)
{
  resize(nRows, nCols, !zero);
}
// sets all elements to zero but preserves sparcity
template<class T>
void SparseVector<T>::zero()
{
  STORE::iterator it;
  for (it=data_.begin(); it!=data_.end(); it++) it->second=T(0);
}


template<class T>
void SparseVector<T>::copy(const T* ptr, INDEX nRows, INDEX nCols)
{

}


template<class T>
void SparseVector<T>::write_restart(FILE *F) const
{

}

// writes a stream to a matlab script to recreate this variable 
template<class T>
void SparseVector<T>::matlab(std::ostream &o, const std::string &s) const
{
  o << s << "=sparse(" << nRows() << ",1);\n";
  o << std::showbase << std::scientific;
  STORE::const_iterator it;
  for (it=data_.begin(); it!=data_.end(); it++)
    o << s << "(" << it->first+1 << ") = " << it->second << ";\n";
}

} // end namespace
