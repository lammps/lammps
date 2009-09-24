// SparseVector-inl.h: provides templated functions for SparseVector in
// separate header

template<class T>
SparseVector<T> sparse_rand(unsigned n, unsigned fill, int seed=1234)
{
  srand(seed);
  const double rmax_inv = 1.0/double(RAND_MAX);
  SparseVector<T> r(n);
  while (r.size()<fill) 
    r(rand()%r.nRows())=double(::rand()*rmax_inv);
  return r;
}

// Multiplies a Matrix by a SparseVector (M*v) and returns a DenseVector.
template<class T>
DenseVector<T> operator*(const Matrix<T> &M, const SparseVector<T> &v)
{
  DenseVector<T> y(M.nRows());
  STORE::const_iterator it=v.data_.begin();
  for (; it!=v.data_.end(); it++) {
    const unsigned j  = it->first;
    const T& vj = it->second;
     for (unsigned i=0; i<M.nRows(); i++) y(i)+=M(i,j)*vj;
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
    const unsigned i  = it->first;
    const T& vi = it->second;
     for (unsigned j=0; j<M.nCols(); j++) y(j)+=vi*M(i,j);
  }
  return y;
}

// Computes the dot product between two SparseVectors of equal length.
template<class T>
T dot(const SparseVector<T> &a, const SparseVector<T> &b)
{
  double v = 0.0;
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
  for (unsigned i=0; i<M.nRows(); i++) {
    double yi=0.0;
    for (unsigned ij=M._ia[i]; ij<M._ia[i+1]; ij++) {
      const unsigned j = M._ja[ij];
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
  for (unsigned i=0; i<M.nRows(); i++) {
    STORE::const_iterator it = v._data.find(i);
    if (it == v._data.end()) continue;
    for (unsigned ij=M._ia[i]; ij<M._ia[i+1]; ij++) 
      y(M._ja[ij]) += it->second * M._v[ij];
  }
  return y;
}

// General constructor - sets the vector length.
template<class T>
SparseVector<T>::SparseVector(unsigned n) : length_(n){}

// Outputs the vector to string
template<class T>
std::string SparseVector<T>::tostring() const
{
  if (size() > nRows()/2) return Vector<T>::tostring();
  STORE::const_iterator it=data_.begin();
  std::string str;
  using ATC_STRING::tostring;
  for (; it!=data_.end(); it++)
    str += tostring(it->first)+": "+tostring(it->second)+'\n';
  return str;
}

// Indexes the ith component of the vector or returns zero if not found. 
template<class T>
T SparseVector<T>::operator()(unsigned i, unsigned j) const
{
  STORE::const_iterator it = data_.find(i);
  if (it == data_.end()) return 0.0;
  return it->second;
}

// Indexes the ith component of the vector or returns zero if not found. 
template<class T>
T& SparseVector<T>::operator()(unsigned i, unsigned j)
{
  return data_[i]; 
}

// Indexes the ith component of the vector or returns zero if not found. 
template<class T> T SparseVector<T>::operator[](unsigned i) const
{
  return (*this)(i);
}

// Indexes the ith component of the vector or returns zero if not found. 
template<class T> T& SparseVector<T>::operator[](unsigned i)
{
  return (*this)[i];
}

// Returns the number of nonzeros in the sparse vector.
template<class T> inline unsigned SparseVector<T>::size() const
{ return data_.size(); }

// Returns the number of nonzeros in the sparse vector.
template<class T> inline unsigned SparseVector<T>::nRows() const
{ return length_; }

// Changes the size of the SparseVector
template<class T>
void SparseVector<T>::resize(unsigned nRows, unsigned nCols, bool copy)
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
void SparseVector<T>::reset(unsigned nRows, unsigned nCols, bool zero)
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

// TODO
template<class T>
void SparseVector<T>::copy(const T* ptr, unsigned nRows, unsigned nCols)
{

}

// TODO
template<class T>
void SparseVector<T>::write_restart(FILE *F) const
{

}

// writes a stream to a matlab script to recreate this variable 
template<class T>
void SparseVector<T>::matlab(ostream &o, const string &s) const
{
  o << s << "=sparse(" << nRows() << ",1);\n";
  o << showbase << scientific;
  STORE::const_iterator it;
  for (it=data_.begin(); it!=data_.end(); it++)
    o << s << "(" << it->first+1 << ") = " << it->second << ";\n";
}
