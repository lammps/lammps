#ifndef MATRIX_H
#define MATRIX_H
#include "MatrixDef.h"

/*******************************************************************************
* abstract class Matrix - base class for linear algebra subsystem
*******************************************************************************/
template<typename T>
class Matrix
{
protected:
  Matrix(const Matrix &c);
public:
  Matrix() {}
  virtual ~Matrix() {}

  //* stream output functions
  void print(ostream &o) const { o << tostring(); }
  void print(ostream &o, const string &name) const;
  friend ostream& operator<<(ostream &o, const Matrix<T> &m){m.print(o); return o;}
  void print() const;
  virtual void print(const string &name) const;
  virtual string tostring() const;

  // element by element operations
  DenseMatrix<T> operator/(const Matrix<T>& B) const;
  DenseMatrix<T> pow(int n) const;
  DenseMatrix<T> pow(double n) const;

  // functions that return a copy 
  DenseMatrix<T> transpose() const;

  // matrix to scalar functions
  T sum()    const;
  T stdev()  const;
  T max()    const;
  T min()    const;
  T maxabs() const;
  T minabs() const;
  T norm()   const;
  T mean()   const;
  T dot(const Matrix<T> &r) const;
  T trace() const;

  // row and column operations
  T row_sum  (INDEX i=0) const { return    row(*this,i).sum();   }
  T row_mean (INDEX i=0) const { return    row(*this,i).mean();  }
  T row_norm (INDEX i=0) const { return    row(*this,i).norm();  }
  T row_stdev(INDEX i=0) const { return    row(*this,i).stdev(); }
  T col_sum  (INDEX i=0) const { return column(*this,i).sum();   }
  T col_mean (INDEX i=0) const { return column(*this,i).mean();  }
  T col_norm (INDEX i=0) const { return column(*this,i).norm();  }
  T col_stdev(INDEX i=0) const { return column(*this,i).stdev(); }
 
  // pure virtual functions (required to implement these) ---------------------
  //* reference index operator
  virtual T& operator()(INDEX i, INDEX j)=0;      
  //* value index operator
  virtual T  operator()(INDEX i, INDEX j)const=0; 
  //* value flat index operator
  virtual T& operator [](INDEX i)=0;              
  //* reference flat index operator
  virtual T operator [](INDEX i) const=0;         
  //* returns the # of rows
  virtual INDEX nRows() const=0;                  
  //* returns the # of columns
  virtual INDEX nCols() const=0;                  
  //* returns a pointer to the data (dangerous)
  virtual T * get_ptr() const=0;                  
  //* resizes the matrix, copy what fits default to OFF
  virtual void resize(INDEX nRows, INDEX nCols=1, bool copy=false)=0;
  //* resizes the matrix, zero it out default to ON
  virtual void  reset(INDEX nRows, INDEX nCols=1, bool zero=true)=0;
  //* resizes and copies data
  virtual void copy(const T * ptr, INDEX nRows, INDEX nCols=1)=0;
  //* create restart file
  virtual void write_restart(FILE *f) const=0;    
  //* writes a matlab command to recreate this in a variable named s
  virtual void matlab(ostream &o, const string &s="M") const;

  // output to matlab, with variable name s
  void matlab(const string &s="M") const;

  Matrix<T>& operator+=(const Matrix &r);
  Matrix<T>& operator-=(const Matrix &r);
  //* NOTE these two functions are scaling operations not A.B
  Matrix<T>& operator*=(const Matrix<T>& R);
  Matrix<T>& operator/=(const Matrix<T>& R);
  Matrix<T>& operator+=(const T v);
  Matrix<T>& operator-=(const T v);
  Matrix<T>& operator*=(const T v);
  Matrix<T>& operator/=(T v);

  Matrix<T>& operator=(const T &v);
  Matrix<T>& operator=(const Matrix<T> &c); 
  virtual void set_all_elements_to(const T &v);
  //* adds a matrix scaled by factor s to this one.
  void add_scaled(const Matrix<T> &A, const T& s);

  //* sets all elements to zero
  Matrix& zero();
  //* returns the total number of elements
  virtual INDEX size() const;
  //* returns true if (i,j) is within the range of the matrix
  bool in_range(INDEX i, INDEX j) const;
  //* returns true if the matrix size is rs x cs
  bool is_size(INDEX rs, INDEX cs) const;
  //* returns true if the matrix is square and not empty
  bool is_square() const;
  //* returns true if Matrix, m, is the same size as this
  bool same_size(const Matrix &m) const;
  //* returns true if Matrix a and Matrix b are the same size
  static bool same_size(const Matrix<T> &a, const Matrix<T> &b);
  //* checks if memory is contiguous, only can be false for clone vector
  virtual bool memory_contiguous() const { return true; }

protected:
  virtual void _set_equal(const Matrix<T> &r);
};

//* Matrix operations 
//@{
//* Sets C as b*C + a*A[tranpose?]*B[transpose?] 
template<typename T>
void MultAB(const Matrix<T> &A, const Matrix<T> &B, DenseMatrix<T> &C, 
            bool At=0, bool Bt=0, T a=1, T b=0);
//* performs a matrix-vector multiply
template<typename T>
void MultMv(const Matrix<T> &A, const Vector<T> &v, DenseVector<T> &c,
            const bool At, T a, T b);
// returns the inverse of a double precision matrix
DenseMatrix<double> inv(const Matrix<double>& A);
//* returns the trace of a matrix
template<typename T> 
T trace(const Matrix<T>& A) { return A.trace(); }
//* computes the determinant of a square matrix
double det(const Matrix<double>& A);
//* Returns the maximum eigenvalue of a matrix.
double max_eigenvalue(const Matrix<double>& A);
//@}

//-----------------------------------------------------------------------------
// computes the sum of the difference squared of each element. 
//-----------------------------------------------------------------------------
template<typename T>
double sum_difference_squared(const Matrix<T>& A, const Matrix<T> &B)
{
  SSCK(A, B, "sum_difference_squared");
  double v=0.0;
  for (unsigned i=0; i<A.size(); i++) {
    double d = A[i]-B[i];
    v += d*d;
  }
  return v;
}

//-----------------------------------------------------------------------------
//* Operator for Matrix-matrix product
//-----------------------------------------------------------------------------
template<typename T>
DenseMatrix<T> operator*(const Matrix<T> &A, const Matrix<T> &B)
{
  DenseMatrix<T> C(0,0,false);
  MultAB(A,B,C);
  return C;
}
//-----------------------------------------------------------------------------
//* Multiply a Matrix by a scalar 
//-----------------------------------------------------------------------------
template<typename T>
DenseMatrix<T> operator*(const Matrix<T> &M, const T s)
{
  DenseMatrix<T> R(M);
  return R*=s; 
}
//-----------------------------------------------------------------------------
//* Multiply a Matrix by a scalar 
template<typename T>
DenseMatrix<T> operator*(const T s, const Matrix<T> &M)
{
  DenseMatrix<T> R(M);
  return R*=s; 
}
//-----------------------------------------------------------------------------
//* inverse scaling operator - must always create memory
template<typename T>
DenseMatrix<T> operator/(const Matrix<T> &M, const T s)
{
  DenseMatrix<T> R(M);
  return R*=(1.0/s); // for integer types this may be worthless
}
//-----------------------------------------------------------------------------
//* Operator for Matrix-matrix sum
template<typename T>
DenseMatrix<T> operator+(const Matrix<T> &A, const Matrix<T> &B)
{
  DenseMatrix<T> C(A);
  return C+=B;
}
//-----------------------------------------------------------------------------
//* Operator for Matrix-matrix subtraction
template<typename T>
DenseMatrix<T> operator-(const Matrix<T> &A, const Matrix<T> &B)
{
  DenseMatrix<T> C(A);
  return C-=B;
}
/******************************************************************************
* Template definitions for class Matrix
******************************************************************************/

//-----------------------------------------------------------------------------
//* performs a matrix-matrix multiply with general type implementation
template<typename T>
void MultAB(const Matrix<T> &A, const Matrix<T> &B, DenseMatrix<T> &C, 
            const bool At, const bool Bt, T a, T b)
{
  const INDEX sA[2] = {A.nRows(), A.nCols()};  // m is sA[At] k is sA[!At]
  const INDEX sB[2] = {B.nRows(), B.nCols()};  // k is sB[Bt] n is sB[!Bt]
  const INDEX M=sA[At], K=sB[Bt], N=sB[!Bt];
  GCK(A, B, sA[!At]!=K, "MultAB<T> shared index not equal size");
  if (!C.is_size(M,N))
  {
    C.resize(M,N);          // set size of C
    C.zero();
  }
  else C *= b;
  for (INDEX p=0; p<M; p++)
    for (INDEX q=0; q<N; q++)
      for (INDEX r=0; r<K; r++)
         C(p,q) += A(p*!At+r*At, p*At+r*!At) * B(r*!Bt+q*Bt, r*Bt+q*!Bt);
}

//-----------------------------------------------------------------------------
//* output operator
template<typename T>
string Matrix<T>::tostring() const
{
  string s;
  for (INDEX i=0; i<nRows(); i++)
  {
    if (i) s += '\n';
    for (INDEX j=0; j<nCols(); j++)
    {
      if (j) s+= '\t';
      s += ATC_STRING::tostring(MIDX(i,j),5)+"   ";
    }
  }
  return s;
}
//-----------------------------------------------------------------------------
//* output operator that wraps the matrix in a nice labeled box
template<typename T>
void Matrix<T>::print(ostream &o, const string &name) const
{
  o << "------- Begin "<<name<<" -----------------\n";
  this->print(o);
  o << "\n------- End "<<name<<" -------------------\n";
}
//-----------------------------------------------------------------------------
//* print operator, use cout by default
template<typename T>
void Matrix<T>::print() const 
{
  print(cout); 
}
//-----------------------------------------------------------------------------
//* named print operator, use cout by default
template<typename T>
void Matrix<T>::print(const string &name) const 
{
  print(cout, name); 
}
//-----------------------------------------------------------------------------
//* element by element division
template<typename T>
DenseMatrix<T> Matrix<T>::operator/ (const Matrix<T>& B) const
{
  SSCK(*this, B, "Matrix<T>::Operator/");
  DenseMatrix<T> R(*this);
  R /= B;
  return R;
}
//-----------------------------------------------------------------------------
//* element-wise raise to a power
template<typename T>
DenseMatrix<T> Matrix<T>::pow(int n) const
{
  DenseMatrix<T> R(*this);
  FORi
  {
    double val = R[i];
    for (int k=1; k<n; k++) val *= R[i];
    for (int k=n; k<1; k++) val /= R[i];
    R[i] = val;
  }
  return R;
}
//-----------------------------------------------------------------------------
//* element-wise raise to a power
template<typename T>
DenseMatrix<T> Matrix<T>::pow(double n) const
{
  DenseMatrix<T> R(*this);
  FORi
  {
    double val = R[i];
    R[i] = pow(val,n);
  }
  return R;
}
//-----------------------------------------------------------------------------
//* returns the transpose of this matrix (makes a copy)
template <typename T>
DenseMatrix<T> Matrix<T>::transpose()                                     const
{
  DenseMatrix<T> t(this->nCols(), this->nRows());
  FORij t(j,i) = (*this)(i,j);
  return t;
}
//-----------------------------------------------------------------------------
//* returns the transpose of a matrix (makes a copy)
template <typename T>
DenseMatrix<T> transpose(const Matrix<T> &A)
{
  return A.transpose();
}
//-----------------------------------------------------------------------------
//* Returns the sum of all matrix elements 
template<typename T>
T Matrix<T>::sum() const
{
  if (!size())  return T(0);
  T v = VIDX(0);
  for (INDEX i=1; i<this->size(); i++) v += VIDX(i);
  return v; 
}
//-----------------------------------------------------------------------------
//* Returns the standard deviation of the matrix
template<typename T>
T Matrix<T>::stdev() const
{
  GCHK(this->size()<2, "Matrix::stdev() size must be > 1");
  T mean = this->mean();
  T diff = VIDX(0)-mean;
  T stdev = diff*diff;
  for (INDEX i=1; i<this->size(); i++) 
  {
    diff = VIDX(i)-mean;
    stdev += diff*diff;
  }
  return sqrt(stdev/T(this->size()-1)); 
}
//-----------------------------------------------------------------------------
//* Returns the maximum of the matrix
template<typename T>
T Matrix<T>::max() const
{
  GCHK(!this->size(), "Matrix::max() size must be > 0");
  T v = VIDX(0);
  for (INDEX i=1; i<this->size(); i++)   v = std::max(v, VIDX(i));
  return v; 
}
//-----------------------------------------------------------------------------
//* Returns the minimum of the matrix
template<typename T>
T Matrix<T>::min() const
{
  GCHK(!this->size(), "Matrix::min() size must be > 0");
  T v = VIDX(0);
  for (INDEX i=1; i<this->size(); i++)   v = std::min(v, VIDX(i));
  return v; 
}
//-----------------------------------------------------------------------------
//* Returns the maximum absolute value of the matrix
template<typename T>
T Matrix<T>::maxabs() const
{
  GCHK(!this->size(), "Matrix::maxabs() size must be > 0");
  T v = VIDX(0);
  for (INDEX i=1; i<this->size(); i++)   v = Utility::MaxAbs(v, VIDX(i));
  return v; 
}
//-----------------------------------------------------------------------------
//* Returns the minimum absoute value of the matrix
template<typename T>
T Matrix<T>::minabs() const
{
  GCHK(!this->size(), "Matrix::minabs() size must be > 0");
  T v = VIDX(0);
  for (INDEX i=1; i<this->size(); i++)   v = Utility::MinAbs(v, VIDX(i));
  return v; 
}
//-----------------------------------------------------------------------------
//* returns the L2 norm of the matrix
template<typename T>
T Matrix<T>::norm()                                                       const    
{
  GCHK(!this->size(), "Matrix::norm() size must be > 0");
  return sqrt(dot(*this)); 
}
//-----------------------------------------------------------------------------
//* returns the average of the matrix
template<typename T>
T Matrix<T>::mean()                                                       const    
{
  GCHK(!this->size(), "Matrix::mean() size must be > 0");
  return sum()/T(this->size());
}
//-----------------------------------------------------------------------------
//* Returns the dot product of two vectors
template<typename T>
T Matrix<T>::dot(const Matrix<T>& r)  const
{
  SSCK(*this, r, "Matrix<T>::dot");
  if (!this->size()) return T(0);
  T v = r[0]*VIDX(0);
  for (INDEX i=1; i<this->size(); i++)  v += r[i]*VIDX(i);
  return v;
}
//-----------------------------------------------------------------------------
// returns the sum of the matrix diagonal
//-----------------------------------------------------------------------------
template<typename T>
T Matrix<T>::trace() const
{
  const INDEX N = std::min(nRows(),nCols());
  if (!N) return T(0);
  T r = MIDX(0,0);
  for (INDEX i=0; i<N; i++)
    r += MIDX(i,i);
  return r;
}
//-----------------------------------------------------------------------------
//* Adds a matrix to this one
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix &r) 
{
  SSCK(*this, r, "operator+= or operator +");
  FORi VIDX(i)+=r[i];
  return *this;
}
//-----------------------------------------------------------------------------
// subtracts a matrix from this one
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix &r) 
{
  SSCK(*this, r, "operator-= or operator -");
  FORi VIDX(i)-=r[i];
  return *this;
}
//-----------------------------------------------------------------------------
// multiplies each element in this by the corresponding element in R
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& R)
{
  SSCK(*this, R, "operator*=");
  FORi 
  {
    VIDX(i) *= R[i]; 
  }
  return *this;
}
//-----------------------------------------------------------------------------
// divides each element in this by the corresponding element in R
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator/=(const Matrix<T>& R)
{
  SSCK(*this, R, "operator/= or operator/");
  FORi 
  {
    GCHK(fabs(R[i])>0,"Operator/: division by zero"); 
    VIDX(i) /= R[i]; 
  }
  return *this;
}
//-----------------------------------------------------------------------------
// scales this matrix by a constant
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T v)       
{
  FORi VIDX(i)*=v;      
  return *this;
}
//-----------------------------------------------------------------------------
// adds a constant to this matrix 
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const T v)       
{
  FORi VIDX(i)+=v;      
  return *this;
}
//-----------------------------------------------------------------------------
// substracts a constant to this matrix 
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const T v)       
{
  FORi VIDX(i)-=v;      
  return *this;
}
//-----------------------------------------------------------------------------
//* scales this matrix by the inverse of a constant
template<typename T>
Matrix<T>& Matrix<T>::operator/=(T v)    
{
  return (*this)*=(1.0/v);
}

//----------------------------------------------------------------------------
// Assigns one matrix to another
//----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &r)
{
  this->_set_equal(r);
  return *this;
}
//----------------------------------------------------------------------------
// general matrix assignment (for densely packed matrices)
//----------------------------------------------------------------------------
template<typename T>
void Matrix<T>::_set_equal(const Matrix<T> &r)
{
  this->resize(r.nRows(), r.nCols());
  const Matrix<T> *pr = &r;
  if (const SparseMatrix<T> *ps = sparse_cast(pr)) 
    copy_sparse_to_matrix(ps, *this);
  
  else if (diag_cast(pr))  // r is Diagonal?
  {
    this->zero();
    for (INDEX i=0; i<r.size(); i++) MIDX(i,i) = r[i];
  }
  else memcpy(this->get_ptr(), r.get_ptr(), r.size()*sizeof(T));
}
//-----------------------------------------------------------------------------
//* sets all elements to a constant
template<typename T>
inline Matrix<T>& Matrix<T>::operator=(const T &v)    
{
  set_all_elements_to(v); 
  return *this;
}
//-----------------------------------------------------------------------------
//* sets all elements to a constant
template<typename T>
void Matrix<T>::set_all_elements_to(const T &v)    
{
  FORi VIDX(i) = v;
}
//----------------------------------------------------------------------------
// adds a matrix scaled by factor s to this one.
//----------------------------------------------------------------------------
template <typename T>
void Matrix<T>::add_scaled(const Matrix<T> &A, const T& s)
{
    SSCK(A, *this, "Matrix::add_scaled");
    FORi VIDX(i) += A[i]*s;
}
//-----------------------------------------------------------------------------
//* writes a matlab command to the console
template<typename T>
void Matrix<T>::matlab(const string &s) const 
{
  this->matlab(cout, s); 
}
//-----------------------------------------------------------------------------
//* Writes a matlab script defining the vector to the stream
template<typename T>
void Matrix<T>::matlab(ostream &o, const string &s) const
{
  o << s <<"=zeros(" << nRows() << ","<<nCols()<<");\n";
  FORij o << s << "("<<i+1<<","<<j+1<<")=" << MIDX(i,j) << ";\n";
}
//-----------------------------------------------------------------------------
//* sets all matrix elements to zero
template<typename T>
inline Matrix<T>& Matrix<T>::zero() 
{ 
  set_all_elements_to(T(0));
  return *this;
}
//-----------------------------------------------------------------------------
//* returns the total number of elements
template<typename T>
inline INDEX Matrix<T>::size() const 
{ 
  return nRows()*nCols(); 
}
//-----------------------------------------------------------------------------
//* returns true if (i,j) is within the range of the matrix
template<typename T>
inline bool Matrix<T>::in_range(INDEX i, INDEX j) const
{ 
  return i<nRows() && j<nCols(); 
}
//-----------------------------------------------------------------------------
//* returns true if the matrix size is rs x cs
template<typename T>
inline bool Matrix<T>::is_size(INDEX rs, INDEX cs) const
{ 
  return nRows()==rs && nCols()==cs; 
}
//-----------------------------------------------------------------------------
//* returns true if the matrix is square and not empty
template<typename T>
inline bool Matrix<T>::is_square() const
{ 
  return nRows()==nCols() && nRows(); 
}
//-----------------------------------------------------------------------------
//* returns true if Matrix, m, is the same size as this
template<typename T>
inline bool Matrix<T>::same_size(const Matrix<T> &m) const
{
  return is_size(m.nRows(), m.nCols());
}
//-----------------------------------------------------------------------------
//* returns true if Matrix a and Matrix b are the same size
template<typename T>
inline bool Matrix<T>::same_size(const Matrix<T> &a, const Matrix<T> &b)
{
  return a.same_size(b); 
}
//-----------------------------------------------------------------------------
//* Displays indexing error message and quits
template<typename T>
void ierror(const Matrix<T> &a, const char *FILE, int LINE, INDEX i, INDEX j)
{
  cout << "Error: Matrix indexing failure ";
  cout << "in file: " << FILE << ", line: "<< LINE <<"\n";
  cout << "Tried accessing index (" << i << ", " << j <<")\n";
  cout << "Matrix size was "<< a.nRows() << "x" << a.nCols() << "\n";
  exit(EXIT_FAILURE);
}
//-----------------------------------------------------------------------------
//* Displays custom message and indexing error and quits
template<typename T>
void ierror(const Matrix<T> &a, INDEX i, INDEX j, const string m)
{
  cout << m << "\n";
  cout << "Tried accessing index (" << i << ", " << j <<")\n";
  cout << "Matrix size was "<< a.nRows() << "x" << a.nCols() << "\n";
  exit(EXIT_FAILURE);
}
//-----------------------------------------------------------------------------
//* Displays matrix compatibility error message
template<typename T>
void merror(const Matrix<T> &a, const Matrix<T> &b, const string m)
{
  cout << "Error: " << m << "\n";
  cout << "Matrix sizes were " << a.nRows() << "x" << a.nCols();
  if (&a != &b) cout << ", and "<< b.nRows() << "x" << b.nCols();
  cout << "\n";
  exit(EXIT_FAILURE);
}

#endif
