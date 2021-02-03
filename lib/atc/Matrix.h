#ifndef MATRIX_H
#define MATRIX_H
#include "MatrixDef.h"

namespace ATC_matrix {
static const int myPrecision = 15;

  /**
   *  @class  Matrix
   *  @brief  Base class for linear algebra subsystem 
   */

template<typename T>
class Matrix
{

protected:
  Matrix(const Matrix &c);
public:
  Matrix() {}
  virtual ~Matrix() {}

  //* stream output functions
  void print(std::ostream &o, int p=myPrecision) const { o << this->to_string(p); }
  void print(std::ostream &o, const std::string &name, int p=myPrecision) const;
  friend std::ostream& operator<<(std::ostream &o, const Matrix<T> &m){m.print(o); return o;}
  void print() const;
  virtual void print(const std::string &name, int p = myPrecision) const;
  virtual std::string to_string(int p) const;
  virtual std::string to_string() const { return to_string(myPrecision); }

  // element by element operations
  DenseMatrix<T> operator/(const Matrix<T>& B) const;
  DenseMatrix<T> pow(int n) const;
  DenseMatrix<T> pow(double n) const;

  // functions that return a copy 
  DenseMatrix<T> transpose() const;
  void row_partition(const std::set<int> & rowsIn, std::set<int> & rows, std::set<int> & colsC,
    DenseMatrix<T> & A1, DenseMatrix<T> & A2, bool complement=true) const;
  std::set<int> row_partition(const std::set<int> & rows, 
    DenseMatrix<T> & A1, DenseMatrix<T> & A2) const;
  void map(const std::set<int>& rows, const std::set<int>& cols, DenseMatrix<T> & A) const;
  void insert(const std::set<int>& rows, const std::set<int>& cols, const DenseMatrix<T> & A);
  void assemble(const std::set<int>& rows, const std::set<int>& cols, const DenseMatrix<T> & A);

  // matrix to scalar functions
  T sum()    const;
  T stdev()  const;
  T max()    const;
  T min()    const;
  T maxabs() const;
  T minabs() const;
  T norm()   const;
  T norm_sq()   const;
  T mean()   const;
  T dot(const Matrix<T> &r) const;
  T trace() const;

  // row and column operations
  T row_sum  (INDEX i=0) const { return    row(*this,i).sum();   }
  T row_mean (INDEX i=0) const { return    row(*this,i).mean();  }
  T row_norm (INDEX i=0) const { return    row(*this,i).norm();  }
  T row_min  (INDEX i=0) const { return    row(*this,i).min();  }
  T row_max  (INDEX i=0) const { return    row(*this,i).max();  }
  T row_stdev(INDEX i=0) const { return    row(*this,i).stdev(); }
  T col_sum  (INDEX i=0) const { return column(*this,i).sum();   }
  T col_mean (INDEX i=0) const { return column(*this,i).mean();  }
  T col_norm (INDEX i=0) const { return column(*this,i).norm();  }
  T col_min  (INDEX i=0) const { return column(*this,i).min();  }
  T col_max  (INDEX i=0) const { return column(*this,i).max();  }
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
  virtual T * ptr() const=0;                  
  //* resizes the matrix, copy what fits default to OFF
  virtual void resize(INDEX nRows, INDEX nCols=1, bool copy=false)=0;
  //* resizes the matrix, zero it out default to ON
  virtual void  reset(INDEX nRows, INDEX nCols=1, bool zero=true)=0;
  //* resizes and copies data
  virtual void copy(const T * ptr, INDEX nRows, INDEX nCols=1)=0;
  //* create restart file
  virtual void write_restart(FILE *f) const=0;    
  //* writes a matlab command to recreate this in a variable named s
  virtual void matlab(std::ostream &o, const std::string &s="M") const;
  //* writes a mathematica command to recreate this in a variable named s
  virtual void mathematica(std::ostream &o, const std::string &s="M") const;

  // output to matlab, with variable name s
  void matlab(const std::string &s="M") const;
  // output to mathematica, with variable name s
  void mathematica(const std::string &s="M") const;

  Matrix<T>& operator+=(const Matrix &r);
  Matrix<T>& operator-=(const Matrix &r);
  
  Matrix<T>& operator*=(const Matrix<T>& R);
  Matrix<T>& operator/=(const Matrix<T>& R);
  Matrix<T>& operator+=(const T v);
  Matrix<T>& operator-=(const T v);
  Matrix<T>& operator*=(const T v);
  Matrix<T>& operator/=(T v);
  Matrix<T>& divide_zero_safe(const Matrix<T>& B);

  Matrix<T>& operator=(const T &v);
  Matrix<T>& operator=(const Matrix<T> &c); 
  virtual void set_all_elements_to(const T &v);
  //* adds a matrix scaled by factor s to this one.
  void add_scaled(const Matrix<T> &A, const T& s);

  //* sets all elements to zero
  Matrix& zero();
  //* sets matrix to the identity
  Matrix& identity(int nrows=0);
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
  //* returns true if Matrix a rows are equal to Matrix b cols
  static bool cols_equals_rows(const Matrix<T> &a, const Matrix<T> &b);
  //* checks if memory is contiguous, only can be false for clone vector
  virtual bool memory_contiguous() const { return true; }
  //* checks if all values are within the prescribed range
  virtual bool check_range(T min, T max) const;

protected:
  virtual void _set_equal(const Matrix<T> &r) = 0;
};

//* Matrix operations 
//@{
//* Sets C as b*C + a*A[transpose?]*B[transpose?]
template<typename T>
void MultAB(const Matrix<T> &A, const Matrix<T> &B, DenseMatrix<T> &C, 
            bool At=0, bool Bt=0, T a=1, T b=0);
//* performs a matrix-vector multiply
template<typename T>
void MultMv(const Matrix<T> &A, const Vector<T> &v, DenseVector<T> &c,
            const bool At, T a, T b);
// returns the inverse of a double precision matrix
DenseMatrix<double> inv(const Matrix<double>& A);
// returns the eigensystem of a pair of double precision matrices
DenseMatrix<double> eigensystem(const Matrix<double>& A, const Matrix<double>& B, DenseMatrix<double> & eVals, bool normalize = true);
// returns the polar decomposition of a double precision matrix
DenseMatrix<double>  polar_decomposition(const Matrix<double>& A, DenseMatrix<double> & rotation, DenseMatrix<double> & stretch, bool leftRotation = true);

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
  for (INDEX i=0; i<A.size(); i++) {
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
            const bool At, const bool Bt, T /* a */, T b)
{
  const INDEX sA[2] = {A.nRows(), A.nCols()};  // m is sA[At] k is sA[!At]
  const INDEX sB[2] = {B.nRows(), B.nCols()};  // k is sB[Bt] n is sB[!Bt]

  const INDEX M=sA[At], K=sB[Bt], N=sB[!Bt]; // M is the number of rows in A or Atrans (sA[At]),
                                    // K is the number of rows in B or Btrans (sB[Bt], sA[!At]),
                                    // N is the number of columns in B or Btrans (sB[!Bt]).

  GCK(A, B, sA[!At]!=K, "MultAB<T> shared index not equal size");
  if (!C.is_size(M,N))
  {
    C.resize(M,N);          // set size of C
    C.zero();
  }
  else C *= b; // Zero C
  for (INDEX p=0; p<M; p++) {
    INDEX p_times_At    = p*At;
    INDEX p_times_notAt = p*!At;
    for (INDEX q=0; q<N; q++) {
        INDEX q_times_Bt    = q*Bt;
        INDEX q_times_notBt = q*!Bt;
      for (INDEX r=0; r<K; r++) {
         INDEX ai = p_times_notAt+r*At;
         INDEX aj = p_times_At+r*!At;
         INDEX bi = r*!Bt+q_times_Bt;
         INDEX bj = r*Bt+q_times_notBt;
         T a_entry = A(ai, aj);
         T b_entry = B(bi, bj);
         T mult    = a_entry * b_entry;
         C(p,q) += mult;
      }
    }
  }
}

//-----------------------------------------------------------------------------
//* output operator
template<typename T>
std::string Matrix<T>::to_string(int p) const
{
  std::string s;
  for (INDEX i=0; i<nRows(); i++) {
    if (i) s += '\n';
    for (INDEX j=0; j<nCols(); j++) {
      //if (j) s+= '\t';
      s += ATC_Utility::to_string((*this)(i,j),p)+" ";
    }
  }
  return s;
}
//-----------------------------------------------------------------------------
//* output operator that wraps the matrix in a nice labeled box
template<typename T>
void Matrix<T>::print(std::ostream &o, const std::string &name, int p) const
{
  o << "------- Begin "<<name<<" -----------------\n";
  this->print(o,p);
  o << "\n------- End "<<name<<" -------------------\n";
}
//-----------------------------------------------------------------------------
//* print operator, use cout by default
template<typename T>
void Matrix<T>::print() const 
{
  print(std::cout); 
}
//-----------------------------------------------------------------------------
//* named print operator, use cout by default
template<typename T>
void Matrix<T>::print(const std::string &name, int p) const 
{
  print(std::cout, name, p); 
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
  int sz=this->size(); for(INDEX i=0; i<sz; i++)
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
  int sz=this->size(); for(INDEX i=0; i<sz; i++)
  {
    double val = R[i];
    R[i] = std::pow(val,n);
  }
  return R;
}
//-----------------------------------------------------------------------------
//* returns the transpose of this matrix (makes a copy)
template <typename T>
DenseMatrix<T> Matrix<T>::transpose()                                     const
{
  DenseMatrix<T> t(this->nCols(), this->nRows());
  int szi = this->nRows();
  int szj = this->nCols(); 
  for (INDEX i = 0; i < szi; i++) 
    for (INDEX j = 0; j < szj; j++)
      t(j,i) = (*this)(i,j);
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
  T v = (*this)[0];
  for (INDEX i=1; i<this->size(); i++) v += (*this)[i];
  return v; 
}
//-----------------------------------------------------------------------------
//* Returns the standard deviation of the matrix
template<typename T>
T Matrix<T>::stdev() const
{
  GCHK(this->size()<2, "Matrix::stdev() size must be > 1");
  T mean = this->mean();
  T diff = (*this)[0]-mean;
  T stdev = diff*diff;
  for (INDEX i=1; i<this->size(); i++) 
  {
    diff = (*this)[i]-mean;
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
  T v = (*this)[0];
  for (INDEX i=1; i<this->size(); i++)   v = std::max(v, (*this)[i]);
  return v; 
}
//-----------------------------------------------------------------------------
//* Returns the minimum of the matrix
template<typename T>
T Matrix<T>::min() const
{
  GCHK(!this->size(), "Matrix::min() size must be > 0");
  T v = (*this)[0];
  for (INDEX i=1; i<this->size(); i++)   v = std::min(v, (*this)[i]);
  return v; 
}
//-----------------------------------------------------------------------------
//* Returns the maximum absolute value of the matrix
template<typename T>
T Matrix<T>::maxabs() const
{
  GCHK(!this->size(), "Matrix::maxabs() size must be > 0");
  T v = (*this)[0];
  for (INDEX i=1; i<this->size(); i++)   v = ATC_Utility::max_abs(v, (*this)[i]);
  return v; 
}
//-----------------------------------------------------------------------------
//* Returns the minimum absoute value of the matrix
template<typename T>
T Matrix<T>::minabs() const
{
  GCHK(!this->size(), "Matrix::minabs() size must be > 0");
  T v = (*this)[0];
  for (INDEX i=1; i<this->size(); i++)   v = ATC_Utility::min_abs(v, (*this)[i]);
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
//* returns the L2 norm of the matrix
template<typename T>
T Matrix<T>::norm_sq()                                                       const    
{
  GCHK(!this->size(), "Matrix::norm() size must be > 0");
  return dot(*this); 
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
  T v = r[0]*(*this)[0];
  for (INDEX i=1; i<this->size(); i++)  v += r[i]*(*this)[i];
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
  T r = (*this)(0,0);
  for (INDEX i=0; i<N; i++)
    r += (*this)(i,i);
  return r;
}
//-----------------------------------------------------------------------------
//* Adds a matrix to this one
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix &r) 
{
  SSCK(*this, r, "operator+= or operator +");
  int sz=this->size(); for(INDEX i=0; i<sz; i++) (*this)[i]+=r[i];
  return *this;
}
//-----------------------------------------------------------------------------
// subtracts a matrix from this one
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix &r) 
{
  SSCK(*this, r, "operator-= or operator -");
  int sz=this->size(); for(INDEX i=0; i<sz; i++) (*this)[i]-=r[i];
  return *this;
}
//-----------------------------------------------------------------------------
// multiplies each element in this by the corresponding element in R
//-----------------------------------------------------------------------------

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& R)
{
  if ((R.nCols()==1) && (this->nCols()>1)) { // multiply every entry in a row by the same value
    int szi = this->nRows();
    int szj = this->nCols(); 
    for (INDEX i = 0; i < szi; i++) 
      for (INDEX j = 0; j < szj; j++)
        {
          (*this)(i,j) *= R[i];
        }
  }
  else if (((R.nCols()==R.size()) && (R.nRows()==R.size())) && !((this->nCols()==this->size()) && (this->nRows()==this->size()))){ 
    int szi = this->nRows();
    int szj = this->nCols(); 
    for (INDEX i = 0; i < szi; i++) 
      for (INDEX j = 0; j < szj; j++)
      {
          (*this)(i,j) *= R[i];
    }
  } 
  else { // multiply each entry by a different value

    int sz = this->size();
    for (INDEX i = 0; i < sz; i++) 
      {
        (*this)[i] *= R[i]; 
      }
  }
  return *this;
}
//-----------------------------------------------------------------------------
// divides each element in this by the corresponding element in R
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator/=(const Matrix<T>& R)
{
  if ((R.nCols()==1) && (this->nCols()>1)) { // divide every entry in a row by the same value
    int szi = this->nRows();
    int szj = this->nCols();
    for (INDEX i = 0; i < szi; i++)
      for (INDEX j = 0; j < szj; j++)
        {
          (*this)(i,j) /= R[i];
        }
  }
  else { // divide each entry by a different value
    SSCK(*this, R, "operator/= or operator/");
    int sz = this->size();
    for(INDEX i = 0; i < sz; i++) 
      {
        GCHK(fabs(R[i])==0,"Operator/: division by zero"); 
        (*this)[i] /= R[i]; 
      }
  }
  return *this;
}
//-----------------------------------------------------------------------------
// divides each element in this by the corresponding element in R unless zero
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::divide_zero_safe(const Matrix<T>& R)
{
  if ((R.nCols()==1) && (this->nCols()>1)) { // divide every entry in a row by the same value
    int szi = this->nRows();
    int szj = this->nCols();
    for (INDEX i = 0; i < szi; i++)
      for (INDEX j = 0; j < szj; j++)
        {
          if(fabs(R[i])!=0) {
            (*this)(i,j) /= R[i];
          }
        }
  }
  else { // divide each entry by a different value
    SSCK(*this, R, "operator/= or operator/");
    int sz = this->size();
    for(INDEX i = 0; i < sz; i++) 
      {
        if(fabs(R[i])!=0) {
          (*this)[i] /= R[i]; 
        }
      }
  }
  return *this;
}
//-----------------------------------------------------------------------------
// scales this matrix by a constant
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T v)       
{
  int sz=this->size(); for(INDEX i=0; i<sz; i++) (*this)[i]*=v;      
  return *this;
}
//-----------------------------------------------------------------------------
// adds a constant to this matrix 
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator+=(const T v)       
{
  int sz=this->size(); for(INDEX i=0; i<sz; i++) (*this)[i]+=v;      
  return *this;
}
//-----------------------------------------------------------------------------
// subtracts a constant to this matrix
//-----------------------------------------------------------------------------
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const T v)       
{
  int sz=this->size(); for(INDEX i=0; i<sz; i++) (*this)[i]-=v;      
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
  int sz=this->size(); for(INDEX i=0; i<sz; i++) (*this)[i] = v;
}
//----------------------------------------------------------------------------
// adds a matrix scaled by factor s to this one.
//----------------------------------------------------------------------------
template <typename T>
void Matrix<T>::add_scaled(const Matrix<T> &A, const T& s)
{
    SSCK(A, *this, "Matrix::add_scaled");
    int sz=this->size(); for(INDEX i=0; i<sz; i++) (*this)[i] += A[i]*s;
}
//-----------------------------------------------------------------------------
//* writes a matlab command to the console
template<typename T>
void Matrix<T>::matlab(const std::string &s) const 
{
  this->matlab(std::cout, s); 
}
//-----------------------------------------------------------------------------
//* Writes a matlab script defining the vector to the stream
template<typename T>
void Matrix<T>::matlab(std::ostream &o, const std::string &s) const
{
  o << s <<"=zeros(" << nRows() << ","<<nCols()<<");\n";
  int szi = this->nRows();
  int szj = this->nCols(); 
  for (INDEX i = 0; i < szi; i++)
    for (INDEX j = 0; j < szj; j++) 
      o << s << "("<<i+1<<","<<j+1<<")=" << (*this)(i,j) << ";\n";
}
//-----------------------------------------------------------------------------
//* writes a mathematica command to the console
template<typename T>
void Matrix<T>::mathematica(const std::string &s) const 
{
  this->mathematica(std::cout, s); 
}
//-----------------------------------------------------------------------------
//* Writes a mathematica script defining the vector to the stream
template<typename T>
void Matrix<T>::mathematica(std::ostream &o, const std::string &s) const
{
  o << s <<" = { \n";
  o.precision(15);
  o << std::fixed;
  for(INDEX i=0; i< nRows(); i++) {
    o <<" { " << (*this)(i,0);
    for(INDEX j=1; j< nCols(); j++) o << ", " << (*this)(i,j);
    if (i+1 == nRows()) { o <<" }  \n"; }
    else                { o <<" }, \n"; }
      
  }
  o << "};\n";
  o << std::scientific;
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
//* sets to identity
template<typename T>
inline Matrix<T>& Matrix<T>::identity(int nrows) 
{ 
  if (nrows == 0) {
    SQCK(*this, "DenseMatrix::inv(), matrix not square"); // check matrix is square
    nrows = nRows();
  }
  reset(nrows,nrows);
  for(INDEX i=0; i< nRows(); i++) (*this)(i,i) = 1;
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
//* returns true if Matrix a rows =  Matrix b cols
template<typename T>
inline bool Matrix<T>::cols_equals_rows(const Matrix<T> &a, const Matrix<T> &b)
{
  return a.nCols() == b.nRows(); 
}
//-----------------------------------------------------------------------------
//* returns true if no value is outside of the range
template<typename T>
inline bool Matrix<T>::check_range(T min, T max) const
{
  for (INDEX i = 0; i < this->nRows(); i++) {
    for (INDEX j = 0; j < this->nCols(); j++) {
      T val = (*this)(i,j);
      if ( (val > max) || (val < min) ) return false;
    }
  }
  return true;
}
//-----------------------------------------------------------------------------
//* Displays indexing error message and quits
template<typename T>
void ierror(const Matrix<T> &a, const char *FILE, int LINE, INDEX i, INDEX j)
{
  std::cout << "Error: Matrix indexing failure ";
  std::cout << "in file: " << FILE << ", line: "<< LINE <<"\n";
  std::cout << "Tried accessing index (" << i << ", " << j <<")\n";
  std::cout << "Matrix size was "<< a.nRows() << "x" << a.nCols() << "\n";
  ERROR_FOR_BACKTRACE
  exit(EXIT_FAILURE);
}
//-----------------------------------------------------------------------------
//* Displays custom message and indexing error and quits
template<typename T>
void ierror(const Matrix<T> &a, INDEX i, INDEX j, const std::string m)
{
  std::cout << m << "\n";
  std::cout << "Tried accessing index (" << i << ", " << j <<")\n";
  std::cout << "Matrix size was "<< a.nRows() << "x" << a.nCols() << "\n";
  ERROR_FOR_BACKTRACE
  exit(EXIT_FAILURE);
}
//-----------------------------------------------------------------------------
//* Displays matrix compatibility error message
template<typename T>
void merror(const Matrix<T> &a, const Matrix<T> &b, const std::string m)
{
  std::cout << "Error: " << m << "\n";
  std::cout << "Matrix sizes were " << a.nRows() << "x" << a.nCols();
  if (&a != &b) std::cout << ", and "<< b.nRows() << "x" << b.nCols();
  std::cout << "\n";
  if (a.size() < 100) a.print("Matrix");
  ERROR_FOR_BACKTRACE
  exit(EXIT_FAILURE);
}

//-----------------------------------------------------------------------------
//* returns upper or lower half of a partitioned matrix
//* A1 is the on-diagonal square matrix, A2 is the off-diagonal matrix
//* rowsIn is the rows to be placed in A1
//* rows is the map for A1, (rows,colsC) is the map for A2

template <typename T>
void Matrix<T>::row_partition(const std::set<int> & rowsIn, 
std::set<int> & rows, std::set<int> & colsC,
DenseMatrix<T> & A1, DenseMatrix<T> & A2, bool complement) const
{
  if (complement) {
    for (INDEX i = 0; i < this->nRows();  i++) {
      if (rowsIn.find(i) == rowsIn.end() ) rows.insert(i);
    }
  }
  else  rows = rowsIn;
  // complement of set "rows" in set of this.cols is "cols"
  for (INDEX i = 0; i < this->nCols();  i++) {
    if (rows.find(i) == rows.end() ) colsC.insert(i);
  }
  // degenerate cases
  if      (int(rows.size()) == this->nCols()) {
    A1 = (*this);
    A2.reset(0,0);
    return;
  }
  else if (rows.size() == 0) {
    A1.reset(0,0);
    A2 = (*this);
    return;
  }
  // non-degenerate case
  int nrows  = rows.size();
  int ncolsC = colsC.size();
  A1.reset(nrows,nrows);
  A2.reset(nrows,ncolsC);
  std::set<int>::const_iterator itrI, itrJ;
  INDEX i =0;
  for (itrI = rows.begin(); itrI != rows.end(); itrI++) {
    INDEX j = 0;
    for (itrJ = rows.begin(); itrJ != rows.end(); itrJ++) {
      A1(i,j) = (*this)(*itrI,*itrJ);
      j++;
    }
    j = 0;
    for (itrJ = colsC.begin(); itrJ != colsC.end(); itrJ++) {
      A2(i,j) = (*this)(*itrI,*itrJ);
      j++;
    }
    i++;
  }
}

template <typename T>
std::set<int> Matrix<T>::row_partition(const std::set<int> & rows, 
DenseMatrix<T> & A1, DenseMatrix<T> & A2) const
{
  // complement of set "rows" in set of this.cols is "cols"
  std::set<int> colsC;
  for (INDEX i = 0; i < this->nCols();  i++) {
    if (rows.find(i) == rows.end() ) colsC.insert(i);
  }
  // degenerate cases
  if      (int(rows.size()) == this->nCols()) {
    A1 = (*this);
    A2.reset(0,0);
    return colsC;
  }
  else if (rows.size() == 0) {
    A1.reset(0,0);
    A2 = (*this);
    return colsC;
  }
  // non-degenerate case
  int nrows  = rows.size();
  int ncolsC = colsC.size();
  A1.reset(nrows,nrows);
  A2.reset(nrows,ncolsC);
  std::set<int>::const_iterator itrI, itrJ;
  INDEX i =0;
  for (itrI = rows.begin(); itrI != rows.end(); itrI++) {
    INDEX j = 0;
    for (itrJ = rows.begin(); itrJ != rows.end(); itrJ++) {
      A1(i,j) = (*this)(*itrI,*itrJ);
      j++;
    }
    j = 0;
    for (itrJ = colsC.begin(); itrJ != colsC.end(); itrJ++) {
      A2(i,j) = (*this)(*itrI,*itrJ);
      j++;
    }
    i++;
  }
  return colsC;
}

//-----------------------------------------------------------------------------
//* returns row & column mapped matrix
template <typename T>
void Matrix<T>::map(const std::set<int> & rows, const std::set<int> & cols, 
DenseMatrix<T> & A ) const
{
  if      (rows.size() == 0 || cols.size() == 0 ) {
    A.reset(0,0);
    return;
  }
  int nrows = rows.size();
  int ncols = cols.size();
  A.reset(nrows,ncols);
  std::set<int>::const_iterator itrI, itrJ;
  INDEX i =0;
  for (itrI = rows.begin(); itrI != rows.end(); itrI++) {
    INDEX j = 0;
    for (itrJ = cols.begin(); itrJ != cols.end(); itrJ++) {
      A(i,j) = (*this)(*itrI,*itrJ);
      j++;
    }
    i++;
  }
}
//-----------------------------------------------------------------------------
//* inserts elements from a smaller matrix 
template <typename T>
void Matrix<T>::insert(const std::set<int> & rows, const std::set<int> & cols, 
const DenseMatrix<T> & A ) 
{
  if      (rows.size() == 0 || cols.size() == 0 )  return;
  std::set<int>::const_iterator itrI, itrJ;
  int i =0;
  for (itrI = rows.begin(); itrI != rows.end(); itrI++) {
    int j = 0;
    for (itrJ = cols.begin(); itrJ != cols.end(); itrJ++) {
      (*this)(*itrI,*itrJ) = A(i,j);
//std::cout << *itrI << " " << *itrJ << " : " << (*this)(*itrI,*itrJ) << "\n";
      j++;
    }
    i++;
  }
}
//-----------------------------------------------------------------------------
//* assemble elements from a smaller matrix 
template <typename T>
void Matrix<T>::assemble(const std::set<int> & rows, const std::set<int> & cols, 
const DenseMatrix<T> & A ) 
{
  if      (rows.size() == 0 || cols.size() == 0 )  return;
  std::set<int>::const_iterator itrI, itrJ;
  int i =0;
  for (itrI = rows.begin(); itrI != rows.end(); itrI++) {
    int j = 0;
    for (itrJ = cols.begin(); itrJ != cols.end(); itrJ++) {
      (*this)(*itrI,*itrJ) += A(i,j);
      j++;
    }
    i++;
  }
}
//-----------------------------------------------------------------------------

} // end namespace

#endif
