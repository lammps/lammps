#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <cstdlib>
#include <string>
#include <iostream>
#include <cstdlib>
#include <stdio.h>

#include "Array.h"

// for macros
#include "MatrixDef.h"

namespace ATC_matrix {

  /**
   *  @class  Array2D
   *  @brief  Base class for creating, sizing and operating on 2-D arrays of data
   */

template<typename T>
class Array2D {
public: 
   Array2D();
   Array2D(int nrows, int ncols);
   Array2D(const Array2D<T>& A); // copy constructor
  ~Array2D();
 
   // Resize and reinitalize matrix
   void reset(int nrows, int ncols);
   // Access method to get the (i,j) element:
   T& operator() (int i, int j);       
   // Access method to get the i-th col 
   AliasArray<T> column(int i) const;       
   // Access method to get the (i,j) element:
   const T& operator() (int i, int j) const;       
   // Copy operator
   Array2D<T>& operator= (const Array2D<T>& other);
   // assignment operator
   Array2D<T>& operator= (const T other);
   // Get size of Array2D
   int nRows() const;
   int nCols() const;
   // Do I have this element?
   bool has_member(T val) const;
   // print
   void print(std::string name ="") const;
   // Dump templated type to disk; operation not safe for all types
   void write_restart(FILE *f) const;

private: 
   int nrows_, ncols_;
   T *data_; 
}; 

template<typename T>
Array2D<T>::Array2D() {
   nrows_ = 0;
   ncols_ = 0;
   data_  = NULL;
}

template<typename T>
Array2D<T>::Array2D(int nrows, int ncols) {
   nrows_ = nrows;
   ncols_ = ncols;
   data_  = new T[nrows_ * ncols_];
}

template<typename T>
Array2D<T>::Array2D(const Array2D<T>& A) {
   nrows_ = A.nrows_;
   ncols_ = A.ncols_;
   if (A.data_==NULL)
      data_ = NULL;
   else {
      data_  = new T[nrows_ * ncols_];
      for(int i=0;i<nrows_*ncols_;i++)
         data_[i] = A.data_[i];
   }
}

template<typename T>
void Array2D<T>::reset(int nrows, int ncols) {
   if (nrows_ == nrows && ncols_ == ncols) { // no size change; don't realloc memory
      return;
   }
   else { // size changed; realloc memory
      nrows_ = nrows;
      ncols_ = ncols;
      if (data_ != NULL)
         delete [] data_;
      if (ncols_ > 0 && nrows_ > 0)
         data_ = new T[nrows_ * ncols_];
      else {
         data_ = NULL;
         nrows_ = 0;
         ncols_ = 0;
      }
   }
}

template<typename T>
T& Array2D<T>::operator() (int row, int col) {
   // Array bounds checking
   return data_[col*nrows_ + row];
}

template<typename T>
const T& Array2D<T>::operator() (int row, int col) const {
   // Array bounds checking
   return data_[col*nrows_ + row];
}

template<typename T>
AliasArray<T> Array2D<T>::column(int col) const {
   // Array bounds checking
   return AliasArray<T>(nrows_,&(data_[col*nrows_]));
}

template<typename T>
Array2D<T>& Array2D<T>::operator= (const Array2D<T>& other) {
   if (data_ == NULL) {  // initialize my internal storage to match LHS
      nrows_ = other.nrows_;
      ncols_ = other.ncols_;
      if (other.data_==NULL)
         data_ = NULL;
      else
         data_  = new T[nrows_ * ncols_];
   }
   for(int i=0;i<nrows_*ncols_;i++)
      data_[i] = other.data_[i];
   return *this;
}

template<typename T>
Array2D<T>& Array2D<T>::operator= (const T other) {
  for(int i=0;i<nrows_*ncols_;i++)
     data_[i] = other;
  return *this;
}

template<typename T>
int Array2D<T>::nRows() const {
   return nrows_;
}

template<typename T>
int Array2D<T>::nCols() const {
   return ncols_;
}

template<typename T>
bool Array2D<T>::has_member(T val) const {
   int i;
   bool retval = false;
   for(i=0;i<nrows_*ncols_;i++)
      if (val == data_[i])
         retval = true;
   return(retval);
}

template<typename T>
void Array2D<T>::write_restart(FILE *f) const {
  fwrite(&nrows_,sizeof(int),1,f);
  fwrite(&ncols_,sizeof(int),1,f);
  if (nrows_*ncols_ > 0)
    fwrite(data_,sizeof(T),nrows_*ncols_,f);
}

template<typename T>
Array2D<T>::~Array2D() {
   if (data_ != NULL)
     delete[] data_;
} 

template<typename T>
void Array2D<T>::print(std::string name) const {
  std::cout << "------- Begin "<<name<<" -----------------\n";
  if (data_ != NULL) {
    for(int col=0;col<ncols_;col++) {
      for(int row=0;row<nrows_;row++) {
       std::cout << data_[col*nrows_ + row] << " ";
      }
      std::cout << "\n";
    }
  }
  std::cout << "\n------- End "<<name<<" -------------------\n";

}

} // end namespace
#endif // Array2D.h
