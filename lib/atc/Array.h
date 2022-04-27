#ifndef ARRAY_H
#define ARRAY_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <cstdio>

// for macros
#include "MatrixDef.h"

namespace ATC_matrix {

  /**
   *  @class  Array
   *  @brief  Base class for creating, sizing and operating on 1-D arrays of data
   */

template<typename T>
class Array {
public:
   Array();
   Array(int len);
   Array(const Array<T>& A);
   virtual ~Array();

   // Resize and reinitialize the array
   virtual void reset(int len);
   //* resizes the matrix, copy what fits default to OFF
   virtual void resize(int len, bool copy=false);
   // Access method to get the element i:
   T& operator() (int i);
   const T&  operator() (int i) const;
   // Assignment
   virtual Array<T>& operator= (const Array<T> &other);
   virtual Array<T>& operator= (const T &value);
   // Get length of array
   int size() const;
   // Do I have this element?
   bool has_member(T val) const;
   // range
   bool check_range(T min, T max) const;
   void range(T & min, T & max) const;
   // search an ordered array
   int index(T& val) const;
   // Return pointer to internal data
   const T* data() const;
   T* ptr() const;
   // print
   void print(std::string name = "") const;
   // Dump templated type to disk; operation not safe for all types
   void write_restart(FILE *f) const;

protected:
   int len_;
   T *data_;
};

template<typename T>
class AliasArray {
public:
   AliasArray();
   AliasArray(const Array<T>& A);
   AliasArray(const AliasArray<T>& A);
   AliasArray(int len, T * A);
   virtual ~AliasArray();
   virtual AliasArray<T>& operator= (const Array<T> &other);
   virtual AliasArray<T>& operator= (const T &value);

   const T&  operator() (int i) const;
   int size() const;
   T* ptr() const;
protected:
   int len_;
   T *data_;
};

template<typename T>
Array<T>::Array(void) {
   len_  = 0;
   data_ = nullptr;
}

template<typename T>
Array<T>::Array(int len) {
   len_  = len;
   data_ = new T[len_];
}

template<typename T>
Array<T>::Array(const Array<T>& A) {
   len_ = A.len_;
   if (A.data_==nullptr)
      data_ = nullptr;
   else {
      data_  = new T[len_];
      for(int i=0;i<len_;i++)
         data_[i] = A.data_[i];
   }
}

template<typename T>
Array<T>::~Array() {
  if (data_ != nullptr) delete[] data_;
}

template<typename T>
void Array<T>::reset(int len) {
   if (len_ == len) { // no size change; don't realloc memory
      return;
   }
   else { // size change, realloc memory
      len_ = len;
      if (data_ != nullptr)
        delete[] data_;
      if (len_ > 0)
        data_ = new T[len_];
      else {
        data_ = nullptr;
        len_  = 0;
      }
   }
}

template<typename T>
void Array<T>::resize(int len, bool copy) {
   if (len_ == len) { // no size change; don't realloc memory
      return;
   }
   else { // size change, realloc memory
      len_ = len;
      if (len_ > 0) {
        if (copy && data_ != nullptr) {
          Array<T> temp(*this);
          delete[] data_;
          data_ = new T[len_];
          for (int i = 0 ; i < len_; i++) {
            if (i < temp.size())
              data_[i] = temp.data_[i];
          }
        }
        else {
          if (data_ != nullptr) delete[] data_;
          data_ = new T[len_];
        }
      }
      else {
        data_ = nullptr;
        len_  = 0;
      }
   }
}

template<typename T>
T& Array<T>::operator() (int i) {
   return data_[i];
}

template<typename T>
Array<T>& Array<T>::operator= (const Array<T> &other) {
  if (data_ == nullptr) { // initialize my internal storage to match LHS
     len_  = other.len_;
     if (other.data_==nullptr)
        data_ = nullptr;
     else
        data_ = new T[len_];
  }
  for(int i=0;i<len_;i++)
     data_[i] = other.data_[i];
  return *this;
}

template<typename T>
Array<T>& Array<T>::operator= (const T &value) {
  for(int i=0;i<len_;i++) data_[i] = value;
  return *this;
}

template<typename T>
const T& Array<T>::operator() (int i) const {
   return data_[i];
}

template<typename T>
int Array<T>::size(void) const {
   return len_;
}

template<typename T>
bool Array<T>::has_member(T val) const {
   int i;
   bool retval = false;
   for(i=0;i<len_;i++)
      if (val == data_[i])
         retval = true;
   return(retval);
}

template<typename T>
bool Array<T>::check_range(T min, T max) const {
   int i;
   for(i=0;i<len_;i++) {
      T val = data_[i];
      if      (val > max) return false;
      else if (val < min) return false;
   }
   return true;
}

template<typename T>
void Array<T>::range(T& min, T& max) const {
   int i;
   min = max = data_[0];
   for(i=1;i<len_;i++) {
      T val = data_[i];
      if      (val > max) max = val;
      else if (val < min) min = val;
   }
}


template<typename T>
int Array<T>::index(T& val) const {
   int idx = -1;
   int i;
   for(i=0;i<len_;i++) {
      T x = data_[i];
      if (val <= x) return idx;
      idx++;
   }
   return idx;
}

template<typename T>
void Array<T>::write_restart(FILE *f) const {
  fwrite(&len_,sizeof(int),1,f);
  if (len_ > 0)
    fwrite(data_,sizeof(T),len_,f);
}

template<typename T>
const T* Array<T>::data() const {
   return data_;
}
template<typename T>
T* Array<T>::ptr() const {
   return data_;
}

template<typename T>
void Array<T>::print(std::string name) const {
  std::cout << "------- Begin "<<name<<" -----------------\n";
  if (data_ != nullptr) {
    for(int i=0;i<len_;i++) std::cout << data_[i] << " ";
    std::cout << "\n";
  }
  std::cout << "\n------- End "<<name<<" -------------------\n";
}

template<typename T>
AliasArray<T>::AliasArray(void) {
}

template<typename T>
AliasArray<T>::AliasArray(const AliasArray<T> & other) {
  len_  = other.size();
  data_ = other.ptr();
}

// for a mem continguous slice
template<typename T>
AliasArray<T>::AliasArray(int len, T * ptr) {
  len_  = len;
  data_ = ptr;
}

template<typename T>
AliasArray<T>::AliasArray(const Array<T>& A) {
  len_  = A.len_;
  data_ = A.ptr();
}

template<typename T>
AliasArray<T>::~AliasArray(void) {
  len_  = 0;
  data_ = nullptr; // trick base class into not deleting parent data
}

template<typename T>
AliasArray<T>& AliasArray<T>::operator= (const Array<T> &other) {
  len_  = other.size();
  data_ = other.ptr();
  return *this;
}

template<typename T>
AliasArray<T>& AliasArray<T>::operator= (const T &value) {
  for(int i=0;i < len_;i++)
    data_[i] = value;
  return *this;
}

template<typename T>
const T& AliasArray<T>::operator() (int i) const {
   return data_[i];
}

template<typename T>
int AliasArray<T>::size(void) const {
   return len_;
}

template<typename T>
T* AliasArray<T>::ptr() const {
   return data_;
}


} // end namespace
#endif // Array.h
