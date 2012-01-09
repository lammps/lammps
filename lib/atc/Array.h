#ifndef ARRAY_H
#define ARRAY_H

//#include<stdlib.h>
//#include<stdio.h>
#include<iostream>
#include<string>
#include<cstdio>

template<typename T>
class Array {
public:
   Array();
   Array(int len);
   Array(const Array<T>& A);
  ~Array();

   // Resize and reinitialize the array 
   void reset(int len);
   // Access method to get the element i:
   T& operator() (int i);
   // Access method to get the element i:
   const T&  operator() (int i) const;
   // Assignment
   Array<T>& operator= (const Array<T> &other);
   Array<T>& operator= (const T &value);
   // Get length of array
   int get_length() const;
   int size() const;
   // Do I have this element?
   bool has_member(T val) const;
   // Return pointer to internal data
   const T* get_data() const;
   T* get_ptr() const;
   // print
   void print(std::string name = "") const;
   // Dump templated type to disk; operation not safe for all types   
   void write_restart(FILE *f) const;

private:
   int len_;
   T *data_;
};

template<typename T>
Array<T>::Array(void) {
   len_  = 0;
   data_ = NULL;
}

template<typename T>
Array<T>::Array(int len) {
   len_  = len;
   data_ = new T[len_];
}

template<typename T>
Array<T>::Array(const Array<T>& A) {
   len_ = A.len_;
   if (A.data_==NULL)
      data_ = NULL;
   else {
      data_  = new T[len_];
      for(int i=0;i<len_;i++)
         data_[i] = A.data_[i];
   }
}

template<typename T>
void Array<T>::reset(int len) {
   if (len_ == len) { // no size change; don't realloc memory
      return;
   }
   else { // size change, realloc memory
      len_ = len;
      if (data_ != NULL)
        delete[] data_;
      if (len_ > 0)
        data_ = new T[len_];
      else {
        data_ = NULL;
        len_  = 0;
      }
   }
}

template<typename T>
T& Array<T>::operator() (int i) {
   // Array bounds checking
   return data_[i];
}

template<typename T>
Array<T>& Array<T>::operator= (const Array<T> &other) {
  if (data_ == NULL) { // initialize my internal storage to match LHS
     len_  = other.len_;
     if (other.data_==NULL)
	data_ = NULL;
     else
     	data_ = new T[len_];
  }
  for(int i=0;i<len_;i++)
     data_[i] = other.data_[i];
  return *this;
}

template<typename T>
Array<T>& Array<T>::operator= (const T &value) {
  for(int i=0;i<len_;i++)
       data_[i] = value;
  return *this;
}

template<typename T>
const T& Array<T>::operator() (int i) const {
   // Array bounds checking
   return data_[i];
}

template<typename T>
int Array<T>::get_length(void) const {
   return len_;
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
void Array<T>::write_restart(FILE *f) const {
  fwrite(&len_,sizeof(int),1,f);
  if (len_ > 0)
    fwrite(data_,sizeof(T),len_,f);
}

template<typename T>
const T* Array<T>::get_data() const {
   return data_;
}
template<typename T>
T* Array<T>::get_ptr() const {
   return data_;
}

template<typename T>
Array<T>::~Array() {
  if (data_ != NULL)
     delete[] data_;
}

template<typename T>
void Array<T>::print(std::string name) const {
  std::cout << "------- Begin "<<name<<" -----------------\n";
  if (data_ != NULL) {
    for(int i=0;i<len_;i++) std::cout << data_[i] << " ";
    std::cout << "\n";
  }
  std::cout << "\n------- End "<<name<<" -------------------\n";
}
#endif // Array.h
