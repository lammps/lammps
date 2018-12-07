#ifndef RESOURCE_PTR_H
#define RESOURCE_PTR_H

#include <stddef.h>
#include <stdlib.h>
#include "memory.h"

namespace LAMMPS_NS {
namespace detail {
template <typename T>
struct default_delete {
  void operator() (T* ptr) const {
    delete ptr;
  }
};

template <typename T>
struct default_delete<T[]> {
  void operator() (T* ptr) const {
    delete[] ptr;
  }
};

template <typename T, size_t dim>
struct multidimensional_delete {};

template <typename T, size_t dim>
struct multidimensional_delete<T[], dim> {
  void operator() (T* ptr) const {
    if (ptr) multidimensional_delete<T[], dim-1>().operator()(ptr[0]);
    delete[] ptr;
  }
};

template <typename T>
struct multidimensional_delete<T[], 1> {
  void operator() (T* ptr) const {
    delete[] ptr;
  }
};

template <typename T>
struct multidimensional_delete<T[], 0> {};

template <typename T, size_t dim>
struct C_multidimensional_delete {};

template <typename T, size_t dim>
struct C_multidimensional_delete<T[], dim> {
  void operator() (T* ptr) const {
    if (ptr) C_multidimensional_delete<T[], dim-1>().operator()(ptr[0]);
    free (ptr);
  }
};

template <typename T>
struct C_multidimensional_delete<T[], 1> {
  void operator() (T* ptr) const {
    free (ptr);
  }
};

template <typename T>
struct C_multidimensional_delete<T[], 0> {};

template <typename T>
struct C_noncontiguous_2d_delete {};

template <typename T>
struct C_noncontiguous_2d_delete<T[]> {
  size_t array_size;
  C_noncontiguous_2d_delete (size_t size = size_t()) : array_size(size) {}
  void operator() (T* ptr) const {
    for (int i=0; i<array_size; ++i) free (ptr[i]);
    free (ptr);
  }
};

template <typename T>
struct LAMMPS_delete {
  void operator() (T* ptr) const {
    Memory::sfree (ptr);
  }
};

template <typename T, size_t dim>
struct LAMMPS_multidimensional_delete {};

template <typename T, size_t dim>
struct LAMMPS_multidimensional_delete<T[], dim> {
  void operator() (T* ptr) const {
    if (ptr) LAMMPS_multidimensional_delete<T[], dim-1>().operator()(ptr[0]);
    Memory::sfree (ptr);
  }
};

template <typename T>
struct LAMMPS_multidimensional_delete<T[], 1> {
  void operator() (T* ptr) const {
    Memory::sfree (ptr);
  }
};

template <typename T>
struct LAMMPS_multidimensional_delete<T[], 0> {};

template<bool B, class T, class F>
struct conditional {typedef T type;};
 
template<class T, class F>
struct conditional<false, T, F> {typedef F type;};

template <class T> struct is_reference     {static const bool value = false;};
template <class T> struct is_reference<T&> {static const bool value = true;};

template<class T>
struct remove_extent {typedef T type;};
 
template<class T>
struct remove_extent<T[]> {typedef T type;};
 
template<class T, std::size_t N>
struct remove_extent<T[N]> {typedef T type;};

template <typename T>
struct C_delete {
  typedef typename remove_extent<T>::type* pointer_type;
  void operator() (pointer_type ptr) const {
    free (ptr);
  }
};

}

template <typename T, typename D = detail::default_delete<T> >
class resource_ptr {
  typedef typename detail::conditional<detail::is_reference<D>::value,
                                       D,
                                       const D&>::type deleter_ref_type;
public:
  typedef T element_type;
  typedef T* pointer_type;
  typedef D deleter_type;
  
  resource_ptr (pointer_type p = pointer_type()) : _ptr(p), _deleter() {}
  
  resource_ptr (pointer_type p, deleter_ref_type d) : _ptr(p), _deleter(d) {}
  
  ~resource_ptr () {
    if (_ptr != pointer_type()) {
      _deleter (_ptr);
      _ptr = pointer_type();
    }
  }
  
  bool null () const {return _ptr==pointer_type();}
  
  pointer_type get () const {return _ptr;}
  
  deleter_type & get_deleter () {return _deleter;}
  const deleter_type & get_deleter () const {return _deleter;}
  
  pointer_type* get_address () {return &_ptr;}
  
  element_type & operator* () const {return *_ptr;}
  pointer_type operator-> () const {return _ptr;}
  
  pointer_type release () {
    pointer_type p = _ptr;
    _ptr = pointer_type();
    return p;
  }
  
  void reset (pointer_type p = pointer_type()) {
    if (_ptr != pointer_type()) _deleter (_ptr);
    _ptr = p;
  }
protected:
  pointer_type _ptr;
  deleter_type _deleter;
  
private:
  resource_ptr (const resource_ptr &);
  resource_ptr & operator= (const resource_ptr &);
};

template <typename T, typename D>
bool operator== (const resource_ptr<T, D> &lhs, typename resource_ptr<T, D>::pointer_type rhs) {
  return lhs.get() == rhs;
}

template <typename T, typename D>
bool operator== (typename resource_ptr<T, D>::pointer_type lhs, const resource_ptr<T, D> &rhs) {
  return lhs == rhs.get();
}

template <typename T, typename D>
bool operator== (const resource_ptr<T, D> &lhs, resource_ptr<T, D> &rhs) {
  return lhs.get() == rhs.get();
}

template <typename T, typename D>
bool operator!= (const resource_ptr<T, D> &lhs, typename resource_ptr<T, D>::pointer_type rhs) {
  return lhs.get() != rhs;
}

template <typename T, typename D>
bool operator!= (typename resource_ptr<T, D>::pointer_type lhs, const resource_ptr<T, D> &rhs) {
  return lhs != rhs.get();
}

template <typename T, typename D>
bool operator!= (const resource_ptr<T, D> &lhs, resource_ptr<T, D> &rhs) {
  return lhs.get() != rhs.get();
}

template <typename T, typename D>
class resource_ptr<T[], D> {
  typedef typename detail::conditional<detail::is_reference<D>::value,
                                       D,
                                       const D&>::type deleter_ref_type;
public:
  typedef T element_type;
  typedef T* pointer_type;
  typedef D deleter_type;
  
  resource_ptr (pointer_type p = pointer_type()) : _ptr(p), _deleter() {}
  
  resource_ptr (pointer_type p, deleter_ref_type d) : _ptr(p), _deleter(d) {}
  
  ~resource_ptr () {
    if (_ptr != pointer_type()) {
      _deleter (_ptr);
      _ptr = pointer_type();
    }
  }
  
  pointer_type get () const {return _ptr;}
  
  deleter_type & get_deleter () {return _deleter;}
  const deleter_type & get_deleter () const {return _deleter;}
  
  pointer_type* get_address () {return &_ptr;}
  
  element_type & operator[] (size_t i) const {return _ptr[i];}
  
  pointer_type release () {
    pointer_type p = _ptr;
    _ptr = pointer_type();
    return p;
  }
  
  void reset (pointer_type p = pointer_type()) {
    if (_ptr != pointer_type()) _deleter (_ptr);
    _ptr = p;
  }
protected:
  pointer_type _ptr;
  deleter_type _deleter;
  
private:
  resource_ptr (const resource_ptr &);
  resource_ptr & operator= (const resource_ptr &);
};

template <typename T>
struct memory_ptr : public resource_ptr<T, detail::default_delete<T> > {
  typedef T* pointer_type;
  memory_ptr () : resource_ptr<T, detail::default_delete<T> > () {}
  memory_ptr (pointer_type p) : resource_ptr<T, detail::default_delete<T> > (p) {}
  void operator= (pointer_type p) {this->reset (p);}
};

template <typename T>
struct memory_ptr<T[]> : public resource_ptr<T[], detail::default_delete<T[]> > {
  typedef T* pointer_type;
  memory_ptr () : resource_ptr<T[], detail::default_delete<T[]> > () {}
  memory_ptr (pointer_type p) : resource_ptr<T[], detail::default_delete<T[]> > (p) {}
  void operator= (pointer_type p) {this->reset (p);}
};

template <typename T>
struct C_memory_ptr : public resource_ptr<T, detail::C_delete<T> > {
  typedef typename detail::remove_extent<T>::type* pointer_type;
  C_memory_ptr () : resource_ptr<T, detail::C_delete<T> > () {}
  C_memory_ptr (pointer_type p) : resource_ptr<T, detail::C_delete<T> > (p) {}
  void operator= (pointer_type p) {this->reset (p);}
};

template <typename T>
struct C_memory_ptr_noncontiguous_2d : public resource_ptr<T, detail::C_noncontiguous_2d_delete<T> > {
  typedef typename detail::remove_extent<T>::type pointer_type;
  typedef typename detail::remove_extent<T>::type* pointer2pointer_type;
  typedef detail::C_noncontiguous_2d_delete<T> deleter_type;
  C_memory_ptr_noncontiguous_2d ()
  : resource_ptr<T, deleter_type> (pointer2pointer_type(), deleter_type()) {}
  C_memory_ptr_noncontiguous_2d (pointer2pointer_type p, size_t array_size)
  : resource_ptr<T, deleter_type> (p, deleter_type(array_size)) {}
  void reset_ptr_and_size (pointer2pointer_type p, size_t array_size) {
    this->reset (p);
    this->get_deleter() = deleter_type(array_size);
  }
};
}

#endif
