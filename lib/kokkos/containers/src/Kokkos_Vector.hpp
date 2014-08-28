/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_VECTOR_HPP
#define KOKKOS_VECTOR_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_DualView.hpp>

/* Drop in replacement for std::vector based on Kokkos::DualView
 * Most functions only work on the host (it will not compile if called from device kernel)
 *
 */
  namespace Kokkos {

template <typename Scalar, class Device = Kokkos::DefaultExecutionSpace >
class vector : public DualView<Scalar*,LayoutLeft,Device> {
public:
  typedef Device device_type;
  typedef Scalar value_type;
  typedef Scalar* pointer;
  typedef const Scalar* const_pointer;
  typedef Scalar* reference;
  typedef const Scalar* const_reference;
  typedef Scalar* iterator;
  typedef const Scalar* const_iterator;

private:
  size_t _size;
  typedef size_t size_type;
  float _extra_storage;
  typedef DualView<Scalar*,LayoutLeft,Device> DV;


public:
  inline Scalar& operator() (int i) const {return DV::h_view(i);};
  inline Scalar& operator[] (int i) const {return DV::h_view(i);};


  /* Member functions which behave like std::vector functions */

  vector():DV() {
    _size = 0;
    _extra_storage = 1.1;
    DV::modified_host = 1;
  };


  vector(int n, Scalar val=Scalar()):DualView<Scalar*,LayoutLeft,Device>("Vector",size_t(n*(1.1))) {
    _size = n;
    _extra_storage = 1.1;
    DV::modified_host = 1;

    assign(n,val);
  }


  void resize(size_t n) {
    if(n>=capacity())
      DV::resize(size_t (n*_extra_storage));
    _size = n;
  }

  void resize(size_t n, const Scalar& val) {
    assign(n,val);
  }

  void assign (size_t n, const Scalar& val) {

    /* Resize if necessary (behavour of std:vector) */

    if(n>capacity())
      DV::resize(size_t (n*_extra_storage));
    _size = n;

          /* Assign value either on host or on device */

    if( DV::modified_host >= DV::modified_device ) {
      set_functor_host f(DV::h_view,val);
      parallel_for(n,f);
      DV::t_host::device_type::fence();
      DV::modified_host++;
    } else {
      set_functor f(DV::d_view,val);
      parallel_for(n,f);
      DV::t_dev::device_type::fence();
      DV::modified_device++;
    }
  }

  void reserve(size_t n) {
    DV::resize(size_t (n*_extra_storage));
  }

  void push_back(Scalar val) {
    DV::modified_host++;
    if(_size == capacity()) {
      size_t new_size = _size*_extra_storage;
      if(new_size == _size) new_size++;
      DV::resize(new_size);
    }

    DV::h_view(_size) = val;
    _size++;

  };

  void pop_back() {
    _size--;
  };

  void clear() {
    _size = 0;
  }

  size_type size() const {return _size;};
  size_type max_size() const {return 2000000000;}
  size_type capacity() const {return DV::capacity();};
  bool empty() const {return _size==0;};

  iterator begin() const {return &DV::h_view(0);};

  iterator end() const {return &DV::h_view(_size);};


  /* std::algorithms wich work originally with iterators, here they are implemented as member functions */

  size_t
  lower_bound (const size_t& start,
               const size_t& theEnd,
               const Scalar& comp_val) const
  {
    int lower = start; // FIXME (mfh 24 Apr 2014) narrowing conversion
    int upper = _size > theEnd? theEnd : _size-1; // FIXME (mfh 24 Apr 2014) narrowing conversion
    if (upper <= lower) {
      return theEnd;
    }

    Scalar lower_val = DV::h_view(lower);
    Scalar upper_val = DV::h_view(upper);
    size_t idx = (upper+lower)/2;
    Scalar val = DV::h_view(idx);
    if(val>upper_val) return upper;
    if(val<lower_val) return start;

    while(upper>lower) {
      if(comp_val>val) {
        lower = ++idx;
      } else {
        upper = idx;
      }
      idx = (upper+lower)/2;
      val = DV::h_view(idx);
    }
    return idx;
  }

  bool is_sorted() {
    for(int i=0;i<_size-1;i++) {
      if(DV::h_view(i)>DV::h_view(i+1)) return false;
    }
    return true;
  }

  iterator find(Scalar val) const {
    if(_size == 0) return end();

    int upper,lower,current;
    current = _size/2;
    upper = _size-1;
    lower = 0;

    if((val<DV::h_view(0)) || (val>DV::h_view(_size-1)) ) return end();

    while(upper>lower)
    {
      if(val>DV::h_view(current)) lower = current+1;
      else upper = current;
      current = (upper+lower)/2;
    }

    if(val==DV::h_view(current)) return &DV::h_view(current);
    else return end();
  }

  /* Additional functions for data management */

  void device_to_host(){
    deep_copy(DV::h_view,DV::d_view);
  }
  void host_to_device() const {
    deep_copy(DV::d_view,DV::h_view);
  }

  void on_host() {
    DV::modified_host = DV::modified_device + 1;
  }
  void on_device() {
    DV::modified_device = DV::modified_host + 1;
  }

  void set_overallocation(float extra) {
    _extra_storage = 1.0 + extra;
  }


public:
  struct set_functor {
    typedef typename DV::t_dev::device_type device_type;
    typename DV::t_dev _data;
    Scalar _val;

    set_functor(typename DV::t_dev data, Scalar val) :
      _data(data),_val(val) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int &i) const {
      _data(i) = _val;
    }
  };

  struct set_functor_host {
    typedef typename DV::t_host::device_type device_type;
    typename DV::t_host _data;
    Scalar _val;

    set_functor_host(typename DV::t_host data, Scalar val) :
      _data(data),_val(val) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int &i) const {
      _data(i) = _val;
    }
  };

};


}
#endif
