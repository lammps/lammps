/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <hc.hpp>
//#include <hsa_atomic.h>

#ifdef KOKKOS_ENABLE_ROCM_ATOMICS
namespace Kokkos {
  //ROCm can do:
  //Types int/unsigned int
  //variants: atomic_exchange/compare_exchange/fetch_add/fetch_sub/fetch_max/fetch_min/fetch_and/fetch_or/fetch_xor/fetch_inc/fetch_dec 


  KOKKOS_INLINE_FUNCTION
  int atomic_exchange(int* dest, const int& val) {
    return hc::atomic_exchange_int(dest, val);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned int atomic_exchange(unsigned int* dest, const unsigned int& val) {
    return hc::atomic_exchange_unsigned(dest, val);
  }

  KOKKOS_INLINE_FUNCTION
  int64_t atomic_exchange(int64_t* dest, const int64_t& val) {
    return (int64_t)hc::atomic_exchange_uint64((uint64_t*)dest, (const uint64_t&)val);
  }

  KOKKOS_INLINE_FUNCTION
  uint64_t atomic_exchange(uint64_t* dest, const uint64_t& val) {
    return hc::atomic_exchange_uint64(dest, val);
  }

  KOKKOS_INLINE_FUNCTION
  long long atomic_exchange(long long* dest, const long long& val) {
    return (long long)hc::atomic_exchange_uint64((uint64_t*)dest, (const uint64_t&)val);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned long long atomic_exchange(unsigned long long* dest, const unsigned long long& val) {
    return (unsigned long long)hc::atomic_exchange_uint64((uint64_t*)dest, (const uint64_t&)val);
  }

  KOKKOS_INLINE_FUNCTION
  float atomic_exchange(float* dest, const float& val) {
    union U {
      int i ;
      float f ;
      KOKKOS_INLINE_FUNCTION U() {};
    } idest,ival;
    idest.f = *dest;
    ival.f = val;
    idest.i = hc::atomic_exchange_int((int*)dest, ival.i);
    return idest.f;
  }

  KOKKOS_INLINE_FUNCTION
  double atomic_exchange(double* dest, const double& val) {
    union U {
      uint64_t i ;
      double d ;
      KOKKOS_INLINE_FUNCTION U() {};
    } idest,ival;
    idest.d = *dest;
    ival.d = val;
    idest.i = hc::atomic_exchange_uint64((uint64_t*)dest, ival.i);
    return idest.d;
  }

  KOKKOS_INLINE_FUNCTION
  int atomic_compare_exchange(int* dest, int compare, const int& val);

  KOKKOS_INLINE_FUNCTION
  int64_t atomic_compare_exchange(int64_t* dest, int64_t compare, const int64_t& val);

  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_exchange(T* dest, typename std::enable_if<sizeof(T) == sizeof(int), const T&>::type val) {
    union U {
      int i ;
      T t ;
      KOKKOS_INLINE_FUNCTION U() {};
    } assume , oldval , newval ;

    oldval.t = *dest ;
    assume.i = oldval.i ;
    newval.t = val ;
    atomic_compare_exchange( (int*)(dest) , assume.i, newval.i );

    return oldval.t ;    
  }

  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_exchange(T* dest, typename std::enable_if<sizeof(T) != sizeof(int) && sizeof(T) == sizeof(int64_t), const T&>::type val) {
    union U {
      uint64_t i ;
      T t ;
      KOKKOS_INLINE_FUNCTION U() {};
    } assume , oldval , newval ;

    oldval.t = *dest ;

    assume.i = oldval.i ;
    newval.t = val ;
    atomic_compare_exchange( (int64_t*)(dest) , assume.i, newval.i );

    return oldval.t ;    
  }
 
  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_exchange(T* dest, typename std::enable_if<sizeof(T) != sizeof(int) && sizeof(T) != sizeof(int64_t), const T&>::type val) {
    return val;
  }

  KOKKOS_INLINE_FUNCTION
  int atomic_compare_exchange(int* dest, int compare, const int& val) {
    return hc::atomic_compare_exchange_int(dest, compare, val);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned int atomic_compare_exchange(unsigned int* dest, unsigned int compare, const unsigned int& val) {
    return hc::atomic_compare_exchange_unsigned(dest, compare, val);
  }

  KOKKOS_INLINE_FUNCTION
  int64_t atomic_compare_exchange(int64_t* dest, int64_t compare, const int64_t& val) {
    return (int64_t) hc::atomic_compare_exchange_uint64((uint64_t*)dest, (uint64_t)compare, (const uint64_t&)val);
  }

  KOKKOS_INLINE_FUNCTION
  uint64_t atomic_compare_exchange(uint64_t* dest, uint64_t compare, const uint64_t& val) {
    return hc::atomic_compare_exchange_uint64(dest, compare, val);
  }

  KOKKOS_INLINE_FUNCTION
  long long atomic_compare_exchange(long long* dest, long long compare, const long long& val) {
    return (long long)hc::atomic_compare_exchange_uint64((uint64_t*)(dest), (uint64_t)(compare), (const uint64_t&)(val));
  }

  KOKKOS_INLINE_FUNCTION
  float atomic_compare_exchange(float* dest, float compare, const float& val) {
    union U {
      int i ;
      float f ;
      KOKKOS_INLINE_FUNCTION U() {};
    } idest,icompare,ival;
    idest.f = *dest;
    icompare.f = compare;
    ival.f = val;
    idest.i = hc::atomic_compare_exchange_int(reinterpret_cast<int*>(dest), icompare.i, ival.i);
    return idest.f;
  }

  KOKKOS_INLINE_FUNCTION
  double atomic_compare_exchange(double* dest, double compare, const double& val) {
    union U {
      uint64_t i ;
      double d ;
      KOKKOS_INLINE_FUNCTION U() {};
    } idest,icompare,ival;
    idest.d = *dest;
    icompare.d = compare;
    ival.d = val;
    idest.i = hc::atomic_compare_exchange_uint64(reinterpret_cast<uint64_t*>(dest), icompare.i, ival.i);
    return idest.d;
  }

  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_compare_exchange(volatile T* dest, T compare, typename std::enable_if<sizeof(T) == sizeof(int), const T&>::type val) {
    union U {
      int i ;
      T f ;
      KOKKOS_INLINE_FUNCTION U() {};
    } idest,icompare,ival;
    idest.f = *dest;
    icompare.f = compare;
    ival.f = val;
    idest.i = hc::atomic_compare_exchange_int((int*)(dest), icompare.i, ival.i);
    return idest.f;
  }

  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_compare_exchange(volatile T* dest, T compare, typename std::enable_if<sizeof(T) == sizeof(int64_t), const T&>::type val) {
    union U {
      uint64_t i ;
      T f ;
      KOKKOS_INLINE_FUNCTION U() {};
    } idest,icompare,ival;
    idest.f = *dest;
    icompare.f = compare;
    ival.f = val;
    idest.i = hc::atomic_compare_exchange_uint64((uint64_t*)(dest), icompare.i, ival.i);
    return idest.f;
  }

  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_compare_exchange(volatile T* dest, T compare, typename std::enable_if<(sizeof(T) != sizeof(int32_t)) && (sizeof(T) != sizeof(int64_t)), const T&>::type val) {
    return val;
  }

  KOKKOS_INLINE_FUNCTION
  int atomic_fetch_add (volatile int * dest, const int& val) {
    return hc::atomic_fetch_add((int *)dest, val);
  }
  
  KOKKOS_INLINE_FUNCTION
  unsigned int atomic_fetch_add(unsigned int* dest, const unsigned int& val) {
    return hc::atomic_fetch_add(dest, val);
  }

  KOKKOS_INLINE_FUNCTION
  unsigned long atomic_fetch_add(volatile unsigned long* dest, const unsigned long& val) {
    return (unsigned long)hc::atomic_fetch_add((uint64_t *)dest, (const uint64_t)val);
  }

  KOKKOS_INLINE_FUNCTION
  int64_t atomic_fetch_add(volatile int64_t* dest, const int64_t& val) {
    return (int64_t)hc::atomic_fetch_add((uint64_t *)dest, (const uint64_t&)val);
  }

  KOKKOS_INLINE_FUNCTION
  char atomic_fetch_add(volatile char * dest, const char& val) {
    unsigned int oldval,newval,assume;
    oldval = *(int *)dest ;

    do {
      assume = oldval ;
      newval = assume&0x7fffff00 + ((assume&0xff)+val)&0xff ;
      oldval = hc::atomic_compare_exchange_unsigned((unsigned int*)dest, assume,newval);
    } while ( assume != oldval );

    return oldval ;    
  }


  KOKKOS_INLINE_FUNCTION
  short atomic_fetch_add(volatile short * dest, const short& val) {
    unsigned int oldval,newval,assume;
    oldval = *(int *)dest ;

    do {
      assume = oldval ;
      newval = assume&0x7fff0000 + ((assume&0xffff)+val)&0xffff ;
      oldval = hc::atomic_compare_exchange_unsigned((unsigned int*)dest, assume,newval);
    } while ( assume != oldval );

    return oldval ;    
  }

  KOKKOS_INLINE_FUNCTION
  long long atomic_fetch_add(volatile long long * dest, const long long& val) {
    return (long long)hc::atomic_fetch_add((uint64_t*)dest, (const uint64_t&)val);
  }



  KOKKOS_INLINE_FUNCTION
  int atomic_fetch_sub (volatile int * dest, const int& val) {
    return hc::atomic_fetch_sub((int *)dest, val);
  }
  
  KOKKOS_INLINE_FUNCTION
  unsigned int atomic_fetch_sub(volatile unsigned int* dest, const unsigned int& val) {
    return hc::atomic_fetch_sub((unsigned int *)dest, val);
  }

  KOKKOS_INLINE_FUNCTION
  int64_t atomic_fetch_sub(int64_t* dest, const int64_t& val) {
    return (int64_t)hc::atomic_fetch_add((uint64_t *)dest, -(const uint64_t&)val);
//    return (int64_t)hc::atomic_fetch_sub_uint64((uint64_t*)dest, (const uint64_t&)val);
  }
  
  KOKKOS_INLINE_FUNCTION
  char atomic_fetch_sub(volatile char * dest, const char& val) {
    unsigned int oldval,newval,assume;
    oldval = *(int *)dest ;

    do {
      assume = oldval ;
      newval = assume&0x7fffff00 + ((assume&0xff)-val)&0xff ;
      oldval = hc::atomic_compare_exchange_unsigned((unsigned int*)dest, assume,newval);
    } while ( assume != oldval );

    return oldval ;    
  }

  KOKKOS_INLINE_FUNCTION
  short atomic_fetch_sub(volatile short * dest, const short& val) {
    unsigned int oldval,newval,assume;
    oldval = *(int *)dest ;

    do {
      assume = oldval ;
      newval = assume&0x7fff0000 + ((assume&0xffff)-val)&0xffff;
      oldval = hc::atomic_compare_exchange_unsigned((unsigned int*)dest, assume,newval);
    } while ( assume != oldval );

    return oldval ;    
  }

  KOKKOS_INLINE_FUNCTION
  long long atomic_fetch_sub(volatile long long * dest, const long long& val) {
    return (long long)hc::atomic_fetch_add((uint64_t*)dest, -(const uint64_t&)val);
  }

  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_fetch_add(volatile T* dest, typename std::enable_if<sizeof(T) == sizeof(int), const T&>::type val) {
    union U {
      unsigned int i ;
      T t ;
      KOKKOS_INLINE_FUNCTION U() {};
    } assume , oldval , newval ;

    oldval.t = *dest ;

    do {
      assume.i = oldval.i ;
      newval.t = assume.t + val ;
      oldval.i = atomic_compare_exchange( (unsigned int*)(dest) , assume.i , newval.i );
    } while ( assume.i != oldval.i );

    return oldval.t ;    
  }

  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_fetch_add(volatile T* dest, typename std::enable_if<sizeof(T) != sizeof(int) && sizeof(T) == sizeof(int64_t), const T&>::type val) {
    union U {
      uint64_t i ;
      T t ;
      KOKKOS_INLINE_FUNCTION U() {};
    } assume , oldval , newval ;

    oldval.t = *dest ;

    do {
      assume.i = oldval.i ;
      newval.t = assume.t + val ;
      oldval.i = atomic_compare_exchange( (uint64_t*)dest , assume.i , newval.i );
    } while ( assume.i != oldval.i );

    return oldval.t ;    
  }


  //WORKAROUND
  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_fetch_add(volatile T* dest, typename std::enable_if<sizeof(T) != sizeof(int) && sizeof(T) != sizeof(int64_t), const T&>::type val) {
    return val ;
  }

  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_fetch_sub(volatile T* dest, typename std::enable_if<sizeof(T) == sizeof(int),T>::type & val) {
    union U {
      int i ;
      T t ;
      KOKKOS_INLINE_FUNCTION U() {};
    } assume , oldval , newval ;

    oldval.t = *dest ;

    do {
      assume.i = oldval.i ;
      newval.t = assume.t - val ;
      oldval.i = Kokkos::atomic_compare_exchange( (int*)dest , assume.i , newval.i );
    } while ( assume.i != oldval.i );

    return oldval.t ;
  }

  template<class T>
  KOKKOS_INLINE_FUNCTION
  T atomic_fetch_sub(volatile T* dest, typename std::enable_if<sizeof(T) != sizeof(int) && sizeof(T) == sizeof(int64_t), const T&>::type val) {
    union U {
      int64_t i ;
      T t ;
      KOKKOS_INLINE_FUNCTION U() {};
    } assume , oldval , newval ;

    oldval.t = *dest ;

    do {
      assume.i = oldval.i ;
      newval.t = assume.t - val ;
      oldval.i = atomic_compare_exchange( (int64_t*)dest , assume.i , newval.i );
    } while ( assume.i != oldval.i );

    return oldval.t ;    
  }
}
#endif
