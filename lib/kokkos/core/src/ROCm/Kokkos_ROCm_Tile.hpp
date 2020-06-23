/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <hc.hpp>
#include <type_traits>
#include <vector>
#include <memory>
#include <ROCm/Kokkos_ROCm_Config.hpp>

#if !defined(KOKKOS_ROCM_TILE_H)
#define KOKKOS_ROCM_TILE_H

// Macro to abstract out the enable_if craziness
#define KOKKOS_ROCM_REQUIRES(...)                                    \
  bool KokkosROCmRequiresBool##__LINE__ = true,                      \
       typename std::enable_if < KokkosROCmRequiresBool##__LINE__ && \
           (__VA_ARGS__),                                            \
       int > ::type = 0

// This number uniquely identifies the 1.5 release build.
#if __hcc_workweek__ > 17160
#define ROCM15 1
#endif

namespace Kokkos {
namespace Impl {

template <class T>

#if defined(ROCM15)
using lds_t = T;
#else
// prior to 1.5, needed to decorate LDS addresses
using lds_t = __attribute__((address_space(3))) T;
#endif

#define KOKKOS_ROCM_TILE_RESTRIC_CPU restrict(cpu, amp)

// a set of routines to the replace the std::routines
// that will operate on address space 3 types

#if defined(ROCM15)
// 1.5 can't use std::copy et al for LDS access, so we define our own
// set of routines
template <class I, class O>
void rcopy(I first, I last, O out) [[hc]] {
  while (first != last) *out++ = *first++;
}
template <class I, class F>
void rfor_each(I first, I last, F f) [[hc]] {
  for (; first != last; ++first) f(*first);
}

template <class I, class O, class F>
void rtransform(I first, I last, O out, F f) [[hc]] {
  while (first != last) *out++ = f(*first++);
}
#endif

inline std::size_t get_max_tile_size() KOKKOS_ROCM_TILE_RESTRIC_CPU {
  return hc::accelerator().get_max_tile_static_size() - 1024;
}

inline std::size_t get_max_tile_thread() KOKKOS_ROCM_TILE_RESTRIC_CPU {
  return 64;
}

inline int next_pow_2(int x) restrict(cpu, amp) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return x + 1;
}

template <class T>
inline std::size_t get_tile_size(std::size_t n = 1, std::size_t team = 64,
                                 std::size_t vector = 1)
    KOKKOS_ROCM_TILE_RESTRIC_CPU {
  return team * vector;
  /*
    const auto size = sizeof(T) * n;
    const auto group_size = get_max_tile_size();
    if (size == 0 || size > group_size) return 0;
    // Assume that thread size is a power of 2
    auto thread_size = std::min(team*vector,4*get_max_tile_thread());
    // ensure that we have enough tile static memory to keep
    // threadsize * size elements for reductions
    while(size > (group_size / thread_size) && thread_size > 2)
    { thread_size /= 2;}
    return thread_size;
  */
}

template <class T>
struct array_view {
  T* x;
  std::size_t n;

  array_view(T* xp, std::size_t np) [[hc]] [[cpu]] : x(xp), n(np) {}

  array_view(T* xp, T* yp) [[hc]] [[cpu]] : x(xp), n(yp - xp) {}

  T& operator[](std::size_t i) const [[hc]] [[cpu]] { return x[i]; }

  std::size_t size() const [[hc]] [[cpu]] { return this->n; }

  T* data() const [[hc]] [[cpu]] { return x; }

  T* begin() const [[hc]] [[cpu]] { return x; }

  T* end() const [[hc]] [[cpu]] { return x + this->size(); }
};

template <class T>
struct rocm_char {
  using type = char;
};

template <class T>
struct rocm_char<const T> : std::add_const<typename rocm_char<T>::type> {};
#if !defined(ROCM15)
// earlier compilers required explicit address space decorations
template <class T>
struct rocm_char<__attribute__((address_space(3))) T> {
  using type = __attribute__((address_space(3))) typename rocm_char<T>::type;
};

template <class T>
struct rocm_char<const __attribute__((address_space(3))) T> {
  using type =
      const __attribute__((address_space(3))) typename rocm_char<T>::type;
};
#endif

template <class T, class Char = typename rocm_char<T>::type>
Char* rocm_byte_cast(T& x) restrict(cpu, amp) {
  return reinterpret_cast<Char*>(&x);
}

template <class T, class U>
void rocm_raw_assign(T& x, const U& y) restrict(cpu, amp) {
  auto* src  = rocm_byte_cast(y);
  auto* dest = rocm_byte_cast(x);
#if defined(ROCM15)
  rcopy(src, src + sizeof(T), dest);
#else
  std::copy(src, src + sizeof(T), dest);
#endif
}

template <class T, class U>
void rocm_assign_impl(T& x, const U& y, std::true_type) restrict(cpu, amp) {
  rocm_raw_assign(x, y);
}

template <class T, class U>
void rocm_assign_impl(T& x, const U& y, std::false_type) restrict(cpu, amp) {
  x = y;
}

// Workaround for assigning in and out of LDS memory
template <class T, class U>
void rocm_assign(T& x, const U& y) restrict(cpu, amp) {
  rocm_assign_impl(x, y,
                   std::integral_constant<bool, (sizeof(T) == sizeof(U))>());
}

// Compute the address space of tile
template <class T>
struct tile_type {
#if defined(ROCM15)
  typedef T type;
#else
  typedef __attribute__((address_space(3))) T type;
#endif
};

#if !defined(ROCM15)
template <class T, class Body>
void lds_for(__attribute__((address_space(3))) T& value, Body b) [[hc]] {
  T state = value;
  b(state);
  value = state;
}
#endif

template <class T, class Body>
void lds_for(T& value, Body b) [[hc]] {
  b(value);
}

constexpr std::size_t get_max_tile_array_size() { return 24; }

template <class Derived, class T>
struct single_action {
  template <class Action>
  void action_at(std::size_t i, Action a) [[hc]] {
    auto& value = static_cast<Derived&>(*this)[i];
#ifdef KOKKOS_IMPL_ROCM_CLANG_WORKAROUND
    T state = value;
    a(state);
    value = state;
#else
    a(value);
#endif
  }

  template <class Action>
  void action_at(std::size_t i, std::size_t j, Action a) [[hc]] {
    static_cast<Derived&>(*this).action_at(i, [&](T& x) {
      static_cast<Derived&>(*this).action_at(j, [&](T& y) { a(x, y); });
    });
  }
};

template <class T>
struct tile_buffer : array_view<typename tile_type<T>::type>,
                     single_action<tile_buffer<T>, T> {
  typedef typename tile_type<T>::type element_type;
  typedef array_view<element_type> base;

  using base::base;

  tile_buffer(element_type* xp, std::size_t np, std::size_t) [[hc]] [[cpu]]
  : base(xp, np) {}

  tile_buffer(T* xp, T* yp, std::size_t) [[hc]] [[cpu]] : base(xp, yp) {}
};

template <class T>
struct tile_buffer<T[]> {
  typedef typename tile_type<T>::type element_type;
  typedef typename tile_type<char>::type tchar_type;
  element_type* element_data;
  std::size_t n, m;

  tile_buffer(element_type* xp, std::size_t np, std::size_t mp) [[hc]] [[cpu]]
  : element_data(xp),
    n(np),
    m(mp) {}

  tile_buffer(element_type* xp, element_type* yp, std::size_t mp) [[hc]] [[cpu]]
  : element_data(xp),
    n(yp - xp),
    m(mp) {}

  element_type* operator[](std::size_t i) const [[hc]] [[cpu]] {
    return element_data + i * m;
  }

  template <class Action, class Q = T>
  typename std::enable_if<(sizeof(Q) <= 8), void>::type action_at(std::size_t i,
                                                                  Action a)
      [[hc]] {
    element_type* value = (*this)[i];
#if defined(ROCM15)
    a(value);
#else
#ifdef KOKKOS_IMPL_ROCM_CLANG_WORKAROUND
    if (m > get_max_tile_array_size()) return;
    T state[get_max_tile_array_size()];
    // std::copy(value, value+m, state);
    // Workaround for assigning from LDS memory
    std::transform(value, value + m, state, [](element_type& x) {
      T result;
      rocm_assign(result, x);
      return result;
    });
    a(state);
    std::copy(state, state + m, value);
#endif
#endif
  }

  template <class Action, class Q = T>
  typename std::enable_if<!(sizeof(Q) <= 8), void>::type action_at(
      std::size_t i, Action a) [[hc]] {
    element_type* value = (*this)[i];
#if defined(ROCM15)
    a(value);
#else
    if (m > get_max_tile_array_size()) return;
    T state[get_max_tile_array_size()];
    // std::copy(value, value+m, state);
    // Workaround for assigning from LDS memory
    std::transform(value, value + m, state, [](element_type& x) {
      T result;
      rocm_assign(result, x);
      return result;
    });
    a(state);
    // this workaround required when T is greater than 8 bytes
    tile_static char tv[64 * sizeof(T)];
    size_t sT = sizeof(T);
    for (int j = 0; j < sT; j++) tv[i * sT + j] = ((char*)state)[j];
    for (int j = 0; j < sT; j++) ((tchar_type*)value)[j] = tv[i * sT + j];
#endif
  }

  template <class Action>
  void action_at(std::size_t i, std::size_t j, Action a) [[hc]] {
    this->action_at(i,
                    [&](T* x) { this->action_at(j, [&](T* y) { a(x, y); }); });
  }

  std::size_t size() const [[hc]] [[cpu]] { return this->n; }

  element_type* data() const [[hc]] [[cpu]] { return element_data; }
};

// Zero initialize LDS memory
struct zero_init_f {
  template <class T>
#if defined(ROCM15)
  void operator()(T& x, std::size_t = 1) const [[hc]] {
    auto* start = reinterpret_cast<char*>(&x);
    for (int i = 0; i < sizeof(T); i++) start[i] = 0;
    rocm_raw_assign(x, T());
  }
#else
  void operator()(__attribute__((address_space(3))) T& x, std::size_t = 1) const
      [[hc]] {
    auto* start = reinterpret_cast<__attribute__((address_space(3))) char*>(&x);
    std::fill(start, start + sizeof(T), 0);
    rocm_raw_assign(x, T());
  }
#endif

  template <class T>
#if defined(ROCM15)
  void operator()(T* x, std::size_t size) const [[hc]] {
    rfor_each(x, x + size, *this);
  }
#else
  void operator()(__attribute__((address_space(3))) T* x,
                  std::size_t size) const [[hc]] {
    std::for_each(x, x + size, *this);
  }
#endif
};

static constexpr zero_init_f zero_init = {};

struct tile_desc {
  // Number of work items, or size of extent
  std::size_t elements;
  // number of threads in team
  std::size_t team_size;
  // vector length of team
  std::size_t vector_length;
  // Size of tile
  std::size_t tile_size;
  // Size of array
  std::size_t array_size;
  // Number of tiles
  std::size_t num_tiles;
  // Per team reserved LDS memory, used for reduction
  std::size_t reduce_size;
  // Per team shared memory in LDS, this in addition to reduce shared mem
  std::size_t shared_size;
  std::size_t size;
};

template <class T>
tile_desc get_tile_desc(std::size_t size, std::size_t array_size = 1,
                        std::size_t team_size = 64, std::size_t vector_size = 1,
                        std::size_t shared_size = 0) {
  tile_desc result;
  result.elements      = size;
  result.array_size    = array_size;
  result.vector_length = vector_size;
  result.team_size     = team_size;
  result.tile_size     = get_tile_size<T>(array_size, team_size, vector_size);
  result.num_tiles     = std::ceil(1.0 * size / result.tile_size);
  result.reduce_size   = result.tile_size * sizeof(T) * array_size;
  result.shared_size   = shared_size;
  result.size          = result.tile_size * result.num_tiles;

  return result;
}

template <class U, class F, class T = typename std::remove_extent<U>::type>
hc::completion_future tile_for(tile_desc td, F f) {
  assert(td.array_size <= get_max_tile_array_size() && "Exceed max array size");
  assert(((td.size % td.tile_size) == 0) &&
         "Tile size must be divisible by extent");
  auto grid = hc::extent<1>(td.size).tile_with_dynamic(
      td.tile_size, td.reduce_size + td.shared_size);
  // grid.set_dynamic_group_segment_size(td.reduce_size + td.shared_size);
  return parallel_for_each(
      grid, [=](hc::tiled_index<1> t_idx) [[hc]] {
#if defined(ROCM15)
        typedef T group_t;
#else
        typedef __attribute__((address_space(3))) T group_t;
#endif
        group_t* buffer =
            (group_t*)hc::get_dynamic_group_segment_base_pointer();
        tile_buffer<U> tb(buffer, td.tile_size, td.array_size);
        zero_init(tb[t_idx.local[0]], td.array_size);
        f(t_idx, tb);
      });
}

}  // namespace Impl
}  // namespace Kokkos

#endif
