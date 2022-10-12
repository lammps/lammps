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

#ifndef KOKKOS_VECTOR_HPP
#define KOKKOS_VECTOR_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_VECTOR
#endif

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_DualView.hpp>

/* Drop in replacement for std::vector based on Kokkos::DualView
 * Most functions only work on the host (it will not compile if called from
 * device kernel)
 *
 */
namespace Kokkos {

template <class Scalar, class Arg1Type = void>
class vector : public DualView<Scalar*, LayoutLeft, Arg1Type> {
 public:
  using value_type      = Scalar;
  using pointer         = Scalar*;
  using const_pointer   = const Scalar*;
  using reference       = Scalar&;
  using const_reference = const Scalar&;
  using iterator        = Scalar*;
  using const_iterator  = const Scalar*;
  using size_type       = size_t;

 private:
  size_t _size;
  float _extra_storage;
  using DV = DualView<Scalar*, LayoutLeft, Arg1Type>;

 public:
#ifdef KOKKOS_ENABLE_CUDA_UVM
  KOKKOS_INLINE_FUNCTION reference operator()(int i) const {
    return DV::h_view(i);
  };
  KOKKOS_INLINE_FUNCTION reference operator[](int i) const {
    return DV::h_view(i);
  };
#else
  inline reference operator()(int i) const { return DV::h_view(i); };
  inline reference operator[](int i) const { return DV::h_view(i); };
#endif

  /* Member functions which behave like std::vector functions */

  vector() : DV() {
    _size          = 0;
    _extra_storage = 1.1;
  }

  vector(int n, Scalar val = Scalar())
      : DualView<Scalar*, LayoutLeft, Arg1Type>("Vector", size_t(n * (1.1))) {
    _size                 = n;
    _extra_storage        = 1.1;
    DV::modified_flags(0) = 1;

    assign(n, val);
  }

  void resize(size_t n) {
    if (n >= span()) DV::resize(size_t(n * _extra_storage));
    _size = n;
  }

  void resize(size_t n, const Scalar& val) { assign(n, val); }

  void assign(size_t n, const Scalar& val) {
    /* Resize if necessary (behavior of std:vector) */

    if (n > span()) DV::resize(size_t(n * _extra_storage));
    _size = n;

    /* Assign value either on host or on device */

    if (DV::template need_sync<typename DV::t_dev::device_type>()) {
      set_functor_host f(DV::h_view, val);
      parallel_for("Kokkos::vector::assign", n, f);
      typename DV::t_host::execution_space().fence(
          "Kokkos::vector::assign: fence after assigning values");
      DV::template modify<typename DV::t_host::device_type>();
    } else {
      set_functor f(DV::d_view, val);
      parallel_for("Kokkos::vector::assign", n, f);
      typename DV::t_dev::execution_space().fence(
          "Kokkos::vector::assign: fence after assigning values");
      DV::template modify<typename DV::t_dev::device_type>();
    }
  }

  void reserve(size_t n) { DV::resize(size_t(n * _extra_storage)); }

  void push_back(Scalar val) {
    DV::template sync<typename DV::t_host::device_type>();
    DV::template modify<typename DV::t_host::device_type>();
    if (_size == span()) {
      size_t new_size = _size * _extra_storage;
      if (new_size == _size) new_size++;
      DV::resize(new_size);
    }

    DV::h_view(_size) = val;
    _size++;
  }

  void pop_back() { _size--; }

  void clear() { _size = 0; }

  iterator insert(iterator it, const value_type& val) {
    return insert(it, 1, val);
  }

  iterator insert(iterator it, size_type count, const value_type& val) {
    if ((size() == 0) && (it == begin())) {
      resize(count, val);
      DV::sync_host();
      return begin();
    }
    DV::sync_host();
    DV::modify_host();
    if (std::less<>()(it, begin()) || std::less<>()(end(), it))
      Kokkos::abort("Kokkos::vector::insert : invalid insert iterator");
    if (count == 0) return it;
    ptrdiff_t start = std::distance(begin(), it);
    auto org_size   = size();
    resize(size() + count);

    std::copy_backward(begin() + start, begin() + org_size,
                       begin() + org_size + count);
    std::fill_n(begin() + start, count, val);

    return begin() + start;
  }

 private:
  template <class T>
  struct impl_is_input_iterator
      : /* TODO replace this */ std::integral_constant<
            bool, !std::is_convertible<T, size_type>::value> {};

 public:
  // TODO: can use detection idiom to generate better error message here later
  template <typename InputIterator>
  std::enable_if_t<impl_is_input_iterator<InputIterator>::value, iterator>
  insert(iterator it, InputIterator b, InputIterator e) {
    ptrdiff_t count = std::distance(b, e);

    DV::sync_host();
    DV::modify_host();
    if (std::less<>()(it, begin()) || std::less<>()(end(), it))
      Kokkos::abort("Kokkos::vector::insert : invalid insert iterator");

    ptrdiff_t start = std::distance(begin(), it);
    auto org_size   = size();

    // Note: resize(...) invalidates it; use begin() + start instead
    resize(size() + count);

    std::copy_backward(begin() + start, begin() + org_size,
                       begin() + org_size + count);
    std::copy(b, e, begin() + start);

    return begin() + start;
  }

  KOKKOS_INLINE_FUNCTION constexpr bool is_allocated() const {
    return DV::is_allocated();
  }

  size_type size() const { return _size; }
  size_type max_size() const { return 2000000000; }
  size_type span() const { return DV::span(); }
  bool empty() const { return _size == 0; }

  pointer data() const { return DV::h_view.data(); }

  iterator begin() const { return DV::h_view.data(); }

  iterator end() const {
    return _size > 0 ? DV::h_view.data() + _size : DV::h_view.data();
  }

  reference front() { return DV::h_view(0); }

  reference back() { return DV::h_view(_size - 1); }

  const_reference front() const { return DV::h_view(0); }

  const_reference back() const { return DV::h_view(_size - 1); }

  /* std::algorithms which work originally with iterators, here they are
   * implemented as member functions */

  size_t lower_bound(const size_t& start, const size_t& theEnd,
                     const Scalar& comp_val) const {
    int lower = start;  // FIXME (mfh 24 Apr 2014) narrowing conversion
    int upper =
        _size > theEnd
            ? theEnd
            : _size - 1;  // FIXME (mfh 24 Apr 2014) narrowing conversion
    if (upper <= lower) {
      return theEnd;
    }

    Scalar lower_val = DV::h_view(lower);
    Scalar upper_val = DV::h_view(upper);
    size_t idx       = (upper + lower) / 2;
    Scalar val       = DV::h_view(idx);
    if (val > upper_val) return upper;
    if (val < lower_val) return start;

    while (upper > lower) {
      if (comp_val > val) {
        lower = ++idx;
      } else {
        upper = idx;
      }
      idx = (upper + lower) / 2;
      val = DV::h_view(idx);
    }
    return idx;
  }

  bool is_sorted() {
    for (int i = 0; i < _size - 1; i++) {
      if (DV::h_view(i) > DV::h_view(i + 1)) return false;
    }
    return true;
  }

  iterator find(Scalar val) const {
    if (_size == 0) return end();

    int upper, lower, current;
    current = _size / 2;
    upper   = _size - 1;
    lower   = 0;

    if ((val < DV::h_view(0)) || (val > DV::h_view(_size - 1))) return end();

    while (upper > lower) {
      if (val > DV::h_view(current))
        lower = current + 1;
      else
        upper = current;
      current = (upper + lower) / 2;
    }

    if (val == DV::h_view(current))
      return &DV::h_view(current);
    else
      return end();
  }

  /* Additional functions for data management */

  void device_to_host() { deep_copy(DV::h_view, DV::d_view); }
  void host_to_device() const { deep_copy(DV::d_view, DV::h_view); }

  void on_host() { DV::template modify<typename DV::t_host::device_type>(); }
  void on_device() { DV::template modify<typename DV::t_dev::device_type>(); }

  void set_overallocation(float extra) { _extra_storage = 1.0 + extra; }

 public:
  struct set_functor {
    using execution_space = typename DV::t_dev::execution_space;
    typename DV::t_dev _data;
    Scalar _val;

    set_functor(typename DV::t_dev data, Scalar val) : _data(data), _val(val) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i) const { _data(i) = _val; }
  };

  struct set_functor_host {
    using execution_space = typename DV::t_host::execution_space;
    typename DV::t_host _data;
    Scalar _val;

    set_functor_host(typename DV::t_host data, Scalar val)
        : _data(data), _val(val) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const int& i) const { _data(i) = _val; }
  };
};

}  // namespace Kokkos
#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_VECTOR
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_VECTOR
#endif
#endif
