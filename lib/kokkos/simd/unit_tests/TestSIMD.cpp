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

#include <gtest/gtest.h>

#include <Kokkos_SIMD.hpp>

class gtest_checker {
 public:
  void truth(bool x) const { EXPECT_TRUE(x); }
  template <class T>
  void equality(T const& a, T const& b) const {
    EXPECT_EQ(a, b);
  }
};

class kokkos_checker {
 public:
  KOKKOS_INLINE_FUNCTION void truth(bool x) const {
    if (!x) Kokkos::abort("SIMD unit test truth condition failed on device");
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION void equality(T const& a, T const& b) const {
    if (a != b)
      Kokkos::abort("SIMD unit test equality condition failed on device");
  }
};

template <class T, class Abi>
inline void host_check_equality(
    Kokkos::Experimental::simd<T, Abi> const& expected_result,
    Kokkos::Experimental::simd<T, Abi> const& computed_result,
    std::size_t nlanes) {
  gtest_checker checker;
  for (std::size_t i = 0; i < nlanes; ++i) {
    checker.equality(expected_result[i], computed_result[i]);
  }
  using mask_type = typename Kokkos::Experimental::simd<T, Abi>::mask_type;
  mask_type mask(false);
  for (std::size_t i = 0; i < nlanes; ++i) {
    mask[i] = true;
  }
  checker.equality((expected_result == computed_result) && mask, mask);
}

template <class T, class Abi>
KOKKOS_INLINE_FUNCTION void device_check_equality(
    Kokkos::Experimental::simd<T, Abi> const& expected_result,
    Kokkos::Experimental::simd<T, Abi> const& computed_result,
    std::size_t nlanes) {
  kokkos_checker checker;
  for (std::size_t i = 0; i < nlanes; ++i) {
    checker.equality(expected_result[i], computed_result[i]);
  }
  using mask_type = typename Kokkos::Experimental::simd<T, Abi>::mask_type;
  mask_type mask(false);
  for (std::size_t i = 0; i < nlanes; ++i) {
    mask[i] = true;
  }
  checker.equality((expected_result == computed_result) && mask, mask);
}

class load_element_aligned {
 public:
  template <class T, class Abi>
  bool host_load(T const* mem, std::size_t n,
                 Kokkos::Experimental::simd<T, Abi>& result) const {
    if (n < result.size()) return false;
    result.copy_from(mem, Kokkos::Experimental::element_aligned_tag());
    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::simd<T, Abi>& result) const {
    if (n < result.size()) return false;
    result.copy_from(mem, Kokkos::Experimental::element_aligned_tag());
    return true;
  }
};

class load_masked {
 public:
  template <class T, class Abi>
  bool host_load(T const* mem, std::size_t n,
                 Kokkos::Experimental::simd<T, Abi>& result) const {
    using mask_type = typename Kokkos::Experimental::simd<T, Abi>::mask_type;
    mask_type mask(false);
    for (std::size_t i = 0; i < n; ++i) {
      mask[i] = true;
    }
    where(mask, result)
        .copy_from(mem, Kokkos::Experimental::element_aligned_tag());
    where(!mask, result) = 0;
    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::simd<T, Abi>& result) const {
    using mask_type = typename Kokkos::Experimental::simd<T, Abi>::mask_type;
    mask_type mask(false);
    for (std::size_t i = 0; i < n; ++i) {
      mask[i] = true;
    }
    where(mask, result)
        .copy_from(mem, Kokkos::Experimental::element_aligned_tag());
    where(!mask, result) = T(0);
    return true;
  }
};

class load_as_scalars {
 public:
  template <class T, class Abi>
  bool host_load(T const* mem, std::size_t n,
                 Kokkos::Experimental::simd<T, Abi>& result) const {
    for (std::size_t i = 0; i < n; ++i) {
      result[i] = mem[i];
    }
    for (std::size_t i = n; i < result.size(); ++i) {
      result[i] = T(0);
    }
    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::simd<T, Abi>& result) const {
    for (std::size_t i = 0; i < n; ++i) {
      result[i] = mem[i];
    }
    for (std::size_t i = n; i < result.size(); ++i) {
      result[i] = T(0);
    }
    return true;
  }
};

template <class Abi, class Loader, class BinaryOp, class T>
void host_check_binary_op_one_loader(BinaryOp binary_op, std::size_t n,
                                     T const* first_args,
                                     T const* second_args) {
  Loader loader;
  using simd_type             = Kokkos::Experimental::simd<T, Abi>;
  std::size_t constexpr width = simd_type::size();
  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
    simd_type first_arg;
    bool const loaded_first_arg =
        loader.host_load(first_args + i, nlanes, first_arg);
    simd_type second_arg;
    bool const loaded_second_arg =
        loader.host_load(second_args + i, nlanes, second_arg);
    if (!(loaded_first_arg && loaded_second_arg)) continue;
    simd_type expected_result;
    for (std::size_t lane = 0; lane < nlanes; ++lane) {
      expected_result[lane] =
          binary_op.on_host(first_arg[lane], second_arg[lane]);
    }
    simd_type const computed_result = binary_op.on_host(first_arg, second_arg);
    host_check_equality(expected_result, computed_result, nlanes);
  }
}

template <class Abi, class Loader, class BinaryOp, class T>
KOKKOS_INLINE_FUNCTION void device_check_binary_op_one_loader(
    BinaryOp binary_op, std::size_t n, T const* first_args,
    T const* second_args) {
  Loader loader;
  using simd_type             = Kokkos::Experimental::simd<T, Abi>;
  std::size_t constexpr width = simd_type::size();
  for (std::size_t i = 0; i < n; i += width) {
    std::size_t const nremaining = n - i;
    std::size_t const nlanes     = Kokkos::min(nremaining, width);
    simd_type first_arg;
    bool const loaded_first_arg =
        loader.device_load(first_args + i, nlanes, first_arg);
    simd_type second_arg;
    bool const loaded_second_arg =
        loader.device_load(second_args + i, nlanes, second_arg);
    if (!(loaded_first_arg && loaded_second_arg)) continue;
    simd_type expected_result;
    for (std::size_t lane = 0; lane < nlanes; ++lane) {
      expected_result[lane] =
          binary_op.on_device(first_arg[lane], second_arg[lane]);
    }
    simd_type const computed_result =
        binary_op.on_device(first_arg, second_arg);
    device_check_equality(expected_result, computed_result, nlanes);
  }
}

template <class Abi, class BinaryOp, class T>
inline void host_check_binary_op_all_loaders(BinaryOp binary_op, std::size_t n,
                                             T const* first_args,
                                             T const* second_args) {
  host_check_binary_op_one_loader<Abi, load_element_aligned>(
      binary_op, n, first_args, second_args);
  host_check_binary_op_one_loader<Abi, load_masked>(binary_op, n, first_args,
                                                    second_args);
  host_check_binary_op_one_loader<Abi, load_as_scalars>(
      binary_op, n, first_args, second_args);
}

template <class Abi, class BinaryOp, class T>
KOKKOS_INLINE_FUNCTION void device_check_binary_op_all_loaders(
    BinaryOp binary_op, std::size_t n, T const* first_args,
    T const* second_args) {
  device_check_binary_op_one_loader<Abi, load_element_aligned>(
      binary_op, n, first_args, second_args);
  device_check_binary_op_one_loader<Abi, load_masked>(binary_op, n, first_args,
                                                      second_args);
  device_check_binary_op_one_loader<Abi, load_as_scalars>(
      binary_op, n, first_args, second_args);
}

class plus {
 public:
  template <class T>
  auto on_host(T const& a, T const& b) const {
    return a + b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return a + b;
  }
};

class minus {
 public:
  template <class T>
  auto on_host(T const& a, T const& b) const {
    return a - b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return a - b;
  }
};

class multiplies {
 public:
  template <class T>
  auto on_host(T const& a, T const& b) const {
    return a * b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return a * b;
  }
};

class divides {
 public:
  template <class T>
  auto on_host(T const& a, T const& b) const {
    return a / b;
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION auto on_device(T const& a, T const& b) const {
    return a / b;
  }
};

template <class Abi>
inline void host_check_math_ops() {
  std::size_t constexpr n     = 11;
  double const first_args[n]  = {1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2};
  double const second_args[n] = {1, 2, 1, 1, 1, -3, -2, 1, 13, -3, -2};
  host_check_binary_op_all_loaders<Abi>(plus(), n, first_args, second_args);
  host_check_binary_op_all_loaders<Abi>(minus(), n, first_args, second_args);
  host_check_binary_op_all_loaders<Abi>(multiplies(), n, first_args,
                                        second_args);
  host_check_binary_op_all_loaders<Abi>(divides(), n, first_args, second_args);
}

template <class Abi>
KOKKOS_INLINE_FUNCTION void device_check_math_ops() {
  std::size_t constexpr n     = 11;
  double const first_args[n]  = {1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2};
  double const second_args[n] = {1, 2, 1, 1, 1, -3, -2, 1, 13, -3, -2};
  device_check_binary_op_all_loaders<Abi>(plus(), n, first_args, second_args);
  device_check_binary_op_all_loaders<Abi>(minus(), n, first_args, second_args);
  device_check_binary_op_all_loaders<Abi>(multiplies(), n, first_args,
                                          second_args);
  device_check_binary_op_all_loaders<Abi>(divides(), n, first_args,
                                          second_args);
}

template <class Abi>
inline void host_check_abi() {
  host_check_math_ops<Abi>();
}

template <class Abi>
KOKKOS_INLINE_FUNCTION void device_check_abi() {
  device_check_math_ops<Abi>();
}

inline void host_check_abis(Kokkos::Experimental::Impl::abi_set<>) {}

KOKKOS_INLINE_FUNCTION void device_check_abis(
    Kokkos::Experimental::Impl::abi_set<>) {}

template <class FirstAbi, class... RestAbis>
inline void host_check_abis(
    Kokkos::Experimental::Impl::abi_set<FirstAbi, RestAbis...>) {
  host_check_abi<FirstAbi>();
  host_check_abis(Kokkos::Experimental::Impl::abi_set<RestAbis...>());
}

template <class FirstAbi, class... RestAbis>
KOKKOS_INLINE_FUNCTION void device_check_abis(
    Kokkos::Experimental::Impl::abi_set<FirstAbi, RestAbis...>) {
  device_check_abi<FirstAbi>();
  device_check_abis(Kokkos::Experimental::Impl::abi_set<RestAbis...>());
}

TEST(simd, host) {
  host_check_abis(Kokkos::Experimental::Impl::host_abi_set());
}

class simd_device_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_abis(Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, device) {
  Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::IndexType<int>>(0, 1),
                       simd_device_functor());
}
