//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

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
    // gcc 8.4.0 warns if using nlanes as upper bound about first_arg and/or
    // second_arg being uninitialized
    for (std::size_t lane = 0; lane < simd_type::size(); ++lane) {
      if (lane < nlanes)
        expected_result[lane] =
            binary_op.on_host(T(first_arg[lane]), T(second_arg[lane]));
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

template <typename Abi, typename DataType, size_t n>
inline void host_check_all_math_ops(const DataType (&first_args)[n],
                                    const DataType (&second_args)[n]) {
  host_check_binary_op_all_loaders<Abi>(plus(), n, first_args, second_args);
  host_check_binary_op_all_loaders<Abi>(minus(), n, first_args, second_args);
  host_check_binary_op_all_loaders<Abi>(multiplies(), n, first_args,
                                        second_args);

  // TODO: Place fallback division implementations for all simd integer types
  if constexpr (std::is_same_v<DataType, double>)
    host_check_binary_op_all_loaders<Abi>(divides(), n, first_args,
                                          second_args);
}

template <typename Abi, typename DataType>
inline void host_check_abi_size() {
  using simd_type = Kokkos::Experimental::simd<DataType, Abi>;
  using mask_type = typename simd_type::mask_type;
  static_assert(simd_type::size() == mask_type::size());
}

template <class Abi, typename DataType>
inline void host_check_math_ops() {
  constexpr size_t n = 11;

  host_check_abi_size<Abi, DataType>();

  if constexpr (std::is_signed_v<DataType>) {
    DataType const first_args[n]  = {1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2};
    DataType const second_args[n] = {1, 2, 1, 1, 1, -3, -2, 1, 13, -3, -2};
    host_check_all_math_ops<Abi>(first_args, second_args);
  } else {
    DataType const first_args[n]  = {1, 2, 1, 10, 0, 1, 2, 10, 0, 1, 2};
    DataType const second_args[n] = {1, 2, 1, 1, 1, 3, 2, 1, 13, 3, 2};
    host_check_all_math_ops<Abi>(first_args, second_args);
  }
}

template <class Abi>
inline void host_check_mask_ops() {
  using mask_type = Kokkos::Experimental::simd_mask<double, Abi>;
  EXPECT_FALSE(none_of(mask_type(true)));
  EXPECT_TRUE(none_of(mask_type(false)));
  EXPECT_TRUE(all_of(mask_type(true)));
  EXPECT_FALSE(all_of(mask_type(false)));
}

template <class Abi>
inline void host_check_conversions() {
  {
    auto a = Kokkos::Experimental::simd<std::uint64_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::int64_t, Abi>(a);
    EXPECT_TRUE(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd<std::int32_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::uint64_t, Abi>(a);
    EXPECT_TRUE(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd<std::uint64_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::int32_t, Abi>(a);
    EXPECT_TRUE(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<double, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(a);
    EXPECT_TRUE(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::uint64_t, Abi>(a);
    EXPECT_TRUE(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::int64_t, Abi>(a);
    EXPECT_TRUE(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<double, Abi>(a);
    EXPECT_TRUE(b == decltype(b)(true));
  }
}

template <class Abi>
inline void host_check_shifts() {
  auto a = Kokkos::Experimental::simd<std::uint64_t, Abi>(8);
  auto b = a >> 1;
  EXPECT_TRUE(all_of(b == decltype(b)(4)));
}

template <class Abi>
inline void host_check_condition() {
  auto a = Kokkos::Experimental::condition(
      Kokkos::Experimental::simd<std::int32_t, Abi>(1) > 0,
      Kokkos::Experimental::simd<std::uint64_t, Abi>(16),
      Kokkos::Experimental::simd<std::uint64_t, Abi>(20));
  EXPECT_TRUE(all_of(a == decltype(a)(16)));
}

template <typename Abi, typename DataType, size_t n>
KOKKOS_INLINE_FUNCTION void device_check_all_math_ops(
    const DataType (&first_args)[n], const DataType (&second_args)[n]) {
  device_check_binary_op_all_loaders<Abi>(plus(), n, first_args, second_args);
  device_check_binary_op_all_loaders<Abi>(minus(), n, first_args, second_args);
  device_check_binary_op_all_loaders<Abi>(multiplies(), n, first_args,
                                          second_args);

  if constexpr (std::is_same_v<DataType, double>)
    device_check_binary_op_all_loaders<Abi>(divides(), n, first_args,
                                            second_args);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_abi_size() {
  using simd_type = Kokkos::Experimental::simd<DataType, Abi>;
  using mask_type = typename simd_type::mask_type;
  static_assert(simd_type::size() == mask_type::size());
}

template <class Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_math_ops() {
  constexpr size_t n = 11;

  device_check_abi_size<Abi, DataType>();

  if constexpr (std::is_signed_v<DataType>) {
    DataType const first_args[n]  = {1, 2, -1, 10, 0, 1, -2, 10, 0, 1, -2};
    DataType const second_args[n] = {1, 2, 1, 1, 1, -3, -2, 1, 13, -3, -2};
    device_check_all_math_ops<Abi>(first_args, second_args);
  } else {
    DataType const first_args[n]  = {1, 2, 1, 10, 0, 1, 2, 10, 0, 1, 2};
    DataType const second_args[n] = {1, 2, 1, 1, 1, 3, 2, 1, 13, 3, 2};
    device_check_all_math_ops<Abi>(first_args, second_args);
  }
}

template <class Abi>
KOKKOS_INLINE_FUNCTION void device_check_mask_ops() {
  using mask_type = Kokkos::Experimental::simd_mask<double, Abi>;
  kokkos_checker checker;
  checker.truth(!none_of(mask_type(true)));
  checker.truth(none_of(mask_type(false)));
  checker.truth(all_of(mask_type(true)));
  checker.truth(!all_of(mask_type(false)));
}

template <class Abi>
KOKKOS_INLINE_FUNCTION void device_check_conversions() {
  kokkos_checker checker;
  {
    auto a = Kokkos::Experimental::simd<std::uint64_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::int64_t, Abi>(a);
    checker.truth(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd<std::int32_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::uint64_t, Abi>(a);
    checker.truth(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd<std::uint64_t, Abi>(1);
    auto b = Kokkos::Experimental::simd<std::int32_t, Abi>(a);
    checker.truth(all_of(b == decltype(b)(1)));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<double, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(a);
    checker.truth(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::uint64_t, Abi>(a);
    checker.truth(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<std::int64_t, Abi>(a);
    checker.truth(b == decltype(b)(true));
  }
  {
    auto a = Kokkos::Experimental::simd_mask<std::int32_t, Abi>(true);
    auto b = Kokkos::Experimental::simd_mask<double, Abi>(a);
    checker.truth(b == decltype(b)(true));
  }
}

template <class Abi>
KOKKOS_INLINE_FUNCTION void device_check_shifts() {
  kokkos_checker checker;
  auto a = Kokkos::Experimental::simd<std::uint64_t, Abi>(8);
  auto b = a >> 1;
  checker.truth(all_of(b == decltype(b)(4)));
}

template <class Abi>
KOKKOS_INLINE_FUNCTION void device_check_condition() {
  kokkos_checker checker;
  auto a = Kokkos::Experimental::condition(
      Kokkos::Experimental::simd<std::int32_t, Abi>(1) > 0,
      Kokkos::Experimental::simd<std::uint64_t, Abi>(16),
      Kokkos::Experimental::simd<std::uint64_t, Abi>(20));
  checker.truth(all_of(a == decltype(a)(16)));
}

template <typename Abi, typename... DataTypes>
inline void host_check_math_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_math_ops<Abi, DataTypes>(), ...);
}

template <class Abi>
inline void host_check_abi() {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;

  host_check_math_ops_all_types<Abi>(DataTypes());
  host_check_mask_ops<Abi>();
  host_check_conversions<Abi>();
  host_check_shifts<Abi>();
  host_check_condition<Abi>();
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_math_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_math_ops<Abi, DataTypes>(), ...);
}

template <class Abi>
KOKKOS_INLINE_FUNCTION void device_check_abi() {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;

  device_check_math_ops_all_types<Abi>(DataTypes());
  device_check_mask_ops<Abi>();
  device_check_conversions<Abi>();
  device_check_shifts<Abi>();
  device_check_condition<Abi>();
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
