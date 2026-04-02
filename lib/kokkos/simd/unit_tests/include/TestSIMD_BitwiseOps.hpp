// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SIMD_BITWISE_OPS_HPP
#define KOKKOS_TEST_SIMD_BITWISE_OPS_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.simd;
#else
#include <Kokkos_SIMD.hpp>
#endif
#include <SIMDTesting_Utilities.hpp>

// FIXME GCC <= 11.2
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 1130)
#include <limits>
#endif

template <class Abi, class Loader, bool is_compound_op, class BinaryOp, class T>
void host_check_bitwise_op_one_loader(BinaryOp binary_op,
                                      Kokkos::Experimental::Impl::simd_size_t n,
                                      T const* first_args,
                                      T const* second_args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  for (size_type i = 0; i < n; i += width) {
    size_type const nremaining = n - i;
    size_type const nlanes     = Kokkos::min(nremaining, width);
    if ((std::is_same_v<BinaryOp, divides> ||
         std::is_same_v<BinaryOp, divides_eq>)&&nremaining < width)
      continue;
    simd_type first_arg;
    bool const loaded_first_arg =
        loader.host_load(first_args + i, nlanes, first_arg);
    simd_type second_arg;
    bool const loaded_second_arg =
        loader.host_load(second_args + i, nlanes, second_arg);
    if (!(loaded_first_arg && loaded_second_arg)) continue;

    T expected_val[width];
    for (size_type lane = 0; lane < width; ++lane) {
      T tmp              = first_arg[lane];
      expected_val[lane] = binary_op.on_host(tmp, T(second_arg[lane]));
    }

    simd_type expected_result =
        Kokkos::Experimental::simd_unchecked_load<simd_type>(
            expected_val, Kokkos::Experimental::simd_flag_default);
    simd_type const computed_result = binary_op.on_host(first_arg, second_arg);
    host_check_equality(expected_result, computed_result, nlanes);
    if constexpr (is_compound_op) {
      host_check_equality(first_arg, expected_result, nlanes);
    }
  }
}

template <class Abi, class Loader, bool, class UnaryOp, class T>
void host_check_bitwise_op_one_loader(UnaryOp unary_op,
                                      Kokkos::Experimental::Impl::simd_size_t n,
                                      T const* args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr auto width = simd_type::size();
  for (size_type i = 0; i < n; i += width) {
    size_type const nremaining = n - i;
    size_type const nlanes     = Kokkos::min(nremaining, width);
    simd_type arg;
    bool const loaded_arg = loader.host_load(args + i, nlanes, arg);
    if (!loaded_arg) continue;

    auto unary_op_result   = unary_op.on_host(arg);
    using result_simd_type = decltype(unary_op_result);

    typename result_simd_type::value_type expected_val[width];
    for (size_type lane = 0; lane < width; ++lane) {
      expected_val[lane] = unary_op.on_host(T(arg[lane]));
    }

    result_simd_type expected_result =
        Kokkos::Experimental::simd_unchecked_load<result_simd_type>(
            expected_val, Kokkos::Experimental::simd_flag_default);
    auto computed_result = unary_op.on_host(arg);
    host_check_equality(expected_result, computed_result, nlanes);
  }
}

template <class Abi, bool is_compound_op, class Op, class... T>
inline void host_check_bitwise_op_all_loaders(
    Op op, Kokkos::Experimental::Impl::simd_size_t n, T const*... args) {
  host_check_bitwise_op_one_loader<Abi, load_element_aligned, is_compound_op>(
      op, n, args...);
  host_check_bitwise_op_one_loader<Abi, load_masked, is_compound_op>(op, n,
                                                                     args...);
  host_check_bitwise_op_one_loader<Abi, load_as_scalars, is_compound_op>(
      op, n, args...);
  host_check_bitwise_op_one_loader<Abi, load_vector_aligned, is_compound_op>(
      op, n, args...);
}

template <typename Abi, typename DataType,
          Kokkos::Experimental::Impl::simd_size_t n>
inline void host_check_all_bitwise_ops(const DataType (&first_args)[n],
                                       const DataType (&second_args)[n]) {
  host_check_bitwise_op_all_loaders<Abi, false>(bitwise_not(), n, first_args);
  host_check_bitwise_op_all_loaders<Abi, false>(bitwise_and(), n, first_args,
                                                second_args);
  host_check_bitwise_op_all_loaders<Abi, false>(bitwise_or(), n, first_args,
                                                second_args);
  host_check_bitwise_op_all_loaders<Abi, false>(bitwise_xor(), n, first_args,
                                                second_args);
  host_check_bitwise_op_all_loaders<Abi, true>(bitwise_and_eq(), n, first_args,
                                               second_args);
  host_check_bitwise_op_all_loaders<Abi, true>(bitwise_or_eq(), n, first_args,
                                               second_args);
  host_check_bitwise_op_all_loaders<Abi, true>(bitwise_xor_eq(), n, first_args,
                                               second_args);
}

template <typename Abi, typename DataType>
inline void host_check_bitwise_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi> &&
                std::is_integral_v<DataType>) {
    constexpr size_t alignment =
        Kokkos::Experimental::basic_simd<DataType, Abi>::size() *
        sizeof(DataType);

    constexpr int half_shift   = (CHAR_BIT * sizeof(DataType)) / 2;
    constexpr DataType zero    = static_cast<DataType>(0);
    constexpr DataType all_set = ~zero;

    // FIXME GCC <= 11.2, GCC warns about left shift of negative values even
    // though it is well defined in C++20, see
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=103826
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 1130)
    constexpr int signed_half_shift =
        std::is_signed_v<DataType> ? half_shift - 1 : half_shift;
    constexpr DataType lo_set =
        std::numeric_limits<DataType>::max() >> signed_half_shift;
    constexpr DataType hi_set = ~lo_set;
#else
    constexpr DataType hi_set = all_set << half_shift;
    constexpr DataType lo_set = ~hi_set;
#endif

    alignas(alignment) DataType const first_args[] = {
        0,         0,          1,        all_set,  all_set,   all_set,
        all_set,   lo_set,     hi_set,   0,        704475968, 1845076239,
        432747131, 1285171335, 17011965, 139561533};
    alignas(alignment) DataType const second_args[] = {
        0,         1,          1,          0,         all_set,    hi_set,
        lo_set,    lo_set,     hi_set,     hi_set,    1841853988, 747428271,
        357605498, 1412297337, 1663131103, 2062867687};
    host_check_all_bitwise_ops<Abi>(first_args, second_args);
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_bitwise_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_bitwise_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_bitwise_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_bitwise_ops_all_types<Abis>(DataTypes()), ...);
}

template <class Abi, class Loader, bool is_compound_op, class BinaryOp, class T>
KOKKOS_INLINE_FUNCTION void device_check_bitwise_op_one_loader(
    BinaryOp binary_op, Kokkos::Experimental::Impl::simd_size_t n,
    T const* first_args, T const* second_args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  for (size_type i = 0; i < n; i += width) {
    size_type const nremaining = n - i;
    size_type const nlanes     = Kokkos::min(nremaining, width);
    if ((std::is_same_v<BinaryOp, divides> ||
         std::is_same_v<BinaryOp, divides_eq>)&&nremaining < width)
      continue;
    simd_type first_arg;
    bool const loaded_first_arg =
        loader.device_load(first_args + i, nlanes, first_arg);
    simd_type second_arg;
    bool const loaded_second_arg =
        loader.device_load(second_args + i, nlanes, second_arg);
    if (!(loaded_first_arg && loaded_second_arg)) continue;

    T expected_val[width];
    for (size_type lane = 0; lane < width; ++lane) {
      T tmp              = first_arg[lane];
      expected_val[lane] = binary_op.on_device(tmp, T(second_arg[lane]));
    }

    simd_type expected_result =
        Kokkos::Experimental::simd_unchecked_load<simd_type>(
            expected_val, Kokkos::Experimental::simd_flag_default);
    simd_type const computed_result =
        binary_op.on_device(first_arg, second_arg);
    device_check_equality(expected_result, computed_result, nlanes);
    if constexpr (is_compound_op) {
      device_check_equality(first_arg, expected_result, nlanes);
    }
  }
}

template <class Abi, class Loader, bool, class UnaryOp, class T>
KOKKOS_INLINE_FUNCTION void device_check_bitwise_op_one_loader(
    UnaryOp unary_op, Kokkos::Experimental::Impl::simd_size_t n,
    T const* args) {
  Loader loader;
  using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
  using size_type = Kokkos::Experimental::Impl::simd_size_t;

  constexpr size_type width = simd_type::size();
  for (size_type i = 0; i < n; i += width) {
    size_type const nremaining = n - i;
    size_type const nlanes     = Kokkos::min(nremaining, width);
    simd_type arg;
    bool const loaded_arg = loader.device_load(args + i, nlanes, arg);
    if (!loaded_arg) continue;

    auto unary_op_result   = unary_op.on_device(arg);
    using result_simd_type = decltype(unary_op_result);

    typename result_simd_type::value_type expected_val[width];
    for (size_type lane = 0; lane < width; ++lane) {
      expected_val[lane] = unary_op.on_device(T(arg[lane]));
    }

    result_simd_type expected_result =
        Kokkos::Experimental::simd_unchecked_load<result_simd_type>(
            expected_val, Kokkos::Experimental::simd_flag_default);
    auto computed_result = unary_op.on_device(arg);
    device_check_equality(expected_result, computed_result, nlanes);
  }
}

template <class Abi, bool is_compound_op, class Op, class... T>
KOKKOS_INLINE_FUNCTION void device_check_bitwise_op_all_loaders(
    Op op, Kokkos::Experimental::Impl::simd_size_t n, T const*... args) {
  device_check_bitwise_op_one_loader<Abi, load_element_aligned, is_compound_op>(
      op, n, args...);
  device_check_bitwise_op_one_loader<Abi, load_masked, is_compound_op>(op, n,
                                                                       args...);
  device_check_bitwise_op_one_loader<Abi, load_as_scalars, is_compound_op>(
      op, n, args...);
  device_check_bitwise_op_one_loader<Abi, load_vector_aligned, is_compound_op>(
      op, n, args...);
}

template <typename Abi, typename DataType,
          Kokkos::Experimental::Impl::simd_size_t n>
KOKKOS_INLINE_FUNCTION void device_check_all_bitwise_ops(
    const DataType (&first_args)[n], const DataType (&second_args)[n]) {
  device_check_bitwise_op_all_loaders<Abi, false>(bitwise_not(), n, first_args);
  device_check_bitwise_op_all_loaders<Abi, false>(bitwise_and(), n, first_args,
                                                  second_args);
  device_check_bitwise_op_all_loaders<Abi, false>(bitwise_or(), n, first_args,
                                                  second_args);
  device_check_bitwise_op_all_loaders<Abi, false>(bitwise_xor(), n, first_args,
                                                  second_args);
  device_check_bitwise_op_all_loaders<Abi, true>(bitwise_and_eq(), n,
                                                 first_args, second_args);
  device_check_bitwise_op_all_loaders<Abi, true>(bitwise_or_eq(), n, first_args,
                                                 second_args);
  device_check_bitwise_op_all_loaders<Abi, true>(bitwise_xor_eq(), n,
                                                 first_args, second_args);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_bitwise_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi> &&
                std::is_integral_v<DataType>) {
    constexpr size_t alignment =
        Kokkos::Experimental::basic_simd<DataType, Abi>::size() *
        sizeof(DataType);

    constexpr int half_shift   = (CHAR_BIT * sizeof(DataType)) / 2;
    constexpr DataType zero    = static_cast<DataType>(0);
    constexpr DataType all_set = ~zero;

    // FIXME GCC <= 11.2, GCC warns about left shift of negative values even
    // though it is well defined in C++20, see
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=103826
#if defined(KOKKOS_COMPILER_GNU) && (KOKKOS_COMPILER_GNU < 1130)
    constexpr int signed_half_shift =
        std::is_signed_v<DataType> ? half_shift - 1 : half_shift;
    constexpr DataType lo_set =
        Kokkos::Experimental::finite_max_v<DataType> >> signed_half_shift;
    constexpr DataType hi_set = ~lo_set;
#else
    constexpr DataType hi_set = all_set << half_shift;
    constexpr DataType lo_set = ~hi_set;
#endif

    alignas(alignment) DataType const first_args[] = {
        0,         0,          1,        all_set,  all_set,   all_set,
        all_set,   lo_set,     hi_set,   0,        704475968, 1845076239,
        432747131, 1285171335, 17011965, 139561533};
    alignas(alignment) DataType const second_args[] = {
        0,         1,          1,          0,         all_set,    hi_set,
        lo_set,    lo_set,     hi_set,     hi_set,    1841853988, 747428271,
        357605498, 1412297337, 1663131103, 2062867687};
    device_check_all_bitwise_ops<Abi>(first_args, second_args);
  }
}

template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_bitwise_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_bitwise_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_bitwise_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_bitwise_ops_all_types<Abis>(DataTypes()), ...);
}

class simd_device_bitwise_ops_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_bitwise_ops_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

template <typename Abi, typename DataType>
inline void host_check_mask_bitwise_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    using mask_type = Kokkos::Experimental::basic_simd_mask<DataType, Abi>;
    using size_type = Kokkos::Experimental::Impl::simd_size_t;

    constexpr size_type width = mask_type::size();

    mask_type const all(true);
    mask_type const none(false);

    host_check_mask_equality(~all, none);
    host_check_mask_equality(~none, all);

    host_check_mask_equality(all & all, all);
    host_check_mask_equality(all & none, none);
    host_check_mask_equality(none & none, none);

    host_check_mask_equality(all | all, all);
    host_check_mask_equality(all | none, all);
    host_check_mask_equality(none | none, none);

    host_check_mask_equality(all ^ all, none);
    host_check_mask_equality(all ^ none, all);
    host_check_mask_equality(none ^ none, none);

    if constexpr (!std::is_same_v<Abi,
                                  Kokkos::Experimental::simd_abi::scalar>) {
      mask_type const hi([=](std::size_t const i) { return i >= width / 2; });
      mask_type const lo([=](std::size_t const i) { return i < width / 2; });
      mask_type const even([=](std::size_t const i) { return i % 2 == 0; });
      mask_type const odd([=](std::size_t const i) { return i % 2 == 1; });

      host_check_mask_equality(~hi, lo);
      host_check_mask_equality(~lo, hi);
      host_check_mask_equality(~even, odd);
      host_check_mask_equality(~odd, even);

      host_check_mask_equality(all & hi, hi);
      host_check_mask_equality(hi & hi, hi);
      host_check_mask_equality(lo & lo, lo);
      host_check_mask_equality(hi & lo, none);
      host_check_mask_equality(even & even, even);
      host_check_mask_equality(even & odd, none);

      host_check_mask_equality(all | hi, all);
      host_check_mask_equality(none | hi, hi);
      host_check_mask_equality(hi | hi, hi);
      host_check_mask_equality(hi | lo, all);
      host_check_mask_equality(even | even, even);
      host_check_mask_equality(even | odd, all);

      host_check_mask_equality(all ^ hi, lo);
      host_check_mask_equality(none ^ hi, hi);
      host_check_mask_equality(hi ^ hi, none);
      host_check_mask_equality(hi ^ lo, all);
      host_check_mask_equality(even ^ even, none);
      host_check_mask_equality(even ^ odd, all);
    }
  }
}

template <typename Abi, typename DataType, typename Op>
inline void host_check_mask_bitwise_assignment_op(
    Op op, Kokkos::Experimental::basic_simd_mask<Abi, DataType> lhs,
    Kokkos::Experimental::basic_simd_mask<Abi, DataType> const& rhs,
    Kokkos::Experimental::basic_simd_mask<Abi, DataType> const& expected) {
  auto const res = op.on_host(lhs, rhs);
  host_check_mask_equality(res, expected);
  host_check_mask_equality(lhs, expected);
}

template <typename Abi, typename DataType>
inline void host_check_mask_bitwise_assignment_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    using mask_type = Kokkos::Experimental::basic_simd_mask<DataType, Abi>;
    using size_type = Kokkos::Experimental::Impl::simd_size_t;

    constexpr size_type width = mask_type::size();

    mask_type const all(true);
    mask_type const none(false);

    host_check_mask_bitwise_assignment_op(bitwise_and_eq(), all, all, all);
    host_check_mask_bitwise_assignment_op(bitwise_and_eq(), all, none, none);
    host_check_mask_bitwise_assignment_op(bitwise_and_eq(), none, none, none);

    host_check_mask_bitwise_assignment_op(bitwise_or_eq(), all, all, all);
    host_check_mask_bitwise_assignment_op(bitwise_or_eq(), all, none, all);
    host_check_mask_bitwise_assignment_op(bitwise_or_eq(), none, none, none);

    host_check_mask_bitwise_assignment_op(bitwise_xor_eq(), all, all, none);
    host_check_mask_bitwise_assignment_op(bitwise_xor_eq(), all, none, all);
    host_check_mask_bitwise_assignment_op(bitwise_xor_eq(), none, none, none);

    if constexpr (!std::is_same_v<Abi,
                                  Kokkos::Experimental::simd_abi::scalar>) {
      mask_type const hi([=](std::size_t const i) { return i >= width / 2; });
      mask_type const lo([=](std::size_t const i) { return i < width / 2; });
      mask_type const even([=](std::size_t const i) { return i % 2 == 0; });
      mask_type const odd([=](std::size_t const i) { return i % 2 == 1; });

      host_check_mask_bitwise_assignment_op(bitwise_and_eq(), all, hi, hi);
      host_check_mask_bitwise_assignment_op(bitwise_and_eq(), hi, hi, hi);
      host_check_mask_bitwise_assignment_op(bitwise_and_eq(), lo, lo, lo);
      host_check_mask_bitwise_assignment_op(bitwise_and_eq(), hi, lo, none);
      host_check_mask_bitwise_assignment_op(bitwise_and_eq(), even, even, even);
      host_check_mask_bitwise_assignment_op(bitwise_and_eq(), even, odd, none);

      host_check_mask_bitwise_assignment_op(bitwise_or_eq(), all, hi, all);
      host_check_mask_bitwise_assignment_op(bitwise_or_eq(), none, hi, hi);
      host_check_mask_bitwise_assignment_op(bitwise_or_eq(), hi, hi, hi);
      host_check_mask_bitwise_assignment_op(bitwise_or_eq(), hi, lo, all);
      host_check_mask_bitwise_assignment_op(bitwise_or_eq(), even, even, even);
      host_check_mask_bitwise_assignment_op(bitwise_or_eq(), even, odd, all);

      host_check_mask_bitwise_assignment_op(bitwise_xor_eq(), all, hi, lo);
      host_check_mask_bitwise_assignment_op(bitwise_xor_eq(), none, hi, hi);
      host_check_mask_bitwise_assignment_op(bitwise_xor_eq(), hi, hi, none);
      host_check_mask_bitwise_assignment_op(bitwise_xor_eq(), hi, lo, all);
      host_check_mask_bitwise_assignment_op(bitwise_xor_eq(), even, even, none);
      host_check_mask_bitwise_assignment_op(bitwise_xor_eq(), even, odd, all);
    }
  }
}

template <typename Abi, typename... DataTypes>
inline void host_check_mask_bitwise_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (host_check_mask_bitwise_ops<Abi, DataTypes>(), ...);
  (host_check_mask_bitwise_assignment_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
inline void host_check_mask_bitwise_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (host_check_mask_bitwise_ops_all_types<Abis>(DataTypes()), ...);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_mask_bitwise_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    using mask_type = Kokkos::Experimental::basic_simd_mask<DataType, Abi>;

    mask_type const all(true);
    mask_type const none(false);

    device_check_mask_equality(~all, none);
    device_check_mask_equality(~none, all);

    device_check_mask_equality(all & all, all);
    device_check_mask_equality(all & none, none);
    device_check_mask_equality(none & none, none);

    device_check_mask_equality(all | all, all);
    device_check_mask_equality(all | none, all);
    device_check_mask_equality(none | none, none);

    device_check_mask_equality(all ^ all, none);
    device_check_mask_equality(all ^ none, all);
    device_check_mask_equality(none ^ none, none);
  }
}

template <typename Abi, typename DataType, typename Op>
KOKKOS_INLINE_FUNCTION void device_check_mask_bitwise_assignment_op(
    Op op, Kokkos::Experimental::basic_simd_mask<Abi, DataType> lhs,
    Kokkos::Experimental::basic_simd_mask<Abi, DataType> const& rhs,
    Kokkos::Experimental::basic_simd_mask<Abi, DataType> const& expected) {
  auto const res = op.on_device(lhs, rhs);
  device_check_mask_equality(res, expected);
  device_check_mask_equality(lhs, expected);
}

template <typename Abi, typename DataType>
KOKKOS_INLINE_FUNCTION void device_check_mask_bitwise_assignment_ops() {
  if constexpr (is_simd_avail_v<DataType, Abi>) {
    using mask_type = Kokkos::Experimental::basic_simd_mask<DataType, Abi>;

    mask_type const all(true);
    mask_type const none(false);

    device_check_mask_bitwise_assignment_op(bitwise_and_eq(), all, all, all);
    device_check_mask_bitwise_assignment_op(bitwise_and_eq(), all, none, none);
    device_check_mask_bitwise_assignment_op(bitwise_and_eq(), none, none, none);

    device_check_mask_bitwise_assignment_op(bitwise_or_eq(), all, all, all);
    device_check_mask_bitwise_assignment_op(bitwise_or_eq(), all, none, all);
    device_check_mask_bitwise_assignment_op(bitwise_or_eq(), none, none, none);

    device_check_mask_bitwise_assignment_op(bitwise_xor_eq(), all, all, none);
    device_check_mask_bitwise_assignment_op(bitwise_xor_eq(), all, none, all);
    device_check_mask_bitwise_assignment_op(bitwise_xor_eq(), none, none, none);
  }
}
template <typename Abi, typename... DataTypes>
KOKKOS_INLINE_FUNCTION void device_check_mask_bitwise_ops_all_types(
    Kokkos::Experimental::Impl::data_types<DataTypes...>) {
  (device_check_mask_bitwise_ops<Abi, DataTypes>(), ...);
  (device_check_mask_bitwise_assignment_ops<Abi, DataTypes>(), ...);
}

template <typename... Abis>
KOKKOS_INLINE_FUNCTION void device_check_mask_bitwise_ops_all_abis(
    Kokkos::Experimental::Impl::abi_set<Abis...>) {
  using DataTypes = Kokkos::Experimental::Impl::data_type_set;
  (device_check_mask_bitwise_ops_all_types<Abis>(DataTypes()), ...);
}

class simd_device_mask_bitwise_ops_functor {
 public:
  KOKKOS_INLINE_FUNCTION void operator()(int) const {
    device_check_mask_bitwise_ops_all_abis(
        Kokkos::Experimental::Impl::device_abi_set());
  }
};

TEST(simd, host_bitwise_ops) {
  host_check_bitwise_ops_all_abis(Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, host_mask_bitwise_ops) {
  host_check_mask_bitwise_ops_all_abis(
      Kokkos::Experimental::Impl::host_abi_set());
}

TEST(simd, device_bitwise_ops) {
  Kokkos::parallel_for(1, simd_device_bitwise_ops_functor());
  Kokkos::fence();
}

TEST(simd, device_mask_bitwise_ops) {
  Kokkos::parallel_for(1, simd_device_mask_bitwise_ops_functor());
  Kokkos::fence();
}

#endif
