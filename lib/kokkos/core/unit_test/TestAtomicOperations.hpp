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

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

namespace TestAtomicOperations {

struct AddAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    Kokkos::atomic_add(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_add(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_add_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return old + update;
  }
  static const char* name() { return "add"; }
};

struct SubAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    Kokkos::atomic_sub(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_sub(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_sub_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return old - update;
  }
  static const char* name() { return "sub"; }
};

struct IncAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T) {
    Kokkos::atomic_inc(ptr_op);
    T old_val = Kokkos::atomic_fetch_inc(ptr_fetch_op);
    T new_val = Kokkos::atomic_inc_fetch(ptr_op_fetch);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T) {
    return old + 1;
  }
  static const char* name() { return "inc"; }
};

struct DecAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T) {
    Kokkos::atomic_dec(ptr_op);
    T old_val = Kokkos::atomic_fetch_dec(ptr_fetch_op);
    T new_val = Kokkos::atomic_dec_fetch(ptr_op_fetch);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T) {
    return old - 1;
  }
  static const char* name() { return "dec"; }
};

struct MaxAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    Kokkos::atomic_max(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_max(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_max_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return update > old ? update : old;
  }
  static const char* name() { return "max"; }
};

struct MinAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    Kokkos::atomic_min(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_min(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_min_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return update < old ? update : old;
  }
  static const char* name() { return "min"; }
};

struct MulAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    Kokkos::atomic_mul(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_mul(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_mul_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return old * update;
  }
  static const char* name() { return "mul"; }
};

struct DivAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    Kokkos::atomic_div(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_div(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_div_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return old / update;
  }
  static const char* name() { return "div"; }
};

struct ModAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    // Kokkos::atomic_mod(ptr_op, update);
    (void)Kokkos::atomic_fetch_mod(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_mod(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_mod_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return old % update;
  }
  static const char* name() { return "mod"; }
};

struct AndAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    Kokkos::atomic_and(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_and(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_and_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return old & update;
  }
  static const char* name() { return "and"; }
};

struct OrAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    Kokkos::atomic_or(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_or(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_or_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return old | update;
  }
  static const char* name() { return "or"; }
};

struct XorAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    // Kokkos::atomic_xor(ptr_op, update);
    (void)Kokkos::atomic_fetch_xor(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_xor(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_xor_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return old ^ update;
  }
  static const char* name() { return "xor"; }
};

struct NandAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    // Kokkos::atomic_nand(ptr_op, update);
    (void)Kokkos::atomic_fetch_nand(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_nand(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_nand_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return ~(old & update);
  }
  static const char* name() { return "nand"; }
};

struct LShiftAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    // Kokkos::atomic_lshift(ptr_op, update);
    (void)Kokkos::atomic_fetch_lshift(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_lshift(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_lshift_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return old << update;
  }
  static const char* name() { return "lshift"; }
};

struct RShiftAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    // Kokkos::atomic_rshift(ptr_op, update); not implemented
    (void)Kokkos::atomic_fetch_rshift(ptr_op, update);
    T old_val = Kokkos::atomic_fetch_rshift(ptr_fetch_op, update);
    T new_val = Kokkos::atomic_rshift_fetch(ptr_op_fetch, update);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T update) {
    return old >> update;
  }
  static const char* name() { return "rshift"; }
};

struct LoadStoreAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T update) {
    T old_val = Kokkos::atomic_load(ptr_op);
    Kokkos::atomic_store(ptr_op, update);
    Kokkos::atomic_store(ptr_op_fetch, update);
    Kokkos::atomic_store(ptr_fetch_op, update);
    return Kokkos::pair<T, T>(old_val, update);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T, T update) {
    return update;
  }
  static const char* name() { return "load/store"; }
};

struct IncModAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T wrap_value) {
    // no atomic_inc_mod in desul
    (void)desul::atomic_fetch_inc_mod(ptr_op, wrap_value,
                                      desul::MemoryOrderRelaxed(),
                                      desul::MemoryScopeDevice());
    T old_val = desul::atomic_fetch_inc_mod(ptr_fetch_op, wrap_value,
                                            desul::MemoryOrderRelaxed(),
                                            desul::MemoryScopeDevice());
    // no atomic_inc_mod_fetch in desul
    (void)desul::atomic_fetch_inc_mod(ptr_op_fetch, wrap_value,
                                      desul::MemoryOrderRelaxed(),
                                      desul::MemoryScopeDevice());
    T new_val = op(old_val, wrap_value);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T wrap_value) {
    return old + 1 > wrap_value ? 0 : old + 1;
  }
  static const char* name() { return "inc_mod"; }
};

struct DecModAtomicTest {
  template <class T>
  KOKKOS_FUNCTION static auto atomic_op(T* ptr_op, T* ptr_fetch_op,
                                        T* ptr_op_fetch, T wrap_value) {
    // no atomic_dec_mod in desul
    (void)desul::atomic_fetch_dec_mod(ptr_op, wrap_value,
                                      desul::MemoryOrderRelaxed(),
                                      desul::MemoryScopeDevice());
    T old_val = desul::atomic_fetch_dec_mod(ptr_fetch_op, wrap_value,
                                            desul::MemoryOrderRelaxed(),
                                            desul::MemoryScopeDevice());
    // no atomic_dec_mod_fetch in desul
    (void)desul::atomic_fetch_dec_mod(ptr_op_fetch, wrap_value,
                                      desul::MemoryOrderRelaxed(),
                                      desul::MemoryScopeDevice());
    T new_val = op(old_val, wrap_value);
    return Kokkos::pair<T, T>(old_val, new_val);
  }
  template <class T>
  KOKKOS_FUNCTION static T op(T old, T wrap_value) {
    return ((old == 0) || (old > wrap_value)) ? wrap_value : old - 1;
  }
  static const char* name() { return "dec_mod"; }
};

template <class Op, class T, class ExecSpace>
bool atomic_op_test(T old_val, T update) {
  Kokkos::View<T[3], ExecSpace> op_data("op_data");
  Kokkos::deep_copy(op_data, old_val);
  int result = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecSpace>(0, 1),
      KOKKOS_LAMBDA(int, int& local_result) {
        auto fetch_result =
            Op::atomic_op(&op_data(0), &op_data(1), &op_data(2), update);
        T expected_val = Op::op(old_val, update);
        Kokkos::memory_fence();
        if (op_data(0) != expected_val) local_result += 1;
        if (op_data(1) != expected_val) local_result += 2;
        if (op_data(2) != expected_val) local_result += 4;
        if (fetch_result.first != old_val) local_result += 8;
        if (fetch_result.second != expected_val) local_result += 16;
      },
      result);
  if ((result & 1) != 0)
    printf("atomic_%s failed with type %s\n", Op::name(), typeid(T).name());
  if ((result & 2) != 0)
    printf("atomic_fetch_%s failed with type %s\n", Op::name(),
           typeid(T).name());
  if ((result & 4) != 0)
    printf("atomic_%s_fetch failed with type %s\n", Op::name(),
           typeid(T).name());
  if ((result & 8) != 0)
    printf("atomic_fetch_%s did not return old value with type %s\n",
           Op::name(), typeid(T).name());
  if ((result & 16) != 0)
    printf("atomic_%s_fetch did not return updated value with type %s\n",
           Op::name(), typeid(T).name());

  return result == 0;
}

template <class T>
constexpr T relative_error_threshold = T(1.0e-15);

template <class Op, class T, class ExecSpace>
bool atomic_op_test_rel(T old_val, T update) {
  Kokkos::View<T[3], ExecSpace> op_data("op_data");
  Kokkos::deep_copy(op_data, old_val);
  int result = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecSpace>(0, 1),
      KOKKOS_LAMBDA(int, int& local_result) {
        auto fetch_result =
            Op::atomic_op(&op_data(0), &op_data(1), &op_data(2), update);
        T expected_val = Op::op(old_val, update);
        Kokkos::memory_fence();
        if (expected_val == T(0)) {
          if (fabs(op_data(0)) > relative_error_threshold<T>) local_result += 1;
          if (fabs(op_data(1)) > relative_error_threshold<T>) local_result += 2;
          if (fabs(op_data(2)) > relative_error_threshold<T>) local_result += 4;
          if (fetch_result.first != old_val) local_result += 8;
          if (fabs(fetch_result.second) > relative_error_threshold<T>)
            local_result += 16;
        } else {
          if (fabs((op_data(0) - expected_val) / expected_val) >
              relative_error_threshold<T>)
            local_result += 1;
          if (fabs((op_data(1) - expected_val) / expected_val) >
              relative_error_threshold<T>)
            local_result += 2;
          if (fabs((op_data(2) - expected_val) / expected_val) >
              relative_error_threshold<T>)
            local_result += 4;
          if (fetch_result.first != old_val) local_result += 8;
          if (fabs((fetch_result.second - expected_val) / expected_val) >
              relative_error_threshold<T>)
            local_result += 16;
        }
      },
      result);
  if ((result & 1) != 0)
    printf("atomic_%s failed with type %s\n", Op::name(), typeid(T).name());
  if ((result & 2) != 0)
    printf("atomic_fetch_%s failed with type %s\n", Op::name(),
           typeid(T).name());
  if ((result & 4) != 0)
    printf("atomic_%s_fetch failed with type %s\n", Op::name(),
           typeid(T).name());
  if ((result & 8) != 0)
    printf("atomic_fetch_%s did not return old value with type %s\n",
           Op::name(), typeid(T).name());
  if ((result & 16) != 0)
    printf("atomic_%s_fetch did not return updated value with type %s\n",
           Op::name(), typeid(T).name());

  return result == 0;
}

//---------------------------------------------------
//--------------atomic_test_control------------------
//---------------------------------------------------

template <class T, class ExecSpace>
bool AtomicOperationsTestIntegralType(int old_val_in, int update_in, int test) {
  T old_val = static_cast<T>(old_val_in);
  T update  = static_cast<T>(update_in);
  switch (test) {
    case 0: return atomic_op_test<AddAtomicTest, T, ExecSpace>(old_val, update);
    case 1: return atomic_op_test<SubAtomicTest, T, ExecSpace>(old_val, update);
    case 2: return atomic_op_test<MaxAtomicTest, T, ExecSpace>(old_val, update);
    case 3: return atomic_op_test<MinAtomicTest, T, ExecSpace>(old_val, update);
    case 4: return atomic_op_test<MulAtomicTest, T, ExecSpace>(old_val, update);
    case 5:
      return update != 0
                 ? atomic_op_test<DivAtomicTest, T, ExecSpace>(old_val, update)
                 : true;
    case 6:
      return update != 0
                 ? atomic_op_test<ModAtomicTest, T, ExecSpace>(old_val, update)
                 : true;
    case 7: return atomic_op_test<AndAtomicTest, T, ExecSpace>(old_val, update);
    case 8: return atomic_op_test<OrAtomicTest, T, ExecSpace>(old_val, update);
    case 9: return atomic_op_test<XorAtomicTest, T, ExecSpace>(old_val, update);
    case 10:
      return atomic_op_test<NandAtomicTest, T, ExecSpace>(old_val, update);
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
    // FIXME_NVHPC: atomic-fetch-shift operation fails due to NVHPC OpenACC
    // compiler bugs, which are reported to NVIDIA.
    case 11: return true;
    case 12: return true;
#else
    case 11:
      return (std::make_signed_t<T>(update_in) >= 0 &&
              std::make_signed_t<T>(old_val) >= 0)
                 ? atomic_op_test<LShiftAtomicTest, T, ExecSpace>(old_val,
                                                                  update)
                 : true;
    case 12:
      return update_in >= 0 ? atomic_op_test<RShiftAtomicTest, T, ExecSpace>(
                                  old_val, update)
                            : true;
#endif
    case 13:
      return atomic_op_test<IncAtomicTest, T, ExecSpace>(old_val, update);
    case 14:
      return atomic_op_test<DecAtomicTest, T, ExecSpace>(old_val, update);
    case 15:
      return atomic_op_test<LoadStoreAtomicTest, T, ExecSpace>(old_val, update);
  }

  return true;
}

template <class T, class ExecSpace>
bool AtomicOperationsTestUnsignedIntegralType(int old_val_in, int update_in,
                                              int test) {
  T old_val = static_cast<T>(old_val_in);
  T update  = static_cast<T>(update_in);
  switch (test) {
    case 1:
      return atomic_op_test<IncModAtomicTest, T, ExecSpace>(old_val, update);
    case 2:
      return atomic_op_test<DecModAtomicTest, T, ExecSpace>(old_val, update);
  }

  return true;
}

template <class T, class ExecSpace>
bool AtomicOperationsTestNonIntegralType(int old_val_in, int update_in,
                                         int test) {
  T old_val = static_cast<T>(old_val_in);
  T update  = static_cast<T>(update_in);
  switch (test) {
    case 0: return atomic_op_test<AddAtomicTest, T, ExecSpace>(old_val, update);
    case 1: return atomic_op_test<SubAtomicTest, T, ExecSpace>(old_val, update);
    case 2: return atomic_op_test<MaxAtomicTest, T, ExecSpace>(old_val, update);
    case 3: return atomic_op_test<MinAtomicTest, T, ExecSpace>(old_val, update);
    case 4: return atomic_op_test<MulAtomicTest, T, ExecSpace>(old_val, update);
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
    // NVHPC may use different internal precisions for the device and host
    // atomic operations. Therefore, relative errors are used to compare the
    // host results and device results.
    case 5:
      return update != 0 ? atomic_op_test_rel<DivAtomicTest, T, ExecSpace>(
                               old_val, update)
                         : true;
#else
    case 5:
      return update != 0
                 ? atomic_op_test<DivAtomicTest, T, ExecSpace>(old_val, update)
                 : true;
#endif
    case 6:
      return atomic_op_test<LoadStoreAtomicTest, T, ExecSpace>(old_val, update);
  }

  return true;
}
}  // namespace TestAtomicOperations
