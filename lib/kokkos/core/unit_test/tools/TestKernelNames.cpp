// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif
#include <Kokkos_TypeInfo.hpp>

#include <gtest/gtest.h>

namespace {

#ifdef KOKKOS_ENABLE_IMPL_TYPEINFO
template <class T>
std::string typeid_name(T const&) {
  return std::string(Kokkos::Impl::TypeInfo<T>::name());
}
#else
template <class T>
std::string typeid_name(T const&) {
  return typeid(T).name();
}
#endif

std::string last_parallel_for;
std::string last_parallel_reduce;
std::string last_parallel_scan;

void get_parallel_for_kernel_name(char const* kernelName, uint32_t /*deviceID*/,
                                  uint64_t* /*kernelID*/) {
  last_parallel_for = kernelName;
}

void get_parallel_reduce_kernel_name(char const* kernelName,
                                     uint32_t /*deviceID*/,
                                     uint64_t* /*kernelID*/) {
  last_parallel_reduce = kernelName;
}

void get_parallel_scan_kernel_name(char const* kernelName,
                                   uint32_t /*deviceID*/,
                                   uint64_t* /*kernelID*/) {
  last_parallel_scan = kernelName;
}

struct WorkTag {};

void test_kernel_name_parallel_for() {
  Kokkos::Tools::Experimental::set_begin_parallel_for_callback(
      get_parallel_for_kernel_name);

  using ExecutionSpace = Kokkos::DefaultExecutionSpace;

  {
    std::string const my_label = "my_parallel_for_range_policy";

    auto const my_lambda = KOKKOS_LAMBDA(int){};
    Kokkos::parallel_for(my_label, Kokkos::RangePolicy<ExecutionSpace>(0, 1),
                         my_lambda);
    ASSERT_EQ(last_parallel_for, my_label);

    Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(0, 1), my_lambda);
    ASSERT_EQ(last_parallel_for, typeid_name(my_lambda));
    ASSERT_FALSE(last_parallel_for.starts_with("const "))
        << last_parallel_for << " is const-qualified";

    auto const my_lambda_with_tag = KOKKOS_LAMBDA(WorkTag, int){};
    Kokkos::parallel_for(my_label,
                         Kokkos::RangePolicy<ExecutionSpace, WorkTag>(0, 1),
                         my_lambda_with_tag);
    ASSERT_EQ(last_parallel_for, my_label);

    Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace, WorkTag>(0, 1),
                         my_lambda_with_tag);
    ASSERT_EQ(last_parallel_for,
              typeid_name(my_lambda_with_tag) + "/" + typeid_name(WorkTag{}));
    ASSERT_FALSE(last_parallel_for.starts_with("const "))
        << last_parallel_for << " is const-qualified";
  }

  Kokkos::Tools::Experimental::set_begin_parallel_for_callback(nullptr);
}

void test_kernel_name_parallel_reduce() {
  Kokkos::Tools::Experimental::set_begin_parallel_reduce_callback(
      get_parallel_reduce_kernel_name);

  using ExecutionSpace = Kokkos::DefaultExecutionSpace;

  {
    std::string const my_label = "my_parallel_reduce_range_policy";
    float my_result;

    auto const my_lambda = KOKKOS_LAMBDA(int, float&){};
    Kokkos::parallel_reduce(my_label, Kokkos::RangePolicy<ExecutionSpace>(0, 1),
                            my_lambda, my_result);
    ASSERT_EQ(last_parallel_reduce, my_label);

    Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(0, 1),
                            my_lambda, my_result);
#ifndef KOKKOS_COMPILER_MSVC
    ASSERT_NE(last_parallel_reduce.find(typeid_name(my_lambda)),
              std::string::npos)
        << last_parallel_reduce << " does not contain "
        << typeid_name(
               my_lambda);  // internally using Impl::CombinedFunctorReducer
                            // but the name should still include the lambda as
                            // template parameter
#endif
    ASSERT_FALSE(last_parallel_reduce.starts_with("const "))
        << last_parallel_reduce << " is const-qualified";

    auto const my_lambda_with_tag = KOKKOS_LAMBDA(WorkTag, int, float&){};
    Kokkos::parallel_reduce(my_label,
                            Kokkos::RangePolicy<ExecutionSpace, WorkTag>(0, 1),
                            my_lambda_with_tag, my_result);
    ASSERT_EQ(last_parallel_reduce, my_label);

    Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace, WorkTag>(0, 1),
                            my_lambda_with_tag, my_result);
    auto const suffix = std::string("/") + typeid_name(WorkTag{});
    ASSERT_EQ(last_parallel_reduce.find(suffix),
              last_parallel_reduce.length() - suffix.length());
    ASSERT_FALSE(last_parallel_reduce.starts_with("const "))
        << last_parallel_reduce << " is const-qualified";
  }

  Kokkos::Tools::Experimental::set_begin_parallel_reduce_callback(nullptr);
}

void test_kernel_name_parallel_scan() {
  Kokkos::Tools::Experimental::set_begin_parallel_scan_callback(
      get_parallel_scan_kernel_name);

  using ExecutionSpace = Kokkos::DefaultExecutionSpace;

  {
    std::string const my_label = "my_parallel_scan_range_policy";

    auto const my_lambda = KOKKOS_LAMBDA(int, float&, bool){};
    Kokkos::parallel_scan(my_label, Kokkos::RangePolicy<ExecutionSpace>(0, 1),
                          my_lambda);
    ASSERT_EQ(last_parallel_scan, my_label);

    Kokkos::parallel_scan(Kokkos::RangePolicy<ExecutionSpace>(0, 1), my_lambda);
    ASSERT_EQ(last_parallel_scan, typeid_name(my_lambda));
    ASSERT_FALSE(last_parallel_scan.starts_with("const "))
        << last_parallel_scan << " is const-qualified";

    auto const my_lambda_with_tag = KOKKOS_LAMBDA(WorkTag, int, float&, bool){};
    Kokkos::parallel_scan(my_label,
                          Kokkos::RangePolicy<ExecutionSpace, WorkTag>(0, 1),
                          my_lambda_with_tag);
    ASSERT_EQ(last_parallel_scan, my_label);

    Kokkos::parallel_scan(Kokkos::RangePolicy<ExecutionSpace, WorkTag>(0, 1),
                          my_lambda_with_tag);
    ASSERT_EQ(last_parallel_scan,
              typeid_name(my_lambda_with_tag) + "/" + typeid_name(WorkTag{}));
    ASSERT_FALSE(last_parallel_scan.starts_with("const "))
        << last_parallel_scan << " is const-qualified";
  }

  Kokkos::Tools::Experimental::set_begin_parallel_scan_callback(nullptr);
}

TEST(kokkosp, kernel_name_parallel_for) { test_kernel_name_parallel_for(); }

TEST(kokkosp, kernel_name_parallel_reduce) {
  test_kernel_name_parallel_reduce();
}

TEST(kokkosp, kernel_name_parallel_scan) { test_kernel_name_parallel_scan(); }

TEST(kokkosp, kernel_name_internal) {
  struct ThisType {};
  {
    std::string const label("my_label");
    Kokkos::Impl::ParallelConstructName<ThisType, void> pcn(label);
    ASSERT_EQ(pcn.get(), label);
    std::string const empty_label("");
    Kokkos::Impl::ParallelConstructName<ThisType, void> empty_pcn(empty_label);
    ASSERT_EQ(empty_pcn.get(), typeid_name(ThisType{}));
  }
  {
    std::string const label("my_label");
    Kokkos::Impl::ParallelConstructName<ThisType, WorkTag> pcn(label);
    ASSERT_EQ(pcn.get(), label);
    std::string const empty_label("");
    Kokkos::Impl::ParallelConstructName<ThisType, WorkTag> empty_pcn(
        empty_label);
    ASSERT_EQ(empty_pcn.get(),
              typeid_name(ThisType{}) + "/" + typeid_name(WorkTag{}));
  }
}

}  // namespace
