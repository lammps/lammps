// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_SIMD.hpp>

export module kokkos.simd_impl;

export {
  namespace Kokkos::Experimental {
  namespace simd_abi::Impl {
  using ::Kokkos::Experimental::simd_abi::Impl::host_fixed_native;
  using ::Kokkos::Experimental::simd_abi::Impl::native_abi;
  using ::Kokkos::Experimental::simd_abi::Impl::native_fixed_abi;
  }  // namespace simd_abi::Impl

  namespace Impl {
  using ::Kokkos::Experimental::Impl::abi_set;
  using ::Kokkos::Experimental::Impl::data_type_set;
  using ::Kokkos::Experimental::Impl::data_types;
  using ::Kokkos::Experimental::Impl::device_abi_set;
  using ::Kokkos::Experimental::Impl::host_abi_set;
  using ::Kokkos::Experimental::Impl::Identity;
  }  // namespace Impl
  }  // namespace Kokkos::Experimental
}
