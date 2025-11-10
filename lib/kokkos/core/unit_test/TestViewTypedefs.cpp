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

namespace {

// clang-format off
template<class DataType>
struct data_analysis {
  using data_type = DataType;
  using const_data_type = const DataType;
  using runtime_data_type = DataType;
  using runtime_const_data_type = const DataType;
  using non_const_data_type = std::remove_const_t<DataType>;
};

template<class DataType>
struct data_analysis<DataType*> {
  using data_type = typename data_analysis<DataType>::data_type*;
  using const_data_type = typename data_analysis<DataType>::const_data_type*;
  using runtime_data_type = typename data_analysis<DataType>::runtime_data_type*;
  using runtime_const_data_type = typename data_analysis<DataType>::runtime_const_data_type*;
  using non_const_data_type = typename data_analysis<DataType>::non_const_data_type*;
};

template<class DataType, size_t N>
struct data_analysis<DataType[N]> {
  using data_type = typename data_analysis<DataType>::data_type[N];
  using const_data_type = typename data_analysis<DataType>::const_data_type[N];
  using runtime_data_type = typename data_analysis<DataType>::runtime_data_type*;
  using runtime_const_data_type = typename data_analysis<DataType>::runtime_const_data_type*;
  using non_const_data_type = typename data_analysis<DataType>::non_const_data_type[N];
};

template<class ViewType, class ViewTraitsType, class DataType, class Layout, class Space, class MemoryTraitsType,
         class HostMirrorSpace, class ValueType, class ReferenceType>
constexpr bool test_view_typedefs_impl() {
  // ========================
  // inherited from ViewTraits
  // ========================
  static_assert(std::is_same_v<typename ViewType::data_type, DataType>);
  static_assert(std::is_same_v<typename ViewType::const_data_type, typename data_analysis<DataType>::const_data_type>);
  static_assert(std::is_same_v<typename ViewType::non_const_data_type, typename data_analysis<DataType>::non_const_data_type>);
  
  // FIXME: these should be deprecated and for proper testing (I.e. where this is different from data_type)
  // we would need ensemble types which use the hidden View dimension facility of View (i.e. which make
  // "specialize" not void)
  static_assert(std::is_same_v<typename ViewType::scalar_array_type, DataType>);
  static_assert(std::is_same_v<typename ViewType::const_scalar_array_type, typename data_analysis<DataType>::const_data_type>);
  static_assert(std::is_same_v<typename ViewType::non_const_scalar_array_type, typename data_analysis<DataType>::non_const_data_type>);
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  static_assert(std::is_same_v<typename ViewType::specialize, void>);
#endif

  // FIXME: value_type definition conflicts with mdspan value_type
  static_assert(std::is_same_v<typename ViewType::value_type, ValueType>);
  static_assert(std::is_same_v<typename ViewType::const_value_type, const ValueType>);
  static_assert(std::is_same_v<typename ViewType::non_const_value_type, std::remove_const_t<ValueType>>);

  // FIXME: should maybe be deprecated
  static_assert(std::is_same_v<typename ViewType::array_layout, Layout>);

  // FIXME: should be deprecated and is some complicated impl type
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  static_assert(!std::is_void_v<typename ViewType::dimension>);
#endif

  static_assert(std::is_same_v<typename ViewType::execution_space, typename Space::execution_space>);
  static_assert(std::is_same_v<typename ViewType::memory_space, typename Space::memory_space>);
  static_assert(std::is_same_v<typename ViewType::device_type, Kokkos::Device<typename ViewType::execution_space, typename ViewType::memory_space>>);
  static_assert(std::is_same_v<typename ViewType::memory_traits, MemoryTraitsType>);
  static_assert(std::is_same_v<typename ViewType::host_mirror_space, HostMirrorSpace>);
  static_assert(std::is_same_v<typename ViewType::size_type, typename ViewType::memory_space::size_type>);
 
  // FIXME: should be deprecated in favor of reference
  static_assert(std::is_same_v<typename ViewType::reference_type, ReferenceType>);
  // FIXME: should be deprecated in favor of data_handle_type
  static_assert(std::is_same_v<typename ViewType::pointer_type, ValueType*>);
 
  // =========================================
  // in Legacy View: some helper View variants
  // =========================================
  static_assert(std::is_same_v<typename ViewType::traits, ViewTraitsType>);
  static_assert(std::is_same_v<typename ViewType::array_type,
                               Kokkos::View<typename ViewType::scalar_array_type, typename ViewType::array_layout,
                                            typename ViewType::device_type, typename ViewTraitsType::hooks_policy,
                                            typename ViewType::memory_traits>>);
  static_assert(std::is_same_v<typename ViewType::const_type,
                               Kokkos::View<typename ViewType::const_data_type, typename ViewType::array_layout,
                                            typename ViewType::device_type, typename ViewTraitsType::hooks_policy,
                                            typename ViewType::memory_traits>>);
  static_assert(std::is_same_v<typename ViewType::non_const_type,
                               Kokkos::View<typename ViewType::non_const_data_type, typename ViewType::array_layout,
                                            typename ViewType::device_type, typename ViewTraitsType::hooks_policy,
                                            typename ViewType::memory_traits>>);
  static_assert(std::is_same_v<typename ViewType::host_mirror_type,
                               Kokkos::View<typename ViewType::non_const_data_type, typename ViewType::array_layout,
                                            Kokkos::Device<Kokkos::DefaultHostExecutionSpace,
                                                           typename ViewType::host_mirror_space::memory_space>,
                                            typename ViewTraitsType::hooks_policy>>);

  using uniform_layout_type = std::conditional_t<ViewType::rank()==0 || (ViewType::rank()==0 &&
                                                 std::is_same_v<Layout, Kokkos::LayoutRight>),
                                                 Kokkos::LayoutLeft, Layout>;

  // FIXME: uniformtype removes all memorytraits?
  static_assert(std::is_same_v<typename ViewType::uniform_type,
                               Kokkos::View<typename ViewType::data_type, uniform_layout_type,
                                            typename ViewType::device_type, Kokkos::MemoryTraits<>>>);
  static_assert(std::is_same_v<typename ViewType::uniform_const_type,
                               Kokkos::View<typename ViewType::const_data_type, uniform_layout_type,
                                            typename ViewType::device_type, Kokkos::MemoryTraits<>>>);
  static_assert(std::is_same_v<typename ViewType::uniform_runtime_type,
                               Kokkos::View<typename data_analysis<DataType>::runtime_data_type, uniform_layout_type,
                                            typename ViewType::device_type, Kokkos::MemoryTraits<>>>);
  static_assert(std::is_same_v<typename ViewType::uniform_runtime_const_type,
                               Kokkos::View<typename data_analysis<DataType>::runtime_const_data_type, uniform_layout_type,
                                            typename ViewType::device_type, Kokkos::MemoryTraits<>>>);

  using anonymous_device_type = Kokkos::Device<typename ViewType::execution_space, Kokkos::AnonymousSpace>;
  static_assert(std::is_same_v<typename ViewType::uniform_nomemspace_type,
                               Kokkos::View<typename ViewType::data_type, uniform_layout_type,
                                            anonymous_device_type, Kokkos::MemoryTraits<>>>);
  static_assert(std::is_same_v<typename ViewType::uniform_const_nomemspace_type,
                               Kokkos::View<typename ViewType::const_data_type, uniform_layout_type,
                                            anonymous_device_type, Kokkos::MemoryTraits<>>>);
  static_assert(std::is_same_v<typename ViewType::uniform_runtime_nomemspace_type,
                               Kokkos::View<typename data_analysis<DataType>::runtime_data_type, uniform_layout_type,
                                            anonymous_device_type, Kokkos::MemoryTraits<>>>);
  static_assert(std::is_same_v<typename ViewType::uniform_runtime_const_nomemspace_type,
                               Kokkos::View<typename data_analysis<DataType>::runtime_const_data_type, uniform_layout_type,
                                            anonymous_device_type, Kokkos::MemoryTraits<>>>);


  // ==================================
  // mdspan compatibility
  // ==================================

  // FIXME: This typedef caused some weird issue with MSVC+NVCC
  // static_assert(std::is_same_v<typename ViewType::layout_type, Layout>);
  // FIXME: Not supported yet
  // static_assert(std::is_same_v<typename ViewType::extents_type, >);
  // static_assert(std::is_same_v<typename ViewType::mapping_type, >);
  // static_assert(std::is_same_v<typename ViewType::accessor_type, >);

  static_assert(std::is_same_v<typename ViewType::element_type, ValueType>);
  // FIXME: should be remove_const_t<element_type>
  static_assert(std::is_same_v<typename ViewType::value_type, ValueType>);
  static_assert(std::is_same_v<typename ViewType::size_type, typename Space::memory_space::size_type>);
  // FIXME: we need to evaluate how we want to proceed with this, as with
  // extents index_type also determines the stride, while LegacyView uses size_t strides
  // So we are doing this now to avoid breakage but it means we may use 64 bit indices on the GPU
  #ifndef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  static_assert(std::is_same_v<typename ViewType::index_type, size_t>);
  #endif
  // FIXME: this isn't given in View since for example SYCL has "int" as its size_type
  // static_assert(std::is_same_v<typename ViewType::size_type, std::make_unsigned_t<typename ViewType::index_type>>);
  static_assert(std::is_same_v<typename ViewType::rank_type, size_t>);

  // FIXME: should come from accessor_type
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  static_assert(std::is_same_v<typename ViewType::data_handle_type, typename ViewType::pointer_type>);
#endif
  static_assert(std::is_same_v<typename ViewType::reference, typename ViewType::reference_type>);
  return true;
}

// Helper function to unpack data type and other args from the View, and pass them on
template<class T, class ... ViewArgs>
struct ViewParams {};

template<class L, class S, class M, class HostMirrorSpace, class ValueType, class ReferenceType, class T, class ... ViewArgs>
constexpr bool test_view_typedefs(ViewParams<T, ViewArgs...>) {
  return test_view_typedefs_impl<Kokkos::View<T, ViewArgs...>, Kokkos::ViewTraits<T, ViewArgs...>,
                                 T, L, S, M, HostMirrorSpace, ValueType, ReferenceType>();
}


constexpr bool is_host_exec = std::is_same_v<Kokkos::DefaultExecutionSpace, Kokkos::DefaultHostExecutionSpace>;

#if defined(KOKKOS_ENABLE_CUDA_UVM) || defined(KOKKOS_ENABLE_IMPL_CUDA_UNIFIED_MEMORY) || defined(KOKKOS_IMPL_HIP_UNIFIED_MEMORY)
constexpr bool has_unified_mem_space = true;
#else
constexpr bool has_unified_mem_space = false;
#endif

// The test take explicit template arguments for: LayoutType, Space, MemoryTraits, HostMirrorSpace, ValueType, ReferenceType
// The ViewParams is just a type pack for the View template arguments

// Kokkos::View<int>
namespace TestInt {
  using layout_type = Kokkos::DefaultExecutionSpace::array_layout;
  using space = Kokkos::DefaultExecutionSpace;
  using memory_traits = Kokkos::MemoryTraits<>;
  // HostMirrorSpace is a mess so: if the default exec is a host exec, that is it
  using host_mirror_space = std::conditional_t<is_host_exec, Kokkos::DefaultExecutionSpace,
  // otherwise if unified memory is not on its HostSpace
                               std::conditional_t<!has_unified_mem_space, Kokkos::HostSpace,
  // otherwise its the following Device type
                               Kokkos::Device<Kokkos::DefaultHostExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory_space>>>;
  static_assert(test_view_typedefs<layout_type, space, memory_traits, host_mirror_space, int, int&>(
                     ViewParams<int>{}));
}

// Kokkos::View<int, DefaultExecutionSpace>
namespace TestIntDefaultExecutionSpace {
  using layout_type = Kokkos::DefaultExecutionSpace::array_layout;
  using space = Kokkos::DefaultExecutionSpace;
  using memory_traits = Kokkos::MemoryTraits<>;
  // HostMirrorSpace is a mess so: if the default exec is a host exec, it is HostSpace (note difference from View<int> ...)
  using host_mirror_space = std::conditional_t<is_host_exec, Kokkos::HostSpace,
  // otherwise if unified memory is not on its also HostSpace!
                               std::conditional_t<!has_unified_mem_space, Kokkos::HostSpace,
  // otherwise its the following memory space ...
                               Kokkos::DefaultExecutionSpace::memory_space>>;
  static_assert(test_view_typedefs<layout_type, space, memory_traits, host_mirror_space, int, int&>(
                     ViewParams<int, Kokkos::DefaultExecutionSpace>{}));
}

// Kokkos::View<const float**, Kokkos::HostSpace>
namespace TestFloatPPHostSpace {
  using layout_type = Kokkos::LayoutRight;
  using space = Kokkos::HostSpace;
  using memory_traits = Kokkos::MemoryTraits<>;
  using host_mirror_space = Kokkos::HostSpace;
  static_assert(test_view_typedefs<layout_type, space, memory_traits, host_mirror_space, const float, const float&>(
                     ViewParams<const float**, Kokkos::HostSpace>{}));
}

// Kokkos::View<float*[3], Kokkos::LayoutLeft>
namespace TestFloatP3LayoutLeft {
  using layout_type = Kokkos::LayoutLeft;
  using space = Kokkos::DefaultExecutionSpace;
  using memory_traits = Kokkos::MemoryTraits<>;
  // HostMirrorSpace is a mess so: if the default exec is a host exec, that is it
  using host_mirror_space = std::conditional_t<is_host_exec, Kokkos::DefaultExecutionSpace,
  // otherwise if unified memory is not on its HostSpace
                               std::conditional_t<!has_unified_mem_space, Kokkos::HostSpace,
  // otherwise its the following Device type
                               Kokkos::Device<Kokkos::DefaultHostExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory_space>>>;
  static_assert(test_view_typedefs<layout_type, space, memory_traits, host_mirror_space, float, float&>(
                     ViewParams<float*[3], Kokkos::LayoutLeft>{}));
}

// Kokkos::View<float[2][3], Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>>
namespace TestFloatPPDeviceDefaultHostExecHostSpace {
  using layout_type = Kokkos::LayoutRight;
  using space = Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;
  using memory_traits = Kokkos::MemoryTraits<>;
  using host_mirror_space = Kokkos::HostSpace;
  static_assert(test_view_typedefs<layout_type, space, memory_traits, host_mirror_space, float, float&>(
                     ViewParams<float[2][3], Kokkos::LayoutRight, Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>>{}));
}

// Kokkos::View<int, Kokkos::MemoryTraits<Kokkos::Atomic>>
namespace TestIntAtomic {
  using layout_type = Kokkos::DefaultExecutionSpace::array_layout;
  using space = Kokkos::DefaultExecutionSpace;
  using memory_traits = Kokkos::MemoryTraits<Kokkos::Atomic>;
  // HostMirrorSpace is a mess so: if the default exec is a host exec, that is it
  using host_mirror_space = std::conditional_t<is_host_exec, Kokkos::DefaultExecutionSpace,
  // otherwise if unified memory is not on its HostSpace
                               std::conditional_t<!has_unified_mem_space, Kokkos::HostSpace,
  // otherwise its the following Device type
                               Kokkos::Device<Kokkos::DefaultHostExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory_space>>>;
// clang-format on
static_assert(test_view_typedefs<
              layout_type, space, memory_traits, host_mirror_space, int,
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
              Kokkos::Impl::AtomicDataElement<
                  Kokkos::ViewTraits<int, Kokkos::MemoryTraits<Kokkos::Atomic>>>
#else
              desul::AtomicRef<int, desul::MemoryOrderRelaxed,
                               desul::MemoryScopeDevice>
#endif
              >(ViewParams<int, Kokkos::MemoryTraits<Kokkos::Atomic>>{}));
// clang-format off
}
// clang-format on
}  // namespace
