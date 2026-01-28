// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_Core.hpp>
#include <KokkosExp_InterOp.hpp>

export module kokkos.core_impl;

export {
  namespace Kokkos {
  // execution/memory spaces
  namespace Impl {
  using ::Kokkos::Impl::ExecutionSpaceTag;
  using ::Kokkos::Impl::make_shared_allocation_record;
  using ::Kokkos::Impl::MemorySpaceAccess;
  }  // namespace Impl

  // View-related
  namespace Impl {
  namespace BV {
  using ::Kokkos::Impl::BV::BasicView;
  }
  using ::Kokkos::Impl::AccessorArg_t;
  using ::Kokkos::Impl::AccessorTypeTag;
  using ::Kokkos::Impl::append_formatted_multidimensional_index;
  using ::Kokkos::Impl::ApplyToViewOfStaticRank;
  using ::Kokkos::Impl::are_integral;
  using ::Kokkos::Impl::as_view_of_rank_n;
  using ::Kokkos::Impl::AtomicAccessorRelaxed;
  using ::Kokkos::Impl::check_view_ctor_args_create_mirror;
  using ::Kokkos::Impl::check_view_ctor_args_create_mirror_view_and_copy;
  using ::Kokkos::Impl::CheckedReferenceCountedAccessor;
  using ::Kokkos::Impl::CheckedReferenceCountedRelaxedAtomicAccessor;
  using ::Kokkos::Impl::CheckedRelaxedAtomicAccessor;
  using ::Kokkos::Impl::choose_create_mirror;
  using ::Kokkos::Impl::CommonSubview;
  using ::Kokkos::Impl::DataTypeFromExtents;
  using ::Kokkos::Impl::DeepCopy;
  using ::Kokkos::Impl::ExtentsFromDataType;
  using ::Kokkos::Impl::get_property;
  using ::Kokkos::Impl::has_type;
  using ::Kokkos::Impl::HostMirror;
  using ::Kokkos::Impl::is_view_ctor_property;
  using ::Kokkos::Impl::is_view_label;
  using ::Kokkos::Impl::LabelTag;
  using ::Kokkos::Impl::MDSpanViewTraits;
  using ::Kokkos::Impl::ParseViewExtents;
  using ::Kokkos::Impl::rank_dynamic;
  using ::Kokkos::Impl::RankDataType;
  using ::Kokkos::Impl::ReferenceCountedAccessor;
  using ::Kokkos::Impl::ReferenceCountedDataHandle;
  using ::Kokkos::Impl::runtime_check_memory_access_violation;
  using ::Kokkos::Impl::SharedAllocationHeader;
  using ::Kokkos::Impl::SharedAllocationRecord;
  using ::Kokkos::Impl::SharedAllocationTracker;
  using ::Kokkos::Impl::size_mismatch;
  using ::Kokkos::Impl::SpaceAwareAccessor;
  using ::Kokkos::Impl::SubviewExtents;
  using ::Kokkos::Impl::SubviewLegalArgsCompileTime;
  using ::Kokkos::Impl::ViewArguments;
  using ::Kokkos::Impl::ViewArrayAnalysis;
  using ::Kokkos::Impl::ViewCopy;
  using ::Kokkos::Impl::ViewCtorProp;
  using ::Kokkos::Impl::ViewCustomArguments;
  using ::Kokkos::Impl::ViewDataAnalysis;
  using ::Kokkos::Impl::ViewDataHandle;
  using ::Kokkos::Impl::ViewDimension;
  using ::Kokkos::Impl::ViewMapping;
  using ::Kokkos::Impl::ViewOffset;
  using ::Kokkos::Impl::ViewRemap;
  using ::Kokkos::Impl::with_properties_if_unset;
  using ::Kokkos::Impl::WithoutInitializing_t;
  }  // namespace Impl

  // execution policies
  namespace Impl {
  using ::Kokkos::Impl::get_tile_size_properties;
  using ::Kokkos::Impl::ParallelConstructName;
  using ::Kokkos::Impl::PolicyTraits;
  using ::Kokkos::Impl::PolicyUpdate;
  using ::Kokkos::Impl::WorkTagTrait;
  }  // namespace Impl

  // miscellaneous
  namespace Impl {
  using ::Kokkos::Impl::FunctorAnalysis;
  using ::Kokkos::Impl::python_view_type_impl_t;
  using ::Kokkos::Impl::throw_runtime_exception;
  }  // namespace Impl

  // initialization/finalization
  namespace Impl {
  using ::Kokkos::Impl::post_finalize;
  using ::Kokkos::Impl::post_initialize;
  using ::Kokkos::Impl::pre_finalize;
  using ::Kokkos::Impl::pre_initialize;
  }  // namespace Impl
  }  // namespace Kokkos
}
