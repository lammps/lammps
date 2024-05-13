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

#ifndef KOKKOS_KOKKOS_TRAITS_FWD_HPP
#define KOKKOS_KOKKOS_TRAITS_FWD_HPP

// Without this the CUDA side does proper EBO while MSVC doesn't
// leading to mismatched sizes of the driver objects (CudaParallel)
// leading to illegal memory accesses etc on device
#if defined(_WIN32) && defined(KOKKOS_ENABLE_CUDA)
#define KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND char dummy;
#else
#define KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
#endif

namespace Kokkos {
namespace Impl {

template <class Enable, class... TraitsList>
struct AnalyzeExecPolicy;

template <class Enable, class TraitSpecList, class... Traits>
struct AnalyzeExecPolicyUseMatcher;

template <class AnalysisResults>
struct ExecPolicyTraitsWithDefaults;

template <class TraitSpec, class Trait, class Enable>
struct PolicyTraitMatcher;

template <class TraitSpec, template <class...> class PolicyTemplate,
          class AlreadyProcessedList, class ToProcessList, class NewTrait,
          class Enable = void>
struct PolicyTraitAdaptorImpl;

template <class TraitSpec, class Policy, class NewTrait>
struct PolicyTraitAdaptor;

// A tag class for dependent defaults that must be handled by the
// ExecPolicyTraitsWithDefaults wrapper, since their defaults depend on other
// traits
struct dependent_policy_trait_default;

//==============================================================================
// <editor-fold desc="Execution policy trait specifications"> {{{1

struct ExecutionSpaceTrait;
struct IndexTypeTrait;
struct ScheduleTrait;
struct IterationPatternTrait;
struct WorkItemPropertyTrait;
struct LaunchBoundsTrait;
struct OccupancyControlTrait;
struct GraphKernelTrait;
struct WorkTagTrait;

// Keep these sorted by frequency of use to reduce compilation time
//
// clang-format off
using execution_policy_trait_specifications =
  type_list<
    ExecutionSpaceTrait,
    IndexTypeTrait,
    ScheduleTrait,
    IterationPatternTrait,
    WorkItemPropertyTrait,
    LaunchBoundsTrait,
    OccupancyControlTrait,
    GraphKernelTrait,
    // This one has to be last, unfortunately:
    WorkTagTrait
  >;
// clang-format on

// </editor-fold> end Execution policy trait specifications }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_TRAITS_FWD_HPP
