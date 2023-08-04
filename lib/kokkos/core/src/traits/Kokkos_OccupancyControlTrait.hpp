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

#ifndef KOKKOS_KOKKOS_OCCUPANCYCONTROLTRAIT_HPP
#define KOKKOS_KOKKOS_OCCUPANCYCONTROLTRAIT_HPP

#include <impl/Kokkos_Error.hpp>  // KOKKOS_EXPECTS macro

#include <traits/Kokkos_PolicyTraitAdaptor.hpp>

#include <traits/Kokkos_Traits_fwd.hpp>

namespace Kokkos {

namespace Experimental {

//==============================================================================
// <editor-fold desc="Occupancy control user interface"> {{{1

struct MaximizeOccupancy;

struct DesiredOccupancy {
  int m_occ = 100;
  explicit constexpr DesiredOccupancy(int occ) : m_occ(occ) {
    KOKKOS_EXPECTS(0 <= occ && occ <= 100);
  }
  explicit constexpr operator int() const { return m_occ; }
  constexpr int value() const { return m_occ; }
  DesiredOccupancy() = default;
  explicit DesiredOccupancy(MaximizeOccupancy const&) : DesiredOccupancy() {}
};

struct MaximizeOccupancy {
  explicit MaximizeOccupancy() = default;
};

// </editor-fold> end Occupancy control user interface }}}1
//==============================================================================

}  // end namespace Experimental

namespace Impl {

template <class Policy, class AnalyzeNextTrait>
struct OccupancyControlPolicyMixin;

//==============================================================================
// <editor-fold desc="Occupancy control trait specification"> {{{1

struct OccupancyControlTrait : TraitSpecificationBase<OccupancyControlTrait> {
  struct base_traits {
    using occupancy_control = Kokkos::Experimental::MaximizeOccupancy;
    static constexpr bool experimental_contains_desired_occupancy = false;
    // Default access occupancy_control, for when it is the (stateless) default
    static constexpr occupancy_control impl_get_occupancy_control() {
      return occupancy_control{};
    }
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class OccControl, class AnalyzeNextTrait>
  using mixin_matching_trait =
      OccupancyControlPolicyMixin<OccControl, AnalyzeNextTrait>;
  template <class T>
  using trait_matches_specification = std::bool_constant<
      std::is_same<T, Kokkos::Experimental::DesiredOccupancy>::value ||
      std::is_same<T, Kokkos::Experimental::MaximizeOccupancy>::value>;
};

// </editor-fold> end Occupancy control trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="OccupancyControlPolicyMixin specializations"> {{{1

template <class AnalyzeNextTrait>
struct OccupancyControlPolicyMixin<Kokkos::Experimental::DesiredOccupancy,
                                   AnalyzeNextTrait> : AnalyzeNextTrait {
  using base_t            = AnalyzeNextTrait;
  using occupancy_control = Kokkos::Experimental::DesiredOccupancy;
  static constexpr bool experimental_contains_desired_occupancy = true;

  // Treat this as private, but make it public so that MSVC will still treat
  // this as a standard layout class and make it the right size: storage for a
  // stateful desired occupancy
  //   private:
  occupancy_control m_desired_occupancy = occupancy_control{};

  OccupancyControlPolicyMixin() = default;
  // Converting constructor
  // Just rely on the convertibility of occupancy_control to transfer the data
  template <class Other>
  OccupancyControlPolicyMixin(ExecPolicyTraitsWithDefaults<Other> const& other)
      : base_t(other),
        m_desired_occupancy(other.impl_get_occupancy_control()) {}

  // Converting assignment operator
  // Just rely on the convertibility of occupancy_control to transfer the data
  template <class Other>
  OccupancyControlPolicyMixin& operator=(
      ExecPolicyTraitsWithDefaults<Other> const& other) {
    *static_cast<base_t*>(this) = other;
    this->impl_set_desired_occupancy(
        occupancy_control{other.impl_get_occupancy_control()});
    return *this;
  }

  // Access to occupancy control instance, usable in generic context
  constexpr occupancy_control impl_get_occupancy_control() const {
    return m_desired_occupancy;
  }

  // Access to desired occupancy (getter and setter)
  Kokkos::Experimental::DesiredOccupancy impl_get_desired_occupancy() const {
    return m_desired_occupancy;
  }

  void impl_set_desired_occupancy(occupancy_control desired_occupancy) {
    m_desired_occupancy = desired_occupancy;
  }
};

template <class AnalyzeNextTrait>
struct OccupancyControlPolicyMixin<Kokkos::Experimental::MaximizeOccupancy,
                                   AnalyzeNextTrait> : AnalyzeNextTrait {
  using base_t = AnalyzeNextTrait;
  using base_t::base_t;
  using occupancy_control = Kokkos::Experimental::MaximizeOccupancy;
  static constexpr bool experimental_contains_desired_occupancy = false;
};

// </editor-fold> end OccupancyControlPolicyMixin specializations }}}1
//==============================================================================

}  // end namespace Impl

namespace Experimental {

//==============================================================================
// <editor-fold desc="User interface"> {{{1

template <typename Policy>
auto prefer(Policy const& p, DesiredOccupancy occ) {
  using new_policy_t =
      Kokkos::Impl::OccupancyControlTrait::policy_with_trait<Policy,
                                                             DesiredOccupancy>;
  new_policy_t pwo{p};
  pwo.impl_set_desired_occupancy(occ);
  return pwo;
}

template <typename Policy>
constexpr auto prefer(Policy const& p, MaximizeOccupancy) {
  static_assert(Kokkos::is_execution_policy<Policy>::value, "");
  using new_policy_t =
      Kokkos::Impl::OccupancyControlTrait::policy_with_trait<Policy,
                                                             MaximizeOccupancy>;
  return new_policy_t{p};
}

// </editor-fold> end User interface }}}1
//==============================================================================

}  // end namespace Experimental

}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_OCCUPANCYCONTROLTRAIT_HPP
