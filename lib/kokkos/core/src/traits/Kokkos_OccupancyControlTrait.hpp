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
  };
  template <class T>
  using trait_matches_specification = std::integral_constant<
      bool,
      std::is_same<T, Kokkos::Experimental::DesiredOccupancy>::value ||
          std::is_same<T, Kokkos::Experimental::MaximizeOccupancy>::value>;
};

// </editor-fold> end Occupancy control trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="AnalyzeExecPolicy specializations"> {{{1

// The DesiredOccupancy case has runtime storage, so we need to handle copies
// and assignments
template <class... Traits>
struct AnalyzeExecPolicy<void, Kokkos::Experimental::DesiredOccupancy,
                         Traits...> : AnalyzeExecPolicy<void, Traits...> {
 public:
  using base_t            = AnalyzeExecPolicy<void, Traits...>;
  using occupancy_control = Kokkos::Experimental::DesiredOccupancy;
  static constexpr bool experimental_contains_desired_occupancy = true;

  template <class OccControl>
  using with_occupancy_control = AnalyzeExecPolicy<void, OccControl, Traits...>;

  // Treat this as private, but make it public so that MSVC will still treat
  // this as a standard layout class and make it the right size: storage for a
  // stateful desired occupancy
  //   private:
  occupancy_control m_desired_occupancy;

  AnalyzeExecPolicy() = default;
  // Converting constructor
  // Just rely on the convertibility of occupancy_control to transfer the data
  template <class Other>
  AnalyzeExecPolicy(ExecPolicyTraitsWithDefaults<Other> const& other)
      : base_t(other),
        m_desired_occupancy(other.impl_get_occupancy_control()) {}

  // Converting assignment operator
  // Just rely on the convertibility of occupancy_control to transfer the data
  template <class Other>
  AnalyzeExecPolicy& operator=(
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

template <class... Traits>
struct AnalyzeExecPolicy<void, Kokkos::Experimental::MaximizeOccupancy,
                         Traits...> : AnalyzeExecPolicy<void, Traits...> {
  using base_t = AnalyzeExecPolicy<void, Traits...>;
  using base_t::base_t;
  using occupancy_control = Kokkos::Experimental::MaximizeOccupancy;
  static constexpr bool experimental_contains_desired_occupancy = false;
};

// </editor-fold> end AnalyzeExecPolicy specializations }}}1
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
