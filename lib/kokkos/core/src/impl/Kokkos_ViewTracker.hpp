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

#ifndef KOKKOS_VIEW_TRACKER_HPP
#define KOKKOS_VIEW_TRACKER_HPP

namespace Kokkos {

template <class DataType, class... Properties>
class View;

namespace Impl {

/*
 * \class ViewTracker
 * \brief template class to wrap the shared allocation tracker
 *
 * \section This class is templated on the View and provides
 * constructors that match the view.  The constructors and assignments
 * from view will externalize the logic needed to enable/disable
 * ref counting to provide a single gate to enable further developments
 * which may hing on the same logic.
 *
 */
template <class ParentView>
struct ViewTracker {
  using track_type  = Kokkos::Impl::SharedAllocationTracker;
  using view_traits = typename ParentView::traits;

  track_type m_tracker;

  KOKKOS_INLINE_FUNCTION
  ViewTracker() : m_tracker() {}

  KOKKOS_INLINE_FUNCTION
  ViewTracker(const ViewTracker& vt) noexcept
      : m_tracker(vt.m_tracker, view_traits::is_managed) {}

  KOKKOS_INLINE_FUNCTION
  explicit ViewTracker(const ParentView& vt) noexcept : m_tracker() {
    assign(vt);
  }

  template <class RT, class... RP>
  KOKKOS_INLINE_FUNCTION explicit ViewTracker(
      const View<RT, RP...>& vt) noexcept
      : m_tracker() {
    assign(vt);
  }

  template <class RT, class... RP>
  KOKKOS_INLINE_FUNCTION void assign(const View<RT, RP...>& vt) noexcept {
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    if (view_traits::is_managed &&
        Kokkos::Impl::SharedAllocationRecord<void, void>::tracking_enabled()) {
      m_tracker.assign_direct(vt.m_track.m_tracker);
    } else {
      m_tracker.assign_force_disable(vt.m_track.m_tracker);
    }
#else
    m_tracker.assign_force_disable(vt.m_track.m_tracker);
#endif
  }

  KOKKOS_INLINE_FUNCTION
  ViewTracker& operator=(const ViewTracker& rhs) noexcept {
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    if (view_traits::is_managed &&
        Kokkos::Impl::SharedAllocationRecord<void, void>::tracking_enabled()) {
      m_tracker.assign_direct(rhs.m_tracker);
    } else {
      m_tracker.assign_force_disable(rhs.m_tracker);
    }
#else
    m_tracker.assign_force_disable(rhs.m_tracker);
#endif
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  explicit ViewTracker(const track_type& tt) noexcept
      : m_tracker(tt, view_traits::is_managed) {}
};

}  // namespace Impl

}  // namespace Kokkos

#endif  // KOKKOS_VIEW_TRACKER_HPP
