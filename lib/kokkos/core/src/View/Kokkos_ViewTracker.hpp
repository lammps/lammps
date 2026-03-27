// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
 * which may hinge on the same logic.
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
      : m_tracker(vt.m_tracker, !view_traits::memory_traits::is_unmanaged) {}

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
  KOKKOS_INLINE_FUNCTION void assign(const View<RT, RP...>& vt) {
    if (this == reinterpret_cast<const ViewTracker*>(&vt.m_track)) return;
    KOKKOS_IF_ON_HOST((
        if (!view_traits::memory_traits::is_unmanaged &&
            Kokkos::Impl::SharedAllocationRecord<void,
                                                 void>::tracking_enabled()) {
          m_tracker.assign_direct(vt.m_track.m_tracker);
        } else { m_tracker.assign_force_disable(vt.m_track.m_tracker); }))

    KOKKOS_IF_ON_DEVICE((m_tracker.assign_force_disable(vt.m_track.m_tracker);))
  }

  KOKKOS_INLINE_FUNCTION ViewTracker& operator=(const ViewTracker& rhs) {
    if (this == &rhs) return *this;
    KOKKOS_IF_ON_HOST((
        if (!view_traits::memory_traits::is_unmanaged &&
            Kokkos::Impl::SharedAllocationRecord<void,
                                                 void>::tracking_enabled()) {
          m_tracker.assign_direct(rhs.m_tracker);
        } else { m_tracker.assign_force_disable(rhs.m_tracker); }))

    KOKKOS_IF_ON_DEVICE((m_tracker.assign_force_disable(rhs.m_tracker);))
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  explicit ViewTracker(const track_type& tt) noexcept
      : m_tracker(tt, !view_traits::memory_traits::is_unmanaged) {}
};

}  // namespace Impl

}  // namespace Kokkos

#endif  // KOKKOS_VIEW_TRACKER_HPP
