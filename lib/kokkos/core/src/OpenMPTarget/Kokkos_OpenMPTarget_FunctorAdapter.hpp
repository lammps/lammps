// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENMPTARGET_FUNCTOR_ADAPTER_HPP
#define KOKKOS_OPENMPTARGET_FUNCTOR_ADAPTER_HPP

#include <OpenMPTarget/Kokkos_OpenMPTarget_Macros.hpp>
#include <type_traits>

namespace Kokkos::Experimental::Impl {

template <class Functor, class Policy>
class FunctorAdapter {
  Functor m_functor;
  using WorkTag = typename Policy::work_tag;

 public:
  FunctorAdapter() = default;
  FunctorAdapter(Functor const &functor) : m_functor(functor) {}

  Functor get_functor() const { return m_functor; }

  template <class... Args>
  KOKKOS_FUNCTION void operator()(Args &&...args) const {
    if constexpr (std::is_void_v<WorkTag>) {
      m_functor(static_cast<Args &&>(args)...);
    } else {
      m_functor(WorkTag(), static_cast<Args &&>(args)...);
    }
  }
};

}  // namespace Kokkos::Experimental::Impl

#endif  // KOKKOS_OPENMPTARGET_FUNCTOR_ADAPTER_HPP
