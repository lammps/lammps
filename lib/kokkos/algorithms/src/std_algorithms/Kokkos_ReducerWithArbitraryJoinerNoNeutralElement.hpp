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

#ifndef KOKKOS_STD_ReducerWithArbitraryJoinerNoNeutralElement_hpp_
#define KOKKOS_STD_ReducerWithArbitraryJoinerNoNeutralElement_hpp_

#include <Kokkos_Core.hpp>
#include "Kokkos_ValueWrapperForNoNeutralElement.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

// This reducer is here and not where all other reducers are
// because it is inside Impl and also because it would not work
// for OpenMPTarget backend. We can move this later.

template <class Scalar, class JoinerType, class Space = HostSpace>
struct ReducerWithArbitraryJoinerNoNeutralElement {
  using scalar_type = typename std::remove_cv<Scalar>::type;

 public:
  // Required
  using reducer =
      ReducerWithArbitraryJoinerNoNeutralElement<Scalar, JoinerType, Space>;
  using value_type = ValueWrapperForNoNeutralElement<scalar_type>;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  JoinerType m_joiner;
  result_view_type m_value;
  bool m_references_scalar_v;

 public:
  KOKKOS_FUNCTION
  ReducerWithArbitraryJoinerNoNeutralElement(value_type& value_,
                                             JoinerType joiner_)
      : m_joiner(joiner_), m_value(&value_), m_references_scalar_v(true) {}

  KOKKOS_FUNCTION
  ReducerWithArbitraryJoinerNoNeutralElement(const result_view_type& value_,
                                             JoinerType joiner_)
      : m_joiner(joiner_), m_value(value_), m_references_scalar_v(false) {}

  // Required
  KOKKOS_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest.val = m_joiner(dest.val, src.val);
  }

  KOKKOS_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest.val = m_joiner(dest.val, src.val);
  }

  KOKKOS_FUNCTION
  void init(value_type& val) const {
    // I cannot call reduction_identity, so need to default this
    val = {};
  }

  KOKKOS_FUNCTION
  value_type& reference() const { return *m_value.data(); }

  KOKKOS_FUNCTION
  result_view_type view() const { return m_value; }

  KOKKOS_FUNCTION
  bool references_scalar() const { return m_references_scalar_v; }
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
