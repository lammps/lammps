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

#ifndef KOKKOS_COPYVIEWS_HPP_
#define KOKKOS_COPYVIEWS_HPP_
#include <string>
#include <Kokkos_Parallel.hpp>
#include <KokkosExp_MDRangePolicy.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template <class Layout>
struct ViewFillLayoutSelector {};

template <>
struct ViewFillLayoutSelector<Kokkos::LayoutLeft> {
  static const Kokkos::Iterate iterate = Kokkos::Iterate::Left;
};

template <>
struct ViewFillLayoutSelector<Kokkos::LayoutRight> {
  static const Kokkos::Iterate iterate = Kokkos::Iterate::Right;
};

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

template <class ViewType, class Layout, class ExecSpace, typename iType>
struct ViewFill<ViewType, Layout, ExecSpace, 0, iType> {
  using ST = typename ViewType::non_const_value_type;
  ViewFill(const ViewType& a, const ST& val, const ExecSpace& space) {
    Kokkos::Impl::DeepCopy<typename ViewType::memory_space, Kokkos::HostSpace,
                           ExecSpace>(space, a.data(), &val, sizeof(ST));
  }
};

template <class ViewType, class Layout, class ExecSpace, typename iType>
struct ViewFill<ViewType, Layout, ExecSpace, 1, iType> {
  ViewType a;
  typename ViewType::const_value_type val;
  using policy_type = Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<iType>>;

  ViewFill(const ViewType& a_, typename ViewType::const_value_type& val_,
           const ExecSpace& space)
      : a(a_), val(val_) {
    Kokkos::parallel_for("Kokkos::ViewFill-1D",
                         policy_type(space, 0, a.extent(0)), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i) const { a(i) = val; };
};

template <class ViewType, class Layout, class ExecSpace, typename iType>
struct ViewFill<ViewType, Layout, ExecSpace, 2, iType> {
  ViewType a;
  typename ViewType::const_value_type val;

  using iterate_type = Kokkos::Rank<2, ViewFillLayoutSelector<Layout>::iterate,
                                    ViewFillLayoutSelector<Layout>::iterate>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewFill(const ViewType& a_, typename ViewType::const_value_type& val_,
           const ExecSpace& space)
      : a(a_), val(val_) {
    Kokkos::parallel_for("Kokkos::ViewFill-2D",
                         policy_type(space, {0, 0}, {a.extent(0), a.extent(1)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1) const { a(i0, i1) = val; };
};

template <class ViewType, class Layout, class ExecSpace, typename iType>
struct ViewFill<ViewType, Layout, ExecSpace, 3, iType> {
  ViewType a;
  typename ViewType::const_value_type val;

  using iterate_type = Kokkos::Rank<3, ViewFillLayoutSelector<Layout>::iterate,
                                    ViewFillLayoutSelector<Layout>::iterate>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewFill(const ViewType& a_, typename ViewType::const_value_type& val_,
           const ExecSpace& space)
      : a(a_), val(val_) {
    Kokkos::parallel_for(
        "Kokkos::ViewFill-3D",
        policy_type(space, {0, 0, 0}, {a.extent(0), a.extent(1), a.extent(2)}),
        *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i2) const {
    a(i0, i1, i2) = val;
  };
};

template <class ViewType, class Layout, class ExecSpace, typename iType>
struct ViewFill<ViewType, Layout, ExecSpace, 4, iType> {
  ViewType a;
  typename ViewType::const_value_type val;

  using iterate_type = Kokkos::Rank<4, ViewFillLayoutSelector<Layout>::iterate,
                                    ViewFillLayoutSelector<Layout>::iterate>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewFill(const ViewType& a_, typename ViewType::const_value_type& val_,
           const ExecSpace& space)
      : a(a_), val(val_) {
    Kokkos::parallel_for(
        "Kokkos::ViewFill-4D",
        policy_type(space, {0, 0, 0, 0},
                    {a.extent(0), a.extent(1), a.extent(2), a.extent(3)}),
        *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i2,
                  const iType& i3) const {
    a(i0, i1, i2, i3) = val;
  };
};

template <class ViewType, class Layout, class ExecSpace, typename iType>
struct ViewFill<ViewType, Layout, ExecSpace, 5, iType> {
  ViewType a;
  typename ViewType::const_value_type val;

  using iterate_type = Kokkos::Rank<5, ViewFillLayoutSelector<Layout>::iterate,
                                    ViewFillLayoutSelector<Layout>::iterate>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewFill(const ViewType& a_, typename ViewType::const_value_type& val_,
           const ExecSpace& space)
      : a(a_), val(val_) {
    Kokkos::parallel_for("Kokkos::ViewFill-5D",
                         policy_type(space, {0, 0, 0, 0, 0},
                                     {a.extent(0), a.extent(1), a.extent(2),
                                      a.extent(3), a.extent(4)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i2,
                  const iType& i3, const iType& i4) const {
    a(i0, i1, i2, i3, i4) = val;
  };
};

template <class ViewType, class Layout, class ExecSpace, typename iType>
struct ViewFill<ViewType, Layout, ExecSpace, 6, iType> {
  ViewType a;
  typename ViewType::const_value_type val;

  using iterate_type = Kokkos::Rank<6, ViewFillLayoutSelector<Layout>::iterate,
                                    ViewFillLayoutSelector<Layout>::iterate>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewFill(const ViewType& a_, typename ViewType::const_value_type& val_,
           const ExecSpace& space)
      : a(a_), val(val_) {
    Kokkos::parallel_for("Kokkos::ViewFill-6D",
                         policy_type(space, {0, 0, 0, 0, 0, 0},
                                     {a.extent(0), a.extent(1), a.extent(2),
                                      a.extent(3), a.extent(4), a.extent(5)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i2,
                  const iType& i3, const iType& i4, const iType& i5) const {
    a(i0, i1, i2, i3, i4, i5) = val;
  };
};

template <class ViewType, class Layout, class ExecSpace, typename iType>
struct ViewFill<ViewType, Layout, ExecSpace, 7, iType> {
  ViewType a;
  typename ViewType::const_value_type val;

  using iterate_type = Kokkos::Rank<6, ViewFillLayoutSelector<Layout>::iterate,
                                    ViewFillLayoutSelector<Layout>::iterate>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewFill(const ViewType& a_, typename ViewType::const_value_type& val_,
           const ExecSpace& space)
      : a(a_), val(val_) {
    Kokkos::parallel_for("Kokkos::ViewFill-7D",
                         policy_type(space, {0, 0, 0, 0, 0, 0},
                                     {a.extent(0), a.extent(1), a.extent(2),
                                      a.extent(3), a.extent(5), a.extent(6)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i3,
                  const iType& i4, const iType& i5, const iType& i6) const {
    for (iType i2 = 0; i2 < iType(a.extent(2)); i2++)
      a(i0, i1, i2, i3, i4, i5, i6) = val;
  };
};

template <class ViewType, class Layout, class ExecSpace, typename iType>
struct ViewFill<ViewType, Layout, ExecSpace, 8, iType> {
  ViewType a;
  typename ViewType::const_value_type val;

  using iterate_type = Kokkos::Rank<6, ViewFillLayoutSelector<Layout>::iterate,
                                    ViewFillLayoutSelector<Layout>::iterate>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewFill(const ViewType& a_, typename ViewType::const_value_type& val_,
           const ExecSpace& space)
      : a(a_), val(val_) {
    Kokkos::parallel_for("Kokkos::ViewFill-8D",
                         policy_type(space, {0, 0, 0, 0, 0, 0},
                                     {a.extent(0), a.extent(1), a.extent(3),
                                      a.extent(5), a.extent(6), a.extent(7)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i3,
                  const iType& i5, const iType& i6, const iType& i7) const {
    for (iType i2 = 0; i2 < iType(a.extent(2)); i2++)
      for (iType i4 = 0; i4 < iType(a.extent(4)); i4++)
        a(i0, i1, i2, i3, i4, i5, i6, i7) = val;
  };
};

template <class ViewTypeA, class ViewTypeB, class Layout, class ExecSpace,
          typename iType>
struct ViewCopy<ViewTypeA, ViewTypeB, Layout, ExecSpace, 1, iType> {
  ViewTypeA a;
  ViewTypeB b;

  using policy_type = Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<iType>>;
  using value_type  = typename ViewTypeA::value_type;

  ViewCopy(const ViewTypeA& a_, const ViewTypeB& b_,
           const ExecSpace space = ExecSpace())
      : a(a_), b(b_) {
    Kokkos::parallel_for("Kokkos::ViewCopy-1D",
                         policy_type(space, 0, a.extent(0)), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0) const {
    a(i0) = static_cast<value_type>(b(i0));
  };
};

template <class ViewTypeA, class ViewTypeB, class Layout, class ExecSpace,
          typename iType>
struct ViewCopy<ViewTypeA, ViewTypeB, Layout, ExecSpace, 2, iType> {
  ViewTypeA a;
  ViewTypeB b;
  static const Kokkos::Iterate outer_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::outer_iteration_pattern;
  static const Kokkos::Iterate inner_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::inner_iteration_pattern;
  using iterate_type =
      Kokkos::Rank<2, outer_iteration_pattern, inner_iteration_pattern>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;
  using value_type = typename ViewTypeA::value_type;

  ViewCopy(const ViewTypeA& a_, const ViewTypeB& b_,
           const ExecSpace space = ExecSpace())
      : a(a_), b(b_) {
    Kokkos::parallel_for("Kokkos::ViewCopy-2D",
                         policy_type(space, {0, 0}, {a.extent(0), a.extent(1)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1) const {
    a(i0, i1) = static_cast<value_type>(b(i0, i1));
  };
};

template <class ViewTypeA, class ViewTypeB, class Layout, class ExecSpace,
          typename iType>
struct ViewCopy<ViewTypeA, ViewTypeB, Layout, ExecSpace, 3, iType> {
  ViewTypeA a;
  ViewTypeB b;

  static const Kokkos::Iterate outer_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::outer_iteration_pattern;
  static const Kokkos::Iterate inner_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::inner_iteration_pattern;
  using iterate_type =
      Kokkos::Rank<3, outer_iteration_pattern, inner_iteration_pattern>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;
  using value_type = typename ViewTypeA::value_type;

  ViewCopy(const ViewTypeA& a_, const ViewTypeB& b_,
           const ExecSpace space = ExecSpace())
      : a(a_), b(b_) {
    Kokkos::parallel_for(
        "Kokkos::ViewCopy-3D",
        policy_type(space, {0, 0, 0}, {a.extent(0), a.extent(1), a.extent(2)}),
        *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i2) const {
    a(i0, i1, i2) = static_cast<value_type>(b(i0, i1, i2));
  };
};

template <class ViewTypeA, class ViewTypeB, class Layout, class ExecSpace,
          typename iType>
struct ViewCopy<ViewTypeA, ViewTypeB, Layout, ExecSpace, 4, iType> {
  ViewTypeA a;
  ViewTypeB b;

  static const Kokkos::Iterate outer_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::outer_iteration_pattern;
  static const Kokkos::Iterate inner_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::inner_iteration_pattern;
  using iterate_type =
      Kokkos::Rank<4, outer_iteration_pattern, inner_iteration_pattern>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewCopy(const ViewTypeA& a_, const ViewTypeB& b_,
           const ExecSpace space = ExecSpace())
      : a(a_), b(b_) {
    Kokkos::parallel_for(
        "Kokkos::ViewCopy-4D",
        policy_type(space, {0, 0, 0, 0},
                    {a.extent(0), a.extent(1), a.extent(2), a.extent(3)}),
        *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i2,
                  const iType& i3) const {
    a(i0, i1, i2, i3) = b(i0, i1, i2, i3);
  };
};

template <class ViewTypeA, class ViewTypeB, class Layout, class ExecSpace,
          typename iType>
struct ViewCopy<ViewTypeA, ViewTypeB, Layout, ExecSpace, 5, iType> {
  ViewTypeA a;
  ViewTypeB b;

  static const Kokkos::Iterate outer_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::outer_iteration_pattern;
  static const Kokkos::Iterate inner_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::inner_iteration_pattern;
  using iterate_type =
      Kokkos::Rank<5, outer_iteration_pattern, inner_iteration_pattern>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewCopy(const ViewTypeA& a_, const ViewTypeB& b_,
           const ExecSpace space = ExecSpace())
      : a(a_), b(b_) {
    Kokkos::parallel_for("Kokkos::ViewCopy-5D",
                         policy_type(space, {0, 0, 0, 0, 0},
                                     {a.extent(0), a.extent(1), a.extent(2),
                                      a.extent(3), a.extent(4)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i2,
                  const iType& i3, const iType& i4) const {
    a(i0, i1, i2, i3, i4) = b(i0, i1, i2, i3, i4);
  };
};

template <class ViewTypeA, class ViewTypeB, class Layout, class ExecSpace,
          typename iType>
struct ViewCopy<ViewTypeA, ViewTypeB, Layout, ExecSpace, 6, iType> {
  ViewTypeA a;
  ViewTypeB b;

  static const Kokkos::Iterate outer_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::outer_iteration_pattern;
  static const Kokkos::Iterate inner_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::inner_iteration_pattern;
  using iterate_type =
      Kokkos::Rank<6, outer_iteration_pattern, inner_iteration_pattern>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewCopy(const ViewTypeA& a_, const ViewTypeB& b_,
           const ExecSpace space = ExecSpace())
      : a(a_), b(b_) {
    Kokkos::parallel_for("Kokkos::ViewCopy-6D",
                         policy_type(space, {0, 0, 0, 0, 0, 0},
                                     {a.extent(0), a.extent(1), a.extent(2),
                                      a.extent(3), a.extent(4), a.extent(5)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i2,
                  const iType& i3, const iType& i4, const iType& i5) const {
    a(i0, i1, i2, i3, i4, i5) = b(i0, i1, i2, i3, i4, i5);
  };
};

template <class ViewTypeA, class ViewTypeB, class Layout, class ExecSpace,
          typename iType>
struct ViewCopy<ViewTypeA, ViewTypeB, Layout, ExecSpace, 7, iType> {
  ViewTypeA a;
  ViewTypeB b;

  static const Kokkos::Iterate outer_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::outer_iteration_pattern;
  static const Kokkos::Iterate inner_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::inner_iteration_pattern;
  using iterate_type =
      Kokkos::Rank<6, outer_iteration_pattern, inner_iteration_pattern>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewCopy(const ViewTypeA& a_, const ViewTypeB& b_,
           const ExecSpace space = ExecSpace())
      : a(a_), b(b_) {
    Kokkos::parallel_for("Kokkos::ViewCopy-7D",
                         policy_type(space, {0, 0, 0, 0, 0, 0},
                                     {a.extent(0), a.extent(1), a.extent(3),
                                      a.extent(4), a.extent(5), a.extent(6)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i3,
                  const iType& i4, const iType& i5, const iType& i6) const {
    for (iType i2 = 0; i2 < iType(a.extent(2)); i2++)
      a(i0, i1, i2, i3, i4, i5, i6) = b(i0, i1, i2, i3, i4, i5, i6);
  };
};

template <class ViewTypeA, class ViewTypeB, class Layout, class ExecSpace,
          typename iType>
struct ViewCopy<ViewTypeA, ViewTypeB, Layout, ExecSpace, 8, iType> {
  ViewTypeA a;
  ViewTypeB b;

  static const Kokkos::Iterate outer_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::outer_iteration_pattern;
  static const Kokkos::Iterate inner_iteration_pattern =
      Kokkos::layout_iterate_type_selector<Layout>::inner_iteration_pattern;
  using iterate_type =
      Kokkos::Rank<6, outer_iteration_pattern, inner_iteration_pattern>;
  using policy_type =
      Kokkos::MDRangePolicy<ExecSpace, iterate_type, Kokkos::IndexType<iType>>;

  ViewCopy(const ViewTypeA& a_, const ViewTypeB& b_,
           const ExecSpace space = ExecSpace())
      : a(a_), b(b_) {
    Kokkos::parallel_for("Kokkos::ViewCopy-8D",
                         policy_type(space, {0, 0, 0, 0, 0, 0},
                                     {a.extent(0), a.extent(1), a.extent(3),
                                      a.extent(5), a.extent(6), a.extent(7)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const iType& i0, const iType& i1, const iType& i3,
                  const iType& i5, const iType& i6, const iType& i7) const {
    for (iType i2 = 0; i2 < iType(a.extent(2)); i2++)
      for (iType i4 = 0; i4 < iType(a.extent(4)); i4++)
        a(i0, i1, i2, i3, i4, i5, i6, i7) = b(i0, i1, i2, i3, i4, i5, i6, i7);
  };
};

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

template <class ExecutionSpace, class DstType, class SrcType>
void view_copy(const ExecutionSpace& space, const DstType& dst,
               const SrcType& src) {
  using dst_memory_space = typename DstType::memory_space;
  using src_memory_space = typename SrcType::memory_space;

  enum {
    ExecCanAccessSrc =
        Kokkos::Impl::SpaceAccessibility<ExecutionSpace,
                                         src_memory_space>::accessible
  };
  enum {
    ExecCanAccessDst =
        Kokkos::Impl::SpaceAccessibility<ExecutionSpace,
                                         dst_memory_space>::accessible
  };

  if (!(ExecCanAccessSrc && ExecCanAccessDst)) {
    Kokkos::Impl::throw_runtime_exception(
        "Kokkos::Impl::view_copy called with invalid execution space");
  } else {
    // Figure out iteration order in case we need it
    int64_t strides[DstType::Rank + 1];
    dst.stride(strides);
    Kokkos::Iterate iterate;
    if (Kokkos::is_layouttiled<typename DstType::array_layout>::value) {
      iterate = Kokkos::layout_iterate_type_selector<
          typename DstType::array_layout>::outer_iteration_pattern;
    } else if (std::is_same<typename DstType::array_layout,
                            Kokkos::LayoutRight>::value) {
      iterate = Kokkos::Iterate::Right;
    } else if (std::is_same<typename DstType::array_layout,
                            Kokkos::LayoutLeft>::value) {
      iterate = Kokkos::Iterate::Left;
    } else if (std::is_same<typename DstType::array_layout,
                            Kokkos::LayoutStride>::value) {
      if (strides[0] > strides[DstType::Rank - 1])
        iterate = Kokkos::Iterate::Right;
      else
        iterate = Kokkos::Iterate::Left;
    } else {
      if (std::is_same<typename DstType::execution_space::array_layout,
                       Kokkos::LayoutRight>::value)
        iterate = Kokkos::Iterate::Right;
      else
        iterate = Kokkos::Iterate::Left;
    }

    if ((dst.span() >= size_t(std::numeric_limits<int>::max())) ||
        (src.span() >= size_t(std::numeric_limits<int>::max()))) {
      if (iterate == Kokkos::Iterate::Right)
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutRight, ExecutionSpace, DstType::Rank, int64_t>(
            dst, src, space);
      else
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutLeft, ExecutionSpace, DstType::Rank, int64_t>(
            dst, src, space);
    } else {
      if (iterate == Kokkos::Iterate::Right)
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutRight, ExecutionSpace, DstType::Rank, int>(dst, src,
                                                                     space);
      else
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutLeft, ExecutionSpace, DstType::Rank, int>(dst, src,
                                                                    space);
    }
  }
}

template <class DstType, class SrcType>
void view_copy(const DstType& dst, const SrcType& src) {
  using dst_execution_space = typename DstType::execution_space;
  using src_execution_space = typename SrcType::execution_space;
  using dst_memory_space    = typename DstType::memory_space;
  using src_memory_space    = typename SrcType::memory_space;

  enum {
    DstExecCanAccessSrc =
        Kokkos::Impl::SpaceAccessibility<dst_execution_space,
                                         src_memory_space>::accessible
  };

  enum {
    SrcExecCanAccessDst =
        Kokkos::Impl::SpaceAccessibility<src_execution_space,
                                         dst_memory_space>::accessible
  };

  if (!DstExecCanAccessSrc && !SrcExecCanAccessDst) {
    std::string message(
        "Error: Kokkos::deep_copy with no available copy mechanism: ");
    message += src.label();
    message += " to ";
    message += dst.label();
    Kokkos::Impl::throw_runtime_exception(message);
  }

  // Figure out iteration order in case we need it
  int64_t strides[DstType::Rank + 1];
  dst.stride(strides);
  Kokkos::Iterate iterate;
  if (Kokkos::is_layouttiled<typename DstType::array_layout>::value) {
    iterate = Kokkos::layout_iterate_type_selector<
        typename DstType::array_layout>::outer_iteration_pattern;
  } else if (std::is_same<typename DstType::array_layout,
                          Kokkos::LayoutRight>::value) {
    iterate = Kokkos::Iterate::Right;
  } else if (std::is_same<typename DstType::array_layout,
                          Kokkos::LayoutLeft>::value) {
    iterate = Kokkos::Iterate::Left;
  } else if (std::is_same<typename DstType::array_layout,
                          Kokkos::LayoutStride>::value) {
    if (strides[0] > strides[DstType::Rank - 1])
      iterate = Kokkos::Iterate::Right;
    else
      iterate = Kokkos::Iterate::Left;
  } else {
    if (std::is_same<typename DstType::execution_space::array_layout,
                     Kokkos::LayoutRight>::value)
      iterate = Kokkos::Iterate::Right;
    else
      iterate = Kokkos::Iterate::Left;
  }

  if ((dst.span() >= size_t(std::numeric_limits<int>::max())) ||
      (src.span() >= size_t(std::numeric_limits<int>::max()))) {
    if (DstExecCanAccessSrc) {
      if (iterate == Kokkos::Iterate::Right)
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutRight, dst_execution_space, DstType::Rank, int64_t>(
            dst, src);
      else
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutLeft, dst_execution_space, DstType::Rank, int64_t>(
            dst, src);
    } else {
      if (iterate == Kokkos::Iterate::Right)
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutRight, src_execution_space, DstType::Rank, int64_t>(
            dst, src);
      else
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutLeft, src_execution_space, DstType::Rank, int64_t>(
            dst, src);
    }
  } else {
    if (DstExecCanAccessSrc) {
      if (iterate == Kokkos::Iterate::Right)
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutRight, dst_execution_space, DstType::Rank, int>(dst,
                                                                          src);
      else
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutLeft, dst_execution_space, DstType::Rank, int>(dst,
                                                                         src);
    } else {
      if (iterate == Kokkos::Iterate::Right)
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutRight, src_execution_space, DstType::Rank, int>(dst,
                                                                          src);
      else
        Kokkos::Impl::ViewCopy<
            typename DstType::uniform_runtime_nomemspace_type,
            typename SrcType::uniform_runtime_const_nomemspace_type,
            Kokkos::LayoutLeft, src_execution_space, DstType::Rank, int>(dst,
                                                                         src);
    }
  }
}

template <class DstType, class SrcType, int Rank, class... Args>
struct CommonSubview;

template <class DstType, class SrcType, class Arg0, class... Args>
struct CommonSubview<DstType, SrcType, 1, Arg0, Args...> {
  using dst_subview_type = typename Kokkos::Subview<DstType, Arg0>;
  using src_subview_type = typename Kokkos::Subview<SrcType, Arg0>;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& arg0,
                Args...)
      : dst_sub(dst, arg0), src_sub(src, arg0) {}
};

template <class DstType, class SrcType, class Arg0, class Arg1, class... Args>
struct CommonSubview<DstType, SrcType, 2, Arg0, Arg1, Args...> {
  using dst_subview_type = typename Kokkos::Subview<DstType, Arg0, Arg1>;
  using src_subview_type = typename Kokkos::Subview<SrcType, Arg0, Arg1>;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& arg0,
                const Arg1& arg1, Args...)
      : dst_sub(dst, arg0, arg1), src_sub(src, arg0, arg1) {}
};

template <class DstType, class SrcType, class Arg0, class Arg1, class Arg2,
          class... Args>
struct CommonSubview<DstType, SrcType, 3, Arg0, Arg1, Arg2, Args...> {
  using dst_subview_type = typename Kokkos::Subview<DstType, Arg0, Arg1, Arg2>;
  using src_subview_type = typename Kokkos::Subview<SrcType, Arg0, Arg1, Arg2>;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& arg0,
                const Arg1& arg1, const Arg2& arg2, Args...)
      : dst_sub(dst, arg0, arg1, arg2), src_sub(src, arg0, arg1, arg2) {}
};

template <class DstType, class SrcType, class Arg0, class Arg1, class Arg2,
          class Arg3, class... Args>
struct CommonSubview<DstType, SrcType, 4, Arg0, Arg1, Arg2, Arg3, Args...> {
  using dst_subview_type =
      typename Kokkos::Subview<DstType, Arg0, Arg1, Arg2, Arg3>;
  using src_subview_type =
      typename Kokkos::Subview<SrcType, Arg0, Arg1, Arg2, Arg3>;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& arg0,
                const Arg1& arg1, const Arg2& arg2, const Arg3& arg3,
                const Args...)
      : dst_sub(dst, arg0, arg1, arg2, arg3),
        src_sub(src, arg0, arg1, arg2, arg3) {}
};

template <class DstType, class SrcType, class Arg0, class Arg1, class Arg2,
          class Arg3, class Arg4, class... Args>
struct CommonSubview<DstType, SrcType, 5, Arg0, Arg1, Arg2, Arg3, Arg4,
                     Args...> {
  using dst_subview_type =
      typename Kokkos::Subview<DstType, Arg0, Arg1, Arg2, Arg3, Arg4>;
  using src_subview_type =
      typename Kokkos::Subview<SrcType, Arg0, Arg1, Arg2, Arg3, Arg4>;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& arg0,
                const Arg1& arg1, const Arg2& arg2, const Arg3& arg3,
                const Arg4& arg4, const Args...)
      : dst_sub(dst, arg0, arg1, arg2, arg3, arg4),
        src_sub(src, arg0, arg1, arg2, arg3, arg4) {}
};

template <class DstType, class SrcType, class Arg0, class Arg1, class Arg2,
          class Arg3, class Arg4, class Arg5, class... Args>
struct CommonSubview<DstType, SrcType, 6, Arg0, Arg1, Arg2, Arg3, Arg4, Arg5,
                     Args...> {
  using dst_subview_type =
      typename Kokkos::Subview<DstType, Arg0, Arg1, Arg2, Arg3, Arg4, Arg5>;
  using src_subview_type =
      typename Kokkos::Subview<SrcType, Arg0, Arg1, Arg2, Arg3, Arg4, Arg5>;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& arg0,
                const Arg1& arg1, const Arg2& arg2, const Arg3& arg3,
                const Arg4& arg4, const Arg5& arg5, const Args...)
      : dst_sub(dst, arg0, arg1, arg2, arg3, arg4, arg5),
        src_sub(src, arg0, arg1, arg2, arg3, arg4, arg5) {}
};

template <class DstType, class SrcType, class Arg0, class Arg1, class Arg2,
          class Arg3, class Arg4, class Arg5, class Arg6, class... Args>
struct CommonSubview<DstType, SrcType, 7, Arg0, Arg1, Arg2, Arg3, Arg4, Arg5,
                     Arg6, Args...> {
  using dst_subview_type = typename Kokkos::Subview<DstType, Arg0, Arg1, Arg2,
                                                    Arg3, Arg4, Arg5, Arg6>;
  using src_subview_type = typename Kokkos::Subview<SrcType, Arg0, Arg1, Arg2,
                                                    Arg3, Arg4, Arg5, Arg6>;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& arg0,
                const Arg1& arg1, const Arg2& arg2, const Arg3& arg3,
                const Arg4& arg4, const Arg5& arg5, const Arg6& arg6, Args...)
      : dst_sub(dst, arg0, arg1, arg2, arg3, arg4, arg5, arg6),
        src_sub(src, arg0, arg1, arg2, arg3, arg4, arg5, arg6) {}
};

template <class DstType, class SrcType, class Arg0, class Arg1, class Arg2,
          class Arg3, class Arg4, class Arg5, class Arg6, class Arg7>
struct CommonSubview<DstType, SrcType, 8, Arg0, Arg1, Arg2, Arg3, Arg4, Arg5,
                     Arg6, Arg7> {
  using dst_subview_type =
      typename Kokkos::Subview<DstType, Arg0, Arg1, Arg2, Arg3, Arg4, Arg5,
                               Arg6, Arg7>;
  using src_subview_type =
      typename Kokkos::Subview<SrcType, Arg0, Arg1, Arg2, Arg3, Arg4, Arg5,
                               Arg6, Arg7>;
  dst_subview_type dst_sub;
  src_subview_type src_sub;
  CommonSubview(const DstType& dst, const SrcType& src, const Arg0& arg0,
                const Arg1& arg1, const Arg2& arg2, const Arg3& arg3,
                const Arg4& arg4, const Arg5& arg5, const Arg6& arg6,
                const Arg7& arg7)
      : dst_sub(dst, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7),
        src_sub(src, arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7) {}
};

template <class DstType, class SrcType,
          class ExecSpace = typename DstType::execution_space,
          int Rank        = DstType::Rank>
struct ViewRemap;

template <class DstType, class SrcType, class ExecSpace>
struct ViewRemap<DstType, SrcType, ExecSpace, 1> {
  using p_type = Kokkos::pair<int64_t, int64_t>;

  ViewRemap(const DstType& dst, const SrcType& src) {
    if (dst.extent(0) == src.extent(0)) {
      view_copy(dst, src);
    } else {
      p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
      using sv_adapter_type = CommonSubview<DstType, SrcType, 1, p_type>;
      sv_adapter_type common_subview(dst, src, ext0);
      view_copy(common_subview.dst_sub, common_subview.src_sub);
    }
  }
};

template <class DstType, class SrcType, class ExecSpace>
struct ViewRemap<DstType, SrcType, ExecSpace, 2> {
  using p_type = Kokkos::pair<int64_t, int64_t>;

  ViewRemap(const DstType& dst, const SrcType& src) {
    if (dst.extent(0) == src.extent(0)) {
      if (dst.extent(1) == src.extent(1)) {
        view_copy(dst, src);
      } else {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 2, Kokkos::Impl::ALL_t, p_type>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    } else {
      if (dst.extent(1) == src.extent(1)) {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 2, p_type, Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, ext0, Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 2, p_type, p_type>;
        sv_adapter_type common_subview(dst, src, ext0, ext1);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    }
  }
};

template <class DstType, class SrcType, class ExecSpace>
struct ViewRemap<DstType, SrcType, ExecSpace, 3> {
  using p_type = Kokkos::pair<int64_t, int64_t>;

  ViewRemap(const DstType& dst, const SrcType& src) {
    if (dst.extent(0) == src.extent(0)) {
      if (dst.extent(2) == src.extent(2)) {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 3, Kokkos::Impl::ALL_t, p_type,
                          Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1,
                                       Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 3, Kokkos::Impl::ALL_t, p_type,
                          p_type>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    } else {
      if (dst.extent(2) == src.extent(2)) {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        using sv_adapter_type = CommonSubview<DstType, SrcType, 3, p_type,
                                              p_type, Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 3, p_type, p_type, p_type>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    }
  }
};

template <class DstType, class SrcType, class ExecSpace>
struct ViewRemap<DstType, SrcType, ExecSpace, 4> {
  using p_type = Kokkos::pair<int64_t, int64_t>;

  ViewRemap(const DstType& dst, const SrcType& src) {
    if (dst.extent(0) == src.extent(0)) {
      if (dst.extent(3) == src.extent(3)) {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 4, Kokkos::Impl::ALL_t, p_type,
                          p_type, Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2,
                                       Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 4, Kokkos::Impl::ALL_t, p_type,
                          p_type, p_type>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2, ext3);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    } else {
      if (dst.extent(7) == src.extent(7)) {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 4, p_type, p_type, p_type,
                          Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2, Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 4, p_type, p_type, p_type, p_type>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2, ext3);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    }
  }
};

template <class DstType, class SrcType, class ExecSpace>
struct ViewRemap<DstType, SrcType, ExecSpace, 5> {
  using p_type = Kokkos::pair<int64_t, int64_t>;

  ViewRemap(const DstType& dst, const SrcType& src) {
    if (dst.extent(0) == src.extent(0)) {
      if (dst.extent(4) == src.extent(4)) {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 5, Kokkos::Impl::ALL_t, p_type,
                          p_type, p_type, Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2, ext3,
                                       Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 5, Kokkos::Impl::ALL_t, p_type,
                          p_type, p_type, p_type>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2, ext3,
                                       ext4);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    } else {
      if (dst.extent(4) == src.extent(4)) {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 5, p_type, p_type, p_type, p_type,
                          Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2, ext3,
                                       Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        using sv_adapter_type = CommonSubview<DstType, SrcType, 5, p_type,
                                              p_type, p_type, p_type, p_type>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2, ext3, ext4);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    }
  }
};
template <class DstType, class SrcType, class ExecSpace>
struct ViewRemap<DstType, SrcType, ExecSpace, 6> {
  using p_type = Kokkos::pair<int64_t, int64_t>;

  ViewRemap(const DstType& dst, const SrcType& src) {
    if (dst.extent(0) == src.extent(0)) {
      if (dst.extent(5) == src.extent(5)) {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 6, Kokkos::Impl::ALL_t, p_type,
                          p_type, p_type, p_type, Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2, ext3,
                                       ext4, Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        p_type ext5(0, std::min(dst.extent(5), src.extent(5)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 6, Kokkos::Impl::ALL_t, p_type,
                          p_type, p_type, p_type, p_type>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2, ext3,
                                       ext4, ext5);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    } else {
      if (dst.extent(5) == src.extent(5)) {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));

        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 6, p_type, p_type, p_type, p_type,
                          p_type, Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2, ext3, ext4,
                                       Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        p_type ext5(0, std::min(dst.extent(5), src.extent(5)));

        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 6, p_type, p_type, p_type, p_type,
                          p_type, p_type>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2, ext3, ext4,
                                       ext5);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    }
  }
};

template <class DstType, class SrcType, class ExecSpace>
struct ViewRemap<DstType, SrcType, ExecSpace, 7> {
  using p_type = Kokkos::pair<int64_t, int64_t>;

  ViewRemap(const DstType& dst, const SrcType& src) {
    if (dst.extent(0) == src.extent(0)) {
      if (dst.extent(6) == src.extent(6)) {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        p_type ext5(0, std::min(dst.extent(5), src.extent(5)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 7, Kokkos::Impl::ALL_t, p_type,
                          p_type, p_type, p_type, p_type, Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2, ext3,
                                       ext4, ext5, Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        p_type ext5(0, std::min(dst.extent(5), src.extent(5)));
        p_type ext6(0, std::min(dst.extent(6), src.extent(6)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 7, Kokkos::Impl::ALL_t, p_type,
                          p_type, p_type, p_type, p_type, p_type>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2, ext3,
                                       ext4, ext5, ext6);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    } else {
      if (dst.extent(6) == src.extent(6)) {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        p_type ext5(0, std::min(dst.extent(5), src.extent(5)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 7, p_type, p_type, p_type, p_type,
                          p_type, p_type, Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2, ext3, ext4,
                                       ext5, Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        p_type ext5(0, std::min(dst.extent(5), src.extent(5)));
        p_type ext6(0, std::min(dst.extent(6), src.extent(6)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 7, p_type, p_type, p_type, p_type,
                          p_type, p_type, p_type>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2, ext3, ext4,
                                       ext5, ext6);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    }
  }
};

template <class DstType, class SrcType, class ExecSpace>
struct ViewRemap<DstType, SrcType, ExecSpace, 8> {
  using p_type = Kokkos::pair<int64_t, int64_t>;

  ViewRemap(const DstType& dst, const SrcType& src) {
    if (dst.extent(0) == src.extent(0)) {
      if (dst.extent(7) == src.extent(7)) {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        p_type ext5(0, std::min(dst.extent(5), src.extent(5)));
        p_type ext6(0, std::min(dst.extent(6), src.extent(6)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 8, Kokkos::Impl::ALL_t, p_type,
                          p_type, p_type, p_type, p_type, p_type,
                          Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2, ext3,
                                       ext4, ext5, ext6, Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        p_type ext5(0, std::min(dst.extent(5), src.extent(5)));
        p_type ext6(0, std::min(dst.extent(6), src.extent(6)));
        p_type ext7(0, std::min(dst.extent(7), src.extent(7)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 8, Kokkos::Impl::ALL_t, p_type,
                          p_type, p_type, p_type, p_type, p_type, p_type>;
        sv_adapter_type common_subview(dst, src, Kokkos::ALL, ext1, ext2, ext3,
                                       ext4, ext5, ext6, ext7);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    } else {
      if (dst.extent(7) == src.extent(7)) {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        p_type ext5(0, std::min(dst.extent(5), src.extent(5)));
        p_type ext6(0, std::min(dst.extent(6), src.extent(6)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 8, p_type, p_type, p_type, p_type,
                          p_type, p_type, p_type, Kokkos::Impl::ALL_t>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2, ext3, ext4,
                                       ext5, ext6, Kokkos::ALL);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      } else {
        p_type ext0(0, std::min(dst.extent(0), src.extent(0)));
        p_type ext1(0, std::min(dst.extent(1), src.extent(1)));
        p_type ext2(0, std::min(dst.extent(2), src.extent(2)));
        p_type ext3(0, std::min(dst.extent(3), src.extent(3)));
        p_type ext4(0, std::min(dst.extent(4), src.extent(4)));
        p_type ext5(0, std::min(dst.extent(5), src.extent(5)));
        p_type ext6(0, std::min(dst.extent(6), src.extent(6)));
        p_type ext7(0, std::min(dst.extent(7), src.extent(7)));
        using sv_adapter_type =
            CommonSubview<DstType, SrcType, 8, p_type, p_type, p_type, p_type,
                          p_type, p_type, p_type, p_type>;
        sv_adapter_type common_subview(dst, src, ext0, ext1, ext2, ext3, ext4,
                                       ext5, ext6, ext7);
        view_copy(common_subview.dst_sub, common_subview.src_sub);
      }
    }
  }
};

}  // namespace Impl

/** \brief  Deep copy a value from Host memory into a view.  */
template <class DT, class... DP>
inline void deep_copy(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<std::is_same<
        typename ViewTraits<DT, DP...>::specialize, void>::value>::type* =
        nullptr) {
  using ViewType        = View<DT, DP...>;
  using exec_space_type = typename ViewType::execution_space;

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginDeepCopy(
        Kokkos::Profiling::make_space_handle(ViewType::memory_space::name()),
        dst.label(), dst.data(),
        Kokkos::Profiling::make_space_handle(Kokkos::HostSpace::name()),
        "Scalar", &value, dst.span() * sizeof(typename ViewType::value_type));
  }

  if (dst.data() == nullptr) {
    Kokkos::fence();
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endDeepCopy();
    }
    return;
  }

  Kokkos::fence();
  static_assert(std::is_same<typename ViewType::non_const_value_type,
                             typename ViewType::value_type>::value,
                "deep_copy requires non-const type");

  // If contiguous we can simply do a 1D flat loop
  if (dst.span_is_contiguous()) {
    using ViewTypeFlat = Kokkos::View<
        typename ViewType::value_type*, Kokkos::LayoutRight,
        Kokkos::Device<typename ViewType::execution_space,
                       typename std::conditional<
                           ViewType::Rank == 0, typename ViewType::memory_space,
                           Kokkos::AnonymousSpace>::type>,
        Kokkos::MemoryTraits<0>>;

    ViewTypeFlat dst_flat(dst.data(), dst.size());
    if (dst.span() < static_cast<size_t>(std::numeric_limits<int>::max())) {
      Kokkos::Impl::ViewFill<ViewTypeFlat, Kokkos::LayoutRight, exec_space_type,
                             ViewTypeFlat::Rank, int>(dst_flat, value,
                                                      exec_space_type());
    } else
      Kokkos::Impl::ViewFill<ViewTypeFlat, Kokkos::LayoutRight, exec_space_type,
                             ViewTypeFlat::Rank, int64_t>(dst_flat, value,
                                                          exec_space_type());
    Kokkos::fence();
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endDeepCopy();
    }
    return;
  }

  // Figure out iteration order to do the ViewFill
  int64_t strides[ViewType::Rank + 1];
  dst.stride(strides);
  Kokkos::Iterate iterate;
  if (std::is_same<typename ViewType::array_layout,
                   Kokkos::LayoutRight>::value) {
    iterate = Kokkos::Iterate::Right;
  } else if (std::is_same<typename ViewType::array_layout,
                          Kokkos::LayoutLeft>::value) {
    iterate = Kokkos::Iterate::Left;
  } else if (std::is_same<typename ViewType::array_layout,
                          Kokkos::LayoutStride>::value) {
    if (strides[0] > strides[ViewType::Rank > 0 ? ViewType::Rank - 1 : 0])
      iterate = Kokkos::Iterate::Right;
    else
      iterate = Kokkos::Iterate::Left;
  } else {
    if (std::is_same<typename ViewType::execution_space::array_layout,
                     Kokkos::LayoutRight>::value)
      iterate = Kokkos::Iterate::Right;
    else
      iterate = Kokkos::Iterate::Left;
  }

  // Lets call the right ViewFill functor based on integer space needed and
  // iteration type
  using ViewTypeUniform = typename std::conditional<
      ViewType::Rank == 0, typename ViewType::uniform_runtime_type,
      typename ViewType::uniform_runtime_nomemspace_type>::type;
  if (dst.span() > static_cast<size_t>(std::numeric_limits<int>::max())) {
    if (iterate == Kokkos::Iterate::Right)
      Kokkos::Impl::ViewFill<ViewTypeUniform, Kokkos::LayoutRight,
                             exec_space_type, ViewType::Rank, int64_t>(
          dst, value, exec_space_type());
    else
      Kokkos::Impl::ViewFill<ViewTypeUniform, Kokkos::LayoutLeft,
                             exec_space_type, ViewType::Rank, int64_t>(
          dst, value, exec_space_type());
  } else {
    if (iterate == Kokkos::Iterate::Right)
      Kokkos::Impl::ViewFill<ViewTypeUniform, Kokkos::LayoutRight,
                             exec_space_type, ViewType::Rank, int>(
          dst, value, exec_space_type());
    else
      Kokkos::Impl::ViewFill<ViewTypeUniform, Kokkos::LayoutLeft,
                             exec_space_type, ViewType::Rank, int>(
          dst, value, exec_space_type());
  }
  Kokkos::fence();

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endDeepCopy();
  }
}

/** \brief  Deep copy into a value in Host memory from a view.  */
template <class ST, class... SP>
inline void deep_copy(
    typename ViewTraits<ST, SP...>::non_const_value_type& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<std::is_same<
        typename ViewTraits<ST, SP...>::specialize, void>::value>::type* =
        nullptr) {
  using src_traits       = ViewTraits<ST, SP...>;
  using src_memory_space = typename src_traits::memory_space;

  static_assert(src_traits::rank == 0,
                "ERROR: Non-rank-zero view in deep_copy( value , View )");

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginDeepCopy(
        Kokkos::Profiling::make_space_handle(Kokkos::HostSpace::name()),
        "Scalar", &dst,
        Kokkos::Profiling::make_space_handle(src_memory_space::name()),
        src.label(), src.data(),
        src.span() * sizeof(typename src_traits::value_type));
  }

  if (src.data() == nullptr) {
    Kokkos::fence();
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endDeepCopy();
    }
    return;
  }

  Kokkos::Impl::DeepCopy<HostSpace, src_memory_space>(&dst, src.data(),
                                                      sizeof(ST));
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endDeepCopy();
  }
}

//----------------------------------------------------------------------------
/** \brief  A deep copy between views of compatible type, and rank zero.  */
template <class DT, class... DP, class ST, class... SP>
inline void deep_copy(
    const View<DT, DP...>& dst, const View<ST, SP...>& src,
    typename std::enable_if<(
        std::is_same<typename ViewTraits<DT, DP...>::specialize, void>::value &&
        std::is_same<typename ViewTraits<ST, SP...>::specialize, void>::value &&
        (unsigned(ViewTraits<DT, DP...>::rank) == unsigned(0) &&
         unsigned(ViewTraits<ST, SP...>::rank) == unsigned(0)))>::type* =
        nullptr) {
  using dst_type = View<DT, DP...>;
  using src_type = View<ST, SP...>;

  using value_type       = typename dst_type::value_type;
  using dst_memory_space = typename dst_type::memory_space;
  using src_memory_space = typename src_type::memory_space;

  static_assert(std::is_same<typename dst_type::value_type,
                             typename src_type::non_const_value_type>::value,
                "deep_copy requires matching non-const destination type");

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginDeepCopy(
        Kokkos::Profiling::make_space_handle(dst_memory_space::name()),
        dst.label(), dst.data(),
        Kokkos::Profiling::make_space_handle(src_memory_space::name()),
        src.label(), src.data(),
        src.span() * sizeof(typename dst_type::value_type));
  }

  if (dst.data() == nullptr && src.data() == nullptr) {
    Kokkos::fence();
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endDeepCopy();
    }
    return;
  }

  Kokkos::fence();
  if (dst.data() != src.data()) {
    Kokkos::Impl::DeepCopy<dst_memory_space, src_memory_space>(
        dst.data(), src.data(), sizeof(value_type));
    Kokkos::fence();
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endDeepCopy();
  }
}

//----------------------------------------------------------------------------
/** \brief  A deep copy between views of the default specialization, compatible
 * type, same non-zero rank, same contiguous layout.
 */
template <class DT, class... DP, class ST, class... SP>
inline void deep_copy(
    const View<DT, DP...>& dst, const View<ST, SP...>& src,
    typename std::enable_if<(
        std::is_same<typename ViewTraits<DT, DP...>::specialize, void>::value &&
        std::is_same<typename ViewTraits<ST, SP...>::specialize, void>::value &&
        (unsigned(ViewTraits<DT, DP...>::rank) != 0 ||
         unsigned(ViewTraits<ST, SP...>::rank) != 0))>::type* = nullptr) {
  using dst_type            = View<DT, DP...>;
  using src_type            = View<ST, SP...>;
  using dst_execution_space = typename dst_type::execution_space;
  using src_execution_space = typename src_type::execution_space;
  using dst_memory_space    = typename dst_type::memory_space;
  using src_memory_space    = typename src_type::memory_space;
  using dst_value_type      = typename dst_type::value_type;
  using src_value_type      = typename src_type::value_type;

  static_assert(std::is_same<typename dst_type::value_type,
                             typename dst_type::non_const_value_type>::value,
                "deep_copy requires non-const destination type");

  static_assert((unsigned(dst_type::rank) == unsigned(src_type::rank)),
                "deep_copy requires Views of equal rank");

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginDeepCopy(
        Kokkos::Profiling::make_space_handle(dst_memory_space::name()),
        dst.label(), dst.data(),
        Kokkos::Profiling::make_space_handle(src_memory_space::name()),
        src.label(), src.data(),
        src.span() * sizeof(typename dst_type::value_type));
  }

  if (dst.data() == nullptr || src.data() == nullptr) {
    // throw if dimension mismatch
    if ((src.extent(0) != dst.extent(0)) || (src.extent(1) != dst.extent(1)) ||
        (src.extent(2) != dst.extent(2)) || (src.extent(3) != dst.extent(3)) ||
        (src.extent(4) != dst.extent(4)) || (src.extent(5) != dst.extent(5)) ||
        (src.extent(6) != dst.extent(6)) || (src.extent(7) != dst.extent(7))) {
      std::string message(
          "Deprecation Error: Kokkos::deep_copy extents of views don't "
          "match: ");
      message += dst.label();
      message += "(";
      for (int r = 0; r < dst_type::Rank - 1; r++) {
        message += std::to_string(dst.extent(r));
        message += ",";
      }
      message += std::to_string(dst.extent(dst_type::Rank - 1));
      message += ") ";
      message += src.label();
      message += "(";
      for (int r = 0; r < src_type::Rank - 1; r++) {
        message += std::to_string(src.extent(r));
        message += ",";
      }
      message += std::to_string(src.extent(src_type::Rank - 1));
      message += ") ";

      Kokkos::Impl::throw_runtime_exception(message);
    }
    Kokkos::fence();
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endDeepCopy();
    }
    return;
  }

  enum {
    DstExecCanAccessSrc =
        Kokkos::Impl::SpaceAccessibility<dst_execution_space,
                                         src_memory_space>::accessible
  };

  enum {
    SrcExecCanAccessDst =
        Kokkos::Impl::SpaceAccessibility<src_execution_space,
                                         dst_memory_space>::accessible
  };

  // Checking for Overlapping Views.
  dst_value_type* dst_start = dst.data();
  dst_value_type* dst_end   = dst.data() + dst.span();
  src_value_type* src_start = src.data();
  src_value_type* src_end   = src.data() + src.span();
  if (((std::ptrdiff_t)dst_start == (std::ptrdiff_t)src_start) &&
      ((std::ptrdiff_t)dst_end == (std::ptrdiff_t)src_end) &&
      (dst.span_is_contiguous() && src.span_is_contiguous())) {
    Kokkos::fence();
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endDeepCopy();
    }
    return;
  }

  if ((((std::ptrdiff_t)dst_start < (std::ptrdiff_t)src_end) &&
       ((std::ptrdiff_t)dst_end > (std::ptrdiff_t)src_start)) &&
      ((dst.span_is_contiguous() && src.span_is_contiguous()))) {
    std::string message("Error: Kokkos::deep_copy of overlapping views: ");
    message += dst.label();
    message += "(";
    message += std::to_string((std::ptrdiff_t)dst_start);
    message += ",";
    message += std::to_string((std::ptrdiff_t)dst_end);
    message += ") ";
    message += src.label();
    message += "(";
    message += std::to_string((std::ptrdiff_t)src_start);
    message += ",";
    message += std::to_string((std::ptrdiff_t)src_end);
    message += ") ";
    Kokkos::Impl::throw_runtime_exception(message);
  }

  // Check for same extents
  if ((src.extent(0) != dst.extent(0)) || (src.extent(1) != dst.extent(1)) ||
      (src.extent(2) != dst.extent(2)) || (src.extent(3) != dst.extent(3)) ||
      (src.extent(4) != dst.extent(4)) || (src.extent(5) != dst.extent(5)) ||
      (src.extent(6) != dst.extent(6)) || (src.extent(7) != dst.extent(7))) {
    std::string message(
        "Deprecation Error: Kokkos::deep_copy extents of views don't match: ");
    message += dst.label();
    message += "(";
    for (int r = 0; r < dst_type::Rank - 1; r++) {
      message += std::to_string(dst.extent(r));
      message += ",";
    }
    message += std::to_string(dst.extent(dst_type::Rank - 1));
    message += ") ";
    message += src.label();
    message += "(";
    for (int r = 0; r < src_type::Rank - 1; r++) {
      message += std::to_string(src.extent(r));
      message += ",";
    }
    message += std::to_string(src.extent(src_type::Rank - 1));
    message += ") ";

    Kokkos::Impl::throw_runtime_exception(message);
  }

  // If same type, equal layout, equal dimensions, equal span, and contiguous
  // memory then can byte-wise copy

  if (std::is_same<typename dst_type::value_type,
                   typename src_type::non_const_value_type>::value &&
      (std::is_same<typename dst_type::array_layout,
                    typename src_type::array_layout>::value ||
       (dst_type::rank == 1 && src_type::rank == 1)) &&
      dst.span_is_contiguous() && src.span_is_contiguous() &&
      ((dst_type::rank < 1) || (dst.stride_0() == src.stride_0())) &&
      ((dst_type::rank < 2) || (dst.stride_1() == src.stride_1())) &&
      ((dst_type::rank < 3) || (dst.stride_2() == src.stride_2())) &&
      ((dst_type::rank < 4) || (dst.stride_3() == src.stride_3())) &&
      ((dst_type::rank < 5) || (dst.stride_4() == src.stride_4())) &&
      ((dst_type::rank < 6) || (dst.stride_5() == src.stride_5())) &&
      ((dst_type::rank < 7) || (dst.stride_6() == src.stride_6())) &&
      ((dst_type::rank < 8) || (dst.stride_7() == src.stride_7()))) {
    const size_t nbytes = sizeof(typename dst_type::value_type) * dst.span();
    Kokkos::fence();
    if ((void*)dst.data() != (void*)src.data()) {
      Kokkos::Impl::DeepCopy<dst_memory_space, src_memory_space>(
          dst.data(), src.data(), nbytes);
      Kokkos::fence();
    }
  } else {
    Kokkos::fence();
    Impl::view_copy(dst, src);
    Kokkos::fence();
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endDeepCopy();
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
namespace Experimental {
/** \brief  A local deep copy between views of the default specialization,
 * compatible type, same non-zero rank.
 */
template <class TeamType, class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION
local_deep_copy_contiguous(const TeamType& team, const View<DT, DP...>& dst,
                           const View<ST, SP...>& src) {
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, src.span()),
                       [&](const int& i) { dst.data()[i] = src.data()[i]; });
}
//----------------------------------------------------------------------------
template <class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const View<DT, DP...>& dst, const View<ST, SP...>& src) {
  for (size_t i = 0; i < src.span(); ++i) {
    dst.data()[i] = src.data()[i];
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 1 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 1)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0);

  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N),
                       [&](const int& i) { dst(i) = src(i); });
  team.team_barrier();
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 2 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 2)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0) * dst.extent(1);

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, src);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0      = i % dst.extent(0);
      int i1      = i / dst.extent(0);
      dst(i0, i1) = src(i0, i1);
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 3 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 3)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0) * dst.extent(1) * dst.extent(2);

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, src);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0          = i % dst.extent(0);
      int itmp        = i / dst.extent(0);
      int i1          = itmp % dst.extent(1);
      int i2          = itmp / dst.extent(1);
      dst(i0, i1, i2) = src(i0, i1, i2);
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 4 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 4)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N =
      dst.extent(0) * dst.extent(1) * dst.extent(2) * dst.extent(3);

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, src);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0              = i % dst.extent(0);
      int itmp            = i / dst.extent(0);
      int i1              = itmp % dst.extent(1);
      itmp                = itmp / dst.extent(1);
      int i2              = itmp % dst.extent(2);
      int i3              = itmp / dst.extent(2);
      dst(i0, i1, i2, i3) = src(i0, i1, i2, i3);
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 5 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 5)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0) * dst.extent(1) * dst.extent(2) *
                   dst.extent(3) * dst.extent(4);

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, src);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0                  = i % dst.extent(0);
      int itmp                = i / dst.extent(0);
      int i1                  = itmp % dst.extent(1);
      itmp                    = itmp / dst.extent(1);
      int i2                  = itmp % dst.extent(2);
      itmp                    = itmp / dst.extent(2);
      int i3                  = itmp % dst.extent(3);
      int i4                  = itmp / dst.extent(3);
      dst(i0, i1, i2, i3, i4) = src(i0, i1, i2, i3, i4);
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 6 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 6)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0) * dst.extent(1) * dst.extent(2) *
                   dst.extent(3) * dst.extent(4) * dst.extent(5);

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, src);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0                      = i % dst.extent(0);
      int itmp                    = i / dst.extent(0);
      int i1                      = itmp % dst.extent(1);
      itmp                        = itmp / dst.extent(1);
      int i2                      = itmp % dst.extent(2);
      itmp                        = itmp / dst.extent(2);
      int i3                      = itmp % dst.extent(3);
      itmp                        = itmp / dst.extent(3);
      int i4                      = itmp % dst.extent(4);
      int i5                      = itmp / dst.extent(4);
      dst(i0, i1, i2, i3, i4, i5) = src(i0, i1, i2, i3, i4, i5);
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 7 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 7)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0) * dst.extent(1) * dst.extent(2) *
                   dst.extent(3) * dst.extent(4) * dst.extent(5) *
                   dst.extent(6);

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, src);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0                          = i % dst.extent(0);
      int itmp                        = i / dst.extent(0);
      int i1                          = itmp % dst.extent(1);
      itmp                            = itmp / dst.extent(1);
      int i2                          = itmp % dst.extent(2);
      itmp                            = itmp / dst.extent(2);
      int i3                          = itmp % dst.extent(3);
      itmp                            = itmp / dst.extent(3);
      int i4                          = itmp % dst.extent(4);
      itmp                            = itmp / dst.extent(4);
      int i5                          = itmp % dst.extent(5);
      int i6                          = itmp / dst.extent(5);
      dst(i0, i1, i2, i3, i4, i5, i6) = src(i0, i1, i2, i3, i4, i5, i6);
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst, const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 1 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 1)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0);

  for (size_t i = 0; i < N; ++i) {
    dst(i) = src(i);
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst, const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 2 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 2)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, src);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1) dst(i0, i1) = src(i0, i1);
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst, const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 3 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 3)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, src);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1)
        for (size_t i2 = 0; i2 < dst.extent(2); ++i2)
          dst(i0, i1, i2) = src(i0, i1, i2);
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst, const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 4 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 4)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, src);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1)
        for (size_t i2 = 0; i2 < dst.extent(2); ++i2)
          for (size_t i3 = 0; i3 < dst.extent(3); ++i3)
            dst(i0, i1, i2, i3) = src(i0, i1, i2, i3);
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst, const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 5 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 5)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, src);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1)
        for (size_t i2 = 0; i2 < dst.extent(2); ++i2)
          for (size_t i3 = 0; i3 < dst.extent(3); ++i3)
            for (size_t i4 = 0; i4 < dst.extent(4); ++i4)
              dst(i0, i1, i2, i3, i4) = src(i0, i1, i2, i3, i4);
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst, const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 6 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 6)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, src);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1)
        for (size_t i2 = 0; i2 < dst.extent(2); ++i2)
          for (size_t i3 = 0; i3 < dst.extent(3); ++i3)
            for (size_t i4 = 0; i4 < dst.extent(4); ++i4)
              for (size_t i5 = 0; i5 < dst.extent(5); ++i5)
                dst(i0, i1, i2, i3, i4, i5) = src(i0, i1, i2, i3, i4, i5);
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP, class ST, class... SP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst, const View<ST, SP...>& src,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) == 7 &&
                             unsigned(ViewTraits<ST, SP...>::rank) ==
                                 7)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous() && src.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, src);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1)
        for (size_t i2 = 0; i2 < dst.extent(2); ++i2)
          for (size_t i3 = 0; i3 < dst.extent(3); ++i3)
            for (size_t i4 = 0; i4 < dst.extent(4); ++i4)
              for (size_t i5 = 0; i5 < dst.extent(5); ++i5)
                for (size_t i6 = 0; i6 < dst.extent(6); ++i6)
                  dst(i0, i1, i2, i3, i4, i5, i6) =
                      src(i0, i1, i2, i3, i4, i5, i6);
  }
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Deep copy a value into a view.  */
template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value) {
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, dst.span()),
                       [&](const int& i) { dst.data()[i] = value; });
}
//----------------------------------------------------------------------------
template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy_contiguous(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value) {
  for (size_t i = 0; i < dst.span(); ++i) {
    dst.data()[i] = value;
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             1)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0);

  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N),
                       [&](const int& i) { dst(i) = value; });
  team.team_barrier();
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             2)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0) * dst.extent(1);

  if (dst.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, value);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0      = i % dst.extent(0);
      int i1      = i / dst.extent(0);
      dst(i0, i1) = value;
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             3)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0) * dst.extent(1) * dst.extent(2);

  if (dst.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, value);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0          = i % dst.extent(0);
      int itmp        = i / dst.extent(0);
      int i1          = itmp % dst.extent(1);
      int i2          = itmp / dst.extent(1);
      dst(i0, i1, i2) = value;
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             4)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N =
      dst.extent(0) * dst.extent(1) * dst.extent(2) * dst.extent(3);

  if (dst.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, value);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0              = i % dst.extent(0);
      int itmp            = i / dst.extent(0);
      int i1              = itmp % dst.extent(1);
      itmp                = itmp / dst.extent(1);
      int i2              = itmp % dst.extent(2);
      int i3              = itmp / dst.extent(2);
      dst(i0, i1, i2, i3) = value;
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             5)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0) * dst.extent(1) * dst.extent(2) *
                   dst.extent(3) * dst.extent(4);

  if (dst.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, value);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0                  = i % dst.extent(0);
      int itmp                = i / dst.extent(0);
      int i1                  = itmp % dst.extent(1);
      itmp                    = itmp / dst.extent(1);
      int i2                  = itmp % dst.extent(2);
      itmp                    = itmp / dst.extent(2);
      int i3                  = itmp % dst.extent(3);
      int i4                  = itmp / dst.extent(3);
      dst(i0, i1, i2, i3, i4) = value;
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             6)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0) * dst.extent(1) * dst.extent(2) *
                   dst.extent(3) * dst.extent(4) * dst.extent(5);

  if (dst.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, value);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0                      = i % dst.extent(0);
      int itmp                    = i / dst.extent(0);
      int i1                      = itmp % dst.extent(1);
      itmp                        = itmp / dst.extent(1);
      int i2                      = itmp % dst.extent(2);
      itmp                        = itmp / dst.extent(2);
      int i3                      = itmp % dst.extent(3);
      itmp                        = itmp / dst.extent(3);
      int i4                      = itmp % dst.extent(4);
      int i5                      = itmp / dst.extent(4);
      dst(i0, i1, i2, i3, i4, i5) = value;
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class TeamType, class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const TeamType& team, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             7)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0) * dst.extent(1) * dst.extent(2) *
                   dst.extent(3) * dst.extent(4) * dst.extent(5) *
                   dst.extent(6);

  if (dst.span_is_contiguous()) {
    team.team_barrier();
    local_deep_copy_contiguous(team, dst, value);
    team.team_barrier();
  } else {
    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, N), [&](const int& i) {
      int i0                          = i % dst.extent(0);
      int itmp                        = i / dst.extent(0);
      int i1                          = itmp % dst.extent(1);
      itmp                            = itmp / dst.extent(1);
      int i2                          = itmp % dst.extent(2);
      itmp                            = itmp / dst.extent(2);
      int i3                          = itmp % dst.extent(3);
      itmp                            = itmp / dst.extent(3);
      int i4                          = itmp % dst.extent(4);
      itmp                            = itmp / dst.extent(4);
      int i5                          = itmp % dst.extent(5);
      int i6                          = itmp / dst.extent(5);
      dst(i0, i1, i2, i3, i4, i5, i6) = value;
    });
    team.team_barrier();
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             1)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  const size_t N = dst.extent(0);

  for (size_t i = 0; i < N; ++i) {
    dst(i) = value;
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             2)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, value);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1) dst(i0, i1) = value;
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             3)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, value);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1)
        for (size_t i2 = 0; i2 < dst.extent(2); ++i2) dst(i0, i1, i2) = value;
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             4)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, value);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1)
        for (size_t i2 = 0; i2 < dst.extent(2); ++i2)
          for (size_t i3 = 0; i3 < dst.extent(3); ++i3)
            dst(i0, i1, i2, i3) = value;
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             5)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, value);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1)
        for (size_t i2 = 0; i2 < dst.extent(2); ++i2)
          for (size_t i3 = 0; i3 < dst.extent(3); ++i3)
            for (size_t i4 = 0; i4 < dst.extent(4); ++i4)
              dst(i0, i1, i2, i3, i4) = value;
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             6)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, value);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1)
        for (size_t i2 = 0; i2 < dst.extent(2); ++i2)
          for (size_t i3 = 0; i3 < dst.extent(3); ++i3)
            for (size_t i4 = 0; i4 < dst.extent(4); ++i4)
              for (size_t i5 = 0; i5 < dst.extent(5); ++i5)
                dst(i0, i1, i2, i3, i4, i5) = value;
  }
}
//----------------------------------------------------------------------------
template <class DT, class... DP>
void KOKKOS_INLINE_FUNCTION local_deep_copy(
    const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<(unsigned(ViewTraits<DT, DP...>::rank) ==
                             7)>::type* = nullptr) {
  if (dst.data() == nullptr) {
    return;
  }

  if (dst.span_is_contiguous()) {
    local_deep_copy_contiguous(dst, value);
  } else {
    for (size_t i0 = 0; i0 < dst.extent(0); ++i0)
      for (size_t i1 = 0; i1 < dst.extent(1); ++i1)
        for (size_t i2 = 0; i2 < dst.extent(2); ++i2)
          for (size_t i3 = 0; i3 < dst.extent(3); ++i3)
            for (size_t i4 = 0; i4 < dst.extent(4); ++i4)
              for (size_t i5 = 0; i5 < dst.extent(5); ++i5)
                for (size_t i6 = 0; i6 < dst.extent(6); ++i6)
                  dst(i0, i1, i2, i3, i4, i5, i6) = value;
  }
}
} /* namespace Experimental */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Deep copy a value from Host memory into a view. ExecSpace can access
 * dst */
template <class ExecSpace, class DT, class... DP>
inline void deep_copy(
    const ExecSpace& space, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<
        Kokkos::Impl::is_execution_space<ExecSpace>::value &&
        std::is_same<typename ViewTraits<DT, DP...>::specialize, void>::value &&
        Kokkos::Impl::SpaceAccessibility<
            ExecSpace,
            typename ViewTraits<DT, DP...>::memory_space>::accessible>::type* =
        nullptr) {
  using dst_traits = ViewTraits<DT, DP...>;
  static_assert(std::is_same<typename dst_traits::non_const_value_type,
                             typename dst_traits::value_type>::value,
                "deep_copy requires non-const type");
  using dst_memory_space = typename dst_traits::memory_space;
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginDeepCopy(
        Kokkos::Profiling::make_space_handle(dst_memory_space::name()),
        dst.label(), dst.data(),
        Kokkos::Profiling::make_space_handle(Kokkos::HostSpace::name()),
        "(none)", &value, dst.span() * sizeof(typename dst_traits::value_type));
  }
  if (dst.data() == nullptr) {
    space.fence();
  } else {
    using ViewTypeUniform = typename std::conditional<
        View<DT, DP...>::Rank == 0,
        typename View<DT, DP...>::uniform_runtime_type,
        typename View<DT, DP...>::uniform_runtime_nomemspace_type>::type;
    Kokkos::Impl::ViewFill<ViewTypeUniform, typename dst_traits::array_layout,
                           ExecSpace>(dst, value, space);
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endDeepCopy();
  }
}

/** \brief  Deep copy a value from Host memory into a view. ExecSpace can not
 * access dst */
template <class ExecSpace, class DT, class... DP>
inline void deep_copy(
    const ExecSpace& space, const View<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<
        Kokkos::Impl::is_execution_space<ExecSpace>::value &&
        std::is_same<typename ViewTraits<DT, DP...>::specialize, void>::value &&
        !Kokkos::Impl::SpaceAccessibility<
            ExecSpace,
            typename ViewTraits<DT, DP...>::memory_space>::accessible>::type* =
        nullptr) {
  using dst_traits = ViewTraits<DT, DP...>;
  static_assert(std::is_same<typename dst_traits::non_const_value_type,
                             typename dst_traits::value_type>::value,
                "deep_copy requires non-const type");
  using dst_memory_space = typename dst_traits::memory_space;
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginDeepCopy(
        Kokkos::Profiling::make_space_handle(dst_memory_space::name()),
        dst.label(), dst.data(),
        Kokkos::Profiling::make_space_handle(Kokkos::HostSpace::name()),
        "(none)", &value, dst.span() * sizeof(typename dst_traits::value_type));
  }
  if (dst.data() == nullptr) {
    space.fence();
  } else {
    space.fence();
    using ViewTypeUniform = typename std::conditional<
        View<DT, DP...>::Rank == 0,
        typename View<DT, DP...>::uniform_runtime_type,
        typename View<DT, DP...>::uniform_runtime_nomemspace_type>::type;
    using fill_exec_space = typename dst_traits::memory_space::execution_space;
    Kokkos::Impl::ViewFill<ViewTypeUniform, typename dst_traits::array_layout,
                           fill_exec_space>(dst, value, fill_exec_space());
    fill_exec_space().fence();
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endDeepCopy();
  }
}

/** \brief  Deep copy into a value in Host memory from a view.  */
template <class ExecSpace, class ST, class... SP>
inline void deep_copy(
    const ExecSpace& exec_space,
    typename ViewTraits<ST, SP...>::non_const_value_type& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<
        Kokkos::Impl::is_execution_space<ExecSpace>::value &&
        std::is_same<typename ViewTraits<ST, SP...>::specialize,
                     void>::value>::type* = nullptr) {
  using src_traits       = ViewTraits<ST, SP...>;
  using src_memory_space = typename src_traits::memory_space;
  static_assert(src_traits::rank == 0,
                "ERROR: Non-rank-zero view in deep_copy( value , View )");
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginDeepCopy(
        Kokkos::Profiling::make_space_handle(Kokkos::HostSpace::name()),
        "(none)", &dst,
        Kokkos::Profiling::make_space_handle(src_memory_space::name()),
        src.label(), src.data(), sizeof(ST));
  }

  if (src.data() == nullptr) {
    exec_space.fence();
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endDeepCopy();
    }
    return;
  }

  Kokkos::Impl::DeepCopy<HostSpace, src_memory_space, ExecSpace>(
      exec_space, &dst, src.data(), sizeof(ST));
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endDeepCopy();
  }
}

//----------------------------------------------------------------------------
/** \brief  A deep copy between views of compatible type, and rank zero.  */
template <class ExecSpace, class DT, class... DP, class ST, class... SP>
inline void deep_copy(
    const ExecSpace& exec_space, const View<DT, DP...>& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<(
        Kokkos::Impl::is_execution_space<ExecSpace>::value &&
        std::is_same<typename ViewTraits<DT, DP...>::specialize, void>::value &&
        std::is_same<typename ViewTraits<ST, SP...>::specialize, void>::value &&
        (unsigned(ViewTraits<DT, DP...>::rank) == unsigned(0) &&
         unsigned(ViewTraits<ST, SP...>::rank) == unsigned(0)))>::type* =
        nullptr) {
  using src_traits = ViewTraits<ST, SP...>;
  using dst_traits = ViewTraits<DT, DP...>;

  using src_memory_space = typename src_traits::memory_space;
  using dst_memory_space = typename dst_traits::memory_space;
  static_assert(std::is_same<typename dst_traits::value_type,
                             typename src_traits::non_const_value_type>::value,
                "deep_copy requires matching non-const destination type");

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginDeepCopy(
        Kokkos::Profiling::make_space_handle(dst_memory_space::name()),
        dst.label(), dst.data(),
        Kokkos::Profiling::make_space_handle(src_memory_space::name()),
        src.label(), src.data(), sizeof(DT));
  }

  if (dst.data() == nullptr && src.data() == nullptr) {
    exec_space.fence();
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endDeepCopy();
    }
    return;
  }

  if (dst.data() != src.data()) {
    Kokkos::Impl::DeepCopy<dst_memory_space, src_memory_space, ExecSpace>(
        exec_space, dst.data(), src.data(),
        sizeof(typename dst_traits::value_type));
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endDeepCopy();
  }
}

//----------------------------------------------------------------------------
/** \brief  A deep copy between views of the default specialization, compatible
 * type, same non-zero rank
 */
template <class ExecSpace, class DT, class... DP, class ST, class... SP>
inline void deep_copy(
    const ExecSpace& exec_space, const View<DT, DP...>& dst,
    const View<ST, SP...>& src,
    typename std::enable_if<(
        Kokkos::Impl::is_execution_space<ExecSpace>::value &&
        std::is_same<typename ViewTraits<DT, DP...>::specialize, void>::value &&
        std::is_same<typename ViewTraits<ST, SP...>::specialize, void>::value &&
        (unsigned(ViewTraits<DT, DP...>::rank) != 0 ||
         unsigned(ViewTraits<ST, SP...>::rank) != 0))>::type* = nullptr) {
  using dst_type = View<DT, DP...>;
  using src_type = View<ST, SP...>;

  static_assert(std::is_same<typename dst_type::value_type,
                             typename dst_type::non_const_value_type>::value,
                "deep_copy requires non-const destination type");

  static_assert((unsigned(dst_type::rank) == unsigned(src_type::rank)),
                "deep_copy requires Views of equal rank");

  using dst_execution_space = typename dst_type::execution_space;
  using src_execution_space = typename src_type::execution_space;
  using dst_memory_space    = typename dst_type::memory_space;
  using src_memory_space    = typename src_type::memory_space;
  using dst_value_type      = typename dst_type::value_type;
  using src_value_type      = typename src_type::value_type;

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginDeepCopy(
        Kokkos::Profiling::make_space_handle(dst_memory_space::name()),
        dst.label(), dst.data(),
        Kokkos::Profiling::make_space_handle(src_memory_space::name()),
        src.label(), src.data(), dst.span() * sizeof(dst_value_type));
  }

  dst_value_type* dst_start = dst.data();
  dst_value_type* dst_end   = dst.data() + dst.span();
  src_value_type* src_start = src.data();
  src_value_type* src_end   = src.data() + src.span();

  // Early dropout if identical range
  if ((dst_start == nullptr || src_start == nullptr) ||
      ((std::ptrdiff_t(dst_start) == std::ptrdiff_t(src_start)) &&
       (std::ptrdiff_t(dst_end) == std::ptrdiff_t(src_end)))) {
    // throw if dimension mismatch
    if ((src.extent(0) != dst.extent(0)) || (src.extent(1) != dst.extent(1)) ||
        (src.extent(2) != dst.extent(2)) || (src.extent(3) != dst.extent(3)) ||
        (src.extent(4) != dst.extent(4)) || (src.extent(5) != dst.extent(5)) ||
        (src.extent(6) != dst.extent(6)) || (src.extent(7) != dst.extent(7))) {
      std::string message(
          "Deprecation Error: Kokkos::deep_copy extents of views don't "
          "match: ");
      message += dst.label();
      message += "(";
      for (int r = 0; r < dst_type::Rank - 1; r++) {
        message += std::to_string(dst.extent(r));
        message += ",";
      }
      message += std::to_string(dst.extent(dst_type::Rank - 1));
      message += ") ";
      message += src.label();
      message += "(";
      for (int r = 0; r < src_type::Rank - 1; r++) {
        message += std::to_string(src.extent(r));
        message += ",";
      }
      message += std::to_string(src.extent(src_type::Rank - 1));
      message += ") ";

      Kokkos::Impl::throw_runtime_exception(message);
    }
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endDeepCopy();
    }
    return;
  }

  enum {
    ExecCanAccessSrcDst =
        Kokkos::Impl::SpaceAccessibility<ExecSpace,
                                         dst_memory_space>::accessible &&
        Kokkos::Impl::SpaceAccessibility<ExecSpace,
                                         src_memory_space>::accessible
  };
  enum {
    DstExecCanAccessSrc =
        Kokkos::Impl::SpaceAccessibility<dst_execution_space,
                                         src_memory_space>::accessible
  };

  enum {
    SrcExecCanAccessDst =
        Kokkos::Impl::SpaceAccessibility<src_execution_space,
                                         dst_memory_space>::accessible
  };

  // Error out for non-identical overlapping views.
  if ((((std::ptrdiff_t)dst_start < (std::ptrdiff_t)src_end) &&
       ((std::ptrdiff_t)dst_end > (std::ptrdiff_t)src_start)) &&
      ((dst.span_is_contiguous() && src.span_is_contiguous()))) {
    std::string message("Error: Kokkos::deep_copy of overlapping views: ");
    message += dst.label();
    message += "(";
    message += std::to_string((std::ptrdiff_t)dst_start);
    message += ",";
    message += std::to_string((std::ptrdiff_t)dst_end);
    message += ") ";
    message += src.label();
    message += "(";
    message += std::to_string((std::ptrdiff_t)src_start);
    message += ",";
    message += std::to_string((std::ptrdiff_t)src_end);
    message += ") ";
    Kokkos::Impl::throw_runtime_exception(message);
  }

  // Check for same extents
  if ((src.extent(0) != dst.extent(0)) || (src.extent(1) != dst.extent(1)) ||
      (src.extent(2) != dst.extent(2)) || (src.extent(3) != dst.extent(3)) ||
      (src.extent(4) != dst.extent(4)) || (src.extent(5) != dst.extent(5)) ||
      (src.extent(6) != dst.extent(6)) || (src.extent(7) != dst.extent(7))) {
    std::string message(
        "Deprecation Error: Kokkos::deep_copy extents of views don't match: ");
    message += dst.label();
    message += "(";
    for (int r = 0; r < dst_type::Rank - 1; r++) {
      message += std::to_string(dst.extent(r));
      message += ",";
    }
    message += std::to_string(dst.extent(dst_type::Rank - 1));
    message += ") ";
    message += src.label();
    message += "(";
    for (int r = 0; r < src_type::Rank - 1; r++) {
      message += std::to_string(src.extent(r));
      message += ",";
    }
    message += std::to_string(src.extent(src_type::Rank - 1));
    message += ") ";

    Kokkos::Impl::throw_runtime_exception(message);
  }

  // If same type, equal layout, equal dimensions, equal span, and contiguous
  // memory then can byte-wise copy

  if (std::is_same<typename dst_type::value_type,
                   typename src_type::non_const_value_type>::value &&
      (std::is_same<typename dst_type::array_layout,
                    typename src_type::array_layout>::value ||
       (dst_type::rank == 1 && src_type::rank == 1)) &&
      dst.span_is_contiguous() && src.span_is_contiguous() &&
      ((dst_type::rank < 1) || (dst.stride_0() == src.stride_0())) &&
      ((dst_type::rank < 2) || (dst.stride_1() == src.stride_1())) &&
      ((dst_type::rank < 3) || (dst.stride_2() == src.stride_2())) &&
      ((dst_type::rank < 4) || (dst.stride_3() == src.stride_3())) &&
      ((dst_type::rank < 5) || (dst.stride_4() == src.stride_4())) &&
      ((dst_type::rank < 6) || (dst.stride_5() == src.stride_5())) &&
      ((dst_type::rank < 7) || (dst.stride_6() == src.stride_6())) &&
      ((dst_type::rank < 8) || (dst.stride_7() == src.stride_7()))) {
    const size_t nbytes = sizeof(typename dst_type::value_type) * dst.span();
    if ((void*)dst.data() != (void*)src.data()) {
      Kokkos::Impl::DeepCopy<dst_memory_space, src_memory_space, ExecSpace>(
          exec_space, dst.data(), src.data(), nbytes);
    }
  } else {
    // Copying data between views in accessible memory spaces and either
    // non-contiguous or incompatible shape.
    if (ExecCanAccessSrcDst) {
      Impl::view_copy(exec_space, dst, src);
    } else if (DstExecCanAccessSrc || SrcExecCanAccessDst) {
      using cpy_exec_space =
          typename std::conditional<DstExecCanAccessSrc, dst_execution_space,
                                    src_execution_space>::type;
      exec_space.fence();
      Impl::view_copy(cpy_exec_space(), dst, src);
      cpy_exec_space().fence();
    } else {
      Kokkos::Impl::throw_runtime_exception(
          "deep_copy given views that would require a temporary allocation");
    }
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endDeepCopy();
  }
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Resize a view with copying old data to new data at the corresponding
 * indices. */
template <class T, class... P>
inline typename std::enable_if<
    std::is_same<typename Kokkos::View<T, P...>::array_layout,
                 Kokkos::LayoutLeft>::value ||
    std::is_same<typename Kokkos::View<T, P...>::array_layout,
                 Kokkos::LayoutRight>::value>::type
resize(Kokkos::View<T, P...>& v, const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
  using view_type = Kokkos::View<T, P...>;

  static_assert(Kokkos::ViewTraits<T, P...>::is_managed,
                "Can only resize managed views");

  // Fix #904 by checking dimensions before actually resizing.
  //
  // Rank is known at compile time, so hopefully the compiler will
  // remove branches that are compile-time false.  The upcoming "if
  // constexpr" language feature would make this certain.
  if (view_type::Rank == 1 && n0 == static_cast<size_t>(v.extent(0))) {
    return;
  }
  if (view_type::Rank == 2 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1))) {
    return;
  }
  if (view_type::Rank == 3 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2))) {
    return;
  }
  if (view_type::Rank == 4 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2)) &&
      n3 == static_cast<size_t>(v.extent(3))) {
    return;
  }
  if (view_type::Rank == 5 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2)) &&
      n3 == static_cast<size_t>(v.extent(3)) &&
      n4 == static_cast<size_t>(v.extent(4))) {
    return;
  }
  if (view_type::Rank == 6 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2)) &&
      n3 == static_cast<size_t>(v.extent(3)) &&
      n4 == static_cast<size_t>(v.extent(4)) &&
      n5 == static_cast<size_t>(v.extent(5))) {
    return;
  }
  if (view_type::Rank == 7 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2)) &&
      n3 == static_cast<size_t>(v.extent(3)) &&
      n4 == static_cast<size_t>(v.extent(4)) &&
      n5 == static_cast<size_t>(v.extent(5)) &&
      n6 == static_cast<size_t>(v.extent(6))) {
    return;
  }
  if (view_type::Rank == 8 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2)) &&
      n3 == static_cast<size_t>(v.extent(3)) &&
      n4 == static_cast<size_t>(v.extent(4)) &&
      n5 == static_cast<size_t>(v.extent(5)) &&
      n6 == static_cast<size_t>(v.extent(6)) &&
      n7 == static_cast<size_t>(v.extent(7))) {
    return;
  }
  // If Kokkos ever supports Views of rank > 8, the above code won't
  // be incorrect, because avoiding reallocation in resize() is just
  // an optimization.

  // TODO (mfh 27 Jun 2017) If the old View has enough space but just
  // different dimensions (e.g., if the product of the dimensions,
  // including extra space for alignment, will not change), then
  // consider just reusing storage.  For now, Kokkos always
  // reallocates if any of the dimensions change, even if the old View
  // has enough space.

  view_type v_resized(v.label(), n0, n1, n2, n3, n4, n5, n6, n7);

  Kokkos::Impl::ViewRemap<view_type, view_type>(v_resized, v);

  v = v_resized;
}

/** \brief  Resize a view with copying old data to new data at the corresponding
 * indices. */
template <class I, class T, class... P>
inline typename std::enable_if<
    std::is_same<typename Kokkos::View<T, P...>::array_layout,
                 Kokkos::LayoutLeft>::value ||
    std::is_same<typename Kokkos::View<T, P...>::array_layout,
                 Kokkos::LayoutRight>::value>::type
resize(const I& arg_prop, Kokkos::View<T, P...>& v,
       const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
       const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
  using view_type = Kokkos::View<T, P...>;

  static_assert(Kokkos::ViewTraits<T, P...>::is_managed,
                "Can only resize managed views");

  // Fix #904 by checking dimensions before actually resizing.
  //
  // Rank is known at compile time, so hopefully the compiler will
  // remove branches that are compile-time false.  The upcoming "if
  // constexpr" language feature would make this certain.
  if (view_type::Rank == 1 && n0 == static_cast<size_t>(v.extent(0))) {
    return;
  }
  if (view_type::Rank == 2 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1))) {
    return;
  }
  if (view_type::Rank == 3 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2))) {
    return;
  }
  if (view_type::Rank == 4 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2)) &&
      n3 == static_cast<size_t>(v.extent(3))) {
    return;
  }
  if (view_type::Rank == 5 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2)) &&
      n3 == static_cast<size_t>(v.extent(3)) &&
      n4 == static_cast<size_t>(v.extent(4))) {
    return;
  }
  if (view_type::Rank == 6 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2)) &&
      n3 == static_cast<size_t>(v.extent(3)) &&
      n4 == static_cast<size_t>(v.extent(4)) &&
      n5 == static_cast<size_t>(v.extent(5))) {
    return;
  }
  if (view_type::Rank == 7 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2)) &&
      n3 == static_cast<size_t>(v.extent(3)) &&
      n4 == static_cast<size_t>(v.extent(4)) &&
      n5 == static_cast<size_t>(v.extent(5)) &&
      n6 == static_cast<size_t>(v.extent(6))) {
    return;
  }
  if (view_type::Rank == 8 && n0 == static_cast<size_t>(v.extent(0)) &&
      n1 == static_cast<size_t>(v.extent(1)) &&
      n2 == static_cast<size_t>(v.extent(2)) &&
      n3 == static_cast<size_t>(v.extent(3)) &&
      n4 == static_cast<size_t>(v.extent(4)) &&
      n5 == static_cast<size_t>(v.extent(5)) &&
      n6 == static_cast<size_t>(v.extent(6)) &&
      n7 == static_cast<size_t>(v.extent(7))) {
    return;
  }
  // If Kokkos ever supports Views of rank > 8, the above code won't
  // be incorrect, because avoiding reallocation in resize() is just
  // an optimization.

  // TODO (mfh 27 Jun 2017) If the old View has enough space but just
  // different dimensions (e.g., if the product of the dimensions,
  // including extra space for alignment, will not change), then
  // consider just reusing storage.  For now, Kokkos always
  // reallocates if any of the dimensions change, even if the old View
  // has enough space.

  view_type v_resized(view_alloc(v.label(), std::forward<const I>(arg_prop)),
                      n0, n1, n2, n3, n4, n5, n6, n7);

  Kokkos::Impl::ViewRemap<view_type, view_type>(v_resized, v);

  v = v_resized;
}

/** \brief  Resize a view with copying old data to new data at the corresponding
 * indices. */
template <class T, class... P>
inline void resize(Kokkos::View<T, P...>& v,
                   const typename Kokkos::View<T, P...>::array_layout& layout) {
  using view_type = Kokkos::View<T, P...>;

  static_assert(Kokkos::ViewTraits<T, P...>::is_managed,
                "Can only resize managed views");

  view_type v_resized(v.label(), layout);

  Kokkos::Impl::ViewRemap<view_type, view_type>(v_resized, v);

  v = v_resized;
}

/** \brief  Resize a view with discarding old data. */
template <class T, class... P>
inline typename std::enable_if<
    std::is_same<typename Kokkos::View<T, P...>::array_layout,
                 Kokkos::LayoutLeft>::value ||
    std::is_same<typename Kokkos::View<T, P...>::array_layout,
                 Kokkos::LayoutRight>::value>::type
realloc(Kokkos::View<T, P...>& v,
        const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
  using view_type = Kokkos::View<T, P...>;

  static_assert(Kokkos::ViewTraits<T, P...>::is_managed,
                "Can only realloc managed views");

  const std::string label = v.label();

  v = view_type();  // Deallocate first, if the only view to allocation
  v = view_type(label, n0, n1, n2, n3, n4, n5, n6, n7);
}

/** \brief  Resize a view with discarding old data. */
template <class T, class... P>
inline void realloc(
    Kokkos::View<T, P...>& v,
    const typename Kokkos::View<T, P...>::array_layout& layout) {
  using view_type = Kokkos::View<T, P...>;

  static_assert(Kokkos::ViewTraits<T, P...>::is_managed,
                "Can only realloc managed views");

  const std::string label = v.label();

  v = view_type();  // Deallocate first, if the only view to allocation
  v = view_type(label, layout);
}
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Deduce Mirror Types
template <class Space, class T, class... P>
struct MirrorViewType {
  // The incoming view_type
  using src_view_type = typename Kokkos::View<T, P...>;
  // The memory space for the mirror view
  using memory_space = typename Space::memory_space;
  // Check whether it is the same memory space
  enum {
    is_same_memspace =
        std::is_same<memory_space, typename src_view_type::memory_space>::value
  };
  // The array_layout
  using array_layout = typename src_view_type::array_layout;
  // The data type (we probably want it non-const since otherwise we can't even
  // deep_copy to it.
  using data_type = typename src_view_type::non_const_data_type;
  // The destination view type if it is not the same memory space
  using dest_view_type = Kokkos::View<data_type, array_layout, Space>;
  // If it is the same memory_space return the existsing view_type
  // This will also keep the unmanaged trait if necessary
  using view_type = typename std::conditional<is_same_memspace, src_view_type,
                                              dest_view_type>::type;
};

template <class Space, class T, class... P>
struct MirrorType {
  // The incoming view_type
  using src_view_type = typename Kokkos::View<T, P...>;
  // The memory space for the mirror view
  using memory_space = typename Space::memory_space;
  // Check whether it is the same memory space
  enum {
    is_same_memspace =
        std::is_same<memory_space, typename src_view_type::memory_space>::value
  };
  // The array_layout
  using array_layout = typename src_view_type::array_layout;
  // The data type (we probably want it non-const since otherwise we can't even
  // deep_copy to it.
  using data_type = typename src_view_type::non_const_data_type;
  // The destination view type if it is not the same memory space
  using view_type = Kokkos::View<data_type, array_layout, Space>;
};

}  // namespace Impl

template <class T, class... P>
inline typename Kokkos::View<T, P...>::HostMirror create_mirror(
    const Kokkos::View<T, P...>& src,
    typename std::enable_if<
        std::is_same<typename ViewTraits<T, P...>::specialize, void>::value &&
        !std::is_same<typename Kokkos::ViewTraits<T, P...>::array_layout,
                      Kokkos::LayoutStride>::value>::type* = nullptr) {
  using src_type = View<T, P...>;
  using dst_type = typename src_type::HostMirror;

  return dst_type(
      std::string(src.label()).append("_mirror"),
      src.rank_dynamic > 0 ? src.extent(0) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      src.rank_dynamic > 1 ? src.extent(1) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      src.rank_dynamic > 2 ? src.extent(2) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      src.rank_dynamic > 3 ? src.extent(3) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      src.rank_dynamic > 4 ? src.extent(4) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      src.rank_dynamic > 5 ? src.extent(5) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      src.rank_dynamic > 6 ? src.extent(6) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      src.rank_dynamic > 7 ? src.extent(7) : KOKKOS_IMPL_CTOR_DEFAULT_ARG);
}

template <class T, class... P>
inline typename Kokkos::View<T, P...>::HostMirror create_mirror(
    const Kokkos::View<T, P...>& src,
    typename std::enable_if<
        std::is_same<typename ViewTraits<T, P...>::specialize, void>::value &&
        std::is_same<typename Kokkos::ViewTraits<T, P...>::array_layout,
                     Kokkos::LayoutStride>::value>::type* = nullptr) {
  using src_type = View<T, P...>;
  using dst_type = typename src_type::HostMirror;

  Kokkos::LayoutStride layout;

  layout.dimension[0] = src.extent(0);
  layout.dimension[1] = src.extent(1);
  layout.dimension[2] = src.extent(2);
  layout.dimension[3] = src.extent(3);
  layout.dimension[4] = src.extent(4);
  layout.dimension[5] = src.extent(5);
  layout.dimension[6] = src.extent(6);
  layout.dimension[7] = src.extent(7);

  layout.stride[0] = src.stride_0();
  layout.stride[1] = src.stride_1();
  layout.stride[2] = src.stride_2();
  layout.stride[3] = src.stride_3();
  layout.stride[4] = src.stride_4();
  layout.stride[5] = src.stride_5();
  layout.stride[6] = src.stride_6();
  layout.stride[7] = src.stride_7();

  return dst_type(std::string(src.label()).append("_mirror"), layout);
}

// Create a mirror in a new space (specialization for different space)
template <class Space, class T, class... P>
typename Impl::MirrorType<Space, T, P...>::view_type create_mirror(
    const Space&, const Kokkos::View<T, P...>& src,
    typename std::enable_if<std::is_same<
        typename ViewTraits<T, P...>::specialize, void>::value>::type* =
        nullptr) {
  return typename Impl::MirrorType<Space, T, P...>::view_type(src.label(),
                                                              src.layout());
}

template <class T, class... P>
inline typename Kokkos::View<T, P...>::HostMirror create_mirror_view(
    const Kokkos::View<T, P...>& src,
    typename std::enable_if<
        (std::is_same<
             typename Kokkos::View<T, P...>::memory_space,
             typename Kokkos::View<T, P...>::HostMirror::memory_space>::value &&
         std::is_same<typename Kokkos::View<T, P...>::data_type,
                      typename Kokkos::View<T, P...>::HostMirror::data_type>::
             value)>::type* = nullptr) {
  return src;
}

template <class T, class... P>
inline typename Kokkos::View<T, P...>::HostMirror create_mirror_view(
    const Kokkos::View<T, P...>& src,
    typename std::enable_if<!(
        std::is_same<
            typename Kokkos::View<T, P...>::memory_space,
            typename Kokkos::View<T, P...>::HostMirror::memory_space>::value &&
        std::is_same<typename Kokkos::View<T, P...>::data_type,
                     typename Kokkos::View<T, P...>::HostMirror::data_type>::
            value)>::type* = nullptr) {
  return Kokkos::create_mirror(src);
}

// Create a mirror view in a new space (specialization for same space)
template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type create_mirror_view(
    const Space&, const Kokkos::View<T, P...>& src,
    typename std::enable_if<
        Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type* =
        nullptr) {
  return src;
}

// Create a mirror view in a new space (specialization for different space)
template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type create_mirror_view(
    const Space&, const Kokkos::View<T, P...>& src,
    typename std::enable_if<
        !Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type* =
        nullptr) {
  return typename Impl::MirrorViewType<Space, T, P...>::view_type(src.label(),
                                                                  src.layout());
}

// Create a mirror view and deep_copy in a new space (specialization for same
// space)
template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type
create_mirror_view_and_copy(
    const Space&, const Kokkos::View<T, P...>& src,
    std::string const& name = "",
    typename std::enable_if<
        Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type* =
        nullptr) {
  (void)name;
  fence();  // same behavior as deep_copy(src, src)
  return src;
}

// Create a mirror view and deep_copy in a new space (specialization for
// different space)
template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type
create_mirror_view_and_copy(
    const Space&, const Kokkos::View<T, P...>& src,
    std::string const& name = "",
    typename std::enable_if<
        !Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type* =
        nullptr) {
  using Mirror      = typename Impl::MirrorViewType<Space, T, P...>::view_type;
  std::string label = name.empty() ? src.label() : name;
  auto mirror       = typename Mirror::non_const_type{
      view_alloc(WithoutInitializing, label), src.layout()};
  deep_copy(mirror, src);
  return mirror;
}

// Create a mirror view in a new space without initializing (specialization for
// same space)
template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type create_mirror_view(
    const Space&, const Kokkos::View<T, P...>& src,
    Kokkos::Impl::WithoutInitializing_t,
    typename std::enable_if<
        Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type* =
        nullptr) {
  return src;
}

// Create a mirror view in a new space without initializing (specialization for
// different space)
template <class Space, class T, class... P>
typename Impl::MirrorViewType<Space, T, P...>::view_type create_mirror_view(
    const Space&, const Kokkos::View<T, P...>& src,
    Kokkos::Impl::WithoutInitializing_t,
    typename std::enable_if<
        !Impl::MirrorViewType<Space, T, P...>::is_same_memspace>::type* =
        nullptr) {
  using Mirror = typename Impl::MirrorViewType<Space, T, P...>::view_type;
  return Mirror(view_alloc(WithoutInitializing, src.label()), src.layout());
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif
