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

#ifndef KOKKOS_EXPERIMENTAL_VIEWLAYOUTTILE_HPP
#define KOKKOS_EXPERIMENTAL_VIEWLAYOUTTILE_HPP

#ifndef KOKKOS_ENABLE_DEPRECATED_CODE

#include <Kokkos_Layout.hpp>
#include <Kokkos_View.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

// View offset and mapping for tiled view's

template <Kokkos::Iterate OuterP, Kokkos::Iterate InnerP, unsigned ArgN0,
          unsigned ArgN1>
struct is_array_layout<Kokkos::Experimental::LayoutTiled<
    OuterP, InnerP, ArgN0, ArgN1, 0, 0, 0, 0, 0, 0, true> >
    : public std::true_type {};

template <Kokkos::Iterate OuterP, Kokkos::Iterate InnerP, unsigned ArgN0,
          unsigned ArgN1, unsigned ArgN2>
struct is_array_layout<Kokkos::Experimental::LayoutTiled<
    OuterP, InnerP, ArgN0, ArgN1, ArgN2, 0, 0, 0, 0, 0, true> >
    : public std::true_type {};

template <Kokkos::Iterate OuterP, Kokkos::Iterate InnerP, unsigned ArgN0,
          unsigned ArgN1, unsigned ArgN2, unsigned ArgN3>
struct is_array_layout<Kokkos::Experimental::LayoutTiled<
    OuterP, InnerP, ArgN0, ArgN1, ArgN2, ArgN3, 0, 0, 0, 0, true> >
    : public std::true_type {};

template <Kokkos::Iterate OuterP, Kokkos::Iterate InnerP, unsigned ArgN0,
          unsigned ArgN1, unsigned ArgN2, unsigned ArgN3, unsigned ArgN4>
struct is_array_layout<Kokkos::Experimental::LayoutTiled<
    OuterP, InnerP, ArgN0, ArgN1, ArgN2, ArgN3, ArgN4, 0, 0, 0, true> >
    : public std::true_type {};

template <Kokkos::Iterate OuterP, Kokkos::Iterate InnerP, unsigned ArgN0,
          unsigned ArgN1, unsigned ArgN2, unsigned ArgN3, unsigned ArgN4,
          unsigned ArgN5>
struct is_array_layout<Kokkos::Experimental::LayoutTiled<
    OuterP, InnerP, ArgN0, ArgN1, ArgN2, ArgN3, ArgN4, ArgN5, 0, 0, true> >
    : public std::true_type {};

template <Kokkos::Iterate OuterP, Kokkos::Iterate InnerP, unsigned ArgN0,
          unsigned ArgN1, unsigned ArgN2, unsigned ArgN3, unsigned ArgN4,
          unsigned ArgN5, unsigned ArgN6>
struct is_array_layout<Kokkos::Experimental::LayoutTiled<
    OuterP, InnerP, ArgN0, ArgN1, ArgN2, ArgN3, ArgN4, ArgN5, ArgN6, 0, true> >
    : public std::true_type {};

template <Kokkos::Iterate OuterP, Kokkos::Iterate InnerP, unsigned ArgN0,
          unsigned ArgN1, unsigned ArgN2, unsigned ArgN3, unsigned ArgN4,
          unsigned ArgN5, unsigned ArgN6, unsigned ArgN7>
struct is_array_layout<
    Kokkos::Experimental::LayoutTiled<OuterP, InnerP, ArgN0, ArgN1, ArgN2,
                                      ArgN3, ArgN4, ArgN5, ArgN6, ArgN7, true> >
    : public std::true_type {};

template <class L>
struct is_array_layout_tiled : public std::false_type {};

template <Kokkos::Iterate OuterP, Kokkos::Iterate InnerP, unsigned ArgN0,
          unsigned ArgN1, unsigned ArgN2, unsigned ArgN3, unsigned ArgN4,
          unsigned ArgN5, unsigned ArgN6, unsigned ArgN7, bool IsPowerTwo>
struct is_array_layout_tiled<Kokkos::Experimental::LayoutTiled<
    OuterP, InnerP, ArgN0, ArgN1, ArgN2, ArgN3, ArgN4, ArgN5, ArgN6, ArgN7,
    IsPowerTwo> > : public std::true_type {
};  // Last template parameter "true" meaning this currently only supports
    // powers-of-two

namespace Impl {

template <class Dimension, class Layout>
struct ViewOffset<
    Dimension, Layout,
    typename std::enable_if<((Dimension::rank <= 8) && (Dimension::rank >= 2) &&
                             is_array_layout<Layout>::value &&
                             is_array_layout_tiled<Layout>::value)>::type> {
 public:
  //  enum { outer_pattern = Layout::outer_pattern };
  //  enum { inner_pattern = Layout::inner_pattern };
  static constexpr Kokkos::Iterate outer_pattern = Layout::outer_pattern;
  static constexpr Kokkos::Iterate inner_pattern = Layout::inner_pattern;

  enum { VORank = Dimension::rank };

  enum { SHIFT_0 = Kokkos::Impl::integral_power_of_two(Layout::N0) };
  enum { SHIFT_1 = Kokkos::Impl::integral_power_of_two(Layout::N1) };
  enum { SHIFT_2 = Kokkos::Impl::integral_power_of_two(Layout::N2) };
  enum { SHIFT_3 = Kokkos::Impl::integral_power_of_two(Layout::N3) };
  enum { SHIFT_4 = Kokkos::Impl::integral_power_of_two(Layout::N4) };
  enum { SHIFT_5 = Kokkos::Impl::integral_power_of_two(Layout::N5) };
  enum { SHIFT_6 = Kokkos::Impl::integral_power_of_two(Layout::N6) };
  enum { SHIFT_7 = Kokkos::Impl::integral_power_of_two(Layout::N7) };
  enum { MASK_0 = Layout::N0 - 1 };
  enum { MASK_1 = Layout::N1 - 1 };
  enum { MASK_2 = Layout::N2 - 1 };
  enum { MASK_3 = Layout::N3 - 1 };
  enum { MASK_4 = Layout::N4 - 1 };
  enum { MASK_5 = Layout::N5 - 1 };
  enum { MASK_6 = Layout::N6 - 1 };
  enum { MASK_7 = Layout::N7 - 1 };

  enum { SHIFT_2T = SHIFT_0 + SHIFT_1 };
  enum { SHIFT_3T = SHIFT_0 + SHIFT_1 + SHIFT_2 };
  enum { SHIFT_4T = SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 };
  enum { SHIFT_5T = SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 + SHIFT_4 };
  enum { SHIFT_6T = SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 + SHIFT_4 + SHIFT_5 };
  enum {
    SHIFT_7T =
        SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 + SHIFT_4 + SHIFT_5 + SHIFT_6
  };
  enum {
    SHIFT_8T = SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 + SHIFT_4 + SHIFT_5 +
               SHIFT_6 + SHIFT_7
  };

  // Is an irregular layout that does not have uniform striding for each index.
  using is_mapping_plugin = std::true_type;
  using is_regular        = std::false_type;

  typedef size_t size_type;
  typedef Dimension dimension_type;
  typedef Layout array_layout;

  dimension_type m_dim;
  size_type m_tile_N0;  // Num tiles dim 0
  size_type m_tile_N1;
  size_type m_tile_N2;
  size_type m_tile_N3;
  size_type m_tile_N4;
  size_type m_tile_N5;
  size_type m_tile_N6;
  size_type m_tile_N7;

  //----------------------------------------

#define DEBUG_OUTPUT_CHECK 0

  // Rank 2
  template <typename I0, typename I1>
  KOKKOS_INLINE_FUNCTION size_type operator()(I0 const& i0,
                                              I1 const& i1) const {
    auto tile_offset =
        (outer_pattern == (Kokkos::Iterate::Left))
            ? (((i0 >> SHIFT_0) + m_tile_N0 * ((i1 >> SHIFT_1))) << SHIFT_2T)
            : (((m_tile_N1 * (i0 >> SHIFT_0) + (i1 >> SHIFT_1))) << SHIFT_2T);
    //                     ( num_tiles[1] * ti0     +  ti1 ) * FTD

    auto local_offset = (inner_pattern == (Kokkos::Iterate::Left))
                            ? ((i0 & MASK_0) + ((i1 & MASK_1) << SHIFT_0))
                            : (((i0 & MASK_0) << SHIFT_1) + (i1 & MASK_1));
    //                     ( tile_dim[1] * li0         +  li1 )

#if DEBUG_OUTPUT_CHECK
    std::cout << "Am I Outer Left? "
              << (outer_pattern == (Kokkos::Iterate::Left)) << std::endl;
    std::cout << "Am I Inner Left? "
              << (inner_pattern == (Kokkos::Iterate::Left)) << std::endl;
    std::cout << "i0 = " << i0 << " i1 = " << i1
              << "\ntilei0 = " << (i0 >> SHIFT_0)
              << " tilei1 = " << (i1 >> SHIFT_1)
              << "locali0 = " << (i0 & MASK_0)
              << "\nlocali1 = " << (i1 & MASK_1) << std::endl;
#endif

    return tile_offset + local_offset;
  }

  // Rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_INLINE_FUNCTION size_type operator()(I0 const& i0, I1 const& i1,
                                              I2 const& i2) const {
    auto tile_offset =
        (outer_pattern == Kokkos::Iterate::Left)
            ? (((i0 >> SHIFT_0) +
                m_tile_N0 * ((i1 >> SHIFT_1) + m_tile_N1 * (i2 >> SHIFT_2)))
               << SHIFT_3T)
            : ((m_tile_N2 * (m_tile_N1 * (i0 >> SHIFT_0) + (i1 >> SHIFT_1)) +
                (i2 >> SHIFT_2))
               << SHIFT_3T);

    auto local_offset = (inner_pattern == Kokkos::Iterate::Left)
                            ? ((i0 & MASK_0) + ((i1 & MASK_1) << SHIFT_0) +
                               ((i2 & MASK_2) << (SHIFT_0 + SHIFT_1)))
                            : (((i0 & MASK_0) << (SHIFT_2 + SHIFT_1)) +
                               ((i1 & MASK_1) << (SHIFT_2)) + (i2 & MASK_2));

#if DEBUG_OUTPUT_CHECK
    std::cout << "Am I Outer Left? "
              << (outer_pattern == (Kokkos::Iterate::Left)) << std::endl;
    std::cout << "Am I Inner Left? "
              << (inner_pattern == (Kokkos::Iterate::Left)) << std::endl;
    std::cout << "i0 = " << i0 << " i1 = " << i1 << " i2 = " << i2
              << "\ntilei0 = " << (i0 >> SHIFT_0)
              << " tilei1 = " << (i1 >> SHIFT_1)
              << " tilei2 = " << (i2 >> SHIFT_2)
              << "\nlocali0 = " << (i0 & MASK_0)
              << "locali1 = " << (i1 & MASK_1) << "locali2 = " << (i2 & MASK_2)
              << std::endl;
#endif

    return tile_offset + local_offset;
  }

  // Rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_INLINE_FUNCTION size_type operator()(I0 const& i0, I1 const& i1,
                                              I2 const& i2,
                                              I3 const& i3) const {
    auto tile_offset =
        (outer_pattern == Kokkos::Iterate::Left)
            ? (((i0 >> SHIFT_0) +
                m_tile_N0 * ((i1 >> SHIFT_1) +
                             m_tile_N1 * ((i2 >> SHIFT_2) +
                                          m_tile_N2 * (i3 >> SHIFT_3))))
               << SHIFT_4T)
            : ((m_tile_N3 * (m_tile_N2 * (m_tile_N1 * (i0 >> SHIFT_0) +
                                          (i1 >> SHIFT_1)) +
                             (i2 >> SHIFT_2)) +
                (i3 >> SHIFT_3))
               << SHIFT_4T);

    auto local_offset =
        (inner_pattern == Kokkos::Iterate::Left)
            ? ((i0 & MASK_0) + ((i1 & MASK_1) << SHIFT_0) +
               ((i2 & MASK_2) << (SHIFT_0 + SHIFT_1)) +
               ((i3 & MASK_3) << (SHIFT_0 + SHIFT_1 + SHIFT_2)))
            : (((i0 & MASK_0) << (SHIFT_3 + SHIFT_2 + SHIFT_1)) +
               ((i1 & MASK_1) << (SHIFT_3 + SHIFT_2)) +
               ((i2 & MASK_2) << (SHIFT_3)) + (i3 & MASK_3));

    return tile_offset + local_offset;
  }

  // Rank 5
  template <typename I0, typename I1, typename I2, typename I3, typename I4>
  KOKKOS_INLINE_FUNCTION size_type operator()(I0 const& i0, I1 const& i1,
                                              I2 const& i2, I3 const& i3,
                                              I4 const& i4) const {
    auto tile_offset =
        (outer_pattern == Kokkos::Iterate::Left)
            ? (((i0 >> SHIFT_0) +
                m_tile_N0 *
                    ((i1 >> SHIFT_1) +
                     m_tile_N1 * ((i2 >> SHIFT_2) +
                                  m_tile_N2 * ((i3 >> SHIFT_3) +
                                               m_tile_N3 * (i4 >> SHIFT_4)))))
               << SHIFT_5T)
            : ((m_tile_N4 *
                    (m_tile_N3 * (m_tile_N2 * (m_tile_N1 * (i0 >> SHIFT_0) +
                                               (i1 >> SHIFT_1)) +
                                  (i2 >> SHIFT_2)) +
                     (i3 >> SHIFT_3)) +
                (i4 >> SHIFT_4))
               << SHIFT_5T);

    auto local_offset =
        (inner_pattern == Kokkos::Iterate::Left)
            ? ((i0 & MASK_0) + ((i1 & MASK_1) << SHIFT_0) +
               ((i2 & MASK_2) << (SHIFT_0 + SHIFT_1)) +
               ((i3 & MASK_3) << (SHIFT_0 + SHIFT_1 + SHIFT_2)) +
               ((i4 & MASK_4) << (SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3)))
            : (((i0 & MASK_0) << (SHIFT_4 + SHIFT_3 + SHIFT_2 + SHIFT_1)) +
               ((i1 & MASK_1) << (SHIFT_4 + SHIFT_3 + SHIFT_2)) +
               ((i2 & MASK_2) << (SHIFT_4 + SHIFT_3)) +
               ((i3 & MASK_3) << (SHIFT_4)) + (i4 & MASK_4));

    return tile_offset + local_offset;
  }

  // Rank 6
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5>
  KOKKOS_INLINE_FUNCTION size_type operator()(I0 const& i0, I1 const& i1,
                                              I2 const& i2, I3 const& i3,
                                              I4 const& i4,
                                              I5 const& i5) const {
    auto tile_offset =
        (outer_pattern == Kokkos::Iterate::Left)
            ? (((i0 >> SHIFT_0) +
                m_tile_N0 *
                    ((i1 >> SHIFT_1) +
                     m_tile_N1 *
                         ((i2 >> SHIFT_2) +
                          m_tile_N2 *
                              ((i3 >> SHIFT_3) +
                               m_tile_N3 * ((i4 >> SHIFT_4) +
                                            m_tile_N4 * (i5 >> SHIFT_5))))))
               << SHIFT_6T)
            : ((m_tile_N5 *
                    (m_tile_N4 *
                         (m_tile_N3 *
                              (m_tile_N2 * (m_tile_N1 * (i0 >> SHIFT_0) +
                                            (i1 >> SHIFT_1)) +
                               (i2 >> SHIFT_2)) +
                          (i3 >> SHIFT_3)) +
                     (i4 >> SHIFT_4)) +
                (i5 >> SHIFT_5))
               << SHIFT_6T);

    auto local_offset =
        (inner_pattern == Kokkos::Iterate::Left)
            ? ((i0 & MASK_0) + ((i1 & MASK_1) << SHIFT_0) +
               ((i2 & MASK_2) << (SHIFT_0 + SHIFT_1)) +
               ((i3 & MASK_3) << (SHIFT_0 + SHIFT_1 + SHIFT_2)) +
               ((i4 & MASK_4) << (SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3)) +
               ((i5 & MASK_5)
                << (SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 + SHIFT_4)))
            : (((i0 & MASK_0)
                << (SHIFT_5 + SHIFT_4 + SHIFT_3 + SHIFT_2 + SHIFT_1)) +
               ((i1 & MASK_1) << (SHIFT_5 + SHIFT_4 + SHIFT_3 + SHIFT_2)) +
               ((i2 & MASK_2) << (SHIFT_5 + SHIFT_4 + SHIFT_3)) +
               ((i3 & MASK_3) << (SHIFT_5 + SHIFT_4)) +
               ((i4 & MASK_4) << (SHIFT_5)) + (i5 & MASK_5));

    return tile_offset + local_offset;
  }

  // Rank 7
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6>
  KOKKOS_INLINE_FUNCTION size_type operator()(I0 const& i0, I1 const& i1,
                                              I2 const& i2, I3 const& i3,
                                              I4 const& i4, I5 const& i5,
                                              I6 const& i6) const {
    auto tile_offset =
        (outer_pattern == Kokkos::Iterate::Left)
            ? (((i0 >> SHIFT_0) +
                m_tile_N0 *
                    ((i1 >> SHIFT_1) +
                     m_tile_N1 *
                         ((i2 >> SHIFT_2) +
                          m_tile_N2 *
                              ((i3 >> SHIFT_3) +
                               m_tile_N3 *
                                   ((i4 >> SHIFT_4) +
                                    m_tile_N4 *
                                        ((i5 >> SHIFT_5) +
                                         m_tile_N5 * (i6 >> SHIFT_6)))))))
               << SHIFT_7T)
            : ((m_tile_N6 *
                    (m_tile_N5 *
                         (m_tile_N4 *
                              (m_tile_N3 *
                                   (m_tile_N2 * (m_tile_N1 * (i0 >> SHIFT_0) +
                                                 (i1 >> SHIFT_1)) +
                                    (i2 >> SHIFT_2)) +
                               (i3 >> SHIFT_3)) +
                          (i4 >> SHIFT_4)) +
                     (i5 >> SHIFT_5)) +
                (i6 >> SHIFT_6))
               << SHIFT_7T);

    auto local_offset =
        (inner_pattern == Kokkos::Iterate::Left)
            ? ((i0 & MASK_0) + ((i1 & MASK_1) << SHIFT_0) +
               ((i2 & MASK_2) << (SHIFT_0 + SHIFT_1)) +
               ((i3 & MASK_3) << (SHIFT_0 + SHIFT_1 + SHIFT_2)) +
               ((i4 & MASK_4) << (SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3)) +
               ((i5 & MASK_5)
                << (SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 + SHIFT_4)) +
               ((i6 & MASK_6)
                << (SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 + SHIFT_4 + SHIFT_5)))
            : (((i0 & MASK_0) << (SHIFT_6 + SHIFT_5 + SHIFT_4 + SHIFT_3 +
                                  SHIFT_2 + SHIFT_1)) +
               ((i1 & MASK_1)
                << (SHIFT_6 + SHIFT_5 + SHIFT_4 + SHIFT_3 + SHIFT_2)) +
               ((i2 & MASK_2) << (SHIFT_6 + SHIFT_5 + SHIFT_4 + SHIFT_3)) +
               ((i3 & MASK_3) << (SHIFT_6 + SHIFT_5 + SHIFT_4)) +
               ((i4 & MASK_4) << (SHIFT_6 + SHIFT_5)) +
               ((i5 & MASK_5) << (SHIFT_6)) + (i6 & MASK_6));

    return tile_offset + local_offset;
  }

  // Rank 8
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename I7>
  KOKKOS_INLINE_FUNCTION size_type operator()(I0 const& i0, I1 const& i1,
                                              I2 const& i2, I3 const& i3,
                                              I4 const& i4, I5 const& i5,
                                              I6 const& i6,
                                              I7 const& i7) const {
    auto tile_offset =
        (outer_pattern == Kokkos::Iterate::Left)
            ? (((i0 >> SHIFT_0) +
                m_tile_N0 *
                    ((i1 >> SHIFT_1) +
                     m_tile_N1 *
                         ((i2 >> SHIFT_2) +
                          m_tile_N2 *
                              ((i3 >> SHIFT_3) +
                               m_tile_N3 *
                                   ((i4 >> SHIFT_4) +
                                    m_tile_N4 *
                                        ((i5 >> SHIFT_5) +
                                         m_tile_N5 *
                                             ((i6 >> SHIFT_6) +
                                              m_tile_N6 * (i7 >> SHIFT_7))))))))
               << SHIFT_8T)
            : ((m_tile_N7 *
                    (m_tile_N6 *
                         (m_tile_N5 *
                              (m_tile_N4 *
                                   (m_tile_N3 *
                                        (m_tile_N2 *
                                             (m_tile_N1 * (i0 >> SHIFT_0) +
                                              (i1 >> SHIFT_1)) +
                                         (i2 >> SHIFT_2)) +
                                    (i3 >> SHIFT_3)) +
                               (i4 >> SHIFT_4)) +
                          (i5 >> SHIFT_5)) +
                     (i6 >> SHIFT_6)) +
                (i7 >> SHIFT_7))
               << SHIFT_8T);

    auto local_offset =
        (inner_pattern == Kokkos::Iterate::Left)
            ? ((i0 & MASK_0) + ((i1 & MASK_1) << SHIFT_0) +
               ((i2 & MASK_2) << (SHIFT_0 + SHIFT_1)) +
               ((i3 & MASK_3) << (SHIFT_0 + SHIFT_1 + SHIFT_2)) +
               ((i4 & MASK_4) << (SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3)) +
               ((i5 & MASK_5)
                << (SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 + SHIFT_4)) +
               ((i6 & MASK_6) << (SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 +
                                  SHIFT_4 + SHIFT_5)) +
               ((i7 & MASK_7) << (SHIFT_0 + SHIFT_1 + SHIFT_2 + SHIFT_3 +
                                  SHIFT_4 + SHIFT_5 + SHIFT_6)))
            : (((i0 & MASK_0) << (SHIFT_7 + SHIFT_6 + SHIFT_5 + SHIFT_4 +
                                  SHIFT_3 + SHIFT_2 + SHIFT_1)) +
               ((i1 & MASK_1) << (SHIFT_7 + SHIFT_6 + SHIFT_5 + SHIFT_4 +
                                  SHIFT_3 + SHIFT_2)) +
               ((i2 & MASK_2)
                << (SHIFT_7 + SHIFT_6 + SHIFT_5 + SHIFT_4 + SHIFT_3)) +
               ((i3 & MASK_3) << (SHIFT_7 + SHIFT_6 + SHIFT_5 + SHIFT_4)) +
               ((i4 & MASK_4) << (SHIFT_7 + SHIFT_6 + SHIFT_5)) +
               ((i5 & MASK_5) << (SHIFT_7 + SHIFT_6)) +
               ((i6 & MASK_6) << (SHIFT_7)) + (i7 & MASK_7));

    return tile_offset + local_offset;
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION constexpr array_layout layout() const {
    return array_layout(m_dim.N0, m_dim.N1, m_dim.N2, m_dim.N2, m_dim.N3,
                        m_dim.N4, m_dim.N5, m_dim.N6, m_dim.N7);
  }

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const {
    return m_dim.N0;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const {
    return m_dim.N1;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const {
    return m_dim.N2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const {
    return m_dim.N3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const {
    return m_dim.N4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const {
    return m_dim.N5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const {
    return m_dim.N6;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const {
    return m_dim.N7;
  }

  KOKKOS_INLINE_FUNCTION constexpr size_type size() const {
    return m_dim.N0 * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 *
           m_dim.N6 * m_dim.N7;
  }

  // Strides are meaningless due to irregularity
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const { return 0; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return 0; }

  // Stride with [ rank ] value is the total length
  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride(iType* const s) const {
    s[0] = 0;
    if (0 < dimension_type::rank) {
      s[1] = 0;
    }
    if (1 < dimension_type::rank) {
      s[2] = 0;
    }
    if (2 < dimension_type::rank) {
      s[3] = 0;
    }
    if (3 < dimension_type::rank) {
      s[4] = 0;
    }
    if (4 < dimension_type::rank) {
      s[5] = 0;
    }
    if (5 < dimension_type::rank) {
      s[6] = 0;
    }
    if (6 < dimension_type::rank) {
      s[7] = 0;
    }
    if (7 < dimension_type::rank) {
      s[8] = 0;
    }
  }

  KOKKOS_INLINE_FUNCTION constexpr size_type span() const {
    // Rank2: ( NumTile0 * ( NumTile1 ) ) * TileSize, etc
    return (VORank == 2)
               ? (m_tile_N0 * m_tile_N1) << SHIFT_2T
               : (VORank == 3)
                     ? (m_tile_N0 * m_tile_N1 * m_tile_N2) << SHIFT_3T
                     : (VORank == 4)
                           ? (m_tile_N0 * m_tile_N1 * m_tile_N2 * m_tile_N3)
                                 << SHIFT_4T
                           : (VORank == 5)
                                 ? (m_tile_N0 * m_tile_N1 * m_tile_N2 *
                                    m_tile_N3 * m_tile_N4)
                                       << SHIFT_5T
                                 : (VORank == 6)
                                       ? (m_tile_N0 * m_tile_N1 * m_tile_N2 *
                                          m_tile_N3 * m_tile_N4 * m_tile_N5)
                                             << SHIFT_6T
                                       : (VORank == 7)
                                             ? (m_tile_N0 * m_tile_N1 *
                                                m_tile_N2 * m_tile_N3 *
                                                m_tile_N4 * m_tile_N5 *
                                                m_tile_N6)
                                                   << SHIFT_7T
                                             : (m_tile_N0 * m_tile_N1 *
                                                m_tile_N2 * m_tile_N3 *
                                                m_tile_N4 * m_tile_N5 *
                                                m_tile_N6 * m_tile_N7)
                                                   << SHIFT_8T;
  }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const {
    return true;
  }

  //----------------------------------------

#ifdef KOKKOS_CUDA_9_DEFAULTED_BUG_WORKAROUND
  KOKKOS_INLINE_FUNCTION ~ViewOffset() {}
  KOKKOS_INLINE_FUNCTION ViewOffset() {}
  KOKKOS_INLINE_FUNCTION ViewOffset(const ViewOffset& rhs)
      : m_dim(rhs.m_dim),
        m_tile_N0(rhs.m_tile_N0),
        m_tile_N1(rhs.m_tile_N1),
        m_tile_N2(rhs.m_tile_N2),
        m_tile_N3(rhs.m_tile_N3),
        m_tile_N4(rhs.m_tile_N4),
        m_tile_N5(rhs.m_tile_N5),
        m_tile_N6(rhs.m_tile_N6),
        m_tile_N7(rhs.m_tile_N7) {}

  KOKKOS_INLINE_FUNCTION ViewOffset& operator=(const ViewOffset& rhs) {
    m_dim     = rhs.m_dim;
    m_tile_N0 = rhs.m_tile_N0;
    m_tile_N1 = rhs.m_tile_N1;
    m_tile_N2 = rhs.m_tile_N2;
    m_tile_N3 = rhs.m_tile_N3;
    m_tile_N4 = rhs.m_tile_N4;
    m_tile_N5 = rhs.m_tile_N5;
    m_tile_N6 = rhs.m_tile_N6;
    m_tile_N7 = rhs.m_tile_N7;
    return *this;
  }

#else
  KOKKOS_INLINE_FUNCTION ~ViewOffset()                 = default;
  KOKKOS_INLINE_FUNCTION ViewOffset()                  = default;
  KOKKOS_INLINE_FUNCTION ViewOffset(const ViewOffset&) = default;
  KOKKOS_INLINE_FUNCTION ViewOffset& operator=(const ViewOffset&) = default;
#endif

  template <unsigned TrivialScalarSize>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      std::integral_constant<unsigned, TrivialScalarSize> const&,
      array_layout const arg_layout)
      : m_dim(arg_layout.dimension[0], arg_layout.dimension[1],
              arg_layout.dimension[2], arg_layout.dimension[3],
              arg_layout.dimension[4], arg_layout.dimension[5],
              arg_layout.dimension[6], arg_layout.dimension[7]),
        m_tile_N0((arg_layout.dimension[0] + MASK_0) >>
                  SHIFT_0 /* number of tiles in first dimension */),
        m_tile_N1((arg_layout.dimension[1] + MASK_1) >> SHIFT_1),
        m_tile_N2((VORank > 2) ? (arg_layout.dimension[2] + MASK_2) >> SHIFT_2
                               : 0),
        m_tile_N3((VORank > 3) ? (arg_layout.dimension[3] + MASK_3) >> SHIFT_3
                               : 0),
        m_tile_N4((VORank > 4) ? (arg_layout.dimension[4] + MASK_4) >> SHIFT_4
                               : 0),
        m_tile_N5((VORank > 5) ? (arg_layout.dimension[5] + MASK_5) >> SHIFT_5
                               : 0),
        m_tile_N6((VORank > 6) ? (arg_layout.dimension[6] + MASK_6) >> SHIFT_6
                               : 0),
        m_tile_N7((VORank > 7) ? (arg_layout.dimension[7] + MASK_7) >> SHIFT_7
                               : 0) {}
};

//----------------------------------------

// ViewMapping assign method needed in order to return a 'subview' tile as a
// proper View The outer iteration pattern determines the mapping of the pointer
// offset to the beginning of requested tile The inner iteration pattern is
// needed for the layout of the tile's View to be returned Rank 2
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P, typename iType0,
          typename iType1>
struct ViewMapping<
    typename std::enable_if<(N2 == 0 && N3 == 0 && N4 == 0 && N5 == 0 &&
                             N6 == 0 && N7 == 0)>::type  // void
    ,
    Kokkos::ViewTraits<
        T**,
        Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                          N5, N6, N7, true>,
        P...>,
    Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4, N5,
                                      N6, N7, true>,
    iType0, iType1> {
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      src_layout;
  typedef Kokkos::ViewTraits<T**, src_layout, P...> src_traits;

  enum { is_outer_left = (OuterP == Kokkos::Iterate::Left) };
  enum { is_inner_left = (InnerP == Kokkos::Iterate::Left) };
  typedef typename std::conditional<is_inner_left, Kokkos::LayoutLeft,
                                    Kokkos::LayoutRight>::type array_layout;
  typedef Kokkos::ViewTraits<T[N0][N1], array_layout, P...> traits;
  typedef Kokkos::View<T[N0][N1], array_layout, P...> type;

  KOKKOS_INLINE_FUNCTION static void assign(
      ViewMapping<traits, void>& dst, const ViewMapping<src_traits, void>& src,
      const src_layout&, const iType0 i_tile0, const iType1 i_tile1) {
    typedef ViewMapping<traits, void> dst_map_type;
    typedef ViewMapping<src_traits, void> src_map_type;
    typedef typename dst_map_type::handle_type dst_handle_type;
    typedef typename dst_map_type::offset_type dst_offset_type;
    typedef typename src_map_type::offset_type src_offset_type;

    dst = dst_map_type(
        dst_handle_type(
            src.m_impl_handle +
            (is_outer_left ? ((i_tile0 + src.m_impl_offset.m_tile_N0 * i_tile1)
                              << src_offset_type::SHIFT_2T)
                           : ((src.m_impl_offset.m_tile_N1 * i_tile0 + i_tile1)
                              << src_offset_type::SHIFT_2T))  // offset to start
                                                              // of the tile
            ),
        dst_offset_type());
  }
};

// Rank 3
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P, typename iType0,
          typename iType1, typename iType2>
struct ViewMapping<typename std::enable_if<(N3 == 0 && N4 == 0 && N5 == 0 &&
                                            N6 == 0 && N7 == 0)>::type  // void
                   ,
                   Kokkos::ViewTraits<T***,
                                      Kokkos::Experimental::LayoutTiled<
                                          OuterP, InnerP, N0, N1, N2, N3, N4,
                                          N5, N6, N7, true>,
                                      P...>,
                   Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                     N3, N4, N5, N6, N7, true>,
                   iType0, iType1, iType2> {
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      src_layout;
  typedef Kokkos::ViewTraits<T***, src_layout, P...> src_traits;

  enum { is_outer_left = (OuterP == Kokkos::Iterate::Left) };
  enum { is_inner_left = (InnerP == Kokkos::Iterate::Left) };
  typedef typename std::conditional<is_inner_left, Kokkos::LayoutLeft,
                                    Kokkos::LayoutRight>::type array_layout;
  typedef Kokkos::ViewTraits<T[N0][N1][N2], array_layout, P...> traits;
  typedef Kokkos::View<T[N0][N1][N2], array_layout, P...> type;

  KOKKOS_INLINE_FUNCTION static void assign(
      ViewMapping<traits, void>& dst, const ViewMapping<src_traits, void>& src,
      const src_layout&, const iType0 i_tile0, const iType1 i_tile1,
      const iType2 i_tile2) {
    typedef ViewMapping<traits, void> dst_map_type;
    typedef ViewMapping<src_traits, void> src_map_type;
    typedef typename dst_map_type::handle_type dst_handle_type;
    typedef typename dst_map_type::offset_type dst_offset_type;
    typedef typename src_map_type::offset_type src_offset_type;

    dst = dst_map_type(
        dst_handle_type(
            src.m_impl_handle +
            (is_outer_left
                 ? ((i_tile0 +
                     src.m_impl_offset.m_tile_N0 *
                         (i_tile1 + src.m_impl_offset.m_tile_N1 * i_tile2))
                    << src_offset_type::SHIFT_3T)
                 : ((src.m_impl_offset.m_tile_N2 *
                         (src.m_impl_offset.m_tile_N1 * i_tile0 + i_tile1) +
                     i_tile2)
                    << src_offset_type::SHIFT_3T)))  // offset to start of the
                                                     // tile
        ,
        dst_offset_type());
  }
};

// Rank 4
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P, typename iType0,
          typename iType1, typename iType2, typename iType3>
struct ViewMapping<typename std::enable_if<(N4 == 0 && N5 == 0 && N6 == 0 &&
                                            N7 == 0)>::type  // void
                   ,
                   Kokkos::ViewTraits<T****,
                                      Kokkos::Experimental::LayoutTiled<
                                          OuterP, InnerP, N0, N1, N2, N3, N4,
                                          N5, N6, N7, true>,
                                      P...>,
                   Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                     N3, N4, N5, N6, N7, true>,
                   iType0, iType1, iType2, iType3> {
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      src_layout;
  typedef Kokkos::ViewTraits<T****, src_layout, P...> src_traits;

  enum { is_outer_left = (OuterP == Kokkos::Iterate::Left) };
  enum { is_inner_left = (InnerP == Kokkos::Iterate::Left) };
  typedef typename std::conditional<is_inner_left, Kokkos::LayoutLeft,
                                    Kokkos::LayoutRight>::type array_layout;
  typedef Kokkos::ViewTraits<T[N0][N1][N2][N3], array_layout, P...> traits;
  typedef Kokkos::View<T[N0][N1][N2][N3], array_layout, P...> type;

  KOKKOS_INLINE_FUNCTION static void assign(
      ViewMapping<traits, void>& dst, const ViewMapping<src_traits, void>& src,
      const src_layout&, const iType0 i_tile0, const iType1 i_tile1,
      const iType2 i_tile2, const iType3 i_tile3) {
    typedef ViewMapping<traits, void> dst_map_type;
    typedef ViewMapping<src_traits, void> src_map_type;
    typedef typename dst_map_type::handle_type dst_handle_type;
    typedef typename dst_map_type::offset_type dst_offset_type;
    typedef typename src_map_type::offset_type src_offset_type;

    dst = dst_map_type(
        dst_handle_type(
            src.m_impl_handle +
            (is_outer_left
                 ? ((i_tile0 +
                     src.m_impl_offset.m_tile_N0 *
                         (i_tile1 + src.m_impl_offset.m_tile_N1 *
                                        (i_tile2 + src.m_impl_offset.m_tile_N2 *
                                                       i_tile3)))
                    << src_offset_type::SHIFT_4T)
                 : ((src.m_impl_offset.m_tile_N3 *
                         (src.m_impl_offset.m_tile_N2 *
                              (src.m_impl_offset.m_tile_N1 * i_tile0 +
                               i_tile1) +
                          i_tile2) +
                     i_tile3)
                    << src_offset_type::SHIFT_4T)))  // offset to start of the
                                                     // tile
        ,
        dst_offset_type());
  }
};

// Rank 5
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P, typename iType0,
          typename iType1, typename iType2, typename iType3, typename iType4>
struct ViewMapping<
    typename std::enable_if<(N5 == 0 && N6 == 0 && N7 == 0)>::type  // void
    ,
    Kokkos::ViewTraits<
        T*****,
        Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                          N5, N6, N7, true>,
        P...>,
    Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4, N5,
                                      N6, N7, true>,
    iType0, iType1, iType2, iType3, iType4> {
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      src_layout;
  typedef Kokkos::ViewTraits<T*****, src_layout, P...> src_traits;

  enum { is_outer_left = (OuterP == Kokkos::Iterate::Left) };
  enum { is_inner_left = (InnerP == Kokkos::Iterate::Left) };
  typedef typename std::conditional<is_inner_left, Kokkos::LayoutLeft,
                                    Kokkos::LayoutRight>::type array_layout;
  typedef Kokkos::ViewTraits<T[N0][N1][N2][N3][N4], array_layout, P...> traits;
  typedef Kokkos::View<T[N0][N1][N2][N3][N4], array_layout, P...> type;

  KOKKOS_INLINE_FUNCTION static void assign(
      ViewMapping<traits, void>& dst, const ViewMapping<src_traits, void>& src,
      const src_layout&, const iType0 i_tile0, const iType1 i_tile1,
      const iType2 i_tile2, const iType3 i_tile3, const iType4 i_tile4) {
    typedef ViewMapping<traits, void> dst_map_type;
    typedef ViewMapping<src_traits, void> src_map_type;
    typedef typename dst_map_type::handle_type dst_handle_type;
    typedef typename dst_map_type::offset_type dst_offset_type;
    typedef typename src_map_type::offset_type src_offset_type;

    dst = dst_map_type(
        dst_handle_type(
            src.m_impl_handle +
            (is_outer_left
                 ? ((i_tile0 +
                     src.m_impl_offset.m_tile_N0 *
                         (i_tile1 +
                          src.m_impl_offset.m_tile_N1 *
                              (i_tile2 +
                               src.m_impl_offset.m_tile_N2 *
                                   (i_tile3 +
                                    src.m_impl_offset.m_tile_N3 * i_tile4))))
                    << src_offset_type::SHIFT_5T)
                 : ((src.m_impl_offset.m_tile_N4 *
                         (src.m_impl_offset.m_tile_N3 *
                              (src.m_impl_offset.m_tile_N2 *
                                   (src.m_impl_offset.m_tile_N1 * i_tile0 +
                                    i_tile1) +
                               i_tile2) +
                          i_tile3) +
                     i_tile4)
                    << src_offset_type::SHIFT_5T)))  // offset to start of the
                                                     // tile
        ,
        dst_offset_type());
  }
};

// Rank 6
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P, typename iType0,
          typename iType1, typename iType2, typename iType3, typename iType4,
          typename iType5>
struct ViewMapping<typename std::enable_if<(N6 == 0 && N7 == 0)>::type  // void
                   ,
                   Kokkos::ViewTraits<T******,
                                      Kokkos::Experimental::LayoutTiled<
                                          OuterP, InnerP, N0, N1, N2, N3, N4,
                                          N5, N6, N7, true>,
                                      P...>,
                   Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                     N3, N4, N5, N6, N7, true>,
                   iType0, iType1, iType2, iType3, iType4, iType5> {
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      src_layout;
  typedef Kokkos::ViewTraits<T******, src_layout, P...> src_traits;

  enum { is_outer_left = (OuterP == Kokkos::Iterate::Left) };
  enum { is_inner_left = (InnerP == Kokkos::Iterate::Left) };
  typedef typename std::conditional<is_inner_left, Kokkos::LayoutLeft,
                                    Kokkos::LayoutRight>::type array_layout;
  typedef Kokkos::ViewTraits<T[N0][N1][N2][N3][N4][N5], array_layout, P...>
      traits;
  typedef Kokkos::View<T[N0][N1][N2][N3][N4][N5], array_layout, P...> type;

  KOKKOS_INLINE_FUNCTION static void assign(
      ViewMapping<traits, void>& dst, const ViewMapping<src_traits, void>& src,
      const src_layout&, const iType0 i_tile0, const iType1 i_tile1,
      const iType2 i_tile2, const iType3 i_tile3, const iType4 i_tile4,
      const iType5 i_tile5) {
    typedef ViewMapping<traits, void> dst_map_type;
    typedef ViewMapping<src_traits, void> src_map_type;
    typedef typename dst_map_type::handle_type dst_handle_type;
    typedef typename dst_map_type::offset_type dst_offset_type;
    typedef typename src_map_type::offset_type src_offset_type;

    dst = dst_map_type(
        dst_handle_type(
            src.m_impl_handle +
            (is_outer_left
                 ? ((i_tile0 +
                     src.m_impl_offset.m_tile_N0 *
                         (i_tile1 +
                          src.m_impl_offset.m_tile_N1 *
                              (i_tile2 +
                               src.m_impl_offset.m_tile_N2 *
                                   (i_tile3 +
                                    src.m_impl_offset.m_tile_N3 *
                                        (i_tile4 + src.m_impl_offset.m_tile_N4 *
                                                       i_tile5)))))
                    << src_offset_type::SHIFT_6T)
                 : ((src.m_impl_offset.m_tile_N5 *
                         (src.m_impl_offset.m_tile_N4 *
                              (src.m_impl_offset.m_tile_N3 *
                                   (src.m_impl_offset.m_tile_N2 *
                                        (src.m_impl_offset.m_tile_N1 * i_tile0 +
                                         i_tile1) +
                                    i_tile2) +
                               i_tile3) +
                          i_tile4) +
                     i_tile5)
                    << src_offset_type::SHIFT_6T)))  // offset to start of the
                                                     // tile
        ,
        dst_offset_type());
  }
};

// Rank 7
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P, typename iType0,
          typename iType1, typename iType2, typename iType3, typename iType4,
          typename iType5, typename iType6>
struct ViewMapping<typename std::enable_if<(N7 == 0)>::type  // void
                   ,
                   Kokkos::ViewTraits<T*******,
                                      Kokkos::Experimental::LayoutTiled<
                                          OuterP, InnerP, N0, N1, N2, N3, N4,
                                          N5, N6, N7, true>,
                                      P...>,
                   Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                     N3, N4, N5, N6, N7, true>,
                   iType0, iType1, iType2, iType3, iType4, iType5, iType6> {
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      src_layout;
  typedef Kokkos::ViewTraits<T*******, src_layout, P...> src_traits;

  enum { is_outer_left = (OuterP == Kokkos::Iterate::Left) };
  enum { is_inner_left = (InnerP == Kokkos::Iterate::Left) };
  typedef typename std::conditional<is_inner_left, Kokkos::LayoutLeft,
                                    Kokkos::LayoutRight>::type array_layout;
  typedef Kokkos::ViewTraits<T[N0][N1][N2][N3][N4][N5][N6], array_layout, P...>
      traits;
  typedef Kokkos::View<T[N0][N1][N2][N3][N4][N5][N6], array_layout, P...> type;

  KOKKOS_INLINE_FUNCTION static void assign(
      ViewMapping<traits, void>& dst, const ViewMapping<src_traits, void>& src,
      const src_layout&, const iType0 i_tile0, const iType1 i_tile1,
      const iType2 i_tile2, const iType3 i_tile3, const iType4 i_tile4,
      const iType5 i_tile5, const iType6 i_tile6) {
    typedef ViewMapping<traits, void> dst_map_type;
    typedef ViewMapping<src_traits, void> src_map_type;
    typedef typename dst_map_type::handle_type dst_handle_type;
    typedef typename dst_map_type::offset_type dst_offset_type;
    typedef typename src_map_type::offset_type src_offset_type;

    dst = dst_map_type(
        dst_handle_type(
            src.m_impl_handle +
            (is_outer_left
                 ? ((i_tile0 +
                     src.m_impl_offset.m_tile_N0 *
                         (i_tile1 +
                          src.m_impl_offset.m_tile_N1 *
                              (i_tile2 +
                               src.m_impl_offset.m_tile_N2 *
                                   (i_tile3 +
                                    src.m_impl_offset.m_tile_N3 *
                                        (i_tile4 +
                                         src.m_impl_offset.m_tile_N4 *
                                             (i_tile5 +
                                              src.m_impl_offset.m_tile_N5 *
                                                  i_tile6))))))
                    << src_offset_type::SHIFT_7T)
                 : ((src.m_impl_offset.m_tile_N6 *
                         (src.m_impl_offset.m_tile_N5 *
                              (src.m_impl_offset.m_tile_N4 *
                                   (src.m_impl_offset.m_tile_N3 *
                                        (src.m_impl_offset.m_tile_N2 *
                                             (src.m_impl_offset.m_tile_N1 *
                                                  i_tile0 +
                                              i_tile1) +
                                         i_tile2) +
                                    i_tile3) +
                               i_tile4) +
                          i_tile5) +
                     i_tile6)
                    << src_offset_type::SHIFT_7T)))  // offset to start of the
                                                     // tile
        ,
        dst_offset_type());
  }
};

// Rank 8
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P, typename iType0,
          typename iType1, typename iType2, typename iType3, typename iType4,
          typename iType5, typename iType6, typename iType7>
struct ViewMapping<
    typename std::enable_if<(N0 != 0 && N1 != 0 && N2 != 0 && N3 != 0 &&
                             N4 != 0 && N5 != 0 && N6 != 0 &&
                             N7 != 0)>::type  // void
    ,
    Kokkos::ViewTraits<
        T********,
        Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                          N5, N6, N7, true>,
        P...>,
    Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4, N5,
                                      N6, N7, true>,
    iType0, iType1, iType2, iType3, iType4, iType5, iType6, iType7> {
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      src_layout;
  typedef Kokkos::ViewTraits<T********, src_layout, P...> src_traits;

  enum { is_outer_left = (OuterP == Kokkos::Iterate::Left) };
  enum { is_inner_left = (InnerP == Kokkos::Iterate::Left) };
  typedef typename std::conditional<is_inner_left, Kokkos::LayoutLeft,
                                    Kokkos::LayoutRight>::type array_layout;
  typedef Kokkos::ViewTraits<T[N0][N1][N2][N3][N4][N5][N6][N7], array_layout,
                             P...>
      traits;
  typedef Kokkos::View<T[N0][N1][N2][N3][N4][N5][N6][N7], array_layout, P...>
      type;

  KOKKOS_INLINE_FUNCTION static void assign(
      ViewMapping<traits, void>& dst, const ViewMapping<src_traits, void>& src,
      const src_layout&, const iType0 i_tile0, const iType1 i_tile1,
      const iType2 i_tile2, const iType3 i_tile3, const iType4 i_tile4,
      const iType5 i_tile5, const iType6 i_tile6, const iType7 i_tile7) {
    typedef ViewMapping<traits, void> dst_map_type;
    typedef ViewMapping<src_traits, void> src_map_type;
    typedef typename dst_map_type::handle_type dst_handle_type;
    typedef typename dst_map_type::offset_type dst_offset_type;
    typedef typename src_map_type::offset_type src_offset_type;

    dst = dst_map_type(
        dst_handle_type(
            src.m_impl_handle +
            (is_outer_left
                 ? ((i_tile0 +
                     src.m_impl_offset.m_tile_N0 *
                         (i_tile1 +
                          src.m_impl_offset.m_tile_N1 *
                              (i_tile2 +
                               src.m_impl_offset.m_tile_N2 *
                                   (i_tile3 +
                                    src.m_impl_offset.m_tile_N3 *
                                        (i_tile4 +
                                         src.m_impl_offset.m_tile_N4 *
                                             (i_tile5 +
                                              src.m_impl_offset.m_tile_N5 *
                                                  (i_tile6 +
                                                   src.m_impl_offset.m_tile_N6 *
                                                       i_tile7)))))))
                    << src_offset_type::SHIFT_8T)
                 : ((src.m_impl_offset.m_tile_N7 *
                         (src.m_impl_offset.m_tile_N6 *
                              (src.m_impl_offset.m_tile_N5 *
                                   (src.m_impl_offset.m_tile_N4 *
                                        (src.m_impl_offset.m_tile_N3 *
                                             (src.m_impl_offset.m_tile_N2 *
                                                  (src.m_impl_offset.m_tile_N1 *
                                                       i_tile0 +
                                                   i_tile1) +
                                              i_tile2) +
                                         i_tile3) +
                                    i_tile4) +
                               i_tile5) +
                          i_tile6) +
                     i_tile7)
                    << src_offset_type::SHIFT_8T)))  // offset to start of the
                                                     // tile
        ,
        dst_offset_type());
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------

namespace Kokkos {

// Rank 2
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P>
KOKKOS_INLINE_FUNCTION Kokkos::View<
    T[N0][N1],
    typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                              Kokkos::LayoutLeft, Kokkos::LayoutRight>::type,
    P...>
tile_subview(const Kokkos::View<
                 T**,
                 Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                   N3, N4, N5, N6, N7, true>,
                 P...>& src,
             const size_t i_tile0, const size_t i_tile1) {
  // Force the specialized ViewMapping for extracting a tile
  // by using the first subview argument as the layout.
  typedef
      typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                                Kokkos::LayoutLeft, Kokkos::LayoutRight>::type
          array_layout;
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      SrcLayout;

  return Kokkos::View<T[N0][N1], array_layout, P...>(src, SrcLayout(), i_tile0,
                                                     i_tile1);
}

// Rank 3
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P>
KOKKOS_INLINE_FUNCTION Kokkos::View<
    T[N0][N1][N2],
    typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                              Kokkos::LayoutLeft, Kokkos::LayoutRight>::type,
    P...>
tile_subview(const Kokkos::View<
                 T***,
                 Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                   N3, N4, N5, N6, N7, true>,
                 P...>& src,
             const size_t i_tile0, const size_t i_tile1, const size_t i_tile2) {
  // Force the specialized ViewMapping for extracting a tile
  // by using the first subview argument as the layout.
  typedef
      typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                                Kokkos::LayoutLeft, Kokkos::LayoutRight>::type
          array_layout;
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      SrcLayout;

  return Kokkos::View<T[N0][N1][N2], array_layout, P...>(
      src, SrcLayout(), i_tile0, i_tile1, i_tile2);
}

// Rank 4
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P>
KOKKOS_INLINE_FUNCTION Kokkos::View<
    T[N0][N1][N2][N3],
    typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                              Kokkos::LayoutLeft, Kokkos::LayoutRight>::type,
    P...>
tile_subview(const Kokkos::View<
                 T****,
                 Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                   N3, N4, N5, N6, N7, true>,
                 P...>& src,
             const size_t i_tile0, const size_t i_tile1, const size_t i_tile2,
             const size_t i_tile3) {
  // Force the specialized ViewMapping for extracting a tile
  // by using the first subview argument as the layout.
  typedef
      typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                                Kokkos::LayoutLeft, Kokkos::LayoutRight>::type
          array_layout;
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      SrcLayout;

  return Kokkos::View<T[N0][N1][N2][N3], array_layout, P...>(
      src, SrcLayout(), i_tile0, i_tile1, i_tile2, i_tile3);
}

// Rank 5
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P>
KOKKOS_INLINE_FUNCTION Kokkos::View<
    T[N0][N1][N2][N3][N4],
    typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                              Kokkos::LayoutLeft, Kokkos::LayoutRight>::type,
    P...>
tile_subview(const Kokkos::View<
                 T*****,
                 Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                   N3, N4, N5, N6, N7, true>,
                 P...>& src,
             const size_t i_tile0, const size_t i_tile1, const size_t i_tile2,
             const size_t i_tile3, const size_t i_tile4) {
  // Force the specialized ViewMapping for extracting a tile
  // by using the first subview argument as the layout.
  typedef
      typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                                Kokkos::LayoutLeft, Kokkos::LayoutRight>::type
          array_layout;
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      SrcLayout;

  return Kokkos::View<T[N0][N1][N2][N3][N4], array_layout, P...>(
      src, SrcLayout(), i_tile0, i_tile1, i_tile2, i_tile3, i_tile4);
}

// Rank 6
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P>
KOKKOS_INLINE_FUNCTION Kokkos::View<
    T[N0][N1][N2][N3][N4][N5],
    typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                              Kokkos::LayoutLeft, Kokkos::LayoutRight>::type,
    P...>
tile_subview(const Kokkos::View<
                 T******,
                 Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                   N3, N4, N5, N6, N7, true>,
                 P...>& src,
             const size_t i_tile0, const size_t i_tile1, const size_t i_tile2,
             const size_t i_tile3, const size_t i_tile4, const size_t i_tile5) {
  // Force the specialized ViewMapping for extracting a tile
  // by using the first subview argument as the layout.
  typedef
      typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                                Kokkos::LayoutLeft, Kokkos::LayoutRight>::type
          array_layout;
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      SrcLayout;

  return Kokkos::View<T[N0][N1][N2][N3][N4][N5], array_layout, P...>(
      src, SrcLayout(), i_tile0, i_tile1, i_tile2, i_tile3, i_tile4, i_tile5);
}

// Rank 7
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P>
KOKKOS_INLINE_FUNCTION Kokkos::View<
    T[N0][N1][N2][N3][N4][N5][N6],
    typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                              Kokkos::LayoutLeft, Kokkos::LayoutRight>::type,
    P...>
tile_subview(const Kokkos::View<
                 T*******,
                 Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                   N3, N4, N5, N6, N7, true>,
                 P...>& src,
             const size_t i_tile0, const size_t i_tile1, const size_t i_tile2,
             const size_t i_tile3, const size_t i_tile4, const size_t i_tile5,
             const size_t i_tile6) {
  // Force the specialized ViewMapping for extracting a tile
  // by using the first subview argument as the layout.
  typedef
      typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                                Kokkos::LayoutLeft, Kokkos::LayoutRight>::type
          array_layout;
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      SrcLayout;

  return Kokkos::View<T[N0][N1][N2][N3][N4][N5][N6], array_layout, P...>(
      src, SrcLayout(), i_tile0, i_tile1, i_tile2, i_tile3, i_tile4, i_tile5,
      i_tile6);
}

// Rank 8
template <typename T, Kokkos::Iterate OuterP, Kokkos::Iterate InnerP,
          unsigned N0, unsigned N1, unsigned N2, unsigned N3, unsigned N4,
          unsigned N5, unsigned N6, unsigned N7, class... P>
KOKKOS_INLINE_FUNCTION Kokkos::View<
    T[N0][N1][N2][N3][N4][N5][N6][N7],
    typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                              Kokkos::LayoutLeft, Kokkos::LayoutRight>::type,
    P...>
tile_subview(const Kokkos::View<
                 T********,
                 Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2,
                                                   N3, N4, N5, N6, N7, true>,
                 P...>& src,
             const size_t i_tile0, const size_t i_tile1, const size_t i_tile2,
             const size_t i_tile3, const size_t i_tile4, const size_t i_tile5,
             const size_t i_tile6, const size_t i_tile7) {
  // Force the specialized ViewMapping for extracting a tile
  // by using the first subview argument as the layout.
  typedef
      typename std::conditional<(InnerP == Kokkos::Iterate::Left),
                                Kokkos::LayoutLeft, Kokkos::LayoutRight>::type
          array_layout;
  typedef Kokkos::Experimental::LayoutTiled<OuterP, InnerP, N0, N1, N2, N3, N4,
                                            N5, N6, N7, true>
      SrcLayout;

  return Kokkos::View<T[N0][N1][N2][N3][N4][N5][N6][N7], array_layout, P...>(
      src, SrcLayout(), i_tile0, i_tile1, i_tile2, i_tile3, i_tile4, i_tile5,
      i_tile6, i_tile7);
}

} /* namespace Kokkos */
#endif  //! defined(KOKKOS_ENABLE_DEPRECATED_CODE
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXPERIENTAL_VIEWLAYOUTTILE_HPP */
