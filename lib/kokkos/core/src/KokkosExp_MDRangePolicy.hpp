/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CORE_EXP_MD_RANGE_POLICY_HPP
#define KOKKOS_CORE_EXP_MD_RANGE_POLICY_HPP

#include <Kokkos_ExecPolicy.hpp>
#include <Kokkos_Parallel.hpp>
#include <initializer_list>

#if defined(KOKKOS_OPT_RANGE_AGGRESSIVE_VECTORIZATION) && defined(KOKKOS_HAVE_PRAGMA_IVDEP) && !defined(__CUDA_ARCH__)
#define KOKKOS_MDRANGE_IVDEP
#endif

namespace Kokkos { namespace Experimental {

enum class Iterate
{
  Default, // Default for the device
  Left,    // Left indices stride fastest
  Right,   // Right indices stride fastest
  Flat,    // Do not tile, only valid for inner direction
};

template <typename ExecSpace>
struct default_outer_direction
{
  using type = Iterate;
  static constexpr Iterate value = Iterate::Right;
};

template <typename ExecSpace>
struct default_inner_direction
{
  using type = Iterate;
  static constexpr Iterate value = Iterate::Right;
};


// Iteration Pattern
template < unsigned N
         , Iterate OuterDir = Iterate::Default
         , Iterate InnerDir = Iterate::Default
         >
struct Rank
{
  static_assert( N != 0u, "Kokkos Error: rank 0 undefined");
  static_assert( N != 1u, "Kokkos Error: rank 1 is not a multi-dimensional range");
  static_assert( N < 4u, "Kokkos Error: Unsupported rank...");

  using iteration_pattern = Rank<N, OuterDir, InnerDir>;

  static constexpr int rank = N;
  static constexpr Iterate outer_direction = OuterDir;
  static constexpr Iterate inner_direction = InnerDir;
};



// multi-dimensional iteration pattern
template <typename... Properties>
struct MDRangePolicy
{
  using range_policy = RangePolicy<Properties...>;

  static_assert( !std::is_same<range_policy,void>::value
               , "Kokkos Error: MD iteration pattern not defined" );

  using iteration_pattern   = typename range_policy::iteration_pattern;
  using work_tag            = typename range_policy::work_tag;

  static constexpr int rank = iteration_pattern::rank;

  static constexpr int outer_direction = static_cast<int> (
      (iteration_pattern::outer_direction != Iterate::Default && iteration_pattern::outer_direction != Iterate::Flat)
    ? iteration_pattern::outer_direction
    : default_outer_direction< typename range_policy::execution_space>::value );

  static constexpr int inner_direction = static_cast<int> (
      iteration_pattern::inner_direction != Iterate::Default
    ? iteration_pattern::inner_direction
    : default_inner_direction< typename range_policy::execution_space>::value ) ;


  // Ugly ugly workaround intel 14 not handling scoped enum correctly
  static constexpr int Flat = static_cast<int>( Iterate::Flat );
  static constexpr int Right = static_cast<int>( Iterate::Right );


  using size_type   = typename range_policy::index_type;
  using index_type  = typename std::make_signed<size_type>::type;


  template <typename I>
  MDRangePolicy( std::initializer_list<I> upper_corner )
  {
    static_assert( std::is_integral<I>::value, "Kokkos Error: corner defined with non-integral type" );

    // TODO check size of lists equal to rank
    // static_asserts on initializer_list.size() require c++14

    //static_assert( upper_corner.size() == rank, "Kokkos Error: upper_corner has incorrect rank" );

    const auto u = upper_corner.begin();

    m_num_tiles = 1;
    for (int i=0; i<rank; ++i) {
      m_offset[i] = static_cast<index_type>(0);
      m_dim[i]    = static_cast<index_type>(u[i]);
      if (inner_direction != Flat) {
        // default tile size to 4
        m_tile[i] = 4;
      } else {
        m_tile[i] = 1;
      }
      m_tile_dim[i] = (m_dim[i] + (m_tile[i] - 1)) / m_tile[i];
      m_num_tiles *= m_tile_dim[i];
    }
  }

  template <typename IA, typename IB>
  MDRangePolicy( std::initializer_list<IA> corner_a
               , std::initializer_list<IB> corner_b
               )
  {
    static_assert( std::is_integral<IA>::value, "Kokkos Error: corner A defined with non-integral type" );
    static_assert( std::is_integral<IB>::value, "Kokkos Error: corner B defined with non-integral type" );

    // TODO check size of lists equal to rank
    // static_asserts on initializer_list.size() require c++14
    //static_assert( corner_a.size() == rank, "Kokkos Error: corner_a has incorrect rank" );
    //static_assert( corner_b.size() == rank, "Kokkos Error: corner_b has incorrect rank" );


    using A = typename std::make_signed<IA>::type;
    using B = typename std::make_signed<IB>::type;

    const auto a = [=](int i) { return static_cast<A>(corner_a.begin()[i]); };
    const auto b = [=](int i) { return static_cast<B>(corner_b.begin()[i]); };

    m_num_tiles = 1;
    for (int i=0; i<rank; ++i) {
      m_offset[i] = static_cast<index_type>(a(i) <= b(i) ? a(i) : b(i));
      m_dim[i]    = static_cast<index_type>(a(i) <= b(i) ? b(i) - a(i) : a(i) - b(i));
      if (inner_direction != Flat) {
        // default tile size to 4
        m_tile[i] = 4;
      } else {
        m_tile[i] = 1;
      }
      m_tile_dim[i] = (m_dim[i] + (m_tile[i] - 1)) / m_tile[i];
      m_num_tiles *= m_tile_dim[i];
    }
  }

  template <typename IA, typename IB, typename T>
  MDRangePolicy( std::initializer_list<IA> corner_a
               , std::initializer_list<IB> corner_b
               , std::initializer_list<T> tile
               )
  {
    static_assert( std::is_integral<IA>::value, "Kokkos Error: corner A defined with non-integral type" );
    static_assert( std::is_integral<IB>::value, "Kokkos Error: corner B defined with non-integral type" );
    static_assert( std::is_integral<T>::value, "Kokkos Error: tile defined with non-integral type" );
    static_assert( inner_direction != Flat, "Kokkos Error: tiling not support with flat iteration" );

    // TODO check size of lists equal to rank
    // static_asserts on initializer_list.size() require c++14
    //static_assert( corner_a.size() == rank, "Kokkos Error: corner_a has incorrect rank" );
    //static_assert( corner_b.size() == rank, "Kokkos Error: corner_b has incorrect rank" );
    //static_assert( tile.size() == rank, "Kokkos Error: tile has incorrect rank" );

    using A = typename std::make_signed<IA>::type;
    using B = typename std::make_signed<IB>::type;

    const auto a = [=](int i) { return static_cast<A>(corner_a.begin()[i]); };
    const auto b = [=](int i) { return static_cast<B>(corner_b.begin()[i]); };
    const auto t = tile.begin();

    m_num_tiles = 1;
    for (int i=0; i<rank; ++i) {
      m_offset[i] = static_cast<index_type>(a(i) <= b(i) ? a(i) : b(i));
      m_dim[i]    = static_cast<index_type>(a(i) <= b(i) ? b(i) - a(i) : a(i) - b(i));
      m_tile[i]   = static_cast<int>(t[i] > (T)0 ? t[i] : (T)1 );
      m_tile_dim[i] = (m_dim[i] + (m_tile[i] - 1)) / m_tile[i];
      m_num_tiles *= m_tile_dim[i];
    }
  }

  index_type   m_offset[rank];
  index_type   m_dim[rank];
  int          m_tile[rank];
  index_type   m_tile_dim[rank];
  size_type    m_num_tiles;       // product of tile dims
};

namespace Impl {

// Serial, Threads, OpenMP
// use enable_if to overload for Cuda
template < typename MDRange, typename Functor, typename Enable = void >
struct MDForFunctor
{
  using work_tag   = typename MDRange::work_tag;
  using index_type = typename MDRange::index_type;
  using size_type  = typename MDRange::size_type;

  MDRange m_range;
  Functor m_func;

  KOKKOS_INLINE_FUNCTION
  MDForFunctor( MDRange const& range, Functor const& f )
    : m_range(range)
    , m_func( f )
  {}

  KOKKOS_INLINE_FUNCTION
  MDForFunctor( MDRange const& range, Functor && f )
    : m_range(range)
    , m_func( std::forward<Functor>(f) )
  {}

  KOKKOS_INLINE_FUNCTION
  MDForFunctor( MDRange && range, Functor const& f )
    : m_range( std::forward<MDRange>(range) )
    , m_func( f )
  {}

  KOKKOS_INLINE_FUNCTION
  MDForFunctor( MDRange && range, Functor && f )
    : m_range( std::forward<MDRange>(range) )
    , m_func( std::forward<Functor>(f) )
  {}


  KOKKOS_INLINE_FUNCTION
  MDForFunctor( MDForFunctor const& ) = default;

  KOKKOS_INLINE_FUNCTION
  MDForFunctor& operator=( MDForFunctor const& ) = default;

  KOKKOS_INLINE_FUNCTION
  MDForFunctor( MDForFunctor && ) = default;

  KOKKOS_INLINE_FUNCTION
  MDForFunctor& operator=( MDForFunctor && ) = default;

  // Rank-2, Flat, No Tag
  template <typename Idx>
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<(  std::is_integral<Idx>::value
                          && std::is_same<void, work_tag>::value
                          && MDRange::rank == 2
                          && MDRange::inner_direction == MDRange::Flat
                          )>::type
  operator()(Idx t) const
  {
    if (  MDRange::outer_direction == MDRange::Right ) {
      m_func( m_range.m_offset[0] + ( t / m_range.m_dim[1] )
            , m_range.m_offset[1] + ( t % m_range.m_dim[1] ) );
    } else {
      m_func( m_range.m_offset[0] + ( t % m_range.m_dim[0] )
            , m_range.m_offset[1] + ( t / m_range.m_dim[0] ) );
    }
  }

  // Rank-2, Flat, Tag
  template <typename Idx>
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<(  std::is_integral<Idx>::value
                          && !std::is_same<void, work_tag>::value
                          && MDRange::rank == 2
                          && MDRange::inner_direction == MDRange::Flat
                          )>::type
  operator()(Idx t) const
  {
    if (  MDRange::outer_direction == MDRange::Right ) {
      m_func( work_tag{}, m_range.m_offset[0] + ( t / m_range.m_dim[1] )
            , m_range.m_offset[1] + ( t % m_range.m_dim[1] ) );
    } else {
      m_func( work_tag{}, m_range.m_offset[0] + ( t % m_range.m_dim[0] )
            , m_range.m_offset[1] + ( t / m_range.m_dim[0] ) );
    }
  }

  // Rank-2, Not Flat, No Tag
  template <typename Idx>
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<(  std::is_integral<Idx>::value
                          && std::is_same<void, work_tag>::value
                          && MDRange::rank == 2
                          && MDRange::inner_direction != MDRange::Flat
                          )>::type
  operator()(Idx t) const
  {
    index_type t0, t1;
    if (  MDRange::outer_direction == MDRange::Right ) {
      t0 = t / m_range.m_tile_dim[1];
      t1 = t % m_range.m_tile_dim[1];
    } else {
      t0 = t % m_range.m_tile_dim[0];
      t1 = t / m_range.m_tile_dim[0];
    }

    const index_type b0 = t0 * m_range.m_tile[0] + m_range.m_offset[0];
    const index_type b1 = t1 * m_range.m_tile[1] + m_range.m_offset[1];

    const index_type e0 = b0 + m_range.m_tile[0] <= (m_range.m_dim[0] + m_range.m_offset[0] ) ? b0 + m_range.m_tile[0] : ( m_range.m_dim[0] + m_range.m_offset[0] );
    const index_type e1 = b1 + m_range.m_tile[1] <= (m_range.m_dim[1] + m_range.m_offset[1] ) ? b1 + m_range.m_tile[1] : ( m_range.m_dim[1] + m_range.m_offset[1] );

    if (  MDRange::inner_direction == MDRange::Right ) {
      for (int i0=b0; i0<e0; ++i0) {
      #if defined(KOKKOS_MDRANGE_IVDEP)
      #pragma ivdep
      #endif
      for (int i1=b1; i1<e1; ++i1) {
        m_func( i0, i1 );
      }}
    } else {
      for (int i1=b1; i1<e1; ++i1) {
      #if defined(KOKKOS_MDRANGE_IVDEP)
      #pragma ivdep
      #endif
      for (int i0=b0; i0<e0; ++i0) {
        m_func( i0, i1 );
      }}
    }
  }

  // Rank-2, Not Flat, Tag
  template <typename Idx>
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<(  std::is_integral<Idx>::value
                          && !std::is_same<void, work_tag>::value
                          && MDRange::rank == 2
                          && MDRange::inner_direction != MDRange::Flat
                          )>::type
  operator()(Idx t) const
  {
    work_tag tag;

    index_type t0, t1;
    if (  MDRange::outer_direction == MDRange::Right ) {
      t0 = t / m_range.m_tile_dim[1];
      t1 = t % m_range.m_tile_dim[1];
    } else {
      t0 = t % m_range.m_tile_dim[0];
      t1 = t / m_range.m_tile_dim[0];
    }

    const index_type b0 = t0 * m_range.m_tile[0] + m_range.m_offset[0];
    const index_type b1 = t1 * m_range.m_tile[1] + m_range.m_offset[1];

    const index_type e0 = b0 + m_range.m_tile[0] <= (m_range.m_dim[0] + m_range.m_offset[0] ) ? b0 + m_range.m_tile[0] : ( m_range.m_dim[0] + m_range.m_offset[0] );
    const index_type e1 = b1 + m_range.m_tile[1] <= (m_range.m_dim[1] + m_range.m_offset[1] ) ? b1 + m_range.m_tile[1] : ( m_range.m_dim[1] + m_range.m_offset[1] );

    if (  MDRange::inner_direction == MDRange::Right ) {
      for (int i0=b0; i0<e0; ++i0) {
      #if defined(KOKKOS_MDRANGE_IVDEP)
      #pragma ivdep
      #endif
      for (int i1=b1; i1<e1; ++i1) {
        m_func( tag, i0, i1 );
      }}
    } else {
      for (int i1=b1; i1<e1; ++i1) {
      #if defined(KOKKOS_MDRANGE_IVDEP)
      #pragma ivdep
      #endif
      for (int i0=b0; i0<e0; ++i0) {
        m_func( tag, i0, i1 );
      }}
    }
  }

  //---------------------------------------------------------------------------

  // Rank-3, Flat, No Tag
  template <typename Idx>
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<(  std::is_integral<Idx>::value
                          && std::is_same<void, work_tag>::value
                          && MDRange::rank == 3
                          && MDRange::inner_direction == MDRange::Flat
                          )>::type
  operator()(Idx t) const
  {
    if (  MDRange::outer_direction == MDRange::Right ) {
    const int64_t tmp_prod = m_range.m_dim[1]*m_range.m_dim[2];
    m_func( m_range.m_offset[0] + (  t / tmp_prod )
          , m_range.m_offset[1] + ( (t % tmp_prod) / m_range.m_dim[2] )
          , m_range.m_offset[2] + ( (t % tmp_prod) % m_range.m_dim[2] )
          );
    } else {
    const int64_t tmp_prod = m_range.m_dim[0]*m_range.m_dim[1];
    m_func( m_range.m_offset[0] + ( (t % tmp_prod) % m_range.m_dim[0] )
          , m_range.m_offset[1] + ( (t % tmp_prod) / m_range.m_dim[0] )
          , m_range.m_offset[2] + (  t / tmp_prod )
          );
    }
  }

  // Rank-3, Flat, Tag
  template <typename Idx>
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<(  std::is_integral<Idx>::value
                          && !std::is_same<void, work_tag>::value
                          && MDRange::rank == 3
                          && MDRange::inner_direction == MDRange::Flat
                          )>::type
  operator()(Idx t) const
  {
    if (  MDRange::outer_direction == MDRange::Right ) {
      const int64_t tmp_prod = m_range.m_dim[1]*m_range.m_dim[2];
      m_func( work_tag{}
            , m_range.m_offset[0] + (  t / tmp_prod )
            , m_range.m_offset[1] + ( (t % tmp_prod) / m_range.m_dim[2] )
            , m_range.m_offset[2] + ( (t % tmp_prod) % m_range.m_dim[2] )
            );
    } else {
      const int64_t tmp_prod = m_range.m_dim[0]*m_range.m_dim[1];
      m_func( work_tag{}
            , m_range.m_offset[0] + ( (t % tmp_prod) % m_range.m_dim[0] )
            , m_range.m_offset[1] + ( (t % tmp_prod) / m_range.m_dim[0] )
            , m_range.m_offset[2] + (  t / tmp_prod )
            );
    }
  }

  // Rank-3, Not Flat, No Tag
  template <typename Idx>
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<(  std::is_integral<Idx>::value
                          && std::is_same<void, work_tag>::value
                          && MDRange::rank == 3
                          && MDRange::inner_direction != MDRange::Flat
                          )>::type
  operator()(Idx t) const
  {
    index_type t0, t1, t2;
    if (  MDRange::outer_direction == MDRange::Right ) {
      const index_type tmp_prod = ( m_range.m_tile_dim[1]*m_range.m_tile_dim[2]);
      t0 = t / tmp_prod;
      t1 = ( t % tmp_prod ) / m_range.m_tile_dim[2];
      t2 = ( t % tmp_prod ) % m_range.m_tile_dim[2];
    } else {
      const index_type tmp_prod = ( m_range.m_tile_dim[0]*m_range.m_tile_dim[1]);
      t0 = ( t % tmp_prod ) % m_range.m_tile_dim[0];
      t1 = ( t % tmp_prod ) / m_range.m_tile_dim[0];
      t2 = t / tmp_prod;
    }

    const index_type b0 = t0 * m_range.m_tile[0] + m_range.m_offset[0];
    const index_type b1 = t1 * m_range.m_tile[1] + m_range.m_offset[1];
    const index_type b2 = t2 * m_range.m_tile[2] + m_range.m_offset[2];

    const index_type e0 = b0 + m_range.m_tile[0] <= (m_range.m_dim[0] + m_range.m_offset[0] ) ? b0 + m_range.m_tile[0] : ( m_range.m_dim[0] + m_range.m_offset[0] );
    const index_type e1 = b1 + m_range.m_tile[1] <= (m_range.m_dim[1] + m_range.m_offset[1] ) ? b1 + m_range.m_tile[1] : ( m_range.m_dim[1] + m_range.m_offset[1] );
    const index_type e2 = b2 + m_range.m_tile[2] <= (m_range.m_dim[2] + m_range.m_offset[2] ) ? b2 + m_range.m_tile[2] : ( m_range.m_dim[2] + m_range.m_offset[2] );

    if (  MDRange::inner_direction == MDRange::Right ) {
      for (int i0=b0; i0<e0; ++i0) {
      for (int i1=b1; i1<e1; ++i1) {
      #if defined(KOKKOS_MDRANGE_IVDEP)
      #pragma ivdep
      #endif
      for (int i2=b2; i2<e2; ++i2) {
        m_func( i0, i1, i2 );
      }}}
    } else {
      for (int i2=b2; i2<e2; ++i2) {
      for (int i1=b1; i1<e1; ++i1) {
      #if defined(KOKKOS_MDRANGE_IVDEP)
      #pragma ivdep
      #endif
      for (int i0=b0; i0<e0; ++i0) {
        m_func( i0, i1, i2 );
      }}}
    }
  }

  // Rank-3, Not Flat, Tag
  template <typename Idx>
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<(  std::is_integral<Idx>::value
                          && !std::is_same<void, work_tag>::value
                          && MDRange::rank == 3
                          && MDRange::inner_direction != MDRange::Flat
                          )>::type
  operator()(Idx t) const
  {
    work_tag tag;

    index_type t0, t1, t2;
    if (  MDRange::outer_direction == MDRange::Right ) {
      const index_type tmp_prod = ( m_range.m_tile_dim[1]*m_range.m_tile_dim[2]);
      t0 = t / tmp_prod;
      t1 = ( t % tmp_prod ) / m_range.m_tile_dim[2];
      t2 = ( t % tmp_prod ) % m_range.m_tile_dim[2];
    } else {
      const index_type tmp_prod = ( m_range.m_tile_dim[0]*m_range.m_tile_dim[1]);
      t0 = ( t % tmp_prod ) % m_range.m_tile_dim[0];
      t1 = ( t % tmp_prod ) / m_range.m_tile_dim[0];
      t2 = t / tmp_prod;
    }

    const index_type b0 = t0 * m_range.m_tile[0] + m_range.m_offset[0];
    const index_type b1 = t1 * m_range.m_tile[1] + m_range.m_offset[1];
    const index_type b2 = t2 * m_range.m_tile[2] + m_range.m_offset[2];

    const index_type e0 = b0 + m_range.m_tile[0] <= (m_range.m_dim[0] + m_range.m_offset[0] ) ? b0 + m_range.m_tile[0] : ( m_range.m_dim[0] + m_range.m_offset[0] );
    const index_type e1 = b1 + m_range.m_tile[1] <= (m_range.m_dim[1] + m_range.m_offset[1] ) ? b1 + m_range.m_tile[1] : ( m_range.m_dim[1] + m_range.m_offset[1] );
    const index_type e2 = b2 + m_range.m_tile[2] <= (m_range.m_dim[2] + m_range.m_offset[2] ) ? b2 + m_range.m_tile[2] : ( m_range.m_dim[2] + m_range.m_offset[2] );

    if (  MDRange::inner_direction == MDRange::Right ) {
      for (int i0=b0; i0<e0; ++i0) {
      for (int i1=b1; i1<e1; ++i1) {
      #if defined(KOKKOS_MDRANGE_IVDEP)
      #pragma ivdep
      #endif
      for (int i2=b2; i2<e2; ++i2) {
        m_func( tag, i0, i1, i2 );
      }}}
    } else {
      for (int i2=b2; i2<e2; ++i2) {
      for (int i1=b1; i1<e1; ++i1) {
      #if defined(KOKKOS_MDRANGE_IVDEP)
      #pragma ivdep
      #endif
      for (int i0=b0; i0<e0; ++i0) {
        m_func( tag, i0, i1, i2 );
      }}}
    }
  }
};



} // namespace Impl


template <typename MDRange, typename Functor>
void md_parallel_for( MDRange const& range
                    , Functor const& f
                    , const std::string& str = ""
                    )
{
  Impl::MDForFunctor<MDRange, Functor> g(range, f);

  using range_policy = typename MDRange::range_policy;

  Kokkos::parallel_for( range_policy(0, range.m_num_tiles).set_chunk_size(1), g, str );
}

template <typename MDRange, typename Functor>
void md_parallel_for( const std::string& str
                    , MDRange const& range
                    , Functor const& f
                    )
{
  Impl::MDForFunctor<MDRange, Functor> g(range, f);

  using range_policy = typename MDRange::range_policy;

  Kokkos::parallel_for( range_policy(0, range.m_num_tiles).set_chunk_size(1), g, str );
}

}} // namespace Kokkos::Experimental

#endif //KOKKOS_CORE_EXP_MD_RANGE_POLICY_HPP

