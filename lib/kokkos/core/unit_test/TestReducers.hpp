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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <limits>

#include <Kokkos_Core.hpp>


//--------------------------------------------------------------------------

namespace Test {

struct ReducerTag {};

template< class Scalar, class ExecSpace = Kokkos::DefaultExecutionSpace >
struct TestReducers {
  struct SumFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value += values( i );
    }
  };

  struct ProdFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value *= values( i );
    }
  };

  struct MinFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      if ( values( i ) < value ) value = values( i );
    }
  };

  struct MaxFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      if ( values( i ) > value ) value = values( i );
    }
  };

  struct MinLocFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, typename Kokkos::MinLoc< Scalar, int >::value_type & value ) const {
      if ( values( i ) < value.val ) {
        value.val = values( i );
        value.loc = i;
      }
    }
  };

  struct MaxLocFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, typename Kokkos::MaxLoc< Scalar, int >::value_type & value ) const {
      if ( values( i ) > value.val ) {
        value.val = values( i );
        value.loc = i;
      }
    }
  };

  struct MinMaxLocFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, typename Kokkos::MinMaxLoc< Scalar, int >::value_type & value ) const {
      if ( values( i ) > value.max_val ) {
        value.max_val = values( i );
        value.max_loc = i;
      }

      if ( values( i ) < value.min_val ) {
        value.min_val = values( i );
        value.min_loc = i;
      }
    }
  };

  struct BAndFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value = value & values( i );
    }
  };

  struct BOrFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value = value | values( i );
    }
  };

  struct LAndFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value = value && values( i );
    }
  };

  struct LOrFunctor {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const int & i, Scalar & value ) const {
      value = value || values( i );
    }
  };

  struct SumFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, Scalar & value ) const {
      value += values( i );
    }
  };

  struct ProdFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, Scalar & value ) const {
      value *= values( i );
    }
  };

  struct MinFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, Scalar & value ) const {
      if ( values( i ) < value ) value = values( i );
    }
  };

  struct MaxFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, Scalar & value ) const {
      if ( values( i ) > value ) value = values( i );
    }
  };

  struct MinLocFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, typename Kokkos::MinLoc< Scalar, int >::value_type & value ) const {
      if ( values( i ) < value.val ) {
        value.val = values( i );
        value.loc = i;
      }
    }
  };

  struct MaxLocFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, typename Kokkos::MaxLoc< Scalar, int >::value_type & value ) const {
      if ( values( i ) > value.val ) {
        value.val = values( i );
        value.loc = i;
      }
    }
  };

  struct MinMaxLocFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, typename Kokkos::MinMaxLoc< Scalar, int >::value_type & value ) const {
      if ( values( i ) > value.max_val ) {
        value.max_val = values( i );
        value.max_loc = i;
      }

      if ( values( i ) < value.min_val ) {
        value.min_val = values( i );
        value.min_loc = i;
      }
    }
  };

  struct BAndFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, Scalar & value ) const {
      value = value & values( i );
    }
  };

  struct BOrFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, Scalar & value ) const {
      value = value | values( i );
    }
  };

  struct LAndFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, Scalar & value ) const {
      value = value && values( i );
    }
  };

  struct LOrFunctorTag {
    Kokkos::View< const Scalar*, ExecSpace > values;

    KOKKOS_INLINE_FUNCTION
    void operator()( const ReducerTag, const int & i, Scalar & value ) const {
      value = value || values( i );
    }
  };
  static void test_sum( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_sum = 0;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100 );
      reference_sum += h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    SumFunctor f;
    f.values = values;
    SumFunctorTag f_tag;
    f_tag.values = values;
    Scalar init = 0;

    {
      Scalar sum_scalar = init;
      Kokkos::Sum< Scalar > reducer_scalar( sum_scalar );
      
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
      ASSERT_EQ( sum_scalar, reference_sum );
     
      sum_scalar = init;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
      ASSERT_EQ( sum_scalar, reference_sum );

      Scalar sum_scalar_view = reducer_scalar.reference();
      ASSERT_EQ( sum_scalar_view, reference_sum );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace> sum_view( "View" );
      sum_view() = init;
      Kokkos::Sum< Scalar > reducer_view( sum_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
      Kokkos::fence();

      Scalar sum_view_scalar = sum_view();
      ASSERT_EQ( sum_view_scalar, reference_sum );

      Scalar sum_view_view = reducer_view.reference();
      ASSERT_EQ( sum_view_view, reference_sum );
    }
  }

  static void test_prod( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_prod = 1;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 4 + 1 );
      reference_prod *= h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    ProdFunctor f;
    f.values = values;
    ProdFunctorTag f_tag;
    f_tag.values = values;
    Scalar init = 1;

    {
      Scalar prod_scalar = init;
      Kokkos::Prod< Scalar > reducer_scalar( prod_scalar );
   
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
      ASSERT_EQ( prod_scalar, reference_prod );
      
      prod_scalar = init;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
      ASSERT_EQ( prod_scalar, reference_prod );

      Scalar prod_scalar_view = reducer_scalar.reference();
      ASSERT_EQ( prod_scalar_view, reference_prod );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > prod_view( "View" );
      prod_view() = init;
      Kokkos::Prod< Scalar > reducer_view( prod_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
      Kokkos::fence();

      Scalar prod_view_scalar = prod_view();
      ASSERT_EQ( prod_view_scalar, reference_prod );

      Scalar prod_view_view = reducer_view.reference();
      ASSERT_EQ( prod_view_view, reference_prod );
    }
  }

  static void test_min( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_min = std::numeric_limits< Scalar >::max();

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100000 );

      if ( h_values( i ) < reference_min ) reference_min = h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    MinFunctor f;
    f.values = values;
    MinFunctorTag f_tag;
    f_tag.values = values;
    Scalar init = std::numeric_limits< Scalar >::max();

    {
      Scalar min_scalar = init;
      Kokkos::Min< Scalar > reducer_scalar( min_scalar );
     
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
      ASSERT_EQ( min_scalar, reference_min );

      min_scalar = init;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
      ASSERT_EQ( min_scalar, reference_min );

      Scalar min_scalar_view = reducer_scalar.reference();
      ASSERT_EQ( min_scalar_view, reference_min );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > min_view( "View" );
      min_view() = init;
      Kokkos::Min< Scalar > reducer_view( min_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
      Kokkos::fence();

      Scalar min_view_scalar = min_view();
      ASSERT_EQ( min_view_scalar, reference_min );

      Scalar min_view_view = reducer_view.reference();
      ASSERT_EQ( min_view_view, reference_min );
    }
  }

  static void test_max( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_max = std::numeric_limits< Scalar >::min();

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100000 + 1 );

      if ( h_values( i ) > reference_max ) reference_max = h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    MaxFunctor f;
    f.values = values;
    MaxFunctorTag f_tag;
    f_tag.values = values;
    Scalar init = std::numeric_limits< Scalar >::min();

    {
      Scalar max_scalar = init;
      Kokkos::Max< Scalar > reducer_scalar( max_scalar );

      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
      ASSERT_EQ( max_scalar, reference_max );

      max_scalar = init;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
      ASSERT_EQ( max_scalar, reference_max );

      Scalar max_scalar_view = reducer_scalar.reference();
      ASSERT_EQ( max_scalar_view, reference_max );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > max_view( "View" );
      max_view() = init;
      Kokkos::Max< Scalar > reducer_view( max_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
      Kokkos::fence();

      Scalar max_view_scalar = max_view();
      ASSERT_EQ( max_view_scalar, reference_max );

      Scalar max_view_view = reducer_view.reference();
      ASSERT_EQ( max_view_view, reference_max );
    }
  }

  static void test_minloc( int N ) {
    typedef typename Kokkos::MinLoc< Scalar, int >::value_type value_type;

    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_min = std::numeric_limits< Scalar >::max();
    int reference_loc = -1;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100000 + 2 );

      if ( h_values( i ) < reference_min ) {
        reference_min = h_values( i );
        reference_loc = i;
      }
      else if ( h_values( i ) == reference_min ) {
        // Make min unique.
        h_values( i ) += Scalar(1);
      }
    }
    Kokkos::deep_copy( values, h_values );

    MinLocFunctor f;
    f.values = values;
    MinLocFunctorTag f_tag;
    f_tag.values = values;

    {
      value_type min_scalar;
      Kokkos::MinLoc< Scalar, int > reducer_scalar( min_scalar );

      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
      ASSERT_EQ( min_scalar.val, reference_min );
      ASSERT_EQ( min_scalar.loc, reference_loc );

      min_scalar = value_type();
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
      ASSERT_EQ( min_scalar.val, reference_min );
      ASSERT_EQ( min_scalar.loc, reference_loc );

      value_type min_scalar_view = reducer_scalar.reference();
      ASSERT_EQ( min_scalar_view.val, reference_min );
      ASSERT_EQ( min_scalar_view.loc, reference_loc );
    }

    {
      Kokkos::View< value_type, Kokkos::HostSpace > min_view( "View" );
      Kokkos::MinLoc< Scalar, int > reducer_view( min_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
      Kokkos::fence();

      value_type min_view_scalar = min_view();
      ASSERT_EQ( min_view_scalar.val, reference_min );
      ASSERT_EQ( min_view_scalar.loc, reference_loc );

      value_type min_view_view = reducer_view.reference();
      ASSERT_EQ( min_view_view.val, reference_min );
      ASSERT_EQ( min_view_view.loc, reference_loc );
    }
  }

  static void test_maxloc( int N ) {
    typedef typename Kokkos::MaxLoc< Scalar, int >::value_type value_type;

    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_max = std::numeric_limits< Scalar >::min();
    int reference_loc = -1;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100000 + 2 );

      if ( h_values( i ) > reference_max ) {
        reference_max = h_values( i );
        reference_loc = i;
      }
      else if ( h_values( i ) == reference_max ) {
        // Make max unique.
        h_values( i ) -= Scalar(1);
      }
    }
    Kokkos::deep_copy( values, h_values );

    MaxLocFunctor f;
    f.values = values;
    MaxLocFunctorTag f_tag;
    f_tag.values = values;

    {
      value_type max_scalar;
      Kokkos::MaxLoc< Scalar, int > reducer_scalar( max_scalar );

      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
      ASSERT_EQ( max_scalar.val, reference_max );
      ASSERT_EQ( max_scalar.loc, reference_loc );

      max_scalar = value_type();
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
      ASSERT_EQ( max_scalar.val, reference_max );
      ASSERT_EQ( max_scalar.loc, reference_loc );

      value_type max_scalar_view = reducer_scalar.reference();
      ASSERT_EQ( max_scalar_view.val, reference_max );
      ASSERT_EQ( max_scalar_view.loc, reference_loc );
    }

    {
      Kokkos::View< value_type, Kokkos::HostSpace > max_view( "View" );
      Kokkos::MaxLoc< Scalar, int > reducer_view( max_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
      Kokkos::fence();

      value_type max_view_scalar = max_view();
      ASSERT_EQ( max_view_scalar.val, reference_max );
      ASSERT_EQ( max_view_scalar.loc, reference_loc );

      value_type max_view_view = reducer_view.reference();
      ASSERT_EQ( max_view_view.val, reference_max );
      ASSERT_EQ( max_view_view.loc, reference_loc );
    }
  }

  static void test_minmaxloc( int N ) {
     typedef typename Kokkos::MinMaxLoc< Scalar, int >::value_type value_type;

     Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
     auto h_values = Kokkos::create_mirror_view( values );
     Scalar reference_max = std::numeric_limits< Scalar >::min();
     Scalar reference_min = std::numeric_limits< Scalar >::max();
     int reference_minloc = -1;
     int reference_maxloc = -1;

     for ( int i = 0; i < N; i++ ) {
       h_values( i ) = (Scalar) ( rand() % 100000 + 2);
     }

     for ( int i = 0; i < N; i++ ) {
       if ( h_values( i ) > reference_max ) {
         reference_max = h_values( i );
         reference_maxloc = i;
       }
       else if ( h_values( i ) == reference_max ) {
         // Make max unique.
         h_values( i ) -= Scalar(1);
       }
     }

     for ( int i = 0; i < N; i++ ) {
       if ( h_values( i ) < reference_min ) {
         reference_min = h_values( i );
         reference_minloc = i;
       }
       else if ( h_values( i ) == reference_min ) {
         // Make min unique.
         h_values( i ) += Scalar(1);
       }
     }

     Kokkos::deep_copy( values, h_values );

     MinMaxLocFunctor f;
     f.values = values;
     MinMaxLocFunctorTag f_tag;
     f_tag.values = values;

     {
       value_type minmax_scalar;
       Kokkos::MinMaxLoc< Scalar, int > reducer_scalar( minmax_scalar );

       Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
       ASSERT_EQ( minmax_scalar.min_val, reference_min );

       for ( int i = 0; i < N; i++ ) {
         if ( ( i == minmax_scalar.min_loc ) && ( h_values( i ) == reference_min ) ) {
           reference_minloc = i;
         }
       }

       ASSERT_EQ( minmax_scalar.min_loc, reference_minloc );
       ASSERT_EQ( minmax_scalar.max_val, reference_max );

       for ( int i = 0; i < N; i++ ) {
         if ( ( i == minmax_scalar.max_loc ) && ( h_values( i ) == reference_max ) ) {
           reference_maxloc = i;
         }
       }

       ASSERT_EQ( minmax_scalar.max_loc, reference_maxloc );

       minmax_scalar = value_type();
       Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
       ASSERT_EQ( minmax_scalar.min_val, reference_min );

       for ( int i = 0; i < N; i++ ) {
         if ( ( i == minmax_scalar.min_loc ) && ( h_values( i ) == reference_min ) ) {
           reference_minloc = i;
         }
       }

       ASSERT_EQ( minmax_scalar.min_loc, reference_minloc );
       ASSERT_EQ( minmax_scalar.max_val, reference_max );

       for ( int i = 0; i < N; i++ ) {
         if ( ( i == minmax_scalar.max_loc ) && ( h_values( i ) == reference_max ) ) {
           reference_maxloc = i;
         }
       }

       ASSERT_EQ( minmax_scalar.max_loc, reference_maxloc );

       value_type minmax_scalar_view = reducer_scalar.reference();
       ASSERT_EQ( minmax_scalar_view.min_val, reference_min );
       ASSERT_EQ( minmax_scalar_view.min_loc, reference_minloc );
       ASSERT_EQ( minmax_scalar_view.max_val, reference_max );
       ASSERT_EQ( minmax_scalar_view.max_loc, reference_maxloc );
     }

     {
       Kokkos::View< value_type, Kokkos::HostSpace > minmax_view( "View" );
       Kokkos::MinMaxLoc< Scalar, int > reducer_view( minmax_view );
       Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
       Kokkos::fence();

       value_type minmax_view_scalar = minmax_view();
       ASSERT_EQ( minmax_view_scalar.min_val, reference_min );
       ASSERT_EQ( minmax_view_scalar.min_loc, reference_minloc );
       ASSERT_EQ( minmax_view_scalar.max_val, reference_max );
       ASSERT_EQ( minmax_view_scalar.max_loc, reference_maxloc );

       value_type minmax_view_view = reducer_view.reference();
       ASSERT_EQ( minmax_view_view.min_val, reference_min );
       ASSERT_EQ( minmax_view_view.min_loc, reference_minloc );
       ASSERT_EQ( minmax_view_view.max_val, reference_max );
       ASSERT_EQ( minmax_view_view.max_loc, reference_maxloc );
     }
   }

  static void test_BAnd( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_band = Scalar() | ( ~Scalar() );

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 100000 + 1 );
      reference_band = reference_band & h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    BAndFunctor f;
    f.values = values;
    BAndFunctorTag f_tag;
    f_tag.values = values;
    Scalar init = Scalar() | ( ~Scalar() );

    {
      Scalar band_scalar = init;
      Kokkos::BAnd< Scalar > reducer_scalar( band_scalar );

      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
      ASSERT_EQ( band_scalar, reference_band );

      band_scalar = init;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
      ASSERT_EQ( band_scalar, reference_band );

      Scalar band_scalar_view = reducer_scalar.reference();

      ASSERT_EQ( band_scalar_view, reference_band );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > band_view( "View" );
      band_view() = init;
      Kokkos::BAnd< Scalar > reducer_view( band_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
      Kokkos::fence();

      Scalar band_view_scalar = band_view();
      ASSERT_EQ( band_view_scalar, reference_band );

      Scalar band_view_view = reducer_view.reference();
      ASSERT_EQ( band_view_view, reference_band );
    }
  }

  static void test_BOr( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_bor = Scalar() & ( ~Scalar() );

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( ( rand() % 100000 + 1 ) * 2 );
      reference_bor = reference_bor | h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    BOrFunctor f;
    f.values = values;
    BOrFunctorTag f_tag;
    f_tag.values = values;
    Scalar init = Scalar() & ( ~Scalar() );

    {
      Scalar bor_scalar = init;
      Kokkos::BOr< Scalar > reducer_scalar( bor_scalar );

      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
      ASSERT_EQ( bor_scalar, reference_bor );

      bor_scalar = init;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
      ASSERT_EQ( bor_scalar, reference_bor );

      Scalar bor_scalar_view = reducer_scalar.reference();
      ASSERT_EQ( bor_scalar_view, reference_bor );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > bor_view( "View" );
      bor_view() = init;
      Kokkos::BOr< Scalar > reducer_view( bor_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
      Kokkos::fence();

      Scalar bor_view_scalar = bor_view();
      ASSERT_EQ( bor_view_scalar, reference_bor );

      Scalar bor_view_view = reducer_view.reference();
      ASSERT_EQ( bor_view_view, reference_bor );
    }
  }

  static void test_LAnd( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_land = 1;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 2 );
      reference_land = reference_land && h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    LAndFunctor f;
    f.values = values;
    LAndFunctorTag f_tag;
    f_tag.values = values;
    Scalar init = 1;

    {
      Scalar land_scalar = init;
      Kokkos::LAnd< Scalar > reducer_scalar( land_scalar );

      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
      ASSERT_EQ( land_scalar, reference_land );

      land_scalar = init;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
      ASSERT_EQ( land_scalar, reference_land );

      Scalar land_scalar_view = reducer_scalar.reference();
      ASSERT_EQ( land_scalar_view, reference_land );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > land_view( "View" );
      land_view() = init;
      Kokkos::LAnd< Scalar > reducer_view( land_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
      Kokkos::fence();

      Scalar land_view_scalar = land_view();
      ASSERT_EQ( land_view_scalar, reference_land );

      Scalar land_view_view = reducer_view.reference();
      ASSERT_EQ( land_view_view, reference_land );
    }
  }

  static void test_LOr( int N ) {
    Kokkos::View< Scalar*, ExecSpace > values( "Values", N );
    auto h_values = Kokkos::create_mirror_view( values );
    Scalar reference_lor = 0;

    for ( int i = 0; i < N; i++ ) {
      h_values( i ) = (Scalar) ( rand() % 2 );
      reference_lor = reference_lor || h_values( i );
    }
    Kokkos::deep_copy( values, h_values );

    LOrFunctor f;
    f.values = values;
    LOrFunctorTag f_tag;
    f_tag.values = values;
    Scalar init = 0;

    {
      Scalar lor_scalar = init;
      Kokkos::LOr< Scalar > reducer_scalar( lor_scalar );

      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_scalar );
      ASSERT_EQ( lor_scalar, reference_lor );

      lor_scalar = init;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , ReducerTag >( 0, N ), f_tag, reducer_scalar );
      ASSERT_EQ( lor_scalar, reference_lor );

      Scalar lor_scalar_view = reducer_scalar.reference();
      ASSERT_EQ( lor_scalar_view, reference_lor );
    }

    {
      Kokkos::View< Scalar, Kokkos::HostSpace > lor_view( "View" );
      lor_view() = init;
      Kokkos::LOr< Scalar > reducer_view( lor_view );
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, N ), f, reducer_view );
      Kokkos::fence();

      Scalar lor_view_scalar = lor_view();
      ASSERT_EQ( lor_view_scalar, reference_lor );

      Scalar lor_view_view = reducer_view.reference();
      ASSERT_EQ( lor_view_view, reference_lor );
    }
  }

  static void execute_float() {
    test_sum( 10001 );
    test_prod( 35 );
    test_min( 10003 );
    test_minloc( 10003 );
    test_max( 10007 );
    test_maxloc( 10007 );
    test_minmaxloc( 10007 );
  }

  static void execute_integer() {
    test_sum( 10001 );
    test_prod( 35 );
    test_min( 10003 );
    test_minloc( 10003 );
    test_max( 10007 );
    test_maxloc( 10007 );
    test_minmaxloc( 10007 );
    test_BAnd( 35 );
    test_BOr( 35 );
    test_LAnd( 35 );
    test_LOr( 35 );
  }

  static void execute_basic() {
    test_sum( 10001 );
    test_prod( 35 );
  }
};

} // namespace Test
