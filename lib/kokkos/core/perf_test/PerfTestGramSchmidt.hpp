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

#include <cmath>
#include <PerfTestBlasKernels.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Test {

// Reduction   : result = dot( Q(:,j) , Q(:,j) );
// PostProcess : R(j,j) = result ; inv = 1 / result ;
template< class VectorView , class ValueView  >
struct InvNorm2 : public Kokkos::DotSingle< VectorView > {

  typedef typename Kokkos::DotSingle< VectorView >::value_type value_type ;

  ValueView  Rjj ;
  ValueView  inv ;

  InvNorm2( const VectorView & argX ,
            const ValueView  & argR ,
            const ValueView  & argInv )
    : Kokkos::DotSingle< VectorView >( argX )
    , Rjj( argR )
    , inv( argInv )
    {}

  KOKKOS_INLINE_FUNCTION
  void final( value_type & result ) const
  {
    result = sqrt( result );
    Rjj() = result ;
    inv() = ( 0 < result ) ? 1.0 / result : 0 ;
  }
};

template< class VectorView , class ValueView >
inline
void invnorm2( const VectorView & x ,
               const ValueView  & r ,
               const ValueView  & r_inv )
{
  Kokkos::parallel_reduce( x.dimension_0() , InvNorm2< VectorView , ValueView >( x , r , r_inv ) );
}

// PostProcess : tmp = - ( R(j,k) = result );
template< class VectorView , class ValueView  >
struct DotM : public Kokkos::Dot< VectorView > {

  typedef typename Kokkos::Dot< VectorView >::value_type value_type ;

  ValueView  Rjk ;
  ValueView  tmp ;

  DotM( const VectorView & argX ,
        const VectorView & argY ,
        const ValueView & argR ,
        const ValueView & argTmp )
    : Kokkos::Dot< VectorView >( argX , argY )
    , Rjk( argR )
    , tmp( argTmp )
    {}

  KOKKOS_INLINE_FUNCTION
  void final( value_type & result ) const
  {
     Rjk()  = result ;
     tmp()  = - result ;
  }
};

template< class VectorView , class ValueView >
inline
void dot_neg( const VectorView & x ,
              const VectorView & y ,
              const ValueView  & r ,
              const ValueView  & r_neg )
{
  Kokkos::parallel_reduce( x.dimension_0() , DotM< VectorView , ValueView >( x , y , r , r_neg ) );
}


template< typename Scalar , class DeviceType >
struct ModifiedGramSchmidt
{
  typedef DeviceType  execution_space ;
  typedef typename execution_space::size_type  size_type ;

  typedef Kokkos::View< Scalar** ,
                        Kokkos::LayoutLeft ,
                        execution_space > multivector_type ;

  typedef Kokkos::View< Scalar* ,
                        Kokkos::LayoutLeft ,
                        execution_space > vector_type ;

  typedef Kokkos::View< Scalar ,
                        Kokkos::LayoutLeft ,
                        execution_space > value_view ;


  multivector_type Q ;
  multivector_type R ;

  static double factorization( const multivector_type Q_ ,
                               const multivector_type R_ )
  {
    const size_type count  = Q_.dimension_1();
    value_view tmp("tmp");
    value_view one("one");

    Kokkos::deep_copy( one , (Scalar) 1 );

    Kokkos::Impl::Timer timer ;

    for ( size_type j = 0 ; j < count ; ++j ) {
      // Reduction   : tmp = dot( Q(:,j) , Q(:,j) );
      // PostProcess : tmp = sqrt( tmp ); R(j,j) = tmp ; tmp = 1 / tmp ;
      const vector_type Qj  = Kokkos::subview( Q_ , Kokkos::ALL() , j );
      const value_view  Rjj = Kokkos::subview( R_ , j , j );

      invnorm2( Qj , Rjj , tmp );

      // Q(:,j) *= ( 1 / R(j,j) ); => Q(:,j) *= tmp ;
      Kokkos::scale( tmp , Qj );

      for ( size_t k = j + 1 ; k < count ; ++k ) {
        const vector_type Qk = Kokkos::subview( Q_ , Kokkos::ALL() , k );
        const value_view  Rjk = Kokkos::subview( R_ , j , k );

        // Reduction   : R(j,k) = dot( Q(:,j) , Q(:,k) );
        // PostProcess : tmp = - R(j,k);
        dot_neg( Qj , Qk , Rjk , tmp );

        // Q(:,k) -= R(j,k) * Q(:,j); => Q(:,k) += tmp * Q(:,j)
        Kokkos::axpby( tmp , Qj , one , Qk );
      }
    }

    execution_space::fence();

    return timer.seconds();
  }

  //--------------------------------------------------------------------------

  static double test( const size_t length ,
                      const size_t count ,
                      const size_t iter = 1 )
  {
    multivector_type Q_( "Q" , length , count );
    multivector_type R_( "R" , count , count );

    typename multivector_type::HostMirror A =
      Kokkos::create_mirror( Q_ );

    // Create and fill A on the host

    for ( size_type j = 0 ; j < count ; ++j ) {
      for ( size_type i = 0 ; i < length ; ++i ) {
        A(i,j) = ( i + 1 ) * ( j + 1 );
      }
    }

    double dt_min = 0 ;

    for ( size_t i = 0 ; i < iter ; ++i ) {

      Kokkos::deep_copy( Q_ , A );

      // A = Q * R

      const double dt = factorization( Q_ , R_ );

      if ( 0 == i ) dt_min = dt ;
      else dt_min = dt < dt_min ? dt : dt_min ;
    }

    return dt_min ;
  }
};

}

