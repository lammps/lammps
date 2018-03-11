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

#ifndef KOKKOS_OPENMP_WORKGRAPHPOLICY_HPP
#define KOKKOS_OPENMP_WORKGRAPHPOLICY_HPP

namespace Kokkos {
namespace Impl {

template< class FunctorType , class ... Traits >
class ParallelFor< FunctorType ,
                   Kokkos::WorkGraphPolicy< Traits ... > ,
                   Kokkos::OpenMP
                 >
{
private:

  typedef Kokkos::WorkGraphPolicy< Traits ... > Policy ;

  Policy       m_policy ;
  FunctorType  m_functor ;

  template< class TagType >
  typename std::enable_if< std::is_same< TagType , void >::value >::type
  exec_one( const std::int32_t w ) const noexcept
    { m_functor( w ); }

  template< class TagType >
  typename std::enable_if< ! std::is_same< TagType , void >::value >::type
  exec_one( const std::int32_t w ) const noexcept
    { const TagType t{} ; m_functor( t , w ); }

public:

  inline
  void execute()
  {
    const int pool_size = OpenMP::thread_pool_size();

    #pragma omp parallel num_threads(pool_size)
    {
      // Spin until COMPLETED_TOKEN.
      // END_TOKEN indicates no work is currently available.

      for ( std::int32_t w = Policy::END_TOKEN ;
            Policy::COMPLETED_TOKEN != ( w = m_policy.pop_work() ) ; ) {
        if ( Policy::END_TOKEN != w ) {
          exec_one< typename Policy::work_tag >( w );
          m_policy.completed_work(w);
        }
      }
    }
  }

  inline
  ParallelFor( const FunctorType & arg_functor
             , const Policy      & arg_policy )
    : m_policy( arg_policy )
    , m_functor( arg_functor )
  {}
};

} // namespace Impl
} // namespace Kokkos

#endif /* #define KOKKOS_OPENMP_WORKGRAPHPOLICY_HPP */

