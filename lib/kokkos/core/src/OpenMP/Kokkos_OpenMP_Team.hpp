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

#ifndef KOKKOS_OPENMP_TEAM_HPP
#define KOKKOS_OPENMP_TEAM_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_OPENMP )

#include <OpenMP/Kokkos_OpenMP_Exec.hpp>

namespace Kokkos { namespace Impl {

template< class ... Properties >
class TeamPolicyInternal< Kokkos::OpenMP, Properties ... >: public PolicyTraits<Properties ...>
{
public:

  //! Tag this class as a kokkos execution policy
  typedef TeamPolicyInternal      execution_policy ;

  typedef PolicyTraits<Properties ... > traits;

  TeamPolicyInternal& operator = (const TeamPolicyInternal& p) {
    m_league_size = p.m_league_size;
    m_team_size = p.m_team_size;
    m_team_alloc = p.m_team_alloc;
    m_team_iter = p.m_team_iter;
    m_team_scratch_size[0] = p.m_team_scratch_size[0];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_team_scratch_size[1] = p.m_team_scratch_size[1];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size = p.m_chunk_size;
    return *this;
  }

  //----------------------------------------

  template< class FunctorType >
  inline static
  int team_size_max( const FunctorType & ) {
      int pool_size = traits::execution_space::thread_pool_size(1);
      int max_host_team_size =  Impl::HostThreadTeamData::max_team_members;
      return pool_size<max_host_team_size?pool_size:max_host_team_size;
    }

  template< class FunctorType >
  inline static
  int team_size_recommended( const FunctorType & )
    { return traits::execution_space::thread_pool_size(2); }

  template< class FunctorType >
  inline static
  int team_size_recommended( const FunctorType &, const int& )
    { return traits::execution_space::thread_pool_size(2); }

  //----------------------------------------

private:

  int m_league_size ;
  int m_team_size ;
  int m_team_alloc ;
  int m_team_iter ;

  size_t m_team_scratch_size[2];
  size_t m_thread_scratch_size[2];

  int m_chunk_size;

  inline void init( const int league_size_request
                  , const int team_size_request )
    {
      const int pool_size  = traits::execution_space::thread_pool_size(0);
      const int max_host_team_size =  Impl::HostThreadTeamData::max_team_members;
      const int team_max   = pool_size<max_host_team_size?pool_size:max_host_team_size;
      const int team_grain = traits::execution_space::thread_pool_size(2);

      m_league_size = league_size_request ;

      m_team_size = team_size_request < team_max ?
                    team_size_request : team_max ;

      // Round team size up to a multiple of 'team_gain'
      const int team_size_grain = team_grain * ( ( m_team_size + team_grain - 1 ) / team_grain );
      const int team_count      = pool_size / team_size_grain ;

      // Constraint : pool_size = m_team_alloc * team_count
      m_team_alloc = pool_size / team_count ;

      // Maxumum number of iterations each team will take:
      m_team_iter  = ( m_league_size + team_count - 1 ) / team_count ;

      set_auto_chunk_size();
    }

public:

  inline int team_size()   const { return m_team_size ; }
  inline int league_size() const { return m_league_size ; }

  inline size_t scratch_size(const int& level, int team_size_ = -1) const {
    if(team_size_ < 0) team_size_ = m_team_size;
    return m_team_scratch_size[level] + team_size_*m_thread_scratch_size[level] ;
  }

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal( typename traits::execution_space &
            , int league_size_request
            , int team_size_request
            , int /* vector_length_request */ = 1 )
            : m_team_scratch_size { 0 , 0 }
            , m_thread_scratch_size { 0 , 0 }
            , m_chunk_size(0)
    { init( league_size_request , team_size_request ); }

  TeamPolicyInternal( typename traits::execution_space &
            , int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , int /* vector_length_request */ = 1)
            : m_team_scratch_size { 0 , 0 }
            , m_thread_scratch_size { 0 , 0 }
            , m_chunk_size(0)
    { init( league_size_request , traits::execution_space::thread_pool_size(2) ); }

  TeamPolicyInternal( int league_size_request
            , int team_size_request
            , int /* vector_length_request */ = 1 )
            : m_team_scratch_size { 0 , 0 }
            , m_thread_scratch_size { 0 , 0 }
            , m_chunk_size(0)
    { init( league_size_request , team_size_request ); }

  TeamPolicyInternal( int league_size_request
            , const Kokkos::AUTO_t & /* team_size_request */
            , int /* vector_length_request */ = 1 )
            : m_team_scratch_size { 0 , 0 }
            , m_thread_scratch_size { 0 , 0 }
            , m_chunk_size(0)
    { init( league_size_request , traits::execution_space::thread_pool_size(2) ); }

  inline int team_alloc() const { return m_team_alloc ; }
  inline int team_iter()  const { return m_team_iter ; }

  inline int chunk_size() const { return m_chunk_size ; }

  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal set_chunk_size(typename traits::index_type chunk_size_) const {
    TeamPolicyInternal p = *this;
    p.m_chunk_size = chunk_size_;
    return p;
  }

  inline TeamPolicyInternal set_scratch_size(const int& level, const PerTeamValue& per_team) const {
    TeamPolicyInternal p = *this;
    p.m_team_scratch_size[level] = per_team.value;
    return p;
  };

  inline TeamPolicyInternal set_scratch_size(const int& level, const PerThreadValue& per_thread) const {
    TeamPolicyInternal p = *this;
    p.m_thread_scratch_size[level] = per_thread.value;
    return p;
  };

  inline TeamPolicyInternal set_scratch_size(const int& level, const PerTeamValue& per_team, const PerThreadValue& per_thread) const {
    TeamPolicyInternal p = *this;
    p.m_team_scratch_size[level] = per_team.value;
    p.m_thread_scratch_size[level] = per_thread.value;
    return p;
  };

private:
  /** \brief finalize chunk_size if it was set to AUTO*/
  inline void set_auto_chunk_size() {

    int concurrency = traits::execution_space::thread_pool_size(0)/m_team_alloc;
    if( concurrency==0 ) concurrency=1;

    if(m_chunk_size > 0) {
      if(!Impl::is_integral_power_of_two( m_chunk_size ))
        Kokkos::abort("TeamPolicy blocking granularity must be power of two" );
    }

    int new_chunk_size = 1;
    while(new_chunk_size*100*concurrency < m_league_size)
      new_chunk_size *= 2;
    if(new_chunk_size < 128) {
      new_chunk_size = 1;
      while( (new_chunk_size*40*concurrency < m_league_size ) && (new_chunk_size<128) )
        new_chunk_size*=2;
    }
    m_chunk_size = new_chunk_size;
  }

public:
  typedef Impl::HostThreadTeamMember< Kokkos::OpenMP > member_type ;
};

}} // namespace Kokkos::Impl

#endif
#endif /* KOKKOS_OPENMP_TEAM_HPP */


