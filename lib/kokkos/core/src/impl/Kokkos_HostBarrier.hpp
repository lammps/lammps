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

#ifndef KOKKOS_HOST_BARRIER_HPP
#define KOKKOS_HOST_BARRIER_HPP

#include <cstddef>
#include <cstdint>

namespace Kokkos { namespace Impl {

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

enum : int { RENDEZVOUS_ALIGNMENT = 128
           , RENDEZVOUS_HEADER    = RENDEZVOUS_ALIGNMENT
           };

inline constexpr int rendezvous_buffer_size( const int nthreads ) noexcept
{
  return RENDEZVOUS_HEADER + ((nthreads-1 + RENDEZVOUS_ALIGNMENT-1) / RENDEZVOUS_ALIGNMENT) * RENDEZVOUS_ALIGNMENT;
}

void rendezvous_initialize( volatile void * buffer
                          , const int       size
                          , const int       rank
                          ) noexcept;


bool rendezvous( volatile void * buffer
               , uint64_t &      step
               , const int       size
               , const int       rank
               , bool            active_wait = true
               ) noexcept;

void rendezvous_release( volatile void * buffer
                       , const uint64_t  step
                       ) noexcept;


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


class HostBarrier
{
public:

  enum : int { ALIGNMENT = RENDEZVOUS_ALIGNMENT };
  enum : int { HEADER    = ALIGNMENT};

  enum Policy : int { ACTIVE, PASSIVE };

  inline static constexpr int buffer_size( const int nthreads ) noexcept
  {
    return rendezvous_buffer_size(nthreads);
  }

  HostBarrier( volatile void * arg_buffer
             , int             arg_size
             , int             arg_rank
             , Policy          arg_policy
             ) noexcept
    : m_buffer{arg_buffer}
    , m_size{arg_size}
    , m_rank{arg_rank}
    , m_policy{arg_policy}
    , m_step{0}
  {
    rendezvous_initialize( m_buffer, m_size, m_rank );
  }

  bool rendezvous() const noexcept
  {
    return Kokkos::Impl::rendezvous( m_buffer
                                   , m_step
                                   , m_size
                                   , m_rank
                                   , m_policy == ACTIVE
                                   );
  }

  void rendezvous_release() const noexcept
  {
    Kokkos::Impl::rendezvous_release( m_buffer, m_step );
  }

private:
  volatile void *   m_buffer ;
  const int         m_size   ;
  const int         m_rank   ;
  const Policy      m_policy ;
  mutable uint64_t  m_step   ;

private:
  HostBarrier( const HostBarrier &  )             = delete;
  HostBarrier(       HostBarrier && )             = delete;
  HostBarrier & operator=( const HostBarrier &  ) = delete;
  HostBarrier & operator=(       HostBarrier && ) = delete;
};

}} // namespace Kokkos::Impl

#endif // KOKKOS_HOST_BARRIER_HPP

