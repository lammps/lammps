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

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ATOMIC_HPP ) && ! defined( KOKKOS_MEMORY_FENCE_HPP )
#define KOKKOS_MEMORY_FENCE_HPP

#include <atomic>

namespace Kokkos {

//----------------------------------------------------------------------------

KOKKOS_FORCEINLINE_FUNCTION
void memory_fence()
{
#if defined( __CUDA_ARCH__ )
  __threadfence();
#else
  std::atomic_thread_fence( std::memory_order_seq_cst );
#endif
}

//////////////////////////////////////////////////////
// store_fence()
//
// If possible use a store fence on the architecture, if not run a full memory fence

KOKKOS_FORCEINLINE_FUNCTION
void store_fence()
{
#if defined( __CUDA_ARCH__ )
  __threadfence();
#else
  std::atomic_thread_fence( std::memory_order_seq_cst );
#endif
}

//////////////////////////////////////////////////////
// load_fence()
//
// If possible use a load fence on the architecture, if not run a full memory fence

KOKKOS_FORCEINLINE_FUNCTION
void load_fence()
{
#if defined( __CUDA_ARCH__ )
  __threadfence();
#else
  std::atomic_thread_fence( std::memory_order_seq_cst );
#endif
}

} // namespace kokkos

#endif

