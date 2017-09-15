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

#ifndef KOKKOS_IMPL_ERROR_HPP
#define KOKKOS_IMPL_ERROR_HPP

#include <string>
#include <iosfwd>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA
#include <Cuda/Kokkos_Cuda_abort.hpp>
#endif

namespace Kokkos {
namespace Impl {

void host_abort( const char * const );

void throw_runtime_exception( const std::string & );

void traceback_callstack( std::ostream & );

std::string human_memory_size(size_t arg_bytes);

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


namespace Kokkos {
KOKKOS_INLINE_FUNCTION
void abort( const char * const message ) {
#ifdef __CUDA_ARCH__
  Kokkos::Impl::cuda_abort(message);
#else
  #if !defined(KOKKOS_ENABLE_OPENMPTARGET) && !defined(__HCC_ACCELERATOR__)
    Kokkos::Impl::host_abort(message);
  #endif
#endif
}

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_ERROR_HPP */

