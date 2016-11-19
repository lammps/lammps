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

#include <Kokkos_HostSpace.hpp>

#include <impl/Kokkos_HBWAllocators.hpp>
#include <impl/Kokkos_Error.hpp>


#include <stdint.h>    // uintptr_t
#include <cstdlib>     // for malloc, realloc, and free
#include <cstring>     // for memcpy

#if defined(KOKKOS_POSIX_MEMALIGN_AVAILABLE)
#include <sys/mman.h>  // for mmap, munmap, MAP_ANON, etc
#include <unistd.h>    // for sysconf, _SC_PAGE_SIZE, _SC_PHYS_PAGES
#endif

#include <sstream>
#include <iostream>

#ifdef KOKKOS_HAVE_HBWSPACE
#include <memkind.h>

namespace Kokkos {
namespace Experimental {
namespace Impl {
#define MEMKIND_TYPE MEMKIND_HBW //hbw_get_kind(HBW_PAGESIZE_4KB)
/*--------------------------------------------------------------------------*/

void* HBWMallocAllocator::allocate( size_t size )
{
  std::cout<< "Allocate HBW: " << 1.0e-6*size << "MB" << std::endl;
  void * ptr = NULL;
  if (size) {
    ptr = memkind_malloc(MEMKIND_TYPE,size);

    if (!ptr)
    {
      std::ostringstream msg ;
      msg << name() << ": allocate(" << size << ") FAILED";
      Kokkos::Impl::throw_runtime_exception( msg.str() );
    }
  }
  return ptr;
}

void HBWMallocAllocator::deallocate( void * ptr, size_t /*size*/ )
{
  if (ptr) {
    memkind_free(MEMKIND_TYPE,ptr);
  }
}

void * HBWMallocAllocator::reallocate(void * old_ptr, size_t /*old_size*/, size_t new_size)
{
  void * ptr = memkind_realloc(MEMKIND_TYPE, old_ptr, new_size);

  if (new_size > 0u && ptr == NULL) {
    Kokkos::Impl::throw_runtime_exception("Error: Malloc Allocator could not reallocate memory");
  }
  return ptr;
}

} // namespace Impl
} // namespace Experimental
} // namespace Kokkos
#endif
