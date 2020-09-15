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

#ifndef KOKKOS_CUDA_ALLOCATION_TRACKING_HPP
#define KOKKOS_CUDA_ALLOCATION_TRACKING_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <impl/Kokkos_Traits.hpp>

namespace Kokkos {
namespace Impl {

template <class DestructFunctor>
SharedAllocationRecord* shared_allocation_record(
    Kokkos::CudaSpace const& arg_space, void* const arg_alloc_ptr,
    DestructFunctor const& arg_destruct) {
  SharedAllocationRecord* const record =
      SharedAllocationRecord::get_record(arg_alloc_ptr);

  // assert: record != 0

  // assert: sizeof(DestructFunctor) <= record->m_destruct_size

  // assert: record->m_destruct_function == 0

  DestructFunctor* const functor = reinterpret_cast<DestructFunctor*>(
      reinterpret_cast<uintptr_t>(record) + sizeof(SharedAllocationRecord));

  new (functor) DestructFunctor(arg_destruct);

  record->m_destruct_functor = &shared_allocation_destroy<DestructFunctor>;

  return record;
}

/// class CudaUnmanagedAllocator
/// does nothing when deallocate(ptr,size) is called
struct CudaUnmanagedAllocator {
  static const char* name() { return "Cuda Unmanaged Allocator"; }

  static void deallocate(void* /*ptr*/, size_t /*size*/) {}

  static bool support_texture_binding() { return true; }
};

/// class CudaUnmanagedAllocator
/// does nothing when deallocate(ptr,size) is called
struct CudaUnmanagedUVMAllocator {
  static const char* name() { return "Cuda Unmanaged UVM Allocator"; }

  static void deallocate(void* /*ptr*/, size_t /*size*/) {}

  static bool support_texture_binding() { return true; }
};

/// class CudaUnmanagedHostAllocator
/// does nothing when deallocate(ptr,size) is called
class CudaUnmanagedHostAllocator {
 public:
  static const char* name() { return "Cuda Unmanaged Host Allocator"; }
  // Unmanaged deallocate does nothing
  static void deallocate(void* /*ptr*/, size_t /*size*/) {}
};

/// class CudaMallocAllocator
class CudaMallocAllocator {
 public:
  static const char* name() { return "Cuda Malloc Allocator"; }

  static void* allocate(size_t size);

  static void deallocate(void* ptr, size_t);

  static void* reallocate(void* old_ptr, size_t old_size, size_t new_size);

  static bool support_texture_binding() { return true; }
};

/// class CudaUVMAllocator
class CudaUVMAllocator {
 public:
  static const char* name() { return "Cuda UVM Allocator"; }

  static void* allocate(size_t size);

  static void deallocate(void* ptr, size_t);

  static void* reallocate(void* old_ptr, size_t old_size, size_t new_size);

  static bool support_texture_binding() { return true; }
};

/// class CudaHostAllocator
class CudaHostAllocator {
 public:
  static const char* name() { return "Cuda Host Allocator"; }

  static void* allocate(size_t size);

  static void deallocate(void* ptr, size_t);

  static void* reallocate(void* old_ptr, size_t old_size, size_t new_size);
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_CUDA

#endif  // #ifndef KOKKOS_CUDA_ALLOCATION_TRACKING_HPP
