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

#include <stdio.h>
#include <limits>
#include <iostream>
#include <vector>
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>
#include <iostream>
#include <impl/Kokkos_CPUDiscovery.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>

#ifdef KOKKOS_ENABLE_OPENMPTARGET

namespace Kokkos {
namespace Impl {
namespace {

KOKKOS_INLINE_FUNCTION
int kokkos_omp_in_parallel();

int kokkos_omp_in_critical_region =
    (Kokkos::HostSpace::register_in_parallel(kokkos_omp_in_parallel), 0);

KOKKOS_INLINE_FUNCTION
int kokkos_omp_in_parallel() {
#ifndef __CUDA_ARCH__
  return omp_in_parallel() && !kokkos_omp_in_critical_region;
#else
  return 0;
#endif
}

bool s_using_hwloc = false;

}  // namespace
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Experimental {
bool OpenMPTarget::m_is_initialized = false;
}
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

// int OpenMPTargetExec::m_map_rank[ OpenMPTargetExec::MAX_THREAD_COUNT ] = { 0
// };

// int OpenMPTargetExec::m_pool_topo[ 4 ] = { 0 };

// OpenMPTargetExec * OpenMPTargetExec::m_pool[
// OpenMPTargetExec::MAX_THREAD_COUNT ] = { 0 };

void OpenMPTargetExec::verify_is_process(const char* const label) {
  if (omp_in_parallel()) {
    std::string msg(label);
    msg.append(" ERROR: in parallel");
    Kokkos::Impl::throw_runtime_exception(msg);
  }
}

void OpenMPTargetExec::verify_initialized(const char* const label) {
  if (0 == Kokkos::Experimental::OpenMPTarget::is_initialized()) {
    std::string msg(label);
    msg.append(" ERROR: not initialized");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  if (omp_get_max_threads() !=
      Kokkos::Experimental::OpenMPTarget::thread_pool_size(0)) {
    std::string msg(label);
    msg.append(" ERROR: Initialized but threads modified inappropriately");
    Kokkos::Impl::throw_runtime_exception(msg);
  }
}

void* OpenMPTargetExec::m_scratch_ptr    = NULL;
int64_t OpenMPTargetExec::m_scratch_size = 0;

void OpenMPTargetExec::clear_scratch() {
  Kokkos::Experimental::OpenMPTargetSpace space;
  space.deallocate(m_scratch_ptr, m_scratch_size);
  m_scratch_ptr  = NULL;
  m_scratch_size = NULL;
}

void* OpenMPTargetExec::get_scratch_ptr() { return m_scratch_ptr; }

void OpenMPTargetExec::resize_scratch(int64_t reduce_bytes,
                                      int64_t team_reduce_bytes,
                                      int64_t team_shared_bytes,
                                      int64_t thread_local_bytes) {
  Kokkos::Experimental::OpenMPTargetSpace space;
  uint64_t total_size =
      MAX_ACTIVE_TEAMS * reduce_bytes +         // Inter Team Reduction
      MAX_ACTIVE_TEAMS * team_reduce_bytes +    // Intra Team Reduction
      MAX_ACTIVE_TEAMS * team_shared_bytes +    // Team Local Scratch
      MAX_ACTIVE_THREADS * thread_local_bytes;  // Thread Private Scratch

  if (total_size > m_scratch_size) {
    space.deallocate(m_scratch_ptr, m_scratch_size);
    m_scratch_size = total_size;
    m_scratch_ptr  = space.allocate(total_size);
  }
}
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
//----------------------------------------------------------------------------

int OpenMPTarget::is_initialized() {
  return m_is_initialized;
}  // != Impl::OpenMPTargetExec::m_pool[0]; }

void OpenMPTarget::initialize(unsigned thread_count, unsigned use_numa_count,
                              unsigned use_cores_per_numa) {
  // Before any other call to OMP query the maximum number of threads
  // and save the value for re-initialization unit testing.

  // Init the array for used for arbitrarily sized atomics
  Kokkos::Impl::init_lock_array_host_space();

#ifdef KOKKOS_ENABLE_PROFILING
  Kokkos::Profiling::initialize();
#endif
  m_is_initialized = true;
}

//----------------------------------------------------------------------------

void OpenMPTarget::finalize() {
  Kokkos::Impl::OpenMPTargetExec::verify_initialized("OpenMPTarget::finalize");
  Kokkos::Impl::OpenMPTargetExec::verify_is_process("OpenMPTarget::finalize");

  m_is_initialized = false;

  omp_set_num_threads(1);

  if (Kokkos::Impl::s_using_hwloc && Kokkos::hwloc::can_bind_threads()) {
    hwloc::unbind_this_thread();
  }

#ifdef KOKKOS_ENABLE_PROFILING
  Kokkos::Profiling::finalize();
#endif
}

//----------------------------------------------------------------------------

void OpenMPTarget::print_configuration(std::ostream& s, const bool detail) {
  Kokkos::Impl::OpenMPTargetExec::verify_is_process(
      "OpenMPTarget::print_configuration");
  /*
    s << "Kokkos::Experimental::OpenMPTarget" ;

  #if defined( KOKKOS_ENABLE_OPENMPTARGET )
    s << " KOKKOS_ENABLE_OPENMPTARGET" ;
  #endif
  #if defined( KOKKOS_ENABLE_HWLOC )

    const unsigned numa_count_       =
  Kokkos::hwloc::get_available_numa_count(); const unsigned cores_per_numa   =
  Kokkos::hwloc::get_available_cores_per_numa(); const unsigned threads_per_core
  = Kokkos::hwloc::get_available_threads_per_core();

    s << " hwloc[" << numa_count_ << "x" << cores_per_numa << "x" <<
  threads_per_core << "]"
      << " hwloc_binding_" << ( Impl::s_using_hwloc ? "enabled" : "disabled" )
      ;
  #endif

    const bool is_initialized = 0 != Impl::OpenMPTargetExec::m_pool[0] ;

    if ( is_initialized ) {
      const int numa_count      = Kokkos::Impl::OpenMPTargetExec::m_pool_topo[0]
  / Kokkos::Impl::OpenMPTargetExec::m_pool_topo[1] ; const int core_per_numa   =
  Kokkos::Impl::OpenMPTargetExec::m_pool_topo[1] /
  Kokkos::Impl::OpenMPTargetExec::m_pool_topo[2] ; const int thread_per_core =
  Kokkos::Impl::OpenMPTargetExec::m_pool_topo[2] ;

      s << " thread_pool_topology[ " << numa_count
        << " x " << core_per_numa
        << " x " << thread_per_core
        << " ]"
        << std::endl ;

      if ( detail ) {
        std::vector< std::pair<unsigned,unsigned> > coord(
  Kokkos::Impl::OpenMPTargetExec::m_pool_topo[0] );

  #pragma omp parallel
        {
  #pragma omp critical
          {
            coord[ omp_get_thread_num() ] = hwloc::get_this_thread_coordinate();
          }
  // END #pragma omp critical
        }
  // END #pragma omp parallel

        for ( unsigned i = 0 ; i < coord.size() ; ++i ) {
          s << "  thread omp_rank[" << i << "]"
            << " kokkos_rank[" << Impl::OpenMPTargetExec::m_map_rank[ i ] << "]"
            << " hwloc_coord[" << coord[i].first << "." << coord[i].second <<
  "]"
            << std::endl ;
        }
      }
    }
    else {
      s << " not initialized" << std::endl ;
    }
  */
}

int OpenMPTarget::concurrency() { return thread_pool_size(0); }

const char* OpenMPTarget::name() { return "OpenMPTarget"; }
}  // namespace Experimental
}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_OPENMPTARGET
