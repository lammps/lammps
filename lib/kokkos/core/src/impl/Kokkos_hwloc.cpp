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

#define DEBUG_PRINT 0

#include <iostream>
#include <sstream>
#include <algorithm>

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_hwloc.hpp>
#include <impl/Kokkos_Error.hpp>

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace hwloc {

/* Return 0 if asynchronous, 1 if synchronous and include process. */
unsigned thread_mapping(const char* const label, const bool allow_async,
                        unsigned& thread_count, unsigned& use_numa_count,
                        unsigned& use_cores_per_numa,
                        std::pair<unsigned, unsigned> threads_coord[]) {
  const bool hwloc_avail = Kokkos::hwloc::available();
  const unsigned avail_numa_count =
      hwloc_avail ? hwloc::get_available_numa_count() : 1;
  const unsigned avail_cores_per_numa =
      hwloc_avail ? hwloc::get_available_cores_per_numa() : thread_count;
  const unsigned avail_threads_per_core =
      hwloc_avail ? hwloc::get_available_threads_per_core() : 1;

  // (numa,core) coordinate of the process:
  const std::pair<unsigned, unsigned> proc_coord =
      Kokkos::hwloc::get_this_thread_coordinate();

  //------------------------------------------------------------------------
  // Defaults for unspecified inputs:

  if (!use_numa_count) {
    // Default to use all NUMA regions
    use_numa_count = !thread_count
                         ? avail_numa_count
                         : (thread_count < avail_numa_count ? thread_count
                                                            : avail_numa_count);
  }

  if (!use_cores_per_numa) {
    // Default to use all but one core if asynchronous, all cores if
    // synchronous.
    const unsigned threads_per_numa = thread_count / use_numa_count;

    use_cores_per_numa =
        !threads_per_numa
            ? avail_cores_per_numa - (allow_async ? 1 : 0)
            : (threads_per_numa < avail_cores_per_numa ? threads_per_numa
                                                       : avail_cores_per_numa);
  }

  if (!thread_count) {
    thread_count = use_numa_count * use_cores_per_numa * avail_threads_per_core;
  }

  //------------------------------------------------------------------------
  // Input verification:

  const bool valid_numa = use_numa_count <= avail_numa_count;
  const bool valid_cores =
      use_cores_per_numa && use_cores_per_numa <= avail_cores_per_numa;
  const bool valid_threads =
      thread_count && thread_count <= use_numa_count * use_cores_per_numa *
                                          avail_threads_per_core;
  const bool balanced_numa = !(thread_count % use_numa_count);
  const bool balanced_cores =
      !(thread_count % (use_numa_count * use_cores_per_numa));

  const bool valid_input = valid_numa && valid_cores && valid_threads &&
                           balanced_numa && balanced_cores;

  if (!valid_input) {
    std::ostringstream msg;

    msg << label << " HWLOC ERROR(s)";

    if (!valid_threads) {
      msg << " : thread_count(" << thread_count << ") exceeds capacity("
          << use_numa_count * use_cores_per_numa * avail_threads_per_core
          << ")";
    }
    if (!valid_numa) {
      msg << " : use_numa_count(" << use_numa_count << ") exceeds capacity("
          << avail_numa_count << ")";
    }
    if (!valid_cores) {
      msg << " : use_cores_per_numa(" << use_cores_per_numa
          << ") exceeds capacity(" << avail_cores_per_numa << ")";
    }
    if (!balanced_numa) {
      msg << " : thread_count(" << thread_count << ") imbalanced among numa("
          << use_numa_count << ")";
    }
    if (!balanced_cores) {
      msg << " : thread_count(" << thread_count << ") imbalanced among cores("
          << use_numa_count * use_cores_per_numa << ")";
    }

    Kokkos::Impl::throw_runtime_exception(msg.str());
  }

  const unsigned thread_spawn_synchronous =
      (allow_async && 1 < thread_count &&
       (use_numa_count < avail_numa_count ||
        use_cores_per_numa < avail_cores_per_numa))
          ? 0 /* asyncronous */
          : 1 /* synchronous, threads_coord[0] is process core */;

  // Determine binding coordinates for to-be-spawned threads so that
  // threads may be bound to cores as they are spawned.

  const unsigned threads_per_core =
      thread_count / (use_numa_count * use_cores_per_numa);

  if (thread_spawn_synchronous) {
    // Working synchronously and include process core as threads_coord[0].
    // Swap the NUMA coordinate of the process core with 0
    // Swap the CORE coordinate of the process core with 0
    for (unsigned i = 0, inuma = avail_numa_count - use_numa_count;
         inuma < avail_numa_count; ++inuma) {
      const unsigned numa_coord = 0 == inuma
                                      ? proc_coord.first
                                      : (proc_coord.first == inuma ? 0 : inuma);
      for (unsigned icore = avail_cores_per_numa - use_cores_per_numa;
           icore < avail_cores_per_numa; ++icore) {
        const unsigned core_coord =
            0 == icore ? proc_coord.second
                       : (proc_coord.second == icore ? 0 : icore);
        for (unsigned ith = 0; ith < threads_per_core; ++ith, ++i) {
          threads_coord[i].first  = numa_coord;
          threads_coord[i].second = core_coord;
        }
      }
    }
  } else if (use_numa_count < avail_numa_count) {
    // Working asynchronously and omit the process' NUMA region from the pool.
    // Swap the NUMA coordinate of the process core with ( ( avail_numa_count -
    // use_numa_count ) - 1 )
    const unsigned numa_coord_swap = (avail_numa_count - use_numa_count) - 1;
    for (unsigned i = 0, inuma = avail_numa_count - use_numa_count;
         inuma < avail_numa_count; ++inuma) {
      const unsigned numa_coord =
          proc_coord.first == inuma ? numa_coord_swap : inuma;
      for (unsigned icore = avail_cores_per_numa - use_cores_per_numa;
           icore < avail_cores_per_numa; ++icore) {
        const unsigned core_coord = icore;
        for (unsigned ith = 0; ith < threads_per_core; ++ith, ++i) {
          threads_coord[i].first  = numa_coord;
          threads_coord[i].second = core_coord;
        }
      }
    }
  } else if (use_cores_per_numa < avail_cores_per_numa) {
    // Working asynchronously and omit the process' core from the pool.
    // Swap the CORE coordinate of the process core with ( (
    // avail_cores_per_numa - use_cores_per_numa ) - 1 )
    const unsigned core_coord_swap =
        (avail_cores_per_numa - use_cores_per_numa) - 1;
    for (unsigned i = 0, inuma = avail_numa_count - use_numa_count;
         inuma < avail_numa_count; ++inuma) {
      const unsigned numa_coord = inuma;
      for (unsigned icore = avail_cores_per_numa - use_cores_per_numa;
           icore < avail_cores_per_numa; ++icore) {
        const unsigned core_coord =
            proc_coord.second == icore ? core_coord_swap : icore;
        for (unsigned ith = 0; ith < threads_per_core; ++ith, ++i) {
          threads_coord[i].first  = numa_coord;
          threads_coord[i].second = core_coord;
        }
      }
    }
  }

  return thread_spawn_synchronous;
}

} /* namespace hwloc */
} /* namespace Kokkos */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#if defined(KOKKOS_ENABLE_HWLOC)

#include <iostream>
#include <sstream>
#include <stdexcept>

/*--------------------------------------------------------------------------*/
/* Third Party Libraries */

/* Hardware locality library: http://www.open-mpi.org/projects/hwloc/ */
#include <hwloc.h>

#define REQUIRED_HWLOC_API_VERSION 0x000010300

#if HWLOC_API_VERSION < REQUIRED_HWLOC_API_VERSION
#error \
    "Requires  http://www.open-mpi.org/projects/hwloc/  Version 1.3 or greater"
#endif

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace hwloc {
namespace {

#if DEBUG_PRINT

inline void print_bitmap(std::ostream& s, const hwloc_const_bitmap_t bitmap) {
  s << "{";
  for (int i = hwloc_bitmap_first(bitmap); - 1 != i;
       i     = hwloc_bitmap_next(bitmap, i)) {
    s << " " << i;
  }
  s << " }";
}

#endif

enum { MAX_CORE = 1024 };

std::pair<unsigned, unsigned> s_core_topology(0, 0);
unsigned s_core_capacity(0);
hwloc_topology_t s_hwloc_topology(0);
hwloc_bitmap_t s_hwloc_location(0);
hwloc_bitmap_t s_process_binding(0);
hwloc_bitmap_t s_core[MAX_CORE];
bool s_can_bind_threads(true);

struct Sentinel {
  ~Sentinel();
  Sentinel();
};

bool sentinel() {
  static Sentinel self;

  if (0 == s_hwloc_topology) {
    std::cerr << "Kokkos::hwloc ERROR : Called after return from main()"
              << std::endl;
    std::cerr.flush();
  }

  return 0 != s_hwloc_topology;
}

Sentinel::~Sentinel() {
  hwloc_topology_destroy(s_hwloc_topology);
  hwloc_bitmap_free(s_process_binding);
  hwloc_bitmap_free(s_hwloc_location);

  s_core_topology.first  = 0;
  s_core_topology.second = 0;
  s_core_capacity        = 0;
  s_hwloc_topology       = 0;
  s_hwloc_location       = 0;
  s_process_binding      = 0;
}

Sentinel::Sentinel() {
#if defined(__MIC__)
  static const bool remove_core_0 = true;
#else
  static const bool remove_core_0 = false;
#endif

  s_core_topology   = std::pair<unsigned, unsigned>(0, 0);
  s_core_capacity   = 0;
  s_hwloc_topology  = 0;
  s_hwloc_location  = 0;
  s_process_binding = 0;

  for (unsigned i = 0; i < MAX_CORE; ++i) s_core[i] = 0;

  hwloc_topology_init(&s_hwloc_topology);
  hwloc_topology_load(s_hwloc_topology);

  s_hwloc_location  = hwloc_bitmap_alloc();
  s_process_binding = hwloc_bitmap_alloc();

  hwloc_get_cpubind(s_hwloc_topology, s_process_binding, HWLOC_CPUBIND_PROCESS);

  if (hwloc_bitmap_iszero(s_process_binding)) {
    if (Kokkos::show_warnings()) {
      std::cerr << "WARNING: Cannot detect process binding -- ASSUMING ALL "
                   "processing units"
                << std::endl;
    }
    const int pu_depth = hwloc_get_type_depth(s_hwloc_topology, HWLOC_OBJ_PU);
    int num_pu         = 1;
    if (pu_depth != HWLOC_TYPE_DEPTH_UNKNOWN) {
      num_pu = hwloc_get_nbobjs_by_depth(s_hwloc_topology, pu_depth);
    } else {
      if (Kokkos::show_warnings()) {
        std::cerr << "WARNING: Cannot detect number of processing units -- "
                     "ASSUMING 1 (serial)."
                  << std::endl;
      }
      num_pu = 1;
    }
    hwloc_bitmap_set_range(s_process_binding, 0, num_pu - 1);
    s_can_bind_threads = false;
  }

  if (remove_core_0) {
    const hwloc_obj_t core =
        hwloc_get_obj_by_type(s_hwloc_topology, HWLOC_OBJ_CORE, 0);

    if (hwloc_bitmap_intersects(s_process_binding, core->cpuset)) {
      hwloc_bitmap_t s_process_no_core_zero = hwloc_bitmap_alloc();

      hwloc_bitmap_andnot(s_process_no_core_zero, s_process_binding,
                          core->cpuset);

      bool ok =
          0 == hwloc_set_cpubind(s_hwloc_topology, s_process_no_core_zero,
                                 HWLOC_CPUBIND_PROCESS | HWLOC_CPUBIND_STRICT);

      if (ok) {
        hwloc_get_cpubind(s_hwloc_topology, s_process_binding,
                          HWLOC_CPUBIND_PROCESS);

        ok = 0 !=
             hwloc_bitmap_isequal(s_process_binding, s_process_no_core_zero);
      }

      hwloc_bitmap_free(s_process_no_core_zero);

      if (Kokkos::show_warnings() && !ok) {
        std::cerr << "WARNING: Kokkos::hwloc attempted and failed to move "
                     "process off of core #0"
                  << std::endl;
      }
    }
  }

  // Choose a hwloc object type for the NUMA level, which may not exist.

  hwloc_obj_type_t root_type = HWLOC_OBJ_TYPE_MAX;

  {
    // Object types to search, in order.
    static const hwloc_obj_type_t candidate_root_type[] = {
        HWLOC_OBJ_NODE /* NUMA region     */
        ,
        HWLOC_OBJ_SOCKET /* hardware socket */
        ,
        HWLOC_OBJ_MACHINE /* local machine   */
    };

    enum {
      CANDIDATE_ROOT_TYPE_COUNT =
          sizeof(candidate_root_type) / sizeof(hwloc_obj_type_t)
    };

    for (int k = 0;
         k < CANDIDATE_ROOT_TYPE_COUNT && HWLOC_OBJ_TYPE_MAX == root_type;
         ++k) {
      if (0 <
          hwloc_get_nbobjs_by_type(s_hwloc_topology, candidate_root_type[k])) {
        root_type = candidate_root_type[k];
      }
    }
  }

  // Determine which of these 'root' types are available to this process.
  // The process may have been bound (e.g., by MPI) to a subset of these root
  // types. Determine current location of the master (calling) process>

  hwloc_bitmap_t proc_cpuset_location = hwloc_bitmap_alloc();

  hwloc_get_last_cpu_location(s_hwloc_topology, proc_cpuset_location,
                              HWLOC_CPUBIND_THREAD);

  const unsigned max_root =
      hwloc_get_nbobjs_by_type(s_hwloc_topology, root_type);

  unsigned root_base     = max_root;
  unsigned root_count    = 0;
  unsigned core_per_root = 0;
  unsigned pu_per_core   = 0;
  bool symmetric         = true;

  for (unsigned i = 0; i < max_root; ++i) {
    const hwloc_obj_t root =
        hwloc_get_obj_by_type(s_hwloc_topology, root_type, i);

    if (hwloc_bitmap_intersects(s_process_binding, root->cpuset)) {
      ++root_count;

      // Remember which root (NUMA) object the master thread is running on.
      // This will be logical NUMA rank #0 for this process.

      if (hwloc_bitmap_intersects(proc_cpuset_location, root->cpuset)) {
        root_base = i;
      }

      // Count available cores:

      const unsigned max_core = hwloc_get_nbobjs_inside_cpuset_by_type(
          s_hwloc_topology, root->cpuset, HWLOC_OBJ_CORE);

      unsigned core_count = 0;

      for (unsigned j = 0; j < max_core; ++j) {
        const hwloc_obj_t core = hwloc_get_obj_inside_cpuset_by_type(
            s_hwloc_topology, root->cpuset, HWLOC_OBJ_CORE, j);

        // If process' cpuset intersects core's cpuset then process can access
        // this core. Must use intersection instead of inclusion because the
        // Intel-Phi MPI may bind the process to only one of the core's
        // hyperthreads.
        //
        // Assumption: if the process can access any hyperthread of the core
        // then it has ownership of the entire core.
        // This assumes that it would be performance-detrimental
        // to spawn more than one MPI process per core and use nested threading.

        if (hwloc_bitmap_intersects(s_process_binding, core->cpuset)) {
          ++core_count;

          const unsigned pu_count = hwloc_get_nbobjs_inside_cpuset_by_type(
              s_hwloc_topology, core->cpuset, HWLOC_OBJ_PU);

          if (pu_per_core == 0) pu_per_core = pu_count;

          // Enforce symmetry by taking the minimum:

          pu_per_core = std::min(pu_per_core, pu_count);

          if (pu_count != pu_per_core) symmetric = false;
        }
      }

      if (0 == core_per_root) core_per_root = core_count;

      // Enforce symmetry by taking the minimum:

      core_per_root = std::min(core_per_root, core_count);

      if (core_count != core_per_root) symmetric = false;
    }
  }

  s_core_topology.first  = root_count;
  s_core_topology.second = core_per_root;
  s_core_capacity        = pu_per_core;

  // Fill the 's_core' array for fast mapping from a core coordinate to the
  // hwloc cpuset object required for thread location querying and binding.

  for (unsigned i = 0; i < max_root; ++i) {
    const unsigned root_rank = (i + root_base) % max_root;

    const hwloc_obj_t root =
        hwloc_get_obj_by_type(s_hwloc_topology, root_type, root_rank);

    if (hwloc_bitmap_intersects(s_process_binding, root->cpuset)) {
      const unsigned max_core = hwloc_get_nbobjs_inside_cpuset_by_type(
          s_hwloc_topology, root->cpuset, HWLOC_OBJ_CORE);

      unsigned core_count = 0;

      for (unsigned j = 0; j < max_core && core_count < core_per_root; ++j) {
        const hwloc_obj_t core = hwloc_get_obj_inside_cpuset_by_type(
            s_hwloc_topology, root->cpuset, HWLOC_OBJ_CORE, j);

        if (hwloc_bitmap_intersects(s_process_binding, core->cpuset)) {
          s_core[core_count + core_per_root * i] = core->cpuset;

          ++core_count;
        }
      }
    }
  }

  hwloc_bitmap_free(proc_cpuset_location);

  if (Kokkos::show_warnings() && !symmetric) {
    std::cerr << "Kokkos::hwloc WARNING: Using a symmetric subset of a "
                 "non-symmetric core topology."
              << std::endl;
  }
}

}  // namespace

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

bool available() { return true; }

unsigned get_available_numa_count() {
  sentinel();
  return s_core_topology.first;
}

unsigned get_available_cores_per_numa() {
  sentinel();
  return s_core_topology.second;
}

unsigned get_available_threads_per_core() {
  sentinel();
  return s_core_capacity;
}

bool can_bind_threads() {
  sentinel();
  return s_can_bind_threads;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

unsigned bind_this_thread(const unsigned coordinate_count,
                          std::pair<unsigned, unsigned> coordinate[]) {
  unsigned i = 0;

  try {
    const std::pair<unsigned, unsigned> current = get_this_thread_coordinate();

    // Match one of the requests:
    for (i = 0; i < coordinate_count && current != coordinate[i]; ++i)
      ;

    if (coordinate_count == i) {
      // Match the first request (typically NUMA):
      for (i = 0; i < coordinate_count && current.first != coordinate[i].first;
           ++i)
        ;
    }

    if (coordinate_count == i) {
      // Match any unclaimed request:
      for (i = 0; i < coordinate_count && ~0u == coordinate[i].first; ++i)
        ;
    }

    if (coordinate_count == i || !bind_this_thread(coordinate[i])) {
      // Failed to bind:
      i = ~0u;
    }

    if (i < coordinate_count) {
#if DEBUG_PRINT
      if (current != coordinate[i]) {
        std::cout << "  bind_this_thread: rebinding from (" << current.first
                  << "," << current.second << ") to (" << coordinate[i].first
                  << "," << coordinate[i].second << ")" << std::endl;
      }
#endif

      coordinate[i].first  = ~0u;
      coordinate[i].second = ~0u;
    }
  } catch (...) {
    i = ~0u;
  }

  return i;
}

bool bind_this_thread(const std::pair<unsigned, unsigned> coord) {
  if (!sentinel()) return false;

#if DEBUG_PRINT

  std::cout << "Kokkos::bind_this_thread() at ";

  hwloc_get_last_cpu_location(s_hwloc_topology, s_hwloc_location,
                              HWLOC_CPUBIND_THREAD);

  print_bitmap(std::cout, s_hwloc_location);

  std::cout << " to ";

  print_bitmap(std::cout,
               s_core[coord.second + coord.first * s_core_topology.second]);

  std::cout << std::endl;

#endif

  // As safe and fast as possible.
  // Fast-lookup by caching the coordinate -> hwloc cpuset mapping in 's_core'.
  return coord.first < s_core_topology.first &&
         coord.second < s_core_topology.second &&
         0 == hwloc_set_cpubind(
                  s_hwloc_topology,
                  s_core[coord.second + coord.first * s_core_topology.second],
                  HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT);
}

bool unbind_this_thread() {
  if (!sentinel()) return false;

#define HWLOC_DEBUG_PRINT 0

#if HWLOC_DEBUG_PRINT

  std::cout << "Kokkos::unbind_this_thread() from ";

  hwloc_get_cpubind(s_hwloc_topology, s_hwloc_location, HWLOC_CPUBIND_THREAD);

  print_bitmap(std::cout, s_hwloc_location);

#endif

  const bool result =
      s_hwloc_topology &&
      0 == hwloc_set_cpubind(s_hwloc_topology, s_process_binding,
                             HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT);

#if HWLOC_DEBUG_PRINT

  std::cout << " to ";

  hwloc_get_cpubind(s_hwloc_topology, s_hwloc_location, HWLOC_CPUBIND_THREAD);

  print_bitmap(std::cout, s_hwloc_location);

  std::cout << std::endl;

#endif

  return result;

#undef HWLOC_DEBUG_PRINT
}

//----------------------------------------------------------------------------

std::pair<unsigned, unsigned> get_this_thread_coordinate() {
  std::pair<unsigned, unsigned> coord(0u, 0u);

  if (!sentinel()) return coord;

  const unsigned n = s_core_topology.first * s_core_topology.second;

  // Using the pre-allocated 's_hwloc_location' to avoid memory
  // allocation by this thread.  This call is NOT thread-safe.
  hwloc_get_last_cpu_location(s_hwloc_topology, s_hwloc_location,
                              HWLOC_CPUBIND_THREAD);

  unsigned i = 0;

  while (i < n && !hwloc_bitmap_intersects(s_hwloc_location, s_core[i])) ++i;

  if (i < n) {
    coord.first  = i / s_core_topology.second;
    coord.second = i % s_core_topology.second;
  }

  return coord;
}

//----------------------------------------------------------------------------

} /* namespace hwloc */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#else /* ! defined( KOKKOS_ENABLE_HWLOC ) */

namespace Kokkos {
namespace hwloc {

bool available() { return false; }
bool can_bind_threads() { return false; }

unsigned get_available_numa_count() { return 1; }
unsigned get_available_cores_per_numa() { return 1; }
unsigned get_available_threads_per_core() { return 1; }

unsigned bind_this_thread(const unsigned, std::pair<unsigned, unsigned>[]) {
  return ~0;
}

bool bind_this_thread(const std::pair<unsigned, unsigned>) { return false; }

bool unbind_this_thread() { return true; }

std::pair<unsigned, unsigned> get_this_thread_coordinate() {
  return std::pair<unsigned, unsigned>(0, 0);
}

}  // namespace hwloc
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif
