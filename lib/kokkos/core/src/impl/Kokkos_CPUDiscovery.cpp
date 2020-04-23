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

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#elif !defined(__APPLE__)
#include <unistd.h>
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>

namespace Kokkos {
namespace Impl {

int processors_per_node() {
#ifdef _SC_NPROCESSORS_ONLN
  int const num_procs     = sysconf(_SC_NPROCESSORS_ONLN);
  int const num_procs_max = sysconf(_SC_NPROCESSORS_CONF);
  if ((num_procs < 1) || (num_procs_max < 1)) {
    return -1;
  }
  return num_procs;
#else
  return -1;
#endif
}

int mpi_ranks_per_node() {
  char *str;
  int ppn = 1;
  // if ((str = getenv("SLURM_TASKS_PER_NODE"))) {
  //  ppn = atoi(str);
  //  if(ppn<=0) ppn = 1;
  //}
  if ((str = getenv("MV2_COMM_WORLD_LOCAL_SIZE"))) {
    ppn = atoi(str);
    if (ppn <= 0) ppn = 1;
  }
  if ((str = getenv("OMPI_COMM_WORLD_LOCAL_SIZE"))) {
    ppn = atoi(str);
    if (ppn <= 0) ppn = 1;
  }
  return ppn;
}

int mpi_local_rank_on_node() {
  char *str;
  int local_rank = 0;
  // if ((str = getenv("SLURM_LOCALID"))) {
  //  local_rank = atoi(str);
  //}
  if ((str = getenv("MV2_COMM_WORLD_LOCAL_RANK"))) {
    local_rank = atoi(str);
  }
  if ((str = getenv("OMPI_COMM_WORLD_LOCAL_RANK"))) {
    local_rank = atoi(str);
  }
  return local_rank;
}

}  // namespace Impl
}  // namespace Kokkos
