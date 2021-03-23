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

/*--------------------------------------------------------------------------*/

#ifndef KOKKOS_HIP_INSTANCE_HPP
#define KOKKOS_HIP_INSTANCE_HPP

#include <Kokkos_HIP_Space.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct HIPTraits {
  static int constexpr WarpSize       = 64;
  static int constexpr WarpIndexMask  = 0x003f; /* hexadecimal for 63 */
  static int constexpr WarpIndexShift = 6;      /* WarpSize == 1 << WarpShift*/
  static int constexpr MaxThreadsPerBlock =
      1024;  // FIXME_HIP -- assumed constant for now

  static int constexpr ConstantMemoryUsage        = 0x008000; /* 32k bytes */
  static int constexpr ConstantMemoryUseThreshold = 0x000200; /* 512 bytes */
};

//----------------------------------------------------------------------------

HIP::size_type hip_internal_maximum_warp_count();
HIP::size_type hip_internal_maximum_grid_count();
HIP::size_type hip_internal_multiprocessor_count();

HIP::size_type *hip_internal_scratch_space(const HIP::size_type size);
HIP::size_type *hip_internal_scratch_flags(const HIP::size_type size);

//----------------------------------------------------------------------------

class HIPInternal {
 private:
  HIPInternal(const HIPInternal &);
  HIPInternal &operator=(const HIPInternal &);

 public:
  using size_type = ::Kokkos::Experimental::HIP::size_type;

  int m_hipDev;
  int m_hipArch;
  unsigned m_multiProcCount;
  unsigned m_maxWarpCount;
  unsigned m_maxBlock;
  unsigned m_maxBlocksPerSM;
  unsigned m_maxSharedWords;
  int m_regsPerSM;
  int m_shmemPerSM;
  int m_maxShmemPerBlock;
  int m_maxThreadsPerSM;

  // Scratch Spaces for Reductions
  size_type m_scratchSpaceCount;
  size_type m_scratchFlagsCount;

  size_type *m_scratchSpace;
  size_type *m_scratchFlags;
  uint32_t *m_scratchConcurrentBitset = nullptr;

  hipDeviceProp_t m_deviceProp;

  hipStream_t m_stream;

  // Team Scratch Level 1 Space
  mutable int64_t m_team_scratch_current_size;
  mutable void *m_team_scratch_ptr;

  bool was_finalized = false;

  static HIPInternal &singleton();

  int verify_is_initialized(const char *const label) const;

  int is_initialized() const {
    return m_hipDev >= 0;
  }  // 0 != m_scratchSpace && 0 != m_scratchFlags ; }

  void initialize(int hip_device_id, hipStream_t stream = nullptr);
  void finalize();

  void print_configuration(std::ostream &) const;

  void fence() const;

  ~HIPInternal();

  HIPInternal()
      : m_hipDev(-1),
        m_hipArch(-1),
        m_multiProcCount(0),
        m_maxWarpCount(0),
        m_maxBlock(0),
        m_maxSharedWords(0),
        m_shmemPerSM(0),
        m_maxShmemPerBlock(0),
        m_maxThreadsPerSM(0),
        m_scratchSpaceCount(0),
        m_scratchFlagsCount(0),
        m_scratchSpace(nullptr),
        m_scratchFlags(nullptr),
        m_stream(nullptr),
        m_team_scratch_current_size(0),
        m_team_scratch_ptr(nullptr) {}

  // Resizing of reduction related scratch spaces
  size_type *scratch_space(const size_type size);
  size_type *scratch_flags(const size_type size);

  // Resizing of team level 1 scratch
  void *resize_team_scratch_space(std::int64_t bytes,
                                  bool force_shrink = false);
};

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
