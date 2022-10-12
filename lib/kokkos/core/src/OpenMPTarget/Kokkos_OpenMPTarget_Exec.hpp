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

#ifndef KOKKOS_OPENMPTARGETEXEC_HPP
#define KOKKOS_OPENMPTARGETEXEC_HPP

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Spinwait.hpp>

#include <Kokkos_Atomic.hpp>
#include "Kokkos_OpenMPTarget_Abort.hpp"

// FIXME_OPENMPTARGET - Using this macro to implement a workaround for
// hierarchical reducers. It avoids hitting the code path which we wanted to
// write but doesn't work. undef'ed at the end.
// Intel compilers prefer the non-workaround version.
#ifndef KOKKOS_ARCH_INTEL_GPU
#define KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND
#endif

// FIXME_OPENMPTARGET - Using this macro to implement a workaround for
// hierarchical scan. It avoids hitting the code path which we wanted to
// write but doesn't work. undef'ed at the end.
#ifndef KOKKOS_ARCH_INTEL_GPU
#define KOKKOS_IMPL_TEAM_SCAN_WORKAROUND
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class Reducer>
struct OpenMPTargetReducerWrapper {
  using value_type = typename Reducer::value_type;

  // Using a generic unknown Reducer for the OpenMPTarget backend is not
  // implemented.
  KOKKOS_INLINE_FUNCTION
  static void join(value_type&, const value_type&) = delete;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type&, const volatile value_type&) = delete;

  KOKKOS_INLINE_FUNCTION
  static void init(value_type&) = delete;
};

template <class Scalar, class Space>
struct OpenMPTargetReducerWrapper<Sum<Scalar, Space>> {
 public:
  // Required
  using value_type = std::remove_cv_t<Scalar>;

  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) { dest += src; }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    dest += src;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val = reduction_identity<value_type>::sum();
  }
};

template <class Scalar, class Space>
struct OpenMPTargetReducerWrapper<Prod<Scalar, Space>> {
 public:
  // Required
  using value_type = std::remove_cv_t<Scalar>;

  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) { dest *= src; }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    dest *= src;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val = reduction_identity<value_type>::prod();
  }
};

template <class Scalar, class Space>
struct OpenMPTargetReducerWrapper<Min<Scalar, Space>> {
 public:
  // Required
  using value_type = std::remove_cv_t<Scalar>;

  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    if (src < dest) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    if (src < dest) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val = reduction_identity<value_type>::min();
  }
};

template <class Scalar, class Space>
struct OpenMPTargetReducerWrapper<Max<Scalar, Space>> {
 public:
  // Required
  using value_type = std::remove_cv_t<Scalar>;

  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    if (src > dest) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    if (src > dest) dest = src;
  }

  // Required
  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val = reduction_identity<value_type>::max();
  }
};

template <class Scalar, class Space>
struct OpenMPTargetReducerWrapper<LAnd<Scalar, Space>> {
 public:
  // Required
  using value_type = std::remove_cv_t<Scalar>;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    dest = dest && src;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    dest = dest && src;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val = reduction_identity<value_type>::land();
  }
};

template <class Scalar, class Space>
struct OpenMPTargetReducerWrapper<LOr<Scalar, Space>> {
 public:
  // Required
  using value_type = std::remove_cv_t<Scalar>;

  using result_view_type = Kokkos::View<value_type, Space>;

  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    dest = dest || src;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    dest = dest || src;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val = reduction_identity<value_type>::lor();
  }
};

template <class Scalar, class Space>
struct OpenMPTargetReducerWrapper<BAnd<Scalar, Space>> {
 public:
  // Required
  using value_type = std::remove_cv_t<Scalar>;

  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    dest = dest & src;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    dest = dest & src;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val = reduction_identity<value_type>::band();
  }
};

template <class Scalar, class Space>
struct OpenMPTargetReducerWrapper<BOr<Scalar, Space>> {
 public:
  // Required
  using value_type = std::remove_cv_t<Scalar>;

  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    dest = dest | src;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    dest = dest | src;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val = reduction_identity<value_type>::bor();
  }
};

template <class Scalar, class Index, class Space>
struct OpenMPTargetReducerWrapper<MinLoc<Scalar, Index, Space>> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;

 public:
  // Required
  using value_type = ValLocScalar<scalar_type, index_type>;

  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    if (src.val < dest.val) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    if (src.val < dest.val) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.val = reduction_identity<scalar_type>::min();
    val.loc = reduction_identity<index_type>::min();
  }
};

template <class Scalar, class Index, class Space>
struct OpenMPTargetReducerWrapper<MaxLoc<Scalar, Index, Space>> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;

 public:
  // Required
  using value_type = ValLocScalar<scalar_type, index_type>;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    if (src.val > dest.val) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    if (src.val > dest.val) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.val = reduction_identity<scalar_type>::max();
    val.loc = reduction_identity<index_type>::min();
  }
};

template <class Scalar, class Space>
struct OpenMPTargetReducerWrapper<MinMax<Scalar, Space>> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;

 public:
  // Required
  using value_type = MinMaxScalar<scalar_type>;

  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
    }
    if (src.max_val > dest.max_val) {
      dest.max_val = src.max_val;
    }
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
    }
    if (src.max_val > dest.max_val) {
      dest.max_val = src.max_val;
    }
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.max_val = reduction_identity<scalar_type>::max();
    val.min_val = reduction_identity<scalar_type>::min();
  }
};

template <class Scalar, class Index, class Space>
struct OpenMPTargetReducerWrapper<MinMaxLoc<Scalar, Index, Space>> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;

 public:
  // Required
  using value_type = MinMaxLocScalar<scalar_type, index_type>;

  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    }
    if (src.max_val > dest.max_val) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    }
    if (src.max_val > dest.max_val) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.max_val = reduction_identity<scalar_type>::max();
    val.min_val = reduction_identity<scalar_type>::min();
    val.max_loc = reduction_identity<index_type>::min();
    val.min_loc = reduction_identity<index_type>::min();
  }
};

//
// specialize for MaxFirstLoc
//
template <class Scalar, class Index, class Space>
struct OpenMPTargetReducerWrapper<MaxFirstLoc<Scalar, Index, Space>> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;

 public:
  // Required
  using value_type = ValLocScalar<scalar_type, index_type>;

// WORKAROUND OPENMPTARGET
// This pragma omp declare target should not be necessary, but Intel compiler
// fails without it
#pragma omp declare target
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    if (dest.val < src.val) {
      dest = src;
    } else if (!(src.val < dest.val)) {
      dest.loc = (src.loc < dest.loc) ? src.loc : dest.loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    if (dest.val < src.val) {
      dest = src;
    } else if (!(src.val < dest.val)) {
      dest.loc = (src.loc < dest.loc) ? src.loc : dest.loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.val = reduction_identity<scalar_type>::max();
    val.loc = reduction_identity<index_type>::min();
  }
#pragma omp end declare target
};

//
// specialize for MinFirstLoc
//
template <class Scalar, class Index, class Space>
struct OpenMPTargetReducerWrapper<MinFirstLoc<Scalar, Index, Space>> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;

 public:
  // Required
  using value_type = ValLocScalar<scalar_type, index_type>;

// WORKAROUND OPENMPTARGET
// This pragma omp declare target should not be necessary, but Intel compiler
// fails without it
#pragma omp declare target
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    if (src.val < dest.val) {
      dest = src;
    } else if (!(dest.val < src.val)) {
      dest.loc = (src.loc < dest.loc) ? src.loc : dest.loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    if (src.val < dest.val) {
      dest = src;
    } else if (!(dest.val < src.val)) {
      dest.loc = (src.loc < dest.loc) ? src.loc : dest.loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.val = reduction_identity<scalar_type>::min();
    val.loc = reduction_identity<index_type>::min();
  }
#pragma omp end declare target
};

//
// specialize for MinMaxFirstLastLoc
//
template <class Scalar, class Index, class Space>
struct OpenMPTargetReducerWrapper<MinMaxFirstLastLoc<Scalar, Index, Space>> {
 private:
  using scalar_type = std::remove_cv_t<Scalar>;
  using index_type  = std::remove_cv_t<Index>;

 public:
  // Required
  using value_type = MinMaxLocScalar<scalar_type, index_type>;

// WORKAROUND OPENMPTARGET
// This pragma omp declare target should not be necessary, but Intel compiler
// fails without it
#pragma omp declare target
  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    } else if (!(dest.min_val < src.min_val)) {
      dest.min_loc = (src.min_loc < dest.min_loc) ? src.min_loc : dest.min_loc;
    }

    if (dest.max_val < src.max_val) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    } else if (!(src.max_val < dest.max_val)) {
      dest.max_loc = (src.max_loc > dest.max_loc) ? src.max_loc : dest.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    } else if (!(dest.min_val < src.min_val)) {
      dest.min_loc = (src.min_loc < dest.min_loc) ? src.min_loc : dest.min_loc;
    }

    if (dest.max_val < src.max_val) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    } else if (!(src.max_val < dest.max_val)) {
      dest.max_loc = (src.max_loc > dest.max_loc) ? src.max_loc : dest.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.max_val = reduction_identity<scalar_type>::max();
    val.min_val = reduction_identity<scalar_type>::min();
    val.max_loc = reduction_identity<index_type>::max();
    val.min_loc = reduction_identity<index_type>::min();
  }
#pragma omp end declare target
};

//
// specialize for FirstLoc
//
template <class Index, class Space>
struct OpenMPTargetReducerWrapper<FirstLoc<Index, Space>> {
 private:
  using index_type = std::remove_cv_t<Index>;

 public:
  // Required
  using value_type = FirstLocScalar<index_type>;

// WORKAROUND OPENMPTARGET
// This pragma omp declare target should not be necessary, but Intel compiler
// fails without it
#pragma omp declare target
  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    dest.min_loc_true = (src.min_loc_true < dest.min_loc_true)
                            ? src.min_loc_true
                            : dest.min_loc_true;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    dest.min_loc_true = (src.min_loc_true < dest.min_loc_true)
                            ? src.min_loc_true
                            : dest.min_loc_true;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.min_loc_true = reduction_identity<index_type>::min();
  }
#pragma omp end declare target
};

//
// specialize for LastLoc
//
template <class Index, class Space>
struct OpenMPTargetReducerWrapper<LastLoc<Index, Space>> {
 private:
  using index_type = std::remove_cv_t<Index>;

 public:
  // Required
  using value_type = LastLocScalar<index_type>;

// WORKAROUND OPENMPTARGET
// This pragma omp declare target should not be necessary, but Intel compiler
// fails without it
#pragma omp declare target
  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    dest.max_loc_true = (src.max_loc_true > dest.max_loc_true)
                            ? src.max_loc_true
                            : dest.max_loc_true;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    dest.max_loc_true = (src.max_loc_true > dest.max_loc_true)
                            ? src.max_loc_true
                            : dest.max_loc_true;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.max_loc_true = reduction_identity<index_type>::max();
  }
#pragma omp end declare target
};

//
// specialize for StdIsPartitioned
//
template <class Index, class Space>
struct OpenMPTargetReducerWrapper<StdIsPartitioned<Index, Space>> {
 private:
  using index_type = std::remove_cv_t<Index>;

 public:
  // Required
  using value_type = StdIsPartScalar<index_type>;

// WORKAROUND OPENMPTARGET
// This pragma omp declare target should not be necessary, but Intel compiler
// fails without it
#pragma omp declare target
  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    dest.max_loc_true = (dest.max_loc_true < src.max_loc_true)
                            ? src.max_loc_true
                            : dest.max_loc_true;

    dest.min_loc_false = (dest.min_loc_false < src.min_loc_false)
                             ? dest.min_loc_false
                             : src.min_loc_false;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    dest.max_loc_true = (dest.max_loc_true < src.max_loc_true)
                            ? src.max_loc_true
                            : dest.max_loc_true;

    dest.min_loc_false = (dest.min_loc_false < src.min_loc_false)
                             ? dest.min_loc_false
                             : src.min_loc_false;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.max_loc_true  = ::Kokkos::reduction_identity<index_type>::max();
    val.min_loc_false = ::Kokkos::reduction_identity<index_type>::min();
  }
#pragma omp end declare target
};

//
// specialize for StdPartitionPoint
//
template <class Index, class Space>
struct OpenMPTargetReducerWrapper<StdPartitionPoint<Index, Space>> {
 private:
  using index_type = std::remove_cv_t<Index>;

 public:
  // Required
  using value_type = StdPartPointScalar<index_type>;

// WORKAROUND OPENMPTARGET
// This pragma omp declare target should not be necessary, but Intel compiler
// fails without it
#pragma omp declare target
  // Required
  KOKKOS_INLINE_FUNCTION
  static void join(value_type& dest, const value_type& src) {
    dest.min_loc_false = (dest.min_loc_false < src.min_loc_false)
                             ? dest.min_loc_false
                             : src.min_loc_false;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& dest, const volatile value_type& src) {
    dest.min_loc_false = (dest.min_loc_false < src.min_loc_false)
                             ? dest.min_loc_false
                             : src.min_loc_false;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& val) {
    val.min_loc_false = ::Kokkos::reduction_identity<index_type>::min();
  }
#pragma omp end declare target
};

/*
template<class ReducerType>
class OpenMPTargetReducerWrapper {
  public:
    const ReducerType& reducer;
    using value_type = typename ReducerType::value_type;
    value_type& value;

    KOKKOS_INLINE_FUNCTION
    void join(const value_type& upd) {
      reducer.join(value,upd);
    }

    KOKKOS_INLINE_FUNCTION
    void init(const value_type& upd) {
      reducer.init(value,upd);
    }
};*/

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  Data for OpenMPTarget thread execution */

class OpenMPTargetExec {
 public:
  // FIXME_OPENMPTARGET - Currently the maximum number of
  // teams possible is calculated based on NVIDIA's Volta GPU. In
  // future this value should be based on the chosen architecture for the
  // OpenMPTarget backend.
  enum { MAX_ACTIVE_THREADS = 2080 * 80 };
  enum { MAX_ACTIVE_TEAMS = MAX_ACTIVE_THREADS / 32 };

 private:
  static void* scratch_ptr;

 public:
  static void verify_is_process(const char* const);
  static void verify_initialized(const char* const);

  static int* get_lock_array(int num_teams);
  static void* get_scratch_ptr();
  static void clear_scratch();
  static void clear_lock_array();
  static void resize_scratch(int64_t team_reduce_bytes,
                             int64_t team_shared_bytes,
                             int64_t thread_local_bytes, int64_t league_size);

  static void* m_scratch_ptr;
  static int64_t m_scratch_size;
  static int* m_lock_array;
  static int64_t m_lock_size;
  static uint32_t* m_uniquetoken_ptr;
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

class OpenMPTargetExecTeamMember {
 public:
  enum { TEAM_REDUCE_SIZE = 512 };

  /** \brief  Thread states for team synchronization */
  enum { Active = 0, Rendezvous = 1 };

  using execution_space      = Kokkos::Experimental::OpenMPTarget;
  using scratch_memory_space = execution_space::scratch_memory_space;

  scratch_memory_space m_team_shared;
  size_t m_team_scratch_size[2];
  int m_team_rank;
  int m_team_size;
  int m_league_rank;
  int m_league_size;
  int m_vector_length;
  int m_vector_lane;
  int m_shmem_block_index;
  void* m_glb_scratch;
  void* m_reduce_scratch;

 public:
  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_shmem() const {
    return m_team_shared.set_team_thread_mode(0, 1, 0);
  }

  // set_team_thread_mode routine parameters for future understanding:
  // first parameter - scratch level.
  // second parameter - size multiplier for advancing scratch ptr after a
  // request was serviced. third parameter - offset size multiplier from current
  // scratch ptr when returning a ptr for a request.
  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_scratch(int level) const {
    return m_team_shared.set_team_thread_mode(level, 1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& thread_scratch(int level) const {
    return m_team_shared.set_team_thread_mode(level, team_size(), team_rank());
  }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return m_team_rank; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return m_team_size; }
  KOKKOS_INLINE_FUNCTION void* impl_reduce_scratch() const {
    return m_reduce_scratch;
  }

  KOKKOS_INLINE_FUNCTION void team_barrier() const {
#pragma omp barrier
  }

  template <class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(ValueType& value,
                                             int thread_id) const {
    // Make sure there is enough scratch space:
    using type = std::conditional_t<(sizeof(ValueType) < TEAM_REDUCE_SIZE),
                                    ValueType, void>;
    type* team_scratch =
        reinterpret_cast<type*>(static_cast<char*>(m_glb_scratch) +
                                TEAM_REDUCE_SIZE * omp_get_team_num());
#pragma omp barrier
    if (team_rank() == thread_id) *team_scratch = value;
#pragma omp barrier
    value = *team_scratch;
  }

  template <class Closure, class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(const Closure& f, ValueType& value,
                                             const int& thread_id) const {
    f(value);
    team_broadcast(value, thread_id);
  }

  // FIXME_OPENMPTARGET this function has the wrong interface and currently
  // ignores the reducer passed.
  template <class ValueType, class JoinOp>
  KOKKOS_INLINE_FUNCTION ValueType team_reduce(const ValueType& value,
                                               const JoinOp&) const {
#pragma omp barrier

    using value_type = ValueType;
    //    const JoinLambdaAdapter<value_type, JoinOp> op(op_in);

    // Make sure there is enough scratch space:
    using type = std::conditional_t<(sizeof(value_type) < TEAM_REDUCE_SIZE),
                                    value_type, void>;

    const int n_values = TEAM_REDUCE_SIZE / sizeof(value_type);
    type* team_scratch =
        reinterpret_cast<type*>(static_cast<char*>(m_glb_scratch) +
                                TEAM_REDUCE_SIZE * omp_get_team_num());
    for (int i = m_team_rank; i < n_values; i += m_team_size) {
      team_scratch[i] = value_type();
    }

#pragma omp barrier

    for (int k = 0; k < m_team_size; k += n_values) {
      if ((k <= m_team_rank) && (k + n_values > m_team_rank))
        team_scratch[m_team_rank % n_values] += value;
#pragma omp barrier
    }

    for (int d = 1; d < n_values; d *= 2) {
      if ((m_team_rank + d < n_values) && (m_team_rank % (2 * d) == 0)) {
        team_scratch[m_team_rank] += team_scratch[m_team_rank + d];
      }
#pragma omp barrier
    }
    return team_scratch[0];
  }
  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template <typename ArgType>
  KOKKOS_INLINE_FUNCTION ArgType
  team_scan(const ArgType& /*value*/, ArgType* const /*global_accum*/) const {
    // FIXME_OPENMPTARGET
    /*  // Make sure there is enough scratch space:
      using type =
        std::conditional_t<(sizeof(ArgType) < TEAM_REDUCE_SIZE), ArgType, void>;

      volatile type * const work_value  = ((type*) m_exec.scratch_thread());

      *work_value = value ;

      memory_fence();

      if ( team_fan_in() ) {
        // The last thread to synchronize returns true, all other threads wait
      for team_fan_out()
        // m_team_base[0]                 == highest ranking team member
        // m_team_base[ m_team_size - 1 ] == lowest ranking team member
        //
        // 1) copy from lower to higher rank, initialize lowest rank to zero
        // 2) prefix sum from lowest to highest rank, skipping lowest rank

        type accum = 0 ;

        if ( global_accum ) {
          for ( int i = m_team_size ; i-- ; ) {
            type & val = *((type*) m_exec.pool_rev( m_team_base_rev + i
      )->scratch_thread()); accum += val ;
          }
          accum = atomic_fetch_add( global_accum , accum );
        }

        for ( int i = m_team_size ; i-- ; ) {
          type & val = *((type*) m_exec.pool_rev( m_team_base_rev + i
      )->scratch_thread()); const type offset = accum ; accum += val ; val =
      offset ;
        }

        memory_fence();
      }

      team_fan_out();

      return *work_value ;*/
    return ArgType();
  }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template <typename Type>
  KOKKOS_INLINE_FUNCTION Type team_scan(const Type& value) const {
    return this->template team_scan<Type>(value, 0);
  }

  //----------------------------------------
  // Private for the driver

 private:
  using space = execution_space::scratch_memory_space;

 public:
  // FIXME_OPENMPTARGET - 512(16*32) bytes at the begining of the scratch space
  // for each league is saved for reduction. It should actually be based on the
  // ValueType of the reduction variable.
  inline OpenMPTargetExecTeamMember(
      const int league_rank, const int league_size, const int team_size,
      const int vector_length  // const TeamPolicyInternal< OpenMPTarget,
                               // Properties ...> & team
      ,
      void* const glb_scratch, const int shmem_block_index,
      const size_t shmem_size_L0, const size_t shmem_size_L1)
      : m_team_scratch_size{shmem_size_L0, shmem_size_L1},
        m_team_rank(0),
        m_team_size(team_size),
        m_league_rank(league_rank),
        m_league_size(league_size),
        m_vector_length(vector_length),
        m_shmem_block_index(shmem_block_index),
        m_glb_scratch(glb_scratch) {
    const int omp_tid = omp_get_thread_num();

    // The scratch memory allocated is a sum of TEAM_REDUCE_SIZE, L0 shmem size
    // and L1 shmem size. TEAM_REDUCE_SIZE = 512 bytes saved per team for
    // hierarchical reduction. There is an additional 10% of the requested
    // scratch memory allocated per team as padding. Hence the product with 0.1.
    const int reduce_offset =
        m_shmem_block_index *
        (shmem_size_L0 + shmem_size_L1 +
         ((shmem_size_L0 + shmem_size_L1) * 0.1) + TEAM_REDUCE_SIZE);
    const int l0_offset = reduce_offset + TEAM_REDUCE_SIZE;
    const int l1_offset = l0_offset + shmem_size_L0;
    m_team_shared       = scratch_memory_space(
        (static_cast<char*>(glb_scratch) + l0_offset), shmem_size_L0,
        static_cast<char*>(glb_scratch) + l1_offset, shmem_size_L1);
    m_reduce_scratch = static_cast<char*>(glb_scratch) + reduce_offset;
    m_league_rank    = league_rank;
    m_team_rank      = omp_tid;
    m_vector_lane    = 0;
  }

  static inline int team_reduce_size() { return TEAM_REDUCE_SIZE; }
};

template <class... Properties>
class TeamPolicyInternal<Kokkos::Experimental::OpenMPTarget, Properties...>
    : public PolicyTraits<Properties...> {
 public:
  //! Tag this class as a kokkos execution policy
  using execution_policy = TeamPolicyInternal;

  using traits = PolicyTraits<Properties...>;

  //----------------------------------------

  template <class FunctorType>
  inline static int team_size_max(const FunctorType&, const ParallelForTag&) {
    return 256;
  }

  template <class FunctorType>
  inline static int team_size_max(const FunctorType&,
                                  const ParallelReduceTag&) {
    return 256;
  }

  template <class FunctorType, class ReducerType>
  inline static int team_size_max(const FunctorType&, const ReducerType&,
                                  const ParallelReduceTag&) {
    return 256;
  }

  template <class FunctorType>
  inline static int team_size_recommended(const FunctorType&,
                                          const ParallelForTag&) {
    return 128;
  }

  template <class FunctorType>
  inline static int team_size_recommended(const FunctorType&,
                                          const ParallelReduceTag&) {
    return 128;
  }

  template <class FunctorType, class ReducerType>
  inline static int team_size_recommended(const FunctorType&,
                                          const ReducerType&,
                                          const ParallelReduceTag&) {
    return 128;
  }

  //----------------------------------------

 private:
  int m_league_size;
  int m_team_size;
  int m_vector_length;
  int m_team_alloc;
  int m_team_iter;
  std::array<size_t, 2> m_team_scratch_size;
  std::array<size_t, 2> m_thread_scratch_size;
  bool m_tune_team_size;
  bool m_tune_vector_length;
  constexpr const static size_t default_team_size = 256;
  int m_chunk_size;

  inline void init(const int league_size_request, const int team_size_request,
                   const int vector_length_request) {
    m_league_size = league_size_request;

    // Minimum team size should be 32 for OpenMPTarget backend.
    if (team_size_request < 32) {
      Kokkos::Impl::OpenMPTarget_abort(
          "OpenMPTarget backend requires a minimum of 32 threads per team.\n");
    } else
      m_team_size = team_size_request;

    m_vector_length = vector_length_request;
    set_auto_chunk_size();
  }

  template <typename ExecSpace, typename... OtherProperties>
  friend class TeamPolicyInternal;

 public:
  // FIXME_OPENMPTARGET : Currently this routine is a copy of the Cuda
  // implementation, but this has to be tailored to be architecture specific.
  inline static int scratch_size_max(int level) {
    return (
        level == 0 ? 1024 * 40 :  // 48kB is the max for CUDA, but we need some
                                  // for team_member.reduce etc.
            20 * 1024 *
                1024);  // arbitrarily setting this to 20MB, for a Volta V100
                        // that would give us about 3.2GB for 2 teams per SM
  }
  inline bool impl_auto_team_size() const { return m_tune_team_size; }
  inline bool impl_auto_vector_length() const { return m_tune_vector_length; }
  inline void impl_set_team_size(const size_t size) { m_team_size = size; }
  inline void impl_set_vector_length(const size_t length) {
    m_tune_vector_length = length;
  }
  inline int impl_vector_length() const { return m_vector_length; }
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  KOKKOS_DEPRECATED inline int vector_length() const {
    return impl_vector_length();
  }
#endif
  inline int team_size() const { return m_team_size; }
  inline int league_size() const { return m_league_size; }
  inline size_t scratch_size(const int& level, int team_size_ = -1) const {
    if (team_size_ < 0) team_size_ = m_team_size;
    return m_team_scratch_size[level] +
           team_size_ * m_thread_scratch_size[level];
  }

  inline Kokkos::Experimental::OpenMPTarget space() const {
    return Kokkos::Experimental::OpenMPTarget();
  }

  template <class... OtherProperties>
  TeamPolicyInternal(const TeamPolicyInternal<OtherProperties...>& p)
      : m_league_size(p.m_league_size),
        m_team_size(p.m_team_size),
        m_vector_length(p.m_vector_length),
        m_team_alloc(p.m_team_alloc),
        m_team_iter(p.m_team_iter),
        m_team_scratch_size(p.m_team_scratch_size),
        m_thread_scratch_size(p.m_thread_scratch_size),
        m_tune_team_size(p.m_tune_team_size),
        m_tune_vector_length(p.m_tune_vector_length),
        m_chunk_size(p.m_chunk_size) {}

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request, int team_size_request,
                     int vector_length_request = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(false),
        m_tune_vector_length(false),
        m_chunk_size(0) {
    init(league_size_request, team_size_request, vector_length_request);
  }

  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     int vector_length_request = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(true),
        m_tune_vector_length(false),
        m_chunk_size(0) {
    init(league_size_request, default_team_size / vector_length_request,
         vector_length_request);
  }

  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(true),
        m_tune_vector_length(true),
        m_chunk_size(0) {
    init(league_size_request, default_team_size, 1);
  }
  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request, int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(false),
        m_tune_vector_length(true),
        m_chunk_size(0) {
    init(league_size_request, team_size_request, 1);
  }

  TeamPolicyInternal(int league_size_request, int team_size_request,
                     int vector_length_request = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(false),
        m_tune_vector_length(false),
        m_chunk_size(0) {
    init(league_size_request, team_size_request, vector_length_request);
  }

  TeamPolicyInternal(int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     int vector_length_request = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(true),
        m_tune_vector_length(false),
        m_chunk_size(0) {
    init(league_size_request, default_team_size / vector_length_request,
         vector_length_request);
  }

  TeamPolicyInternal(int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(true),
        m_tune_vector_length(true),
        m_chunk_size(0) {
    init(league_size_request, default_team_size, 1);
  }
  TeamPolicyInternal(int league_size_request, int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(false),
        m_tune_vector_length(true),
        m_chunk_size(0) {
    init(league_size_request, team_size_request, 1);
  }
  inline static size_t vector_length_max() {
    return 32; /* TODO: this is bad. Need logic that is compiler and backend
                  aware */
  }
  inline int team_alloc() const { return m_team_alloc; }
  inline int team_iter() const { return m_team_iter; }

  inline int chunk_size() const { return m_chunk_size; }

  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal& set_chunk_size(
      typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  /** \brief set per team scratch size for a specific level of the scratch
   * hierarchy */
  inline TeamPolicyInternal& set_scratch_size(const int& level,
                                              const PerTeamValue& per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  /** \brief set per thread scratch size for a specific level of the scratch
   * hierarchy */
  inline TeamPolicyInternal& set_scratch_size(
      const int& level, const PerThreadValue& per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  /** \brief set per thread and per team scratch size for a specific level of
   * the scratch hierarchy */
  inline TeamPolicyInternal& set_scratch_size(
      const int& level, const PerTeamValue& per_team,
      const PerThreadValue& per_thread) {
    m_team_scratch_size[level]   = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

 private:
  /** \brief finalize chunk_size if it was set to AUTO*/
  inline void set_auto_chunk_size() {
    int concurrency = 2048 * 128;

    if (concurrency == 0) concurrency = 1;

    if (m_chunk_size > 0) {
      if (!Impl::is_integral_power_of_two(m_chunk_size))
        Kokkos::abort("TeamPolicy blocking granularity must be power of two");
    }

    int new_chunk_size = 1;
    while (new_chunk_size * 100 * concurrency < m_league_size)
      new_chunk_size *= 2;
    if (new_chunk_size < 128) {
      new_chunk_size = 1;
      while ((new_chunk_size * 40 * concurrency < m_league_size) &&
             (new_chunk_size < 128))
        new_chunk_size *= 2;
    }
    m_chunk_size = new_chunk_size;
  }

 public:
  using member_type = Impl::OpenMPTargetExecTeamMember;
};
}  // namespace Impl

}  // namespace Kokkos

namespace Kokkos {

template <typename iType>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    iType, Impl::OpenMPTargetExecTeamMember>
TeamThreadRange(const Impl::OpenMPTargetExecTeamMember& thread,
                const iType& count) {
  return Impl::TeamThreadRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::OpenMPTargetExecTeamMember>
TeamThreadRange(const Impl::OpenMPTargetExecTeamMember& thread,
                const iType1& begin, const iType2& end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::TeamThreadRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, iType(begin),
                                               iType(end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION Impl::ThreadVectorRangeBoundariesStruct<
    iType, Impl::OpenMPTargetExecTeamMember>
ThreadVectorRange(const Impl::OpenMPTargetExecTeamMember& thread,
                  const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::ThreadVectorRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::OpenMPTargetExecTeamMember>
ThreadVectorRange(const Impl::OpenMPTargetExecTeamMember& thread,
                  const iType1& arg_begin, const iType2& arg_end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::ThreadVectorRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, iType(arg_begin),
                                               iType(arg_end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION Impl::TeamVectorRangeBoundariesStruct<
    iType, Impl::OpenMPTargetExecTeamMember>
TeamVectorRange(const Impl::OpenMPTargetExecTeamMember& thread,
                const iType& count) {
  return Impl::TeamVectorRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamVectorRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::OpenMPTargetExecTeamMember>
TeamVectorRange(const Impl::OpenMPTargetExecTeamMember& thread,
                const iType1& arg_begin, const iType2& arg_end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::TeamVectorRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, iType(arg_begin),
                                               iType(arg_end));
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember> PerTeam(
    const Impl::OpenMPTargetExecTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember> PerThread(
    const Impl::OpenMPTargetExecTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember>(thread);
}
}  // namespace Kokkos

namespace Kokkos {

/** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team.
 */
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda) {
#pragma omp for nowait schedule(static, 1)
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) lambda(i);
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team
 * and a summation of val is performed and put into result.
 */

template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<!Kokkos::is_reducer<ValueType>::value>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ValueType& result) {
  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  ValueType* TeamThread_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp barrier
  TeamThread_scratch[0] = ValueType();
#pragma omp barrier

  if constexpr (std::is_arithmetic<ValueType>::value) {
#pragma omp for reduction(+ : TeamThread_scratch[:1])
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      ValueType tmp = ValueType();
      lambda(i, tmp);
      TeamThread_scratch[0] += tmp;
    }
  } else {
#pragma omp declare reduction(custom:ValueType : omp_out += omp_in)

#pragma omp for reduction(custom : TeamThread_scratch[:1])
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      ValueType tmp = ValueType();
      lambda(i, tmp);
      TeamThread_scratch[0] += tmp;
    }
  }

  result = TeamThread_scratch[0];
}

#if !defined(KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND)
// For some reason the actual version we wanted to write doesn't work
// and crashes. We should try this with every new compiler
// This is the variant we actually wanted to write
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ReducerType result) {
  using ValueType = typename ReducerType::value_type;

#pragma omp declare reduction(                                               \
    custominner:ValueType                                                    \
    : Impl::OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(                                                             \
        Impl::OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  ValueType* TeamThread_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp barrier
  Impl::OpenMPTargetReducerWrapper<ReducerType>::init(TeamThread_scratch[0]);
#pragma omp barrier

#pragma omp for reduction(custominner : TeamThread_scratch[:1])
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i, TeamThread_scratch[0]);
  }
  result.reference() = TeamThread_scratch[0];
}
#else
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ReducerType result) {
  using ValueType = typename ReducerType::value_type;

  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  ValueType* TeamThread_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp declare reduction(                                               \
    omp_red_teamthread_reducer:ValueType                                     \
    : Impl::OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(                                                             \
        Impl::OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

#pragma omp barrier
  ValueType tmp;
  result.init(tmp);
  TeamThread_scratch[0] = tmp;
#pragma omp barrier

  iType team_size = iType(omp_get_num_threads());
#pragma omp for reduction(omp_red_teamthread_reducer \
                          : TeamThread_scratch[:1]) schedule(static, 1)
  for (iType t = 0; t < team_size; t++) {
    ValueType tmp2;
    result.init(tmp2);

    for (iType i = loop_boundaries.start + t; i < loop_boundaries.end;
         i += team_size) {
      lambda(i, tmp2);
    }

    // FIXME_OPENMPTARGET: Join should work but doesn't. Every threads gets a
    // private TeamThread_scratch[0] and at the end of the for-loop the `join`
    // operation is performed by OpenMP itself and hence the simple assignment
    // works.
    //    result.join(TeamThread_scratch[0], tmp2);
    TeamThread_scratch[0] = tmp2;
  }

  result.reference() = TeamThread_scratch[0];
}
#endif  // KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a reduction of val is performed using JoinType(ValueType& val, const
 * ValueType& update) and put into init_result. The input value of init_result
 * is used as initializer for temporary variables of ValueType. Therefore the
 * input value should be the neutral element with respect to the join operation
 * (e.g. '0 for +-' or '1 for *').
 */
template <typename iType, class Lambda, typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda, const JoinType& join, ValueType& init_result) {
  ValueType* TeamThread_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  // FIXME_OPENMPTARGET: Still need to figure out how to get value_count here.
  const int value_count = 1;

#pragma omp barrier
  TeamThread_scratch[0] = init_result;
#pragma omp barrier

#pragma omp for
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i, TeamThread_scratch[omp_get_num_threads() * value_count]);
  }

  // Reduce all partial results within a team.
  const int team_size      = omp_get_num_threads();
  int tree_neighbor_offset = 1;
  do {
#pragma omp for
    for (int i = 0; i < team_size - tree_neighbor_offset;
         i += 2 * tree_neighbor_offset) {
      const int neighbor = i + tree_neighbor_offset;
      join(lambda, &TeamThread_scratch[i * value_count],
           &TeamThread_scratch[neighbor * value_count]);
    }
    tree_neighbor_offset *= 2;
  } while (tree_neighbor_offset < team_size);
  init_result = TeamThread_scratch[0];
}

// This is largely the same code as in HIP and CUDA except for the member name
template <typename iType, class FunctorType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_bounds,
    const FunctorType& lambda) {
  using Analysis   = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                                         TeamPolicy<Experimental::OpenMPTarget>,
                                         FunctorType>;
  using value_type = typename Analysis::value_type;

  const auto start = loop_bounds.start;
  const auto end   = loop_bounds.end;
  //   Note this thing is called .member in the CUDA specialization of
  //   TeamThreadRangeBoundariesStruct
  auto& member         = loop_bounds.team;
  const auto team_size = member.team_size();
  const auto team_rank = member.team_rank();

#if defined(KOKKOS_IMPL_TEAM_SCAN_WORKAROUND)
  value_type scan_val = value_type();

  if (team_rank == 0) {
    for (iType i = start; i < end; ++i) {
      lambda(i, scan_val, true);
    }
  }
#pragma omp barrier
#else
  const auto nchunk = (end - start + team_size - 1) / team_size;
  value_type accum  = 0;
  // each team has to process one or
  //      more chunks of the prefix scan
  for (iType i = 0; i < nchunk; ++i) {
    auto ii = start + i * team_size + team_rank;
    // local accumulation for this chunk
    value_type local_accum = 0;
    // user updates value with prefix value
    if (ii < loop_bounds.end) lambda(ii, local_accum, false);
    // perform team scan
    local_accum = member.team_scan(local_accum);
    // add this blocks accum to total accumulation
    auto val = accum + local_accum;
    // user updates their data with total accumulation
    if (ii < loop_bounds.end) lambda(ii, val, true);
    // the last value needs to be propogated to next chunk
    if (team_rank == team_size - 1) accum = val;
    // broadcast last value to rest of the team
    member.team_broadcast(accum, team_size - 1);
  }
#endif
}

}  // namespace Kokkos

namespace Kokkos {
/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 */
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda) {
#pragma omp simd
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) lambda(i);
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a summation of val is performed and put into result.
 */
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda, ValueType& result) {
  ValueType vector_reduce = ValueType();

  if constexpr (std::is_arithmetic<ValueType>::value) {
#pragma omp simd reduction(+ : vector_reduce)
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      ValueType tmp = ValueType();
      lambda(i, tmp);
      vector_reduce += tmp;
    }
  } else {
#pragma omp declare reduction(custom:ValueType : omp_out += omp_in)

#pragma omp simd reduction(custom : vector_reduce)
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      lambda(i, vector_reduce);
    }
  }

  result = vector_reduce;
}

template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ReducerType const& result) {
  using ValueType = typename ReducerType::value_type;

#pragma omp declare reduction(                                               \
    custom:ValueType                                                         \
    : Impl::OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(                                                             \
        Impl::OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

  ValueType vector_reduce;
  Impl::OpenMPTargetReducerWrapper<ReducerType>::init(vector_reduce);

#pragma omp simd reduction(custom : vector_reduce)
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i, vector_reduce);
  }

  result.reference() = vector_reduce;
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a reduction of val is performed using JoinType(ValueType& val, const
 * ValueType& update) and put into init_result. The input value of init_result
 * is used as initializer for temporary variables of ValueType. Therefore the
 * input value should be the neutral element with respect to the join operation
 * (e.g. '0 for +-' or '1 for *').
 */
template <typename iType, class Lambda, typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda, const JoinType& join, ValueType& init_result) {
  ValueType result = init_result;

  // FIXME_OPENMPTARGET think about omp simd
  // join does not work with omp reduction clause
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    ValueType tmp = ValueType();
    lambda(i, tmp);
    join(result, tmp);
  }

  init_result = result;
}

/** \brief  Intra-thread vector parallel exclusive prefix sum. Executes
 * lambda(iType i, ValueType & val, bool final) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes in the thread and a scan
 * operation is performed. Depending on the target execution space the operator
 * might be called twice: once with final=false and once with final=true. When
 * final==true val contains the prefix sum value. The contribution of this "i"
 * needs to be added to val no matter whether final==true or not. In a serial
 * execution (i.e. team_size==1) the operator is only called once with
 * final==true. Scan_val will be set to the final sum value over all vector
 * lanes.
 */
template <typename iType, class FunctorType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const FunctorType& lambda) {
  using Analysis   = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                                         TeamPolicy<Experimental::OpenMPTarget>,
                                         FunctorType>;
  using value_type = typename Analysis::value_type;

  value_type scan_val = value_type();

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; ++i) {
    lambda(i, scan_val, true);
  }
}

}  // namespace Kokkos

#ifdef KOKKOS_IMPL_TEAM_SCAN_WORKAROUND
#undef KOKKOS_IMPL_TEAM_SCAN_WORKAROUND
#endif

namespace Kokkos {
/** \brief  Intra-team vector parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling team.
 */
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamVectorRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda) {
#pragma omp for simd nowait schedule(static, 1)
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) lambda(i);
}

/** \brief  Intra-team vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling team
 * and a summation of val is performed and put into result.
 */
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamVectorRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda, ValueType& result) {
  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  ValueType* TeamVector_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp barrier
  TeamVector_scratch[0] = ValueType();
#pragma omp barrier

  if constexpr (std::is_arithmetic<ValueType>::value) {
#pragma omp for simd reduction(+ : TeamVector_scratch[:1])
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      ValueType tmp = ValueType();
      lambda(i, tmp);
      TeamVector_scratch[0] += tmp;
    }
  } else {
#pragma omp declare reduction(custom:ValueType : omp_out += omp_in)

#pragma omp for simd reduction(custom : TeamVector_scratch[:1])
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      ValueType tmp = ValueType();
      lambda(i, tmp);
      TeamVector_scratch[0] += tmp;
    }
  }

  result = TeamVector_scratch[0];
}

#if !defined(KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND)
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ReducerType const& result) {
  using ValueType = typename ReducerType::value_type;

  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

#pragma omp declare reduction(                                               \
    custom:ValueType                                                         \
    : Impl::OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(                                                             \
        Impl::OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

  ValueType* TeamVector_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp barrier
  Impl::OpenMPTargetReducerWrapper<ReducerType>::init(TeamVector_scratch[0]);
#pragma omp barrier

#pragma omp for simd reduction(custom : TeamVector_scratch[:1])
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i, TeamVector_scratch[0]);
  }

  result.reference() = TeamVector_scratch[0];
}
#else
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ReducerType const& result) {
  using ValueType = typename ReducerType::value_type;

  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  ValueType* TeamVector_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp declare reduction(                                               \
    omp_red_teamthread_reducer:ValueType                                     \
    : Impl::OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(                                                             \
        Impl::OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

#pragma omp barrier
  ValueType tmp;
  result.init(tmp);
  TeamVector_scratch[0] = tmp;
#pragma omp barrier

  iType team_size = iType(omp_get_num_threads());
#pragma omp for simd reduction(omp_red_teamthread_reducer \
                               : TeamVector_scratch[:1]) schedule(static, 1)
  for (iType t = 0; t < team_size; t++) {
    ValueType tmp2;
    result.init(tmp2);

    for (iType i = loop_boundaries.start + t; i < loop_boundaries.end;
         i += team_size) {
      lambda(i, tmp2);
    }
    TeamVector_scratch[0] = tmp2;
  }

  result.reference() = TeamVector_scratch[0];
}
#endif  // KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND
}  // namespace Kokkos

#ifdef KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND
#undef KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND
#endif

namespace Kokkos {

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember>&
    /*single_struct*/,
    const FunctorType& lambda) {
  lambda();
}

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember>&
        single_struct,
    const FunctorType& lambda) {
  if (single_struct.team_member.team_rank() == 0) lambda();
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember>&
    /*single_struct*/,
    const FunctorType& lambda, ValueType& val) {
  lambda(val);
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember>&
        single_struct,
    const FunctorType& lambda, ValueType& val) {
  if (single_struct.team_member.team_rank() == 0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val, 0);
}
}  // namespace Kokkos

#endif /* #ifndef KOKKOS_OPENMPTARGETEXEC_HPP */
