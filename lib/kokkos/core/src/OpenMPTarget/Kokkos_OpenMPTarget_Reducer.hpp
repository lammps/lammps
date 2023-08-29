//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_OPENMPTARGETREDUCER_HPP
#define KOKKOS_OPENMPTARGETREDUCER_HPP

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Spinwait.hpp>

#include <Kokkos_Atomic.hpp>
#include "Kokkos_OpenMPTarget_Abort.hpp"

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

#endif
