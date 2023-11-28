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

#ifndef KOKKOS_BIN_OPS_PUBLIC_API_HPP_
#define KOKKOS_BIN_OPS_PUBLIC_API_HPP_

#include <Kokkos_Macros.hpp>
#include <type_traits>

namespace Kokkos {

template <class KeyViewType>
struct BinOp1D {
  int max_bins_ = {};
  double mul_   = {};
  double min_   = {};

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED BinOp1D() = default;
#else
  BinOp1D() = delete;
#endif

  // Construct BinOp with number of bins, minimum value and maximum value
  BinOp1D(int max_bins__, typename KeyViewType::const_value_type min,
          typename KeyViewType::const_value_type max)
      : max_bins_(max_bins__ + 1),
        // Cast to double to avoid possible overflow when using integer
        mul_(static_cast<double>(max_bins__) /
             (static_cast<double>(max) - static_cast<double>(min))),
        min_(static_cast<double>(min)) {
    // For integral types the number of bins may be larger than the range
    // in which case we can exactly have one unique value per bin
    // and then don't need to sort bins.
    if (std::is_integral<typename KeyViewType::const_value_type>::value &&
        (static_cast<double>(max) - static_cast<double>(min)) <=
            static_cast<double>(max_bins__)) {
      mul_ = 1.;
    }
  }

  // Determine bin index from key value
  template <class ViewType>
  KOKKOS_INLINE_FUNCTION int bin(ViewType& keys, const int& i) const {
    return static_cast<int>(mul_ * (static_cast<double>(keys(i)) - min_));
  }

  // Return maximum bin index + 1
  KOKKOS_INLINE_FUNCTION
  int max_bins() const { return max_bins_; }

  // Compare to keys within a bin if true new_val will be put before old_val
  template <class ViewType, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION bool operator()(ViewType& keys, iType1& i1,
                                         iType2& i2) const {
    return keys(i1) < keys(i2);
  }
};

template <class KeyViewType>
struct BinOp3D {
  int max_bins_[3] = {};
  double mul_[3]   = {};
  double min_[3]   = {};

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED BinOp3D() = default;
#else
  BinOp3D() = delete;
#endif

  BinOp3D(int max_bins__[], typename KeyViewType::const_value_type min[],
          typename KeyViewType::const_value_type max[]) {
    max_bins_[0] = max_bins__[0];
    max_bins_[1] = max_bins__[1];
    max_bins_[2] = max_bins__[2];
    mul_[0]      = static_cast<double>(max_bins__[0]) /
              (static_cast<double>(max[0]) - static_cast<double>(min[0]));
    mul_[1] = static_cast<double>(max_bins__[1]) /
              (static_cast<double>(max[1]) - static_cast<double>(min[1]));
    mul_[2] = static_cast<double>(max_bins__[2]) /
              (static_cast<double>(max[2]) - static_cast<double>(min[2]));
    min_[0] = static_cast<double>(min[0]);
    min_[1] = static_cast<double>(min[1]);
    min_[2] = static_cast<double>(min[2]);
  }

  template <class ViewType>
  KOKKOS_INLINE_FUNCTION int bin(ViewType& keys, const int& i) const {
    return int((((int(mul_[0] * (keys(i, 0) - min_[0])) * max_bins_[1]) +
                 int(mul_[1] * (keys(i, 1) - min_[1]))) *
                max_bins_[2]) +
               int(mul_[2] * (keys(i, 2) - min_[2])));
  }

  KOKKOS_INLINE_FUNCTION
  int max_bins() const { return max_bins_[0] * max_bins_[1] * max_bins_[2]; }

  template <class ViewType, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION bool operator()(ViewType& keys, iType1& i1,
                                         iType2& i2) const {
    if (keys(i1, 0) > keys(i2, 0))
      return true;
    else if (keys(i1, 0) == keys(i2, 0)) {
      if (keys(i1, 1) > keys(i2, 1))
        return true;
      else if (keys(i1, 1) == keys(i2, 1)) {
        if (keys(i1, 2) > keys(i2, 2)) return true;
      }
    }
    return false;
  }
};

}  // namespace Kokkos
#endif
