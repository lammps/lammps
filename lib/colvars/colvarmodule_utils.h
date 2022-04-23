// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#ifndef COLVARMODULE_UTILS_H
#define COLVARMODULE_UTILS_H


#include "colvarmodule.h"


template <typename T>
cvm::real get_force_norm2(T const &x)
{
  return x.norm2();
}


template <>
inline cvm::real get_force_norm2(cvm::real const &x)
{
  return x*x;
}


template <typename T, int flag, bool get_index>
cvm::real compute_norm2_stats(std::vector<T> const &v,
                              int *minmax_index = NULL)
{
  cvm::real result = 0.0;
  if (flag == -1) {
    // Initialize for minimum search, using approx. largest float32 value
    result = 1.0e38;
  }

  typename std::vector<T>::const_iterator xi = v.begin();
  size_t i = 0;

  if (get_index) *minmax_index = -1; // Let's not assume minmax_index is initialized to -1

  for ( ; xi != v.end(); xi++, i++) {
    cvm::real const norm2 = get_force_norm2<T>(*xi);
    if (flag == 0) {
      result += norm2;
    }
    if (flag == 1) {
      if (norm2 > result) {
        result = norm2;
        if (get_index) *minmax_index = i;
      }
    }
    if (flag == -1) {
      if (norm2 < result) {
        result = norm2;
        if (get_index) *minmax_index = i;
      }
    }
  }

  size_t const n = v.size();

  if (flag == 0) {
    if (n > 0) {
      result /= cvm::real(n);
    }
  }

  result = cvm::sqrt(result);

  return result;
}


#endif
