/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with the U.S. Government retains
   certain rights in this software.
------------------------------------------------------------------------- */

#ifndef LMP_UNITTEST_PIMD_TEST_UTILS_H
#define LMP_UNITTEST_PIMD_TEST_UTILS_H

#include "library.h"

#include <string>

namespace LAMMPS_NS {
namespace pimd_test {

constexpr int kDefaultTChain = 3;
constexpr int kNuclearPrefixScalars = 10;

struct UVTVectorIndices {
  int ne;
  int ne_dot;
  int dedn;
  int mu;
  int ne_ke;
  int electronic_potential;
};

inline int nuclear_vector_size(int tchain = kDefaultTChain)
{
  return kNuclearPrefixScalars + 4 * tchain;
}

inline UVTVectorIndices uvt_vector_indices(int tchain = kDefaultTChain)
{
  const int prefix = nuclear_vector_size(tchain);
  return {prefix, prefix + 1, prefix + 2, prefix + 3, prefix + 4, prefix + 5};
}

inline int lammps_fix_index(int zero_based_index)
{
  return zero_based_index + 1;
}

inline double fix_value(void *lmp, const char *id, int index)
{
  void *ptr = lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, index, 0);
  double value = *(double *) ptr;
  lammps_free(ptr);
  return value;
}

inline std::string quadratic_dedn_variable(const char *fix_id = "cp", int tchain = kDefaultTChain)
{
  const auto indices = uvt_vector_indices(tchain);
  return "variable dEdN equal v_k_quad*(f_" + std::string(fix_id) + "[" +
      std::to_string(lammps_fix_index(indices.ne)) + "]-v_N0_quad)";
}

}    // namespace pimd_test
}    // namespace LAMMPS_NS

#endif
