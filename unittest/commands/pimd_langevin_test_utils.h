/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with the U.S. Government retains
   certain rights in this software.
------------------------------------------------------------------------- */

#ifndef LMP_UNITTEST_PIMD_LANGEVIN_TEST_UTILS_H
#define LMP_UNITTEST_PIMD_LANGEVIN_TEST_UTILS_H

#include "library.h"

#include <string>
#include <vector>

namespace LAMMPS_NS {
namespace pimd_langevin_test {

constexpr int kNuclearPrefixScalars = 10;
constexpr int kBosonicVectorSize = 6;
constexpr int kIsoBarostatVectorSize = 15;
constexpr int kAnisoBarostatVectorSize = 17;

enum NuclearVectorIndex {
  KE_BEAD = 0,
  SE_BEAD = 1,
  PE_BEAD = 2,
  TOTAL_ENERGY = 3,
  T_PRIM = 4,
  T_VIR = 5,
  T_CV = 6,
  P_PRIM = 7,
  P_MD = 8,
  P_CV = 9
};

struct IsoBarostatIndices {
  int vw0;
  int barostat_ke;
  int pext_volume;
  int log_volume_term;
  int total_enthalpy;
};

inline IsoBarostatIndices iso_barostat_indices()
{
  return {10, 11, 12, 13, 14};
}

inline int fix_vector_size(void *lmp, const char *id)
{
  auto *ptr = (int *) lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_SIZE_VECTOR, 0, 0);
  return *ptr;
}

inline double fix_value(void *lmp, const char *id, int index)
{
  void *ptr = lammps_extract_fix(lmp, id, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR, index, 0);
  double value = *(double *) ptr;
  lammps_free(ptr);
  return value;
}

inline std::vector<double> fix_vector(void *lmp, const char *id)
{
  const int size = fix_vector_size(lmp, id);
  std::vector<double> values(size);
  for (int i = 0; i < size; ++i) values[i] = fix_value(lmp, id, i);
  return values;
}

inline std::string last_error(void *lmp)
{
  char buffer[2048];
  buffer[0] = '\0';
  if (lammps_get_last_error_message(lmp, buffer, sizeof(buffer)) == 0) return "";
  return buffer;
}

inline int lammps_fix_index(int zero_based_index)
{
  return zero_based_index + 1;
}

}    // namespace pimd_langevin_test
}    // namespace LAMMPS_NS

#endif
