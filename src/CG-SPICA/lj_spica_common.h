/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Common data for the Shinoda, DeVane, Klein (SDK) coarse grain model
   Contributing author: Axel Kohlmeyer (Temple U)
   Contributing author (added LJ12_5) : Yusuke Miyazaki (Okayama U)
------------------------------------------------------------------------- */

#ifndef LMP_LJ_SPICA_COMMON_H
#define LMP_LJ_SPICA_COMMON_H

#include <cstring>

namespace LAMMPS_NS {
namespace LJSPICAParms {

  // LJ type flags. list of supported LJ exponent combinations
  enum { LJ_NOT_SET = 0, LJ9_6, LJ12_4, LJ12_6, LJ12_5, NUM_LJ_TYPES };

#if defined(LMP_NEED_SPICA_FIND_LJ_TYPE)
  static int find_lj_type(const char *label, const char *const *const list)
  {
    for (int i = 0; i < NUM_LJ_TYPES; ++i)
      if (strcmp(label, list[i]) == 0) return i;

    return LJ_NOT_SET;
  }
#endif

  // clang-format off
  static const char *const lj_type_list[] = {"none", "lj9_6", "lj12_4",          "lj12_6", "lj12_5"};
  static constexpr double lj_prefact[]    = {0.0,     6.75,    2.59807621135332,  4.0,      3.20377984125109};
  static constexpr double lj_pow1[]       = {0.0,     9.00,   12.0,              12.0,     12.0};
  static constexpr double lj_pow2[]       = {0.0,     6.00,    4.0,               6.0,      5.0};
  // clang-format on
}    // namespace LJSPICAParms
}    // namespace LAMMPS_NS
#endif
