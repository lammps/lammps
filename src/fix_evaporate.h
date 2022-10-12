/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(evaporate,FixEvaporate);
// clang-format on
#else

#ifndef LMP_FIX_EVAPORATE_H
#define LMP_FIX_EVAPORATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEvaporate : public Fix {
 public:
  FixEvaporate(class LAMMPS *, int, char **);
  ~FixEvaporate() override;
  int setmask() override;
  void init() override;
  void pre_exchange() override;
  double compute_scalar() override;
  double memory_usage() override;

 private:
  int nevery, nflux;
  int molflag;
  int ndeleted;
  char *idregion;
  class Region *region;

  int nmax;
  int *list, *mark;

  class RanPark *random;
};

}    // namespace LAMMPS_NS

#endif
#endif
