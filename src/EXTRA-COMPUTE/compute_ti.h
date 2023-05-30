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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(ti,ComputeTI);
// clang-format on
#else

#ifndef COMPUTE_TI_H
#define COMPUTE_TI_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTI : public Compute {
 public:
  ComputeTI(class LAMMPS *, int, char **);
  ~ComputeTI() override;
  void init() override;
  double compute_scalar() override;

 private:
  int nterms;
  int *which;
  int *ivar1, *ivar2;
  int *ilo, *ihi;
  char **var1, **var2;
  class Pair **pptr;
  char **pstyle;
};

}    // namespace LAMMPS_NS

#endif
#endif
