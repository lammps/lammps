// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nve/intel,FixNVEIntel);
// clang-format on
#else

#ifndef LMP_FIX_NVE_INTEL_H
#define LMP_FIX_NVE_INTEL_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVEIntel : public FixNVE {
 public:
  FixNVEIntel(class LAMMPS *, int, char **);
  ~FixNVEIntel() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void reset_dt() override;
  double memory_usage() override;

 protected:
  double *_dtfm;
  int _nlocal3, _nlocal_max;
};

}    // namespace LAMMPS_NS

#endif
#endif
