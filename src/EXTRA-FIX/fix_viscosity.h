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

#ifdef FIX_CLASS
// clang-format off
FixStyle(viscosity,FixViscosity);
// clang-format on
#else

#ifndef LMP_FIX_VISCOSITY_H
#define LMP_FIX_VISCOSITY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixViscosity : public Fix {
 public:
  FixViscosity(class LAMMPS *, int, char **);
  ~FixViscosity() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  double compute_scalar() override;

 private:
  int me;
  int vdim, pdim, nbin, periodicity;
  int nswap;
  double vtarget;
  double prd, boxlo, boxhi;
  double slablo_lo, slablo_hi, slabhi_lo, slabhi_hi;
  double p_exchange;

  int npositive, nnegative;
  int *pos_index, *neg_index;
  double *pos_delta, *neg_delta;
};

}    // namespace LAMMPS_NS

#endif
#endif
