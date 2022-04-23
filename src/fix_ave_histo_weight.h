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
FixStyle(ave/histo/weight,FixAveHistoWeight);
// clang-format on
#else

#ifndef LMP_FIX_AVE_HISTO_WEIGHT_H
#define LMP_FIX_AVE_HISTO_WEIGHT_H

#include "fix_ave_histo.h"

namespace LAMMPS_NS {

class FixAveHistoWeight : public FixAveHisto {
 public:
  FixAveHistoWeight(class LAMMPS *, int, char **);

  void end_of_step() override;

 private:
  void bin_one_weights(double, double);
  void bin_vector_weights(int, double *, int, double *, int);
  void bin_atoms_weights(double *, int, double *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
