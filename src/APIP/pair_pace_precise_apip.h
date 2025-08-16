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
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pace/precise/apip,PairPACEPreciseAPIP);
// clang-format on
#else

#ifndef LMP_PAIR_PACE_PRECISE_APIP_H
#define LMP_PAIR_PACE_PRECISE_APIP_H

#include "pair_pace_apip.h"

namespace LAMMPS_NS {

class PairPACEPreciseAPIP : public PairPACEAPIP {
 public:
  PairPACEPreciseAPIP(class LAMMPS *);
  void *extract(const char *, int &) override;

 protected:
  double *get_e_ref_ptr() override;
  double compute_factor_lambda(double) override;
  int check_abort_condition(double *, double *, int *, int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
