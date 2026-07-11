/* -*- c++ -*- ---------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(nb3b/screened,PairNb3bScreened);
// clang-format on
#else

#ifndef LMP_PAIR_NB3B_SCREENED_H
#define LMP_PAIR_NB3B_SCREENED_H

#include "pair_nb3b_harmonic.h"

namespace LAMMPS_NS {

class PairNb3bScreened : public PairNb3bHarmonic {
 public:
  PairNb3bScreened(class LAMMPS *);

 protected:
  void threebody(Param *, Param *, Param *, double, double, double *, double *, double *, double *,
                 int, double &) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
