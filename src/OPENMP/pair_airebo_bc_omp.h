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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(airebo/bc/omp,PairAIREBObcOMP);
// clang-format on
#else

#ifndef LMP_PAIR_AIREBO_BC_OMP_H
#define LMP_PAIR_AIREBO_BC_OMP_H

#include "pair_airebo_omp.h"

namespace LAMMPS_NS {

// AIREBO-BC (bond-centric P_CC) with OpenMP threading.  Like airebo/morse/omp,
// it is the threaded base class with one parameter flag set: bcflag = 1
// selects the bond-centric P_CC machinery that lives in PairAIREBO.

class PairAIREBObcOMP : public PairAIREBOOMP {
 public:
  PairAIREBObcOMP(class LAMMPS *);
};

}    // namespace LAMMPS_NS

#endif
#endif
