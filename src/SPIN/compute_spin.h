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
ComputeStyle(spin,ComputeSpin);
// clang-format on
#else

#ifndef LMP_COMPUTE_SPIN_H
#define LMP_COMPUTE_SPIN_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSpin : public Compute {
 public:
  ComputeSpin(class LAMMPS *, int, char **);
  ~ComputeSpin() override;
  void init() override;
  void compute_vector() override;

 private:
  int pair_spin_flag;          // magnetic pair flags
  int long_spin_flag;          // magnetic long-range flag
  int precession_spin_flag;    // magnetic precession flags

  double kb, hbar;

  // pointers to magnetic fixes

  int nprecspin;
  class FixPrecessionSpin **lockprecessionspin;

  // pointers to magnetic pair styles

  int npairs, npairspin;    // # of pairs, and # of spin pairs
  class Pair *pair;
  class PairSpin **spin_pairs;    // vector of spin pairs

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
