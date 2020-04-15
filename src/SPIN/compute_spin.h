/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(spin,ComputeSpin)

#else

#ifndef LMP_COMPUTE_SPIN_H
#define LMP_COMPUTE_SPIN_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSpin : public Compute {
 public:
  ComputeSpin(class LAMMPS *, int, char **);
  ~ComputeSpin();
  void init();
  void compute_vector();

 private:
  int pair_spin_flag;                   // magnetic pair flags
  int long_spin_flag;                   // magnetic long-range flag
  int precession_spin_flag;             // magnetic precession flags

  double kb,hbar;

  // pointers to magnetic fixes

  class FixPrecessionSpin *lockprecessionspin;

  // pointers to magnetic pair styles

  int npairs, npairspin;                // # of pairs, and # of spin pairs
  class Pair *pair;
  class PairSpin **spin_pairs;          // vector of spin pairs

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Chunk/atom compute does not exist for compute compute/spin

Self-explanatory.

E: Compute compute/spin does not use chunk/atom compute

The style of the specified compute is not chunk/atom.

*/
