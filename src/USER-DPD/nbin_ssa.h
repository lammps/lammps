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

#ifdef NBIN_CLASS

NBinStyle(ssa,
          NBinSSA,
          NB_SSA)

#else

#ifndef LMP_NBIN_SSA_H
#define LMP_NBIN_SSA_H

#include "nbin_standard.h"

namespace LAMMPS_NS {

class NBinSSA : public NBinStandard {
 public:

  int *bins_ssa;             // index of next atom in each bin
  int maxbin_ssa;            // size of bins_ssa array
  int *binhead_ssa;          // index of 1st local atom in each bin
  int gairhead_ssa[9];       // index of 1st ghost atom in each AIR
  int maxhead_ssa;           // size of binhead_ssa and gbinhead_ssa arrays

  NBinSSA(class LAMMPS *);
  ~NBinSSA();

  void bin_atoms_setup(int);
  void bin_atoms();

  bigint memory_usage();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
