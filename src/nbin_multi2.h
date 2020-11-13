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

NBinStyle(multi2,
          NBinMulti2,
          NB_MULTI2)

#else

#ifndef LMP_NBIN_MULTI2_H
#define LMP_NBIN_MULTI2_H

#include "nbin.h"

namespace LAMMPS_NS {

class NBinMulti2 : public NBin {
 public:

  NBinMulti2(class LAMMPS *);
  ~NBinMulti2() {}
  void bin_atoms_setup(int);
  void setup_bins(int);
  void bin_atoms();
  double memory_usage();

 private:

  int itype_min();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Domain too large for neighbor bins

UNDOCUMENTED

E: Cannot use neighbor bins - box size << cutoff

UNDOCUMENTED

E: Too many neighbor bins

UNDOCUMENTED

E Non-numeric positions - simulation unstable

UNDOCUMENTED

*/
