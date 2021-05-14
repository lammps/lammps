/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NBIN_CLASS

NBinStyle(standard,
          NBinStandard,
          NB_STANDARD)

#else

#ifndef LMP_NBIN_STANDARD_H
#define LMP_NBIN_STANDARD_H

#include "nbin.h"

namespace LAMMPS_NS {

class NBinStandard : public NBin {
 public:
  NBinStandard(class LAMMPS *);
  ~NBinStandard() {}
  void bin_atoms_setup(int);
  void setup_bins(int);
  void bin_atoms();
  double memory_usage();
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

*/
