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

#ifdef NBIN_CLASS
// clang-format off
NBinStyle(multi,
          NBinMulti,
          NB_MULTI);
// clang-format on
#else

#ifndef LMP_NBIN_MULTI_H
#define LMP_NBIN_MULTI_H

#include "nbin.h"

namespace LAMMPS_NS {

class NBinMulti : public NBin {
 public:
  NBinMulti(class LAMMPS *);
  ~NBinMulti() {}
  void bin_atoms_setup(int);
  void setup_bins(int);
  void bin_atoms();
  double memory_usage();
};

}    // namespace LAMMPS_NS

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
