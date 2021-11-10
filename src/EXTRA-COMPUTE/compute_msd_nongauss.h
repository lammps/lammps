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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(msd/nongauss,ComputeMSDNonGauss);
// clang-format on
#else

#ifndef LMP_COMPUTE_MSD_NONGAUSS_H
#define LMP_COMPUTE_MSD_NONGAUSS_H

#include "compute_msd.h"

namespace LAMMPS_NS {

class ComputeMSDNonGauss : public ComputeMSD {
 public:
  ComputeMSDNonGauss(class LAMMPS *, int, char **);
  ~ComputeMSDNonGauss() {}
  void compute_vector();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find compute msd fix ID

Self-explanatory.

*/
