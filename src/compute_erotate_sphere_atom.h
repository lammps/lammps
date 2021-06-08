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
ComputeStyle(erotate/sphere/atom,ComputeErotateSphereAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_EROTATE_SPHERE_ATOM_H
#define LMP_COMPUTE_EROTATE_SPHERE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeErotateSphereAtom : public Compute {
 public:
  ComputeErotateSphereAtom(class LAMMPS *, int, char **);
  ~ComputeErotateSphereAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double pfactor;
  double *erot;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute erotate/sphere/atom requires atom style sphere

Self-explanatory.

W: More than one compute erotate/sphere/atom

It is not efficient to use compute erorate/sphere/atom more than once.

*/
