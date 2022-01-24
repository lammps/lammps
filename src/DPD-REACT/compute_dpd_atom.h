/* -*- c++ -*- ----------------------------------------------------------------------
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
ComputeStyle(dpd/atom,ComputeDpdAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_DPD_ATOM_H
#define LMP_COMPUTE_DPD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDpdAtom : public Compute {
 public:
  ComputeDpdAtom(class LAMMPS *, int, char **);
  ~ComputeDpdAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double **dpdAtom;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: compute dpd requires atom_style with internal temperature and energies (e.g. dpd)

Self-explanatory

W:  More than one compute dpd/atom command

Self-explanatory

*/
