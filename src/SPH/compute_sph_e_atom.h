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
ComputeStyle(sph/e/atom,ComputeSPHEAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_SPH_E_ATOM_H
#define LMP_COMPUTE_SPH_E_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSPHEAtom : public Compute {
 public:
  ComputeSPHEAtom(class LAMMPS *, int, char **);
  ~ComputeSPHEAtom();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double *evector;
};

}    // namespace LAMMPS_NS

#endif
#endif
