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
ComputeStyle(cna/atom,ComputeCNAAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_CNA_ATOM_H
#define LMP_COMPUTE_CNA_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCNAAtom : public Compute {
 public:
  ComputeCNAAtom(class LAMMPS *, int, char **);
  ~ComputeCNAAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nmax;
  double cutsq;
  class NeighList *list;
  int **nearest;
  int *nnearest;
  double *pattern;
};

}    // namespace LAMMPS_NS

#endif
#endif
