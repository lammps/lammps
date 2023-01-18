/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(entropy/atom,ComputeEntropyAtom);
// clang-format on
#else

#ifndef COMPUTE_ENTROPY_ATOM_H
#define COMPUTE_ENTROPY_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEntropyAtom : public Compute {
 public:
  ComputeEntropyAtom(class LAMMPS *, int, char **);
  ~ComputeEntropyAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nmax, maxneigh, nbin;
  class NeighList *list;
  double *pair_entropy, *pair_entropy_avg;
  double sigma, cutoff, cutoff2;
  double cutsq, cutsq2;
  double deltar;
  int deltabin;
  int avg_flag;
  int local_flag;
};

}    // namespace LAMMPS_NS

#endif
#endif
