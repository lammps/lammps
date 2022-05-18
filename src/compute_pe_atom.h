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
ComputeStyle(pe/atom,ComputePEAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_PE_ATOM_H
#define LMP_COMPUTE_PE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePEAtom : public Compute {
 public:
  ComputePEAtom(class LAMMPS *, int, char **);
  ~ComputePEAtom() override;
  void init() override {}
  void compute_peratom() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

 private:
  int pairflag, bondflag, angleflag, dihedralflag, improperflag;
  int kspaceflag, fixflag;
  int nmax;
  double *energy;
};

}    // namespace LAMMPS_NS

#endif
#endif
