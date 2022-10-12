/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(efield/atom,ComputeEfieldAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_EFIELD_ATOM_H
#define LMP_COMPUTE_EFIELD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEfieldAtom : public Compute {
 public:
  ComputeEfieldAtom(class LAMMPS *, int, char **);
  ~ComputeEfieldAtom() override;
  void init() override;
  void setup() override;
  void compute_peratom() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

 private:
  int pairflag;
  int kspaceflag;
  double **efield_pair, **efield_kspace;

  int nmax;
  double **efield;
};

}    // namespace LAMMPS_NS

#endif
#endif
