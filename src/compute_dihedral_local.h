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
ComputeStyle(dihedral/local,ComputeDihedralLocal);
// clang-format on
#else

#ifndef LMP_COMPUTE_DIHEDRAL_LOCAL_H
#define LMP_COMPUTE_DIHEDRAL_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDihedralLocal : public Compute {
 public:
  ComputeDihedralLocal(class LAMMPS *, int, char **);
  ~ComputeDihedralLocal() override;
  void init() override;
  void compute_local() override;
  double memory_usage() override;

 private:
  int nvalues, nvar, ncount, setflag;

  int pvar;
  int *bstyle, *vvar;
  char *pstr;
  char **vstr;

  int nmax;
  double *vlocal;
  double **alocal;

  int compute_dihedrals(int);
  void reallocate(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
