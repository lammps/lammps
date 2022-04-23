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
ComputeStyle(coord/atom,ComputeCoordAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_COORD_ATOM_H
#define LMP_COMPUTE_COORD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCoordAtom : public Compute {
 public:
  ComputeCoordAtom(class LAMMPS *, int, char **);
  ~ComputeCoordAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double memory_usage() override;
  enum { NONE, CUTOFF, ORIENT };

 protected:
  int nmax, ncol;
  double cutsq;
  class NeighList *list;

  int *typelo, *typehi;
  double *cvec;
  double **carray;

  char *group2;
  int jgroup, jgroupbit;

  class ComputeOrientOrderAtom *c_orientorder;
  char *id_orientorder;
  double threshold;
  double **normv;
  int cstyle, nqlist, l;
};

}    // namespace LAMMPS_NS

#endif
#endif
