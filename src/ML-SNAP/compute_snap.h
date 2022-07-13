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
ComputeStyle(snap,ComputeSnap);
// clang-format on
#else

#ifndef LMP_COMPUTE_SNAP_H
#define LMP_COMPUTE_SNAP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSnap : public Compute {
 public:
  ComputeSnap(class LAMMPS *, int, char **);
  ~ComputeSnap() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  double memory_usage() override;

 private:
  int natoms, nmax, size_peratom, lastcol;
  int ncoeff, nvalues, yoffset, zoffset;
  int ndims_peratom, ndims_force, ndims_virial;
  double **cutsq;
  class NeighList *list;
  double **snap, **snapall;
  double **snap_peratom;
  double rcutfac;
  double *radelem;
  double *wjelem;
  int *map;    // map types to [0,nelements)
  int nelements, chemflag;
  int switchinnerflag;
  double *sinnerelem;
  double *dinnerelem;
  class SNA *snaptr;
  double cutmax;
  int quadraticflag;
  int bikflag, bik_rows, dgradflag, dgrad_rows;

  Compute *c_pe;
  Compute *c_virial;

  void dbdotr_compute();
};

}    // namespace LAMMPS_NS

#endif
#endif
