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
  double compute_scalar() override;
  double memory_usage() override;

 private:
  FILE * fh_d;
  int natoms, nmax, size_peratom, lastcol;
  int ncoeff, nperdim, yoffset, zoffset;
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
  //int bikflag;
  //int bik_rows;
  int bikflag, bik_rows, dbirjflag, dbirj_rows;
  double **dbirj;
  double **dbiri; // dBi/dRi = sum(-dBi/dRj) over neighbors j
  int *nneighs; // number of neighs inside the snap cutoff.
  int *neighsum;
  int *icounter; // counting atoms i for each j.

  Compute *c_pe;
  Compute *c_virial;

  void dbdotr_compute();
  void get_dbirj_length();
};

}    // namespace LAMMPS_NS

#endif
#endif
