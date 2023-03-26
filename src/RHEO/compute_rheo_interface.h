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
ComputeStyle(rheo/interface,ComputeRHEOInterface)
// clang-format on
#else

#ifndef LMP_COMPUTE_RHEO_INTERFACE_H
#define LMP_COMPUTE_RHEO_INTERFACE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRHEOInterface : public Compute {
 public:
  ComputeRHEOInterface(class LAMMPS *, int, char **);
  ~ComputeRHEOInterface();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void correct_v(double *, double *, double *, int, int);
  double correct_rho(int, int);
  double memory_usage();
  void store_forces();
  double *chi;

 private:
  int nmax;
  double cut, cutsq, cs2;
  class NeighList *list;
  double *norm, *normwf, wall_max;

  class ComputeRHEOKernel *compute_kernel;
  char *id_fix_chi;
  class FixStore *fix_chi;

  int index_fx, index_fy, index_fz;
  int comm_stage;
};

}

#endif
#endif
