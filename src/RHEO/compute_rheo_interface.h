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
ComputeStyle(RHEO/INTERFACE,ComputeRHEOInterface)
// clang-format on
#else

#ifndef LMP_COMPUTE_RHEO_INTERFACE_H
#define LMP_COMPUTE_RHEO_INTERFACE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRHEOInterface : public Compute {
 public:
  ComputeRHEOInterface(class LAMMPS *, int, char **);
  ~ComputeRHEOInterface() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
  void correct_v(double *, double *, int, int);
  double correct_rho(int);
  void store_forces();

  double *chi, **fp_store;
  class FixRHEO *fix_rheo;

 private:
  int nmax_store, comm_stage;
  double *rho0, cut, cutsq, wall_max;
  double *norm, *normwf;

  char *id_fix_pa;

  class NeighList *list;
  class ComputeRHEOKernel *compute_kernel;
  class FixRHEOPressure *fix_pressure;
};

}    // namespace LAMMPS_NS

#endif
#endif
