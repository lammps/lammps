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
ComputeStyle(RHEO/SURFACE,ComputeRHEOSurface)
// clang-format on
#else

#ifndef LMP_COMPUTE_RHEO_SURFACE_H
#define LMP_COMPUTE_RHEO_SURFACE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRHEOSurface : public Compute {
 public:
  ComputeRHEOSurface(class LAMMPS *, int, char **);
  ~ComputeRHEOSurface() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  double **nsurface, *rsurface, *divr;
  class FixRHEO *fix_rheo;

 private:
  int surface_style, nmax_store, threshold_z, threshold_splash, interface_flag;
  int threshold_style, comm_stage;

  double cut, cutsq, *rho0, threshold_divr;
  double **B, **gradC;

  class NeighList *list;
  class ComputeRHEOKernel *compute_kernel;
  class ComputeRHEOInterface *compute_interface;

  void grow_arrays(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
