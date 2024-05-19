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
ComputeStyle(RHEO/VSHIFT,ComputeRHEOVShift)
// clang-format on
#else

#ifndef LMP_COMPUTE_RHEO_VSHIFT_H
#define LMP_COMPUTE_RHEO_VSHIFT_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRHEOVShift : public Compute {
 public:
  ComputeRHEOVShift(class LAMMPS *, int, char **);
  ~ComputeRHEOVShift() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
  void correct_surfaces();
  double **vshift;

  class FixRHEO *fix_rheo;

 private:
  int nmax_store;
  double dtv, cut, cutsq, cutthird;
  int surface_flag, interface_flag;
  double *rho0;

  class NeighList *list;
  class ComputeRHEOInterface *compute_interface;
  class ComputeRHEOKernel *compute_kernel;
  class ComputeRHEOSurface *compute_surface;
};

}    // namespace LAMMPS_NS

#endif
#endif
