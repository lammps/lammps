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

#ifndef LMP_COMPUTE_RHEO_INTERFACE_H
#define LMP_COMPUTE_RHEO_INTERFACE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRHEOSurface : public Compute {
 public:
  ComputeRHEOSurface(class LAMMPS *, int, char **);
  ~ComputeRHEOSurface();
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  double **gradC, **n_surface;

 private:
  double cut, cutsq, threshold;
  int surface_style, nmax_old;
  double **B, *divr;
  int comm_stage;

  int index_divr;
  int index_rsurf;

  double divR_limit;
  int coord_limit;

  class NeighList *list;
  class FixRHEO *fix_rheo;
  class ComputeRHEOKernel *compute_kernel;
  class ComputeRHEOSolids *compute_solids;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
