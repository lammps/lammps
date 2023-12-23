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

#ifdef FIX_CLASS
// clang-format off
FixStyle(rheo/tension,FixRHEOTension)
// clang-format on
#else

#ifndef LMP_FIX_RHEO_TENSION_H
#define LMP_FIX_RHEO_TENSION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRHEOTension : public Fix {
 public:
  FixRHEOTension(class LAMMPS *, int, char **);
  ~FixRHEOTension() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void pre_force(int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  void grow_arrays(int) override;

 private:
  int nmax_store, comm_stage, interface_flag, shift_flag;
  int index_ct, index_nt, index_cgradt, index_divnt, index_ft, index_wsame;

  double *ct, **nt, **cgradt, *divnt, *norm, **ft, *wsame;
  double alpha, beta, wmin, cmin, vshift_strength, h, hsq, hinv, hinv3, *rho0;

  class ComputeRHEOKernel *compute_kernel;
  class ComputeRHEOInterface *compute_interface;
  class ComputeRHEOVShift *compute_vshift;
  class FixRHEO *fix_rheo;
  class NeighList *list;
};

}    // namespace LAMMPS_NS

#endif
#endif
