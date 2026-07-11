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
ComputeStyle(RHEO/GRAD,ComputeRHEOGrad)
// clang-format on
#else

#ifndef LMP_COMPUTE_RHEO_GRAD_H
#define LMP_COMPUTE_RHEO_GRAD_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRHEOGrad : public Compute {
 public:
  ComputeRHEOGrad(class LAMMPS *, int, char **);
  ~ComputeRHEOGrad() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;
  void forward_gradients();
  void forward_fields();
  double **gradv;
  double **gradr;
  double **grade;
  double **gradn;
  class FixRHEO *fix_rheo;

 private:
  int comm_stage, ncomm_grad, ncomm_field, nmax_store;
  double cut, cutsq, *rho0;

  int velocity_flag, energy_flag, rho_flag, eta_flag;
  int interface_flag, remap_v_flag;

  class ComputeRHEOKernel *compute_kernel;
  class ComputeRHEOInterface *compute_interface;
  class NeighList *list;

  void grow_arrays(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
