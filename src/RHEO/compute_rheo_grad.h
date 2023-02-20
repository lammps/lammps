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
ComputeStyle(rheo/grad,ComputeRHEOGrad)
// clang-format on
#else

#ifndef LMP_COMPUTE_RHEO_GRAD_H
#define LMP_COMPUTE_RHEO_GRAD_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRHEOGrad : public Compute {
 public:
  ComputeRHEOGrad(class LAMMPS *, int, char **);
  ~ComputeRHEOGrad();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void forward_gradients();
  void forward_fields();
  double **gradv;
  double **gradr;
  double **gradt;
  double **gradn;
  int stage;

 private:
  int dim, comm_stage;
  int ncomm_grad, ncomm_field;
  double cut, cutsq, rho0;
  class NeighList *list;

  class FixRHEO *fix_rheo;
  class ComputeRHEOKernel *compute_kernel;
  class ComputeRHEOInterface *compute_interface;

  int velocity_flag, temperature_flag, rho_flag, eta_flag;
};

}    // namespace LAMMPS_NS

#endif
#endif
