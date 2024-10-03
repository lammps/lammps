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
ComputeStyle(RHEO/KERNEL,ComputeRHEOKernel)
// clang-format on
#else

#ifndef LMP_COMPUTE_RHEO_KERNEL_H
#define LMP_COMPUTE_RHEO_KERNEL_H

#include "compute.h"
#include <unordered_set>

namespace LAMMPS_NS {

class ComputeRHEOKernel : public Compute {
 public:
  ComputeRHEOKernel(class LAMMPS *, int, char **);
  ~ComputeRHEOKernel() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double memory_usage() override;
  void compute_coordination();
  double calc_w_self();
  double calc_w(int, int, double, double, double, double);
  double calc_dw(int, int, double, double, double, double);
  double calc_w_quintic(double);
  double calc_dw_quintic(double, double, double, double, double *, double *);
  double calc_w_wendlandc4(double);
  double calc_dw_wendlandc4(double, double, double, double, double *, double *);
  void grow_arrays(int);

  double dWij[3], dWji[3], Wij, Wji;
  int correction_order;
  int *coordination;
  class FixRHEO *fix_rheo;

 private:
  int comm_stage, comm_forward_save;
  int interface_flag;
  int lapack_error_flag;
  std::unordered_set<tagint> lapack_error_tags;

  int corrections_calculated;
  int kernel_style, zmin, dim, Mdim, ncor;
  int nmax_store;
  double cut, cutsq, cutinv, cutsqinv, pre_w, pre_wp;
  double ***C;
  double *C0;

  class NeighList *list;
  class ComputeRHEOInterface *compute_interface;

  int check_corrections(int);

  double calc_w_rk0(int, int, double);
  double calc_w_rk1(int, int, double, double, double, double);
  double calc_w_rk2(int, int, double, double, double, double);
  void calc_dw_rk1(int, double, double, double, double, double *);
  void calc_dw_rk2(int, double, double, double, double, double *);
};
}    // namespace LAMMPS_NS
#endif
#endif
